"""Campagne de validation (SPEC phase 8).

Quatre familles de controles :

1. geometries synthetiques a reponse analytique connue ;
2. cas de reference exterieur (courbe constructeur publiee) ;
3. invariance par rotation et translation ;
4. robustesse a la decimation du maillage.

Chaque controle renvoie des ecarts chiffres, comparables a une tolerance de
`config` ; rien n'est cache derriere un booleen.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field

from . import config, synthetic
from .geometry import axis as axis_module
from .geometry import blade_angles as blade_module
from .geometry import occupancy as occupancy_module
from .geometry import sections as sections_module
from .geometry import topology as topology_module
from .hydraulics import cavitation as cavitation_module
from .hydraulics import meanline as meanline_module
from .mesh import TriMesh, rotation_matrix


@dataclass
class Check:
    """Un ecart mesure, avec sa tolerance et son verdict."""

    name: str
    measured: float
    expected: float
    tolerance: float
    unit: str = ""
    relative: bool = True

    @property
    def deviation(self) -> float:
        """Ecart, relatif ou absolu selon le controle."""
        if not self.relative:
            return abs(self.measured - self.expected)
        if self.expected == 0.0:
            return 0.0 if self.measured == 0.0 else math.inf
        return abs(self.measured - self.expected) / abs(self.expected)

    @property
    def passed(self) -> bool:
        """Vrai si l'ecart tient dans la tolerance."""
        return self.deviation <= self.tolerance

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "grandeur": self.name,
            "mesure": self.measured,
            "attendu": self.expected,
            "unite": self.unit,
            "ecart": self.deviation,
            "tolerance": self.tolerance,
            "relatif": self.relative,
            "conforme": self.passed,
        }

    def line(self) -> str:
        """Ligne lisible pour la console."""
        mark = "ok  " if self.passed else "ECHEC"
        scale = "%" if self.relative else self.unit
        value = self.deviation * (100.0 if self.relative else 1.0)
        limit = self.tolerance * (100.0 if self.relative else 1.0)
        return (
            f"  {mark} {self.name:<28} mesure {self.measured:>10.4g} {self.unit:<5} "
            f"attendu {self.expected:>10.4g}   ecart {value:6.2f} {scale} (max {limit:g})"
        )


@dataclass
class Report:
    """Resultat d'un controle de validation."""

    title: str
    checks: list[Check] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        """Vrai si tous les ecarts tiennent dans leur tolerance."""
        return all(check.passed for check in self.checks)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "controle": self.title,
            "conforme": self.passed,
            "ecarts": [check.to_dict() for check in self.checks],
            "remarques": list(self.notes),
        }

    def render(self) -> str:
        """Rendu console."""
        lines = [f"{'OK ' if self.passed else 'ECHEC'} - {self.title}"]
        lines.extend(check.line() for check in self.checks)
        lines.extend(f"  note : {note}" for note in self.notes)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Chainage commun
# ---------------------------------------------------------------------------
def extract(mesh: TriMesh, nr: int = config.GRID_NR, nz: int = config.GRID_NZ):
    """Phases 2 a 4 sur un maillage deja en metres."""
    aligned, _ = axis_module.align_to_z(mesh)
    occupancy = occupancy_module.build_occupancy(aligned, nr=nr, nz=nz)
    topology = topology_module.analyse(occupancy, None)
    sections = sections_module.extract_sections(aligned, occupancy, topology)
    geometry = blade_module.analyse(sections, topology)
    return topology, geometry


# ---------------------------------------------------------------------------
# 8.1 Geometries synthetiques
# ---------------------------------------------------------------------------
def check_synthetic(title: str, mesh: TriMesh, expected: dict, nr: int = config.GRID_NR, nz: int = config.GRID_NZ) -> Report:
    """Compare l'extraction aux valeurs imposees a la construction (SPEC 8.1)."""
    report = Report(title=title)
    topology, geometry = extract(mesh, nr, nz)
    report.checks.append(
        Check("nombre de pales", float(topology.blades.n_blades), float(expected["n_blades"]), 0.0)
    )
    for key, measured, unit in (
        ("r_1s", topology.r_1s, "m"),
        ("r_2", topology.r_2, "m"),
    ):
        if key in expected:
            report.checks.append(
                Check(key, measured, expected[key], config.VALID_GEOM_TOL, unit)
            )
    for key, measured in (("beta1_deg", geometry.beta1_deg), ("beta2_deg", geometry.beta2_deg)):
        if key in expected:
            report.checks.append(
                Check(key, measured, expected[key], config.VALID_BETA_DEG, "deg", relative=False)
            )
    if "machine_type" in expected:
        report.notes.append(
            f"type de roue : {topology.machine_type} (attendu {expected['machine_type']})"
        )
        report.checks.append(
            Check(
                "type de roue",
                1.0 if topology.machine_type == expected["machine_type"] else 0.0,
                1.0,
                0.0,
            )
        )
    if "rotation_sign" in expected:
        report.checks.append(
            # La campagne controle la **lecture geometrique** : le sens retenu,
            # lui, vient de l'utilisateur et n'a rien a valider.
            Check(
                "sens de rotation suggere",
                float(geometry.observed_rotation_sign),
                float(expected["rotation_sign"]),
                0.0,
            )
        )
    report.notes.append(
        "l'ecart sur beta est exprime en degres : une tolerance relative de 2 % n'a pas de "
        "sens sur un angle qui varie fortement le long de la pale"
    )
    return report


# ---------------------------------------------------------------------------
# 8.2 Cas de reference exterieur
# ---------------------------------------------------------------------------
def load_reference_case(path: str) -> dict:
    """Charge un cas de reference (geometrie + point de fonctionnement publie)."""
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def check_reference_case(case: dict) -> Report:
    """Compare la hauteur calculee au BEP a la valeur publiee (SPEC 8.2)."""
    report = Report(title=f"cas de reference : {case.get('designation', 'sans nom')}")
    geometry = case["geometrie"]
    data = meanline_module.MeanlineInput(
        r_1=geometry["r_1_m"],
        r_1s=geometry["r_1s_m"],
        r_1h=geometry["r_1h_m"],
        r_2=geometry["r_2_m"],
        area_1=math.pi * (geometry["r_1s_m"] ** 2 - geometry["r_1h_m"] ** 2) * config.TAU_1,
        area_2=2.0 * math.pi * geometry["r_2_m"] * geometry["b_2_m"] * config.TAU_2,
        beta1_deg=geometry["beta1_deg"],
        beta2_deg=geometry["beta2_deg"],
        n_blades=geometry["n_blades"],
    )
    duty = case["point_publie"]
    curve = cavitation_module.apply_to_curve(
        meanline_module.build_curve(data, duty["regime_tr_min"])
    )
    best = curve.best_efficiency_point()
    if best is None:
        report.notes.append("le modele ne produit pas de point exploitable pour cette geometrie")
        report.checks.append(Check("hauteur au BEP", 0.0, duty["H_m"], config.VALID_REFERENCE_TOL, "m"))
        return report
    report.checks.append(
        Check("hauteur au BEP", best.head, duty["H_m"], config.VALID_REFERENCE_TOL, "m")
    )
    report.notes.append(
        f"debit calcule au BEP {best.flow * config.SECONDS_PER_HOUR:.1f} m3/h, "
        f"publie {duty['Q_m3_h']:.1f} m3/h"
    )
    if "source" in case:
        report.notes.append(f"source des donnees publiees : {case['source']}")
    return report


# ---------------------------------------------------------------------------
# 8.3 Invariance
# ---------------------------------------------------------------------------
def check_invariance(
    mesh: TriMesh,
    angle_deg: float = 37.0,
    offset=(0.13, -0.07, 0.21),
    nr: int = config.GRID_NR,
    nz: int = config.GRID_NZ,
) -> Report:
    """Le meme maillage tourne et translate doit donner le meme resultat (SPEC 8.3)."""
    report = Report(title=f"invariance par rotation de {angle_deg:g} deg et translation")
    moved = mesh.transformed(
        matrix=rotation_matrix((0.0, 0.0, 1.0), math.radians(angle_deg)), translation=offset
    )
    reference_topology, reference_geometry = extract(mesh, nr, nz)
    moved_topology, moved_geometry = extract(moved, nr, nz)
    for name, first, second, unit in (
        ("nombre de pales", float(reference_topology.blades.n_blades), float(moved_topology.blades.n_blades), ""),
        ("r_1s", reference_topology.r_1s, moved_topology.r_1s, "m"),
        ("r_2", reference_topology.r_2, moved_topology.r_2, "m"),
        ("beta1", reference_geometry.beta1_deg, moved_geometry.beta1_deg, "deg"),
        ("beta2", reference_geometry.beta2_deg, moved_geometry.beta2_deg, "deg"),
    ):
        tolerance = 0.0 if name == "nombre de pales" else config.VALID_INVARIANCE_TOL
        report.checks.append(Check(name, second, first, tolerance, unit))
    report.checks.append(
        Check(
            "sens de rotation",
            float(moved_geometry.observed_rotation_sign),
            float(reference_geometry.observed_rotation_sign),
            0.0,
        )
    )
    return report


# ---------------------------------------------------------------------------
# 8.4 Robustesse a la decimation
# ---------------------------------------------------------------------------
def check_robustness(
    mesh: TriMesh,
    ratio: float = config.VALID_DECIMATION_RATIO,
    nr: int = config.GRID_NR,
    nz: int = config.GRID_NZ,
) -> Report:
    """Le maillage decime doit donner le meme beta2 a quelques degres pres (SPEC 8.4)."""
    report = Report(title=f"robustesse a une decimation a {ratio * 100.0:.0f} % des triangles")
    coarse = mesh.decimate(ratio)
    report.notes.append(f"{len(mesh.faces)} triangles reduits a {len(coarse.faces)}")
    reference_topology, reference_geometry = extract(mesh, nr, nz)
    coarse_topology, coarse_geometry = extract(coarse, nr, nz)
    report.checks.append(
        Check(
            "nombre de pales",
            float(coarse_topology.blades.n_blades),
            float(reference_topology.blades.n_blades),
            0.0,
        )
    )
    report.checks.append(
        Check(
            "beta2",
            coarse_geometry.beta2_deg,
            reference_geometry.beta2_deg,
            config.VALID_DECIMATION_DEG,
            "deg",
            relative=False,
        )
    )
    return report


# ---------------------------------------------------------------------------
# Campagne complete
# ---------------------------------------------------------------------------
def default_cases(nr: int = config.GRID_NR, nz: int = config.GRID_NZ) -> list[Report]:
    """Campagne par defaut : geometries synthetiques, invariance, robustesse."""
    axial = synthetic.axial_impeller(
        n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0, r_hub=0.030, r_tip=0.100
    )
    centrifugal = synthetic.centrifugal_impeller(
        n_blades=6, beta1_deg=22.0, beta2_deg=25.0, r1=0.035, r2=0.090
    )
    reports = [
        check_synthetic(
            "helice a pales planes",
            synthetic.axial_impeller(n_blades=4, beta_deg=30.0, r_hub=0.030, r_tip=0.100),
            {
                "n_blades": 4,
                "r_1s": 0.100,
                "r_2": None,
                "beta1_deg": 30.0,
                "beta2_deg": 30.0,
                "machine_type": topology_module.AXIAL,
                "rotation_sign": 1,
            },
            nr,
            nz,
        ),
        check_synthetic(
            "helice a angle variable",
            axial,
            {
                "n_blades": 4,
                "r_1s": 0.100,
                "beta1_deg": 20.0,
                "beta2_deg": 35.0,
                "machine_type": topology_module.AXIAL,
                "rotation_sign": 1,
            },
            nr,
            nz,
        ),
        check_synthetic(
            "roue centrifuge a aubes en arc",
            centrifugal,
            {
                "n_blades": 6,
                "r_1s": 0.035,
                "r_2": 0.090,
                "beta1_deg": 22.0,
                "beta2_deg": 25.0,
                "machine_type": topology_module.CENTRIFUGAL,
                "rotation_sign": -1,
            },
            nr,
            nz,
        ),
        check_invariance(centrifugal, nr=nr, nz=nz),
        check_robustness(centrifugal, nr=nr, nz=nz),
    ]
    # Le champ r_2 vaut None pour l'helice a pales planes : on retire le controle.
    reports[0].checks = [check for check in reports[0].checks if check.expected is not None]
    return reports


def main(argv: list[str] | None = None) -> int:  # pragma: no cover - point d'entree
    """Lance la campagne et affiche le compte rendu."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m impeller_analyzer.validation",
        description="Campagne de validation de impeller-analyzer (SPEC phase 8).",
    )
    parser.add_argument("--cas-de-reference", default=None, metavar="JSON",
                        help="fichier de cas de reference (geometrie + point publie)")
    parser.add_argument("--grille", type=int, nargs=2, default=[config.GRID_NR, config.GRID_NZ],
                        metavar=("NR", "NZ"), help="taille de la carte d'occupation")
    parser.add_argument("--json", default=None, metavar="FICHIER", help="ecrire le compte rendu en JSON")
    args = parser.parse_args(argv)

    reports = default_cases(args.grille[0], args.grille[1])
    if args.cas_de_reference and os.path.isfile(args.cas_de_reference):
        reports.append(check_reference_case(load_reference_case(args.cas_de_reference)))
    for report in reports:
        print(report.render())
        print()
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump([report.to_dict() for report in reports], handle, ensure_ascii=False, indent=2)
    failed = [report for report in reports if not report.passed]
    print(f"{len(reports) - len(failed)} / {len(reports)} controles conformes")
    return 0 if not failed else 1


if __name__ == "__main__":  # pragma: no cover - point d'entree
    import sys

    sys.exit(main())
