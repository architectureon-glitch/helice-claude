"""Enchainement complet des phases 1 a 6 sur un fichier de geometrie.

Ce module ne fait aucun calcul propre : il appelle les phases dans l'ordre,
propage les niveaux de confiance et rassemble les avertissements.  La CLI et les
rapports ne connaissent que sa sortie.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from . import config
from .confidence import LOW, ConfidenceMap, worst
from .geometry import axis as axis_module
from .geometry import blade_angles as blade_module
from .geometry import blade_loops as loops_module
from .geometry import blade_normals as normals_module
from .geometry import occupancy as occupancy_module
from .geometry import sections as sections_module
from .geometry import topology as topology_module
from .hydraulics import cavitation as cavitation_module
from .hydraulics import meanline as meanline_module
from .hydraulics import similarity as similarity_module
from .io import loader
from .mesh import TriMesh


@dataclass
class Options:
    """Toutes les entrees utilisateur de la CLI (SPEC phase 7)."""

    unit: str | float | None = None
    r_aspiration_cm: float | None = None
    speeds: tuple[float, ...] = config.DEFAULT_RPM
    blades: int | None = None
    beta1_deg: float | None = None
    beta2_deg: float | None = None
    altitude: float = config.ALTITUDE
    temperature_c: float = config.TEMPERATURE
    suction_height: float = config.HAUTEUR_ASPIRATION
    suction_losses: float = config.PERTES_ASPIRATION
    suction: str = axis_module.SUCTION_AUTO
    grid_nr: int = config.GRID_NR
    grid_nz: int = config.GRID_NZ
    n_theta: int = config.N_THETA
    repair: bool = True
    prefer_trimesh: bool = True
    symmetry_check: bool = True

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "unite": self.unit if self.unit is not None else "cm",
            "r_aspiration_cm": self.r_aspiration_cm,
            "regimes_tr_min": list(self.speeds),
            "pales_imposees": self.blades,
            "beta1_impose_deg": self.beta1_deg,
            "beta2_impose_deg": self.beta2_deg,
            "cote_aspiration": self.suction,
            "altitude_m": self.altitude,
            "temperature_C": self.temperature_c,
            "hauteur_aspiration_m": self.suction_height,
            "pertes_aspiration_m": self.suction_losses,
            "grille": {"nr": self.grid_nr, "nz": self.grid_nz, "n_theta": self.n_theta},
            "reparation": self.repair,
            "controle_de_symetrie": self.symmetry_check,
        }


@dataclass
class AnalysisResult:
    """Resultat complet d'une analyse, tout en SI."""

    options: Options = field(default_factory=Options)
    import_report: loader.ImportReport | None = None
    axis: axis_module.AxisResult | None = None
    suction: axis_module.SuctionResult | None = None
    occupancy: occupancy_module.OccupancyMap | None = None
    topology: topology_module.Topology | None = None
    blades: blade_module.BladeGeometry | None = None
    blade_loops: loops_module.LoopResult | None = None
    blade_normals: normals_module.NormalAngles | None = None
    meanline_input: meanline_module.MeanlineInput | None = None
    curves: list[meanline_module.PerformanceCurve] = field(default_factory=list)
    installation: cavitation_module.Installation | None = None
    speed_limit: cavitation_module.SpeedLimit | None = None
    similarity: similarity_module.SimilarityCheck | None = None
    discharge: dict = field(default_factory=dict)
    mesh: TriMesh | None = None
    warnings: list[str] = field(default_factory=list)
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    elapsed_s: float = 0.0

    def overall_confidence(self) -> str:
        """Niveau de confiance global de l'analyse."""
        return self.confidence.overall()

    def to_dict(self) -> dict:
        """Vue serialisable en JSON, toutes grandeurs en SI (SPEC 7.1)."""
        return {
            "version": __import__("impeller_analyzer").__version__,
            "entrees": self.options.to_dict(),
            "import": self.import_report.to_dict() if self.import_report else None,
            "axe": self.axis.to_dict() if self.axis else None,
            "cote_aspiration": self.suction.to_dict() if self.suction else None,
            "carte_d_occupation": self.occupancy.summary() if self.occupancy else None,
            "topologie": self.topology.to_dict() if self.topology else None,
            "forme_des_aubes": self.blade_loops.to_dict() if self.blade_loops else None,
            "angles_par_normales": self.blade_normals.to_dict() if self.blade_normals else None,
            "pales": self.blades.to_dict() if self.blades else None,
            "sens_de_sortie_du_liquide": self.discharge,
            "installation": self.installation.to_dict() if self.installation else None,
            "regimes": [curve.to_dict() for curve in self.curves],
            "similitude": self.similarity.to_dict() if self.similarity else None,
            "vitesse_maximale": self.speed_limit.to_dict() if self.speed_limit else None,
            "confiance": dict(self.confidence),
            "confiance_globale": self.overall_confidence(),
            "avertissements": list(self.warnings),
            "incertitude_du_modele": {
                "hauteur": config.UNCERTAINTY_H,
                "debit": config.UNCERTAINTY_Q,
                "npshr": config.UNCERTAINTY_NPSH,
                "remarque": (
                    "modele 1D ligne moyenne : a verifier par essai sur banc avant toute "
                    "decision d'achat ou de dimensionnement"
                ),
            },
            "duree_s": self.elapsed_s,
        }


def run(path: str, options: Options | None = None) -> AnalysisResult:
    """Analyse complete d'un fichier de geometrie."""
    options = options or Options()
    started = time.time()
    result = AnalysisResult(options=options)

    # Phase 1 : import.
    mesh, import_report = loader.load_mesh(
        path, unit=options.unit, do_repair=options.repair, prefer_trimesh=options.prefer_trimesh
    )
    result.import_report = import_report
    result.warnings.extend(import_report.warnings)
    result.confidence.update(import_report.confidence)

    # Phase 2 : axe, recentrage, carte d'occupation.
    aligned, axis_result = axis_module.align_to_z(mesh)
    result.mesh = aligned
    result.axis = axis_result
    result.warnings.extend(axis_result.warnings)
    result.confidence.set("axe", axis_result.confidence)
    occupancy = occupancy_module.build_occupancy(
        aligned, nr=options.grid_nr, nz=options.grid_nz, n_theta=options.n_theta
    )

    # La convention impose l'aspiration vers +Z, mais un fichier de CAO n'a
    # aucune raison de la respecter : une roue exportee a l'envers serait lue
    # depuis son cote refoulement, passerait pour axiale et recevrait un sens de
    # sortie faux. Sur une roue a composante radiale la geometrie tranche seule.
    suction = axis_module.detect_suction_side(occupancy, options.suction)
    result.suction = suction
    if suction.flipped:
        aligned = axis_module.flip_axis(aligned)
        result.mesh = aligned
        occupancy = occupancy_module.build_occupancy(
            aligned, nr=options.grid_nr, nz=options.grid_nz, n_theta=options.n_theta
        )
    result.warnings.extend(suction.warnings)
    result.confidence.set("cote_aspiration", suction.confidence)
    result.occupancy = occupancy

    # Phase 3 : topologie.
    topology = topology_module.analyse(
        occupancy,
        aligned if options.symmetry_check else None,
        forced_blades=options.blades,
        r_aspiration_cm=options.r_aspiration_cm,
    )
    result.topology = topology
    result.warnings.extend(topology.warnings)
    result.warnings.extend(topology.blades.warnings)
    result.confidence.update(topology.confidence)

    # Forme des aubes : une aube qui se referme sur elle-meme invalide tout ce
    # que la phase 4 en tirerait, la coupe la traversant deux fois.
    loops = loops_module.detect_looped_blades(occupancy, topology.blades.n_blades)
    result.blade_loops = loops
    result.warnings.extend(loops.warnings)

    # Phase 4 : coupes, angles de pale, sens de rotation.
    sections = sections_module.extract_sections(aligned, occupancy, topology)
    geometry = blade_module.analyse(
        sections, topology, forced_beta1_deg=options.beta1_deg, forced_beta2_deg=options.beta2_deg
    )
    result.blades = geometry
    result.confidence.update(geometry.confidence)
    camber_warnings = list(geometry.warnings)

    # Aube en boucle : la cambrure n'a pas de reponse stable, les normales si.
    if loops.looped and options.beta1_deg is None and options.beta2_deg is None:
        normals = normals_module.analyse(aligned, occupancy, topology)
        result.blade_normals = normals
        if normals.beta1_deg > 0.0 and normals.beta2_deg > 0.0:
            geometry.beta1_deg = normals.beta1_deg
            geometry.beta2_deg = normals.beta2_deg
            geometry.rotation_sign = blade_module.rotation_sense(
                topology.machine_type, normals.slope_sign
            )
            geometry.rotation_label = blade_module.rotation_label(geometry.rotation_sign)
            geometry.notes.extend(normals.notes)
            geometry.warnings.append(
                f"aubes en boucle : beta1 = {normals.beta1_deg:.1f} deg et beta2 = "
                f"{normals.beta2_deg:.1f} deg sont lus sur les normales de la surface d'aube, "
                "la cambrure n'ayant pas de reponse stable sur cette forme. Methode basse de 2 a "
                "5 degres sur des roues d'angles connus, biais non corrige ; imposez --beta1 et "
                "--beta2 si vous les connaissez."
            )
            for quantity in ("angles_de_pale", "sens_de_rotation"):
                result.confidence.set(quantity, normals.confidence)
        else:
            result.warnings.append(
                "aubes en boucle et lecture par les normales infructueuse : les angles de pale "
                "restent ceux de la cambrure, qui ne s'applique pas a cette forme. Imposez "
                "--beta1 et --beta2."
            )
            for quantity in loops_module.INVALIDATED_BY_LOOP:
                result.confidence.set(quantity, LOW)

    if result.blade_normals is not None and result.blade_normals.beta2_deg > 0.0 and camber_warnings:
        # Les reserves de la cambrure portent sur une lecture qui n'a pas ete
        # retenue : les garder sans le dire ferait croire a un doute sur les
        # angles publies.
        result.warnings.append(
            "les reserves qui suivent portent sur la lecture par la cambrure, mise de cote au "
            "profit des normales ; elles n'entament pas les angles publies, elles expliquent "
            "pourquoi la cambrure a ete ecartee"
        )
    result.warnings.extend(camber_warnings)

    # Phases 5 et 6 : hydraulique et cavitation.
    data = meanline_module.MeanlineInput.from_geometry(topology, geometry)
    result.meanline_input = data
    hydraulic_confidence = meanline_module.confidence_of(topology, geometry)
    result.confidence.update(hydraulic_confidence)

    site = cavitation_module.installation(
        altitude=options.altitude,
        temperature_c=options.temperature_c,
        suction_height=options.suction_height,
        suction_losses=options.suction_losses,
    )
    result.installation = site

    curves = []
    for rpm in options.speeds:
        curve = meanline_module.build_curve(data, rpm, hydraulic_confidence)
        cavitation_module.apply_to_curve(curve)
        curves.append(curve)
        for message in curve.warnings:
            if message not in result.warnings:
                result.warnings.append(message)
    result.curves = curves
    result.confidence.set("npshr", worst(hydraulic_confidence.overall(), "medium"))

    if curves:
        reference = curves[0]
        result.speed_limit = cavitation_module.maximum_speed(reference, site, hydraulic_confidence)
        result.warnings.extend(result.speed_limit.warnings)
        if len(curves) > 1:
            result.similarity = similarity_module.check(curves[0], curves[-1])
            result.warnings.extend(result.similarity.warnings)
        # Controle croise du type de roue par la vitesse specifique (SPEC 3.3).
        message = topology_module.cross_check_type(topology, reference.specific_speed)
        if message:
            result.warnings.insert(0, message)
        best = reference.best_efficiency_point()
        if best is not None:
            result.discharge = blade_module.discharge_direction(
                topology, geometry, cm2=best.cm2, cu2=best.cu2
            )
    if not result.discharge:
        result.discharge = blade_module.discharge_direction(topology, geometry)


    result.elapsed_s = time.time() - started
    return result
