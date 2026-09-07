"""Sorties de l'outil : `resultats.json`, `rapport.md`, `courbes.png` (SPEC phase 7)."""

from __future__ import annotations

import json
import os

from .. import config
from ..analysis import AnalysisResult
from . import plot, viewer

JSON_NAME = "resultats.json"
MARKDOWN_NAME = "rapport.md"
CURVES_NAME = "courbes.png"
OCCUPANCY_NAME = "carte_occupation.png"
VIEWER_NAME = "vue3d.html"

#: Traduction des niveaux de confiance pour le rapport lisible.
CONFIDENCE_LABELS = {"high": "haute", "medium": "moyenne", "low": "faible"}


def _mm(value: float | None) -> str:
    """Longueur en millimetres, deux decimales."""
    return "-" if value is None else f"{value * config.MM_PER_M:.2f}"


def _deg(value: float | None) -> str:
    """Angle en degres, une decimale."""
    return "-" if value is None else f"{value:.1f}"


def _confidence(level: str) -> str:
    """Niveau de confiance en francais."""
    return CONFIDENCE_LABELS.get(level, level)


def write_json(result: AnalysisResult, directory: str) -> str:
    """Ecrit `resultats.json` : toutes les grandeurs, en SI, avec les confiances."""
    path = os.path.join(directory, JSON_NAME)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(result.to_dict(), handle, ensure_ascii=False, indent=2, allow_nan=False, default=_fallback)
    return path


def _fallback(value):
    """Serialisation des valeurs non standard (infinis, objets residuels)."""
    return str(value)


def geometry_table(result: AnalysisResult) -> list[tuple[str, str, str]]:
    """Tableau 1 de la SPEC : geometrie extraite."""
    topology = result.topology
    blades = result.blades
    if topology is None or blades is None:
        return []
    confidence = result.confidence
    discharge = result.discharge or {}
    alpha = discharge.get("alpha2_deg")
    outlet = discharge.get("composante_meridienne", "-")
    if alpha is not None:
        outlet = f"{outlet}, alpha2 = {alpha:.1f} deg"
    suffix = " (impose)" if topology.r_aspiration_source == "utilisateur" else ""
    return [
        ("Type de roue", topology.machine_type, _confidence(confidence.get_level("type_de_roue"))),
        (
            "Nombre de pales",
            str(topology.blades.n_blades) + (" (impose)" if topology.blades.forced else ""),
            _confidence(confidence.get_level("nombre_de_pales")),
        ),
        ("Rayon d'aspiration r1s (mm)", _mm(topology.r_1s) + suffix, _confidence(confidence.get_level("rayons"))),
        ("Rayon moyeu r1h (mm)", _mm(topology.r_1h), _confidence(confidence.get_level("rayons"))),
        ("Rayon de sortie r2 (mm)", _mm(topology.r_2), _confidence(confidence.get_level("rayons"))),
        (
            "beta1 / beta2 au rayon moyen (deg)",
            f"{_deg(blades.beta1_deg)} / {_deg(blades.beta2_deg)}"
            + (" (imposes)" if blades.forced_beta else ""),
            _confidence(confidence.get_level("angles_de_pale")),
        ),
        (
            "Sens de rotation",
            blades.rotation_label or "indetermine",
            _confidence(confidence.get_level("sens_de_rotation")),
        ),
        ("Sens de sortie du liquide", outlet, _confidence(confidence.get_level("sens_de_rotation"))),
    ]


def performance_table(result: AnalysisResult) -> tuple[list[str], list[tuple[str, list[str]]]]:
    """Tableau 2 de la SPEC : performances, une colonne par regime."""
    header = [f"{curve.rpm:.0f} tr/min" for curve in result.curves]
    npsha = result.installation.npsha if result.installation else 0.0
    rows: list[tuple[str, list[str]]] = []

    def collect(label: str, formatter) -> None:
        values = []
        for curve in result.curves:
            point = curve.nominal_point()
            values.append("-" if point is None else formatter(curve, point))
        rows.append((label, values))

    collect("Debit nominal (m3/h)", lambda c, p: f"{p.flow * config.SECONDS_PER_HOUR:.1f}")
    collect("Hauteur (m)", lambda c, p: f"{p.head:.2f}")
    collect("Puissance arbre (kW)", lambda c, p: f"{p.shaft_power / config.W_PER_KW:.3f}")
    collect("Couple (N.m)", lambda c, p: f"{p.torque:.2f}")
    collect("Rendement estime (%)", lambda c, p: f"{p.efficiency * 100.0:.1f}")
    collect("Vitesse specifique n_q", lambda c, p: f"{c.specific_speed:.1f}")
    collect("NPSH requis (m)", lambda c, p: f"{p.npshr:.2f}")
    collect("NPSH disponible (m)", lambda c, p: f"{npsha:.2f}")
    collect(
        "Marge NPSHa / NPSHr",
        lambda c, p: f"{npsha / p.npshr:.2f}" if p.npshr > 0.0 else "-",
    )
    return header, rows


def _markdown_table(header: list[str], rows: list[list[str]]) -> list[str]:
    """Tableau Markdown a partir d'un en-tete et de lignes deja formatees."""
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return lines


def write_markdown(result: AnalysisResult, directory: str, source: str = "") -> str:
    """Ecrit `rapport.md` : rapport lisible avec les deux tableaux et les encadres."""
    lines: list[str] = []
    name = os.path.basename(source or (result.import_report.path if result.import_report else ""))
    lines.append(f"# Analyse hydraulique de la roue `{name}`")
    lines.append("")
    lines.append(
        f"Convention : axe de rotation Z, aspiration vers +Z, refoulement vers -Z ou radial. "
        f"Confiance globale de l'analyse : **{_confidence(result.overall_confidence())}**."
    )
    lines.append("")

    if result.warnings:
        lines.append("## Avertissements")
        lines.append("")
        for message in result.warnings:
            lines.append(f"- {message}")
        lines.append("")

    lines.append("## Tableau 1 - Geometrie extraite")
    lines.append("")
    table = geometry_table(result)
    if table:
        lines.extend(_markdown_table(["Grandeur", "Valeur", "Confiance"], [list(row) for row in table]))
    else:
        lines.append("_Geometrie non exploitable._")
    lines.append("")

    lines.append("## Tableau 2 - Performances")
    lines.append("")
    if result.blade_loops is not None and result.blade_loops.looped:
        # Une aube en boucle sort du domaine de la SPEC : le tableau est conserve
        # pour ce qu'il vaut, mais il ne doit pas etre lu comme une performance.
        lines.append(
            "> **Non applicable.** Les aubes se referment sur elles-memes (type toroidal) : "
            "le modele de ligne moyenne suppose un bord d'attaque et un bord de fuite uniques, "
            "et les angles de pale dont derive tout ce tableau ont ete lus entre les deux brins "
            "d'une meme boucle. Les valeurs ci-dessous sont celles qu'aurait donnees une aube "
            "simple de meme trace ; **elles ne decrivent pas cette roue** et ne doivent servir "
            "ni au dimensionnement ni a la comparaison."
        )
        lines.append("")
    header, rows = performance_table(result)
    if header:
        lines.extend(
            _markdown_table([""] + header, [[label] + values for label, values in rows])
        )
    else:
        lines.append("_Aucun regime calcule._")
    lines.append("")

    lines.append("## Niveaux de confiance")
    lines.append("")
    lines.append(
        "Chaque grandeur estimee porte son niveau ; le niveau global est le plus faible "
        "de la table."
    )
    lines.append("")
    lines.extend(
        _markdown_table(
            ["Grandeur", "Confiance"],
            [[name, _confidence(level)] for name, level in sorted(result.confidence.items())],
        )
    )
    lines.append("")

    if result.speed_limit is not None:
        limit = result.speed_limit
        lines.append("## Vitesse maximale sans cavitation")
        lines.append("")
        lines.append(f"> **{limit.rpm_max} tr/min**")
        lines.append(">")
        lines.append(f"> Limite active : **{limit.active_limit}**.")
        lines.append(
            f"> Par le NPSH : {limit.rpm_max_npsh:.0f} tr/min (marge imposee {config.MARGE_NPSH:g}). "
            f"Par la vitesse relative d'entree : {limit.rpm_max_w1s:.0f} tr/min "
            f"(w1s <= {config.W1S_MAX:g} m/s)."
        )
        lines.append(f"> Hypotheses d'installation : {limit.assumptions}.")
        lines.append("")

    if result.similarity is not None:
        check = result.similarity
        verdict = "conforme" if check.passed else "NON CONFORME"
        lines.append("## Controle de similitude")
        lines.append("")
        lines.append(
            f"Ecart entre le calcul direct a {check.rpm_target:.0f} tr/min et l'extrapolation "
            f"depuis {check.rpm_reference:.0f} tr/min : **{check.worst_deviation * 100.0:.2f} %** "
            f"(seuil {config.SIMILARITY_TOL * 100.0:.0f} %) - {verdict}."
        )
        lines.append("")

    if result.blades is not None and result.blades.notes:
        lines.append("## Remarques sur le sens de rotation")
        lines.append("")
        for note in result.blades.notes:
            lines.append(f"- {note}")
        lines.append("")

    lines.append("## Incertitude")
    lines.append("")
    lines.append(
        f"> Modele 1D ligne moyenne : hauteur +/- {config.UNCERTAINTY_H * 100.0:.0f} %, "
        f"debit +/- {config.UNCERTAINTY_Q * 100.0:.0f} %, "
        f"NPSHr +/- {config.UNCERTAINTY_NPSH * 100.0:.0f} %."
    )
    lines.append(">")
    lines.append(
        "> **A verifier par essai sur banc avant toute decision d'achat ou de dimensionnement.**"
    )
    lines.append("")

    if result.import_report is not None:
        report = result.import_report
        lines.append("## Import")
        lines.append("")
        lines.append(
            f"- {report.n_faces} triangles, {report.n_vertices} sommets "
            f"(lecteur `{report.backend}`, unite `{report.unit}`)"
        )
        lines.append(
            f"- boite englobante : {report.extents_mm[0]:.1f} x {report.extents_mm[1]:.1f} x "
            f"{report.extents_mm[2]:.1f} mm"
        )
        lines.append(f"- volume : {report.volume_m3 * config.CM3_PER_M3:.1f} cm3, etanche : {report.watertight}")
        if result.suction is not None:
            side = {1: "vers +Z, conforme a la convention", -1: "vers -Z : le maillage a ete retourne"}.get(
                result.suction.sign, "indetermine, convention +Z conservee"
            )
            lines.append(
                f"- cote aspiration : {side} (asymetrie radiale de la veine "
                f"{result.suction.asymmetry:+.3f})"
            )
        lines.append(f"- duree d'analyse : {result.elapsed_s:.1f} s")
        lines.append("")

    path = os.path.join(directory, MARKDOWN_NAME)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    return path


def write_curves(result: AnalysisResult, directory: str) -> str:
    """Ecrit `courbes.png` : H-Q, rendement-Q et NPSHr-Q, une couleur par regime."""
    path = os.path.join(directory, CURVES_NAME)
    panels = []
    for title, y_label, extract in (
        ("Hauteur", "H (m)", lambda p: p.head),
        ("Rendement", "eta (%)", lambda p: p.efficiency * 100.0),
        ("NPSH requis", "NPSHr (m)", lambda p: p.npshr),
    ):
        series = []
        markers = []
        for index, curve in enumerate(result.curves):
            usable = [point for point in curve.points if point.head > 0.0]
            if not usable:
                continue
            color = plot.SERIES_COLORS[index % len(plot.SERIES_COLORS)]
            series.append(
                {
                    "label": f"{curve.rpm:.0f} tr/min",
                    "x": [point.flow * config.SECONDS_PER_HOUR for point in usable],
                    "y": [extract(point) for point in usable],
                    "color": color,
                }
            )
            best = curve.best_efficiency_point()
            if best is not None and best.head > 0.0:
                markers.append((best.flow * config.SECONDS_PER_HOUR, extract(best), color))
        panels.append(
            {
                "titre": f"{title} - le point marque est le meilleur rendement",
                "x_label": "Q (m3/h)",
                "y_label": y_label,
                "series": series,
                "markers": markers,
            }
        )
    if result.installation is not None:
        for panel in panels:
            if panel["y_label"].startswith("NPSHr") and panel["series"]:
                flows = [value for serie in panel["series"] for value in serie["x"]]
                panel["series"].append(
                    {
                        "label": "NPSH disponible",
                        "x": [min(flows), max(flows)],
                        "y": [result.installation.npsha] * 2,
                        "color": plot.GREY,
                    }
                )
    return plot.curves_png(path, panels, title="Courbes caracteristiques")


def write_occupancy(result: AnalysisResult, directory: str) -> str | None:
    """Ecrit la carte `f(r, z)`, indispensable au controle visuel (SPEC 10)."""
    if result.occupancy is None:
        return None
    path = os.path.join(directory, OCCUPANCY_NAME)
    title = "Carte d'occupation f(r, z)"
    if result.topology is not None:
        title += f" - roue {result.topology.machine_type}, {result.topology.blades.n_blades} pales"
    return plot.occupancy_png(result.occupancy, path, title)


def write_all(
    result: AnalysisResult, directory: str, source: str = "", viewer_page: bool = True
) -> dict[str, str]:
    """Produit les trois sorties de la SPEC, la carte d'occupation et la vue 3D."""
    os.makedirs(directory, exist_ok=True)
    produced = {
        "json": write_json(result, directory),
        "markdown": write_markdown(result, directory, source),
        "curves": write_curves(result, directory),
    }
    occupancy = write_occupancy(result, directory)
    if occupancy:
        produced["occupancy"] = occupancy
    if viewer_page and result.mesh is not None:
        produced["viewer"] = viewer.write_page(result, directory, VIEWER_NAME)
    return produced
