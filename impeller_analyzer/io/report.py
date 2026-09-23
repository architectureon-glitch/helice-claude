"""Sorties de l'outil : `resultats.json`, `rapport.md`, `courbes.png` (SPEC phase 7)."""

from __future__ import annotations

import json
import math
import os

from .. import config
from ..analysis import AnalysisResult
from .. import components as components_module
from ..geometry import topology as topology_module
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


def _rotation_cell(blades) -> str:
    """Case du sens de rotation : ce qui est retenu, et ce que la geometrie suggere.

    Le « (impose) » d'autrefois a disparu de la valeur : c'est la colonne de
    provenance qui porte desormais cette information, pour toutes les grandeurs
    et pas seulement celle-ci.
    """
    if blades is None:
        return "-"
    if blades.forced_rotation:
        return blades.rotation_label
    if blades.observed_rotation_sign:
        return f"**{blades.rotation_label}** -- suggere : {blades.observed_rotation_label}"
    return blades.rotation_label or "indetermine"


def geometry_table(result: AnalysisResult) -> list[tuple[str, str, str, str]]:
    """Tableau 1 de la SPEC : geometrie extraite, avec provenance et confiance.

    La provenance repond a une autre question que la confiance : celle-ci dit a
    quel point la valeur est sure, celle-la d'ou elle vient.  Une valeur imposee
    en ligne de commande peut etre parfaitement sure sans rien devoir au
    maillage, et la lire comme une mesure serait se tromper sur ce que l'outil a
    fait.  Le « (impose) » du sens de rotation faisait deja cela pour lui seul.
    """
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
    source = result.provenance
    hub_label = {
        topology_module.HUB_SOLID: "Rayon de moyeu r1h (mm)",
        topology_module.HUB_BORE: "Rayon interieur de matiere r1h (mm)",
    }.get(topology.hub_kind, "Rayon interieur r1h (mm)")

    def row(label: str, value: str, quantity: str) -> tuple[str, str, str, str]:
        """Une ligne du tableau : libelle, valeur, provenance, confiance."""
        return (label, value, source.label(quantity), _confidence(confidence.get_level(quantity)))

    return [
        row("Type de roue", topology.machine_type, "type_de_roue"),
        row("Nombre de pales", str(topology.blades.n_blades), "nombre_de_pales"),
        row("Rayon d'aspiration r1s (mm)", _mm(topology.r_1s), "rayon_d_aspiration"
            if topology.r_aspiration_source == "utilisateur" else "rayons"),
        row(hub_label, _mm(topology.r_1h), "rayon_de_moyeu"),
        row("Nature du centre", topology.hub_kind, "rayon_de_moyeu"),
        row("Rayon de sortie r2 (mm)", _mm(topology.r_2), "rayons"),
        row(
            "beta1 / beta2 au rayon moyen (deg)",
            (f"{_deg(blades.beta1_deg)} / {_deg(blades.beta2_deg)}"
             if blades.beta2_deg > 0.0 else "non lus"),
            "angles_de_pale",
        ),
        row("Sens de rotation", _rotation_cell(blades), "sens_de_rotation"),
        row("Sens de sortie du liquide", outlet, "sens_de_rotation"),
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
        lines.extend(_markdown_table(
            ["Grandeur", "Valeur", "Provenance", "Confiance"],
            [list(row) for row in table],
        ))
    else:
        lines.append("_Geometrie non exploitable._")
    lines.append("")

    lines.append("## Tableau 2 - Performances")
    lines.append("")
    if result.blade_loops is not None and result.blade_loops.looped:
        # Une aube en boucle sort du domaine de la SPEC : le tableau est conserve
        # pour ce qu'il vaut, mais il ne doit pas etre lu comme une performance.
        lines.append(
            "> **A prendre avec des reserves plus larges que d'ordinaire.** Les aubes se "
            "referment sur elles-memes (type toroidal). La cambrure, qui suppose un bord "
            "d'attaque et un bord de fuite uniques, n'a pas de reponse stable sur cette forme : "
            "beta1 et beta2 sont donc lus sur les **normales de la surface d'aube**, methode "
            "basse de 2 a 5 degres sur des roues d'angles connus. beta2 est le point faible : "
            "il est mesure la ou les deux brins fusionnent en un bout massif, qui bloque plus "
            "qu'il ne guide. Si vous connaissez les angles de dessin, imposez-les par --beta1 "
            "et --beta2 : toute la geometrie autour est juste, seuls les angles sont fragiles."
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

    budget = result.energy_budget
    if budget is not None and budget.total != 0.0:
        centri, diffus, cinet = budget.shares()
        lines.append("## D'ou vient la hauteur")
        lines.append("")
        lines.append(
            "La hauteur theorique vaut `H = u2 cu2 / g`. Ce n'est pas un modele mais un "
            "theoreme : il sort de la conservation du moment cinetique et vaut pour n'importe "
            "quelle forme d'aube. La meme quantite se repartit en trois termes, qui disent par "
            "quel mecanisme l'energie passe au fluide -- elle ne cree rien, elle repartit."
        )
        lines.append("")
        lines.extend(_markdown_table(
            ["Terme", "Hauteur (m)", "Part", "Ce qu'il represente"],
            [
                ["centrifuge `(u2^2-u1^2)/2g`", f"{budget.centrifugal:+.2f}", f"{centri:+.0%}",
                 "ne depend que des rayons et du regime, **pas de la forme des aubes**"],
                ["diffusion `(w1^2-w2^2)/2g`", f"{budget.diffusion:+.2f}", f"{diffus:+.0%}",
                 "vitesse relative convertie en pression ; negative, elle **detruit** de la hauteur"],
                ["cinetique `(c2^2-c1^2)/2g`", f"{budget.kinetic:+.2f}", f"{cinet:+.0%}",
                 "sort en vitesse absolue, a recuperer dans la volute"],
                ["**somme**", f"**{budget.total:.2f}**", "100 %",
                 f"controle : `u2 cu2 / g` = {budget.euler:.2f} m, ecart {budget.residual:.1e}"],
            ],
        ))
        lines.append("")
        if budget.diffusion < 0.0:
            lines.append(
                f"> **Le canal accelere l'ecoulement relatif** au lieu de le ralentir : "
                f"w2 = {budget.w2:.1f} m/s contre w1 = {budget.w1:.1f}. Un canal de pompe bien "
                f"dessine fait l'inverse, et convertit cette vitesse en pression statique. Ici le "
                f"terme de diffusion retire {abs(budget.diffusion):.1f} m des "
                f"{budget.centrifugal:.1f} m que l'effet centrifuge apporte."
            )
            lines.append("")

    losses = result.channel_losses
    if losses is not None and losses.efficiency > 0.0:
        lines.append("## Comparaison de conception")
        lines.append("")
        lines.append(
            "Les pertes de la SPEC sont des fractions du point nominal : le rendement du "
            "tableau 2 ne depend donc pratiquement pas de la forme des aubes, et ne sert pas a "
            "comparer deux roues. Les grandeurs ci-dessous, elles, sortent de la geometrie du "
            "canal, et c'est **l'ecart entre deux roues** qui a un sens, non la valeur absolue."
        )
        lines.append("")
        lines.extend(_markdown_table(
            ["Grandeur", "Valeur"],
            [
                ["Longueur developpee du canal (mm)", _mm(losses.length)],
                ["Diametre hydraulique (mm)", _mm(losses.hydraulic_diameter)],
                ["Elancement L / Dh", f"{losses.slenderness:.1f}"],
                ["Surface mouillee (cm2)", f"{losses.wetted_area * 1e4:.0f}"],
                ["Deceleration relative w2 / w1", f"{losses.de_haller:.2f}"],
                ["Perte de frottement de canal (m)", f"{losses.head_friction:.2f}"],
                ["Perte de diffusion (m)", f"{losses.head_diffusion:.2f}"],
                ["**Rendement de comparaison**", f"**{losses.efficiency * 100.0:.1f} %**"],
            ],
        ))
        lines.append("")
        lines.append(
            f"> Une roue centrifuge ordinaire -- six aubes, beta1/beta2 de 22 et 25 degres -- "
            f"donne **{config.ETA_H * config.ETA_VOL * config.ETA_MEC * 100.0:.1f} %** sur la "
            "meme echelle, par construction : c'est sur elle que le calage est fait. Le chiffre "
            "est independant du diametre et du regime, il ne juge que la forme."
        )
        lines.append("")

    if result.speed_limit is not None and not result.speed_limit.computed:
        limit = result.speed_limit
        lines.append("## Vitesse maximale sans cavitation")
        lines.append("")
        motif = (
            limit.warnings[0]
            if limit.warnings and "aucun point nominal" not in limit.warnings[0]
            else "pas de point de fonctionnement nominal au regime de reference"
        )
        lines.append(f"> **Non calculable** : {motif}.")
        lines.append("")
    elif result.speed_limit is not None:
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

    assembly = result.assembly
    if assembly is not None:
        lines.append("## Import par composants declares")
        lines.append("")
        lines.append(
            "Les pieces ont ete importees separement, et l'outil n'a rien infere de ce qui "
            "suit : le modele hydraulique, la topologie de pale et le nombre de pales sont "
            "**declares**, l'axe vient du contrat d'export, les sections viennent des solides "
            "fluide. Il les a verifies."
        )
        lines.append("")
        lines.extend(_markdown_table(
            ["Declaration", "Valeur"],
            [
                ["Modele hydraulique", assembly.declarations.mode],
                ["Topologie de pale", assembly.declarations.blade_topology],
                ["Nombre de pales", str(assembly.declarations.n_blades)],
            ],
        ))
        lines.append("")
        lines.extend(_markdown_table(
            ["Piece", "Triangles", "Volume (cm3)", "Rayons int./ext. (mm)"],
            [
                [
                    piece.slot,
                    str(len(piece.mesh.faces) if piece.mesh else 0),
                    f"{piece.volume * 1e6:.1f}",
                    "{:.1f} / {:.1f}".format(
                        *[r * config.MM_PER_M for r in piece.radial_extent(assembly.axis_origin)]
                    ),
                ]
                for piece in assembly.components.values()
            ],
        ))
        lines.append("")
        if assembly.inlet and assembly.outlet:
            lines.append(
                f"Sections **mesurees** sur les tranches fluide, `volume / epaisseur` : "
                f"A1 = {assembly.inlet.area * 1e4:.2f} cm2, "
                f"A2 = {assembly.outlet.area * 1e4:.2f} cm2. Elles ne dependent plus de la "
                "detection de r1s et r1h, dont la fragilite est la cause premiere des ecarts "
                "de debit. Sens debitant mesure : "
                f"({', '.join(f'{c:+.2f}' for c in assembly.flow_direction)})."
            )
            lines.append("")
        lines.append("### Ce que l'outil a verifie")
        lines.append("")
        lines.extend(_markdown_table(
            ["Controle", "Resultat", "Detail"],
            [
                [check.name,
                 "conforme" if check.passed else ("**BLOQUANT**" if check.blocking else "reserve"),
                 check.detail]
                for check in assembly.checks
            ],
        ))
        lines.append("")
        lines.append(f"> {components_module.EXPORT_CONTRACT}")
        lines.append("")

    sensitivity = result.sensitivity
    if sensitivity is not None and sensitivity.rows:
        lines.append("## Sensibilite des resultats a la geometrie")
        lines.append("")
        lines.append(
            "L'incertitude annoncee par le modele suppose la geometrie juste. Ce tableau mesure "
            "l'autre moitie de la question : ce qu'une erreur de lecture ferait vraiment. Les "
            "deux premieres colonnes se lisent **par degre** d'erreur sur l'angle de pale. La "
            "troisieme est une elasticite sans dimension : les lois de similitude en donnent la "
            "valeur attendue -- 2 pour la hauteur, 3 pour le debit, 2 pour le NPSHr -- et elle "
            f"sert donc aussi de controle. Au-dela de {config.SENSITIVITY_LOW_PER_DEG:.0%} par "
            "degre, la grandeur passe automatiquement en confiance faible."
        )
        lines.append("")

        def _cell(value: float, suffix: str = " %") -> str:
            if not math.isfinite(value):
                return "**non definie**"
            return f"{value * 100.0:.1f}{suffix}" if suffix else f"{value:.2f}"

        lines.extend(_markdown_table(
            ["Grandeur", "par degre de beta1", "par degre de beta2", "elasticite au diametre"],
            [
                [
                    row.quantity + (" **(declassee)**" if row.alarming() else ""),
                    _cell(row.beta1),
                    _cell(row.beta2),
                    _cell(row.diameter, ""),
                ]
                for row in sensitivity.rows
            ],
        ))
        lines.append("")

    propulsion = result.propulsion
    if propulsion is not None and propulsion.point is not None:
        point = propulsion.point
        lines.append("## En helice libre")
        lines.append("")
        lines.append(
            f"La meme piece tournant **non carenee** dans de l'{propulsion.fluid} a "
            f"{propulsion.rpm:.0f} tr/min, avancant a {point.speed:.1f} m/s "
            f"({point.speed * 3.6:.0f} km/h). Ce n'est pas le regime de pompe des tableaux "
            "ci-dessus : une pompe refoule dans une volute, une helice pousse dans un fluide "
            "libre. Poussee et rendement propulsif viennent d'un bilan par element de pale, "
            "dont le calage, la corde, la cambrure et l'epaisseur sont **lus sur la piece**."
        )
        lines.append("")
        lines.extend(_markdown_table(
            ["Grandeur", "Valeur"],
            [
                ["Parametre d'avance J", f"{point.advance_ratio:.3f}"],
                ["Poussee (N)", f"{point.thrust:.0f}"],
                ["Puissance absorbee (kW)", f"{point.power / config.W_PER_KW:.2f}"],
                ["Couple (N.m)", f"{point.torque:.1f}"],
                ["CT", f"{point.thrust_coefficient:.4f}"],
                ["CP", f"{point.power_coefficient:.4f}"],
                ["Vitesse en bout de pale (m/s)", f"{point.tip_speed:.1f}"],
            ] + ([["Mach en bout de pale", f"{point.tip_mach:.2f}"]] if point.tip_mach else []),
        ))
        lines.append("")
        lines.append("### Le rendement, et son plafond")
        lines.append("")
        lines.append(
            "Un rendement propulsif ne se juge pas dans l'absolu : il se juge contre le maximum "
            "que la conservation de la quantite de mouvement autorise pour **cette poussee-la**. "
            "C'est le disque actif ideal de Froude, pertes de profil, de bout de pale et de "
            "giration toutes mises a zero. Aucune helice ne peut le depasser ; l'ecart dit ce "
            "qu'il reste a gagner par le dessin, et rien de plus."
        )
        lines.append("")
        lines.extend(_markdown_table(
            ["", "Valeur"],
            [
                ["Rendement propulsif calcule", f"{point.efficiency * 100.0:.1f} %"],
                ["**Plafond ideal (Froude)**", f"**{point.froude_efficiency * 100.0:.1f} %**"],
                ["Ecart au plafond", f"{point.ceiling_gap() * 100.0:.1f} points"],
            ],
        ))
        lines.append("")
        for note in propulsion.notes:
            lines.append(f"> {note}")
            lines.append(">")
        if propulsion.notes:
            lines.pop()
        lines.append("")

    if result.blades is not None and result.blades.notes:
        lines.append("## Remarques sur les aubes et le sens de rotation")
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
