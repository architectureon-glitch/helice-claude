"""Vue 3D interactive de la roue analysee, en un seul fichier HTML autonome.

La carte d'occupation est un controle en coupe ; elle ne dit rien du sens de
rotation, qui est pourtant le resultat le plus engageant de l'analyse.  Cette
vue rend le maillage en WebGL, colorie chaque facette selon ce que l'outil en a
compris (moyeu ou pale), trace les rayons caracteristiques a leurs vrais plans
et fait tourner la roue dans le sens calcule.

Le fichier produit ne depend de rien : le maillage y est embarque, le rendu
utilise le WebGL du navigateur sans bibliotheque, et la carte meridienne y est
incluse en image.  Il s'ouvre hors ligne et s'envoie par courriel tel quel.
"""

from __future__ import annotations

import base64
import json
import math
import os
import struct

from .. import config
from ..geometry.occupancy import OccupancyMap
from . import style
from ..mesh import TriMesh, cross, normalize, sub
from . import plot

#: Classes de facette rendues dans la vue.
CLASS_SOLID = 0  # moyeu, flasque : f >= F_SOLIDE juste sous la surface
CLASS_BLADE = 1  # pale : F_MATIERE < f < F_SOLIDE, comme la zone de pales de l'analyse
CLASS_OTHER = 2  # non classee

#: Niveaux de confiance en francais, pour le bandeau de la page.
_CONFIDENCE_FR = {"high": "haute", "medium": "moyenne", "low": "faible"}


def classify_faces(mesh: TriMesh, occupancy: OccupancyMap) -> list[int]:
    """Classe chaque facette en sondant la matiere juste sous elle.

    Le centre de la facette est **sur** la surface, donc ambigu ; on descend
    d'une faible profondeur suivant la normale rentrante avant de lire la carte
    d'occupation.  Une facette de moyeu tombe alors dans du plein, une facette
    de pale dans la bande de pales.
    """
    depth = config.VIEWER_CLASS_DEPTH * max(occupancy.dr, occupancy.dz)
    classes = []
    for a, b, c in mesh.triangles():
        normal = normalize(cross(sub(b, a), sub(c, a)))
        point = tuple((a[k] + b[k] + c[k]) / 3.0 - depth * normal[k] for k in range(3))
        value = occupancy.value(math.hypot(point[0], point[1]), point[2])
        if value >= config.F_SOLIDE:
            classes.append(CLASS_SOLID)
        elif value > config.F_MATIERE:
            classes.append(CLASS_BLADE)
        else:
            classes.append(CLASS_OTHER)
    return classes


def _b64_floats(values) -> str:
    """Encode une suite de flottants en Float32 base64."""
    return base64.b64encode(struct.pack(f"<{len(values)}f", *values)).decode("ascii")


def _b64_bytes(values) -> str:
    """Encode une suite d'octets en base64."""
    return base64.b64encode(bytes(values)).decode("ascii")


def _circle(radius: float, height: float, segments: int = config.VIEWER_CIRCLE_SEGMENTS) -> list[float]:
    """Segments de ligne formant un cercle horizontal."""
    points: list[float] = []
    for index in range(segments):
        a = 2.0 * math.pi * index / segments
        b = 2.0 * math.pi * (index + 1) / segments
        points.extend([radius * math.cos(a), radius * math.sin(a), height])
        points.extend([radius * math.cos(b), radius * math.sin(b), height])
    return points


def _arc_arrow(radius: float, height: float, sign: int, span_deg: float = 220.0) -> list[float]:
    """Arc flechee autour de l'axe, dans le sens de rotation demande."""
    if sign == 0:
        return []
    span = math.radians(span_deg)
    steps = 48
    points: list[float] = []
    for index in range(steps):
        a = sign * (span * index / steps)
        b = sign * (span * (index + 1) / steps)
        points.extend([radius * math.cos(a), radius * math.sin(a), height])
        points.extend([radius * math.cos(b), radius * math.sin(b), height])
    # Pointe de la fleche, dans le plan de l'arc.
    tip = sign * span
    tangent = (-sign * math.sin(tip), sign * math.cos(tip))
    radial = (math.cos(tip), math.sin(tip))
    head = radius * 0.18
    end = (radius * radial[0], radius * radial[1], height)
    for side in (1.0, -1.0):
        barb = (
            end[0] - head * tangent[0] + side * head * 0.5 * radial[0],
            end[1] - head * tangent[1] + side * head * 0.5 * radial[1],
            height,
        )
        points.extend([end[0], end[1], end[2], barb[0], barb[1], barb[2]])
    return points


def _straight_arrow(start, direction, length: float) -> list[float]:
    """Fleche droite, tige et deux barbes."""
    unit = normalize(direction)
    end = tuple(start[k] + unit[k] * length for k in range(3))
    points = [start[0], start[1], start[2], end[0], end[1], end[2]]
    # Une perpendiculaire quelconque pour ouvrir les barbes.
    reference = (0.0, 0.0, 1.0) if abs(unit[2]) < 0.9 else (1.0, 0.0, 0.0)
    side = normalize(cross(unit, reference))
    head = length * 0.25
    for sign in (1.0, -1.0):
        barb = tuple(end[k] - head * unit[k] + sign * head * 0.5 * side[k] for k in range(3))
        points.extend([end[0], end[1], end[2], barb[0], barb[1], barb[2]])
    return points


def build_overlays(result) -> list[dict]:
    """Reperes 3D : axe, rayons caracteristiques, sens de rotation et de sortie."""
    topology = result.topology
    blades = result.blades
    occupancy = result.occupancy
    if topology is None or occupancy is None:
        return []

    z_low, z_high = occupancy.z_min, occupancy.z_max
    height = z_high - z_low
    span = max(topology.r_tip, height) or 1.0
    groups: list[dict] = []

    groups.append(
        {
            "id": "axe",
            "label": "Axe de rotation Z",
            "color": "#5A6472",
            "points": [0.0, 0.0, z_low - 0.15 * height, 0.0, 0.0, z_high + 0.25 * height],
        }
    )
    groups.append(
        {
            "id": "aspiration",
            "label": f"Aspiration : r1s = {topology.r_1s * config.MM_PER_M:.1f} mm",
            "color": "#4FD0DE",
            "points": _circle(topology.r_1s, topology.z_1)
            + (_circle(topology.r_1h, topology.z_1) if topology.r_1h > 0.0 else []),
        }
    )
    groups.append(
        {
            "id": "sortie",
            "label": f"Sortie : r2 = {topology.r_2 * config.MM_PER_M:.1f} mm",
            "color": "#D97136",
            "points": _circle(topology.r_2, topology.z_2),
        }
    )
    # La fleche montre le sens retenu, ou a defaut celui que suggere la geometrie
    # -- c'est justement ce qu'on demande a l'utilisateur de confirmer sur la piece.
    rotation_sign = 0 if blades is None else (
        blades.rotation_sign or blades.observed_rotation_sign
    )
    if rotation_sign:
        groups.append(
            {
                "id": "rotation",
                "label": (
                    blades.rotation_label
                    if blades.rotation_sign
                    else f"{blades.observed_rotation_label} -- suggere, a confirmer"
                ),
                "color": "#F0A202" if blades.rotation_sign else "#9AA3AE",
                "points": _arc_arrow(
                    topology.r_tip * 1.18, z_high + 0.12 * height, rotation_sign
                ),
            }
        )
    if blades is not None:
        radial = result.discharge.get("composante_meridienne", "").startswith("radial")
        if radial:
            start = (topology.r_2 * 1.05, 0.0, topology.z_2)
            direction = (1.0, 0.0, 0.0)
        else:
            start = (topology.r_2 * 0.6, 0.0, topology.z_2)
            direction = (0.0, 0.0, -1.0)
        groups.append(
            {
                "id": "refoulement",
                "label": "Sens de sortie du liquide",
                "color": "#B4531A",
                "points": _straight_arrow(start, direction, span * 0.35),
            }
        )
    return [group for group in groups if group["points"]]


def _rows(result) -> dict:
    """Grandeurs affichees dans le rail de lecture."""
    from .report import CONFIDENCE_LABELS, geometry_table, performance_table

    header, performance = performance_table(result)
    return {
        "geometry": [list(row) for row in geometry_table(result)],
        "speeds": header,
        "performance": [[label] + values for label, values in performance],
        "confidence": [
            [name, CONFIDENCE_LABELS.get(level, level)]
            for name, level in sorted(result.confidence.items())
        ],
    }


def build_payload(result, mesh: TriMesh | None = None) -> dict:
    """Rassemble tout ce que la page a besoin de connaitre."""
    display = mesh if mesh is not None else result.mesh
    if display is None:
        raise ValueError("aucun maillage a afficher")
    if len(display.faces) > config.VIEWER_MAX_FACES:
        display = display.decimate(config.VIEWER_MAX_FACES / len(display.faces))

    positions: list[float] = []
    for a, b, c in display.triangles():
        positions.extend(a)
        positions.extend(b)
        positions.extend(c)
    classes = classify_faces(display, result.occupancy) if result.occupancy else [CLASS_OTHER] * len(display.faces)

    low, high = display.bounds()
    topology = result.topology
    blades = result.blades
    limit = result.speed_limit
    return {
        "nom": os.path.basename(result.import_report.path) if result.import_report else "roue",
        "positions": _b64_floats(positions),
        "classes": _b64_bytes(classes),
        "faces": len(display.faces),
        "bbox": {"min": list(low), "max": list(high)},
        "overlays": build_overlays(result),
        "resume": {
            "type": topology.machine_type if topology else "-",
            "pales": topology.blades.n_blades if topology else 0,
            "rotation": blades.rotation_label if blades else "-",
            "rotation_signe": blades.rotation_sign if blades else 0,
            # Le sens suggere ne sert qu'a faire tourner la vue, jamais a conclure.
            "rotation_suggeree": (blades.observed_rotation_sign if blades else 0),
            "beta1": blades.beta1_deg if blades else 0.0,
            "beta2": blades.beta2_deg if blades else 0.0,
            "confiance": _CONFIDENCE_FR.get(result.overall_confidence(), result.overall_confidence()),
            "vitesse_max": limit.rpm_max if limit else 0,
            "limite": limit.active_limit if limit else "-",
            "hypotheses": limit.assumptions if limit else "",
            "duree": result.elapsed_s,
        },
        "tables": _rows(result),
        "avertissements": list(result.warnings),
        "carte": _occupancy_data_uri(result),
        "grille": occupancy_grid(result),
        "vif": live_curves(result),
        "provenance": dict(result.provenance),
    }


def occupancy_grid(result) -> dict:
    """La carte `f(r, z)` elle-meme, et les reperes qu'on peut y superposer.

    L'image seule ne permet pas de lire une valeur : elle donne une impression.
    La grille, elle, se survole -- on pointe et on lit `r`, `z` et la fraction
    angulaire occupee.  C'est la lecture dont tout le reste decoule, elle merite
    d'etre interrogeable et pas seulement regardable.

    `f` est quantifiee sur un octet et transmise en base64 : une grille 200x200
    tient en 53 ko, la ou le meme tableau en JSON en ferait dix fois plus.
    """
    occupancy = result.occupancy
    if occupancy is None:
        return {}
    topology = result.topology
    loops = result.blade_loops
    valeurs = bytearray()
    for row in occupancy.f:
        for value in row:
            valeurs.append(max(0, min(255, round(value * 255.0))))

    reperes = {}
    if topology is not None:
        reperes = {
            "r1h": topology.r_1h * config.MM_PER_M,
            "r1s": topology.r_1s * config.MM_PER_M,
            "r2": topology.r_2 * config.MM_PER_M,
            "z1": topology.z_1 * config.MM_PER_M,
            "z2": topology.z_2 * config.MM_PER_M,
        }
    if loops is not None and loops.looped:
        reperes["boucle"] = {
            "r_interieur": loops.r_inner * config.MM_PER_M,
            "r_exterieur": loops.r_outer * config.MM_PER_M,
            "r_fusion": (loops.merge_radius * config.MM_PER_M) or None,
        }
    return {
        "nr": occupancy.nr,
        "nz": occupancy.nz,
        "r": [v * config.MM_PER_M for v in occupancy.r_centres],
        "z": [v * config.MM_PER_M for v in occupancy.z_centres],
        "f": base64.b64encode(bytes(valeurs)).decode("ascii"),
        "seuils": {"vide": config.F_VIDE, "matiere": config.F_MATIERE,
                   "solide": config.F_SOLIDE},
        "reperes": reperes,
        # Le sens debitant sur l'axe : +1 si le fluide va vers +Z, -1 sinon.
        "sens_debitant": -1,
    }


def live_curves(result) -> dict:
    """Ce qu'il faut pour deplacer les deux curseurs sans relancer l'analyse.

    Deux mecanismes differents, et il importe de ne pas les confondre :

    * le **regime** se deplace par les lois de similitude, que l'outil verifie
      deja a 0.00 % -- `Q` comme `n`, `H` comme `n^2`, `P` comme `n^3`, `NPSHr`
      comme `n^2`. C'est le meme modele, pas une approximation ;
    * **beta2** ne s'extrapole pas : les courbes sont **precalculees** en Python
      a `beta2 - 1`, `beta2` et `beta2 + 1` degre, et la page interpole entre
      trois courbes reelles.

    Reimplementer le modele de ligne moyenne en JavaScript aurait donne un
    second modele, libre de diverger du premier -- exactement l'incoherence que
    cet outil passe son temps a retirer.
    """
    import dataclasses

    from ..hydraulics import cavitation as cavitation_module
    from ..hydraulics import meanline as meanline_module

    if not result.curves or result.topology is None or result.blades is None:
        return {}
    reference = result.curves[0]
    point = reference.best_efficiency_point()
    if point is None:
        return {}

    limite = result.speed_limit
    data = meanline_module.MeanlineInput.from_geometry(result.topology, result.blades)
    familles = []
    for delta in (-config.BETA_SENSITIVITY_DEG, 0.0, config.BETA_SENSITIVITY_DEG):
        beta2 = data.beta2_deg + delta
        if not 0.0 < beta2 < 90.0:
            continue
        courbe = meanline_module.build_curve(
            dataclasses.replace(data, beta2_deg=beta2), reference.rpm
        )
        cavitation_module.apply_to_curve(courbe)
        familles.append({
            "beta2": beta2,
            "ecart": delta,
            "Q": [p.flow * config.SECONDS_PER_HOUR for p in courbe.points],
            "H": [p.head for p in courbe.points],
            "rendement": [p.efficiency * 100.0 for p in courbe.points],
            "npshr": [p.npshr for p in courbe.points],
        })
    return {
        "rpm_reference": reference.rpm,
        "rpm_min": min(config.DEFAULT_RPM),
        # La borne haute doit permettre d'atteindre la limite de cavitation :
        # c'est la question que l'utilisateur se pose vraiment -- jusqu'ou
        # puis-je monter -- et un curseur qui s'arrete avant la reponse ne la
        # pose meme pas.
        "rpm_max_affiche": max(
            max(config.DEFAULT_RPM),
            reference.rpm * 2.0,
            (limite.rpm_max if limite else 0.0) * config.RPM_SLIDER_MARGIN,
        ),
        "rpm_limite": limite.rpm_max if limite else 0,
        "limite_active": limite.active_limit if limite else "",
        "npsh_disponible": result.installation.npsha if result.installation else 0.0,
        "point": {
            "Q": point.flow * config.SECONDS_PER_HOUR,
            "H": point.head,
            "P": point.shaft_power / config.W_PER_KW,
            "couple": point.torque,
            "npshr": point.npshr,
            "rendement": point.efficiency * 100.0,
        },
        "familles": familles,
        "beta2_provenance": (
            "declare" if result.blades.forced_beta else "mesure"
        ),
    }


def _occupancy_data_uri(result) -> str:
    """Carte meridienne encodee en image, pour l'inclure dans la page."""
    if result.occupancy is None:
        return ""
    canvas = plot.occupancy_canvas(result.occupancy, "")
    return "data:image/png;base64," + base64.b64encode(canvas.encode_png()).decode("ascii")


def build_page(
    result,
    mesh: TriMesh | None = None,
    standalone: bool = True,
    title: str | None = None,
) -> str:
    """Construit la page complete (avec ou sans l'enveloppe HTML)."""
    data = build_payload(result, mesh)
    if title is None:
        title = f"Roue {data['resume']['type']} a {data['resume']['pales']} pales"
    body = _fill(json.dumps(data, ensure_ascii=False, allow_nan=False), server=False)
    return _wrap(body, title, standalone)


def _palette_css(page: str) -> str:
    """Injecte la palette et les piles de polices du style partage.

    Le gabarit porte des marques plutot que des valeurs : la page et les deux
    figures tirent ainsi leurs couleurs du meme endroit, et une page claire ne
    peut plus cotoyer des figures restees aux reglages d'origine.
    """
    remplacements = {
        "__FOND__": config.COULEUR_FOND,
        "__ENCRE__": config.COULEUR_ENCRE,
        "__GRILLE__": config.COULEUR_GRILLE,
        "__MESURE__": config.COULEUR_MESURE,
        "__DECLARE__": config.COULEUR_DECLARE,
        "__LIMITE__": config.COULEUR_LIMITE,
        "__FOND_SOMBRE__": config.COULEUR_FOND_SOMBRE,
        "__ENCRE_SOMBRE__": config.COULEUR_ENCRE_SOMBRE,
        "__GRILLE_SOMBRE__": config.COULEUR_GRILLE_SOMBRE,
        "__PILE_TEXTE__": style.PILE_TEXTE,
        "__PILE_NOMBRES__": style.PILE_NOMBRES,
    }
    for marque, valeur in remplacements.items():
        page = page.replace(marque, valeur)
    return page


def _fill(payload_json: str, server: bool) -> str:
    """Injecte la charge utile, le mode et le port par defaut dans le gabarit."""
    return (
        _palette_css(_TEMPLATE).replace("__DONNEES__", payload_json)
        .replace("__SERVEUR__", "true" if server else "false")
        .replace("__PORT__", str(config.SERVER_PORT))
    )


def _wrap(body: str, title: str, standalone: bool) -> str:
    """Ajoute ou non l'enveloppe HTML complete."""
    if not standalone:
        return f"<title>{title}</title>\n" + body
    return (
        '<!doctype html>\n<html lang="fr">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<style>html,body{margin:0; background:var(--fond); color:var(--encre)}img{max-width:100%}</style>\n'
        f"<title>{title}</title>\n</head>\n<body>\n" + body + "\n</body>\n</html>\n"
    )


def build_app_page(title: str = "Inspecteur de roue") -> str:
    """Coquille de l'application locale : la page s'ouvre vide, en attente d'un fichier."""
    return _wrap(_fill("null", server=True), title, standalone=True)


def write_page(
    result,
    directory: str,
    name: str = "vue3d.html",
    mesh: TriMesh | None = None,
    standalone: bool = True,
    title: str | None = None,
) -> str:
    """Ecrit la vue 3D dans le dossier de sortie."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(build_page(result, mesh, standalone=standalone, title=title))
    return path


#: Page complete. `__DONNEES__` est remplace par la charge utile JSON.
_TEMPLATE = r"""
<!-- Aucune ressource externe : la page s'ouvre depuis le disque, hors ligne,
     et rend la meme chose. Elle chargeait IBM Plex depuis Google Fonts, ce qui
     contredisait cette garantie et degradait en silence des que le reseau
     manquait -- c'est-a-dire sur le poste d'atelier ou elle sert. Les piles de
     polices systeme la remplacent. -->
<style>
/* Palette : une seule source, partagee avec les figures (io/style.py). Teal et
   violet ne decorent pas, ils codent mesure contre declare. Le theme clair est
   le defaut ; le sombre reste en bascule. */
:root{
  --fond:__FOND__; --encre:__ENCRE__; --grille:__GRILLE__;
  --mesure:__MESURE__; --declare:__DECLARE__; --limite:__LIMITE__;

  --paper:var(--fond); --panel:var(--fond); --ink:var(--encre);
  --ink-soft:color-mix(in srgb, var(--encre) 62%, var(--fond));
  --ink-faint:color-mix(in srgb, var(--encre) 42%, var(--fond));
  --line:var(--grille); --line-strong:color-mix(in srgb, var(--grille) 70%, var(--encre));
  --accent:var(--mesure); --accent-wash:color-mix(in srgb, var(--mesure) 10%, var(--fond));
  --oxide:var(--limite); --amber:var(--declare); --steel:var(--ink-faint);
  --stage:color-mix(in srgb, var(--encre) 6%, var(--fond));
  --stage-grid:var(--grille); --stage-ink:var(--encre);
  --ok:var(--mesure); --warn:var(--declare); --erreur:var(--limite);
  --pile-texte:__PILE_TEXTE__;
  --pile-nombres:__PILE_NOMBRES__;
  --rail:22.5rem;
}
[data-theme="sombre"]{
  --fond:__FOND_SOMBRE__; --encre:__ENCRE_SOMBRE__; --grille:__GRILLE_SOMBRE__;
  --mesure:#5FBDB6; --declare:#A192D6; --limite:#E8705E;
  --stage:color-mix(in srgb, var(--encre) 5%, var(--fond));
}


*{box-sizing:border-box}
body{margin:0}
.app{
  height:100vh; overflow:hidden; background:var(--paper); color:var(--ink);
  font-family:var(--pile-texte);
  font-size:14px; line-height:1.5; display:flex; flex-direction:column;
}
@media (max-width:960px){ .app{height:auto; overflow:visible} }
.eyebrow{
  font-family:"IBM Plex Sans Condensed", "IBM Plex Sans", ui-sans-serif, system-ui, sans-serif;
  font-weight:600; text-transform:uppercase; letter-spacing:.11em; font-size:11px;
  color:var(--ink-faint);
}

/* ---- barre de titre ---- */
.bar{
  display:flex; flex-wrap:wrap; align-items:baseline; gap:.35rem 1.5rem;
  padding:.9rem 1.25rem; border-bottom:1px solid var(--line); background:var(--panel); flex:none;
}
.bar h1{margin:0; font-size:1.05rem; font-weight:600; letter-spacing:-.01em; text-wrap:balance}
.bar .file{font-family:var(--pile-nombres); font-size:.78rem; color:var(--ink-soft)}
.verdicts{display:flex; flex-wrap:wrap; gap:.4rem .5rem; margin-left:auto}
.chip{
  display:inline-flex; align-items:center; gap:.4rem; padding:.2rem .55rem;
  border:1px solid var(--line-strong); border-radius:2px; font-size:.78rem; color:var(--ink-soft);
  background:var(--paper);
}
.chip b{color:var(--ink); font-weight:600}
.chip .dot{width:.5rem; height:.5rem; border-radius:50%; flex:none}

/* ---- corps ---- */
/* La carte d'occupation est l'element principal : elle occupe toute la largeur,
   au-dessus de la vue 3D et du rail. C'est la lecture dont tout le reste
   decoule, elle ne se met pas dans un coin. */
main{
  flex:1; min-height:0; overflow-y:auto; display:grid;
  grid-template-columns:1fr var(--rail); grid-template-rows:auto 1fr;
}
.instrument{grid-column:1 / -1}
.stage{grid-row:2}
@media (max-width:960px){ main{grid-template-columns:1fr} }

/* ---- scene ---- */
.stage{position:relative; background:var(--stage); min-height:26rem}
@media (max-width:960px){ .stage{height:62vh} }
.stage canvas{display:block; width:100%; height:100%; touch-action:none; cursor:grab}
.stage canvas:active{cursor:grabbing}
.stage .fallback{
  position:absolute; inset:0; display:none; place-items:center; padding:2rem; text-align:center;
  color:var(--stage-ink); font-size:.9rem;
}
.stage.no-webgl .fallback{display:grid}

.overlay{position:absolute; color:var(--stage-ink); font-size:.75rem}
.overlay.tools{top:.75rem; left:.75rem; display:flex; flex-wrap:wrap; gap:.3rem}
.overlay.legend{bottom:.75rem; left:.75rem; display:grid; gap:.3rem}
.overlay.hint{bottom:.75rem; right:.75rem; text-align:right; color:#767D87; line-height:1.7}
.btn{
  font:inherit; font-size:.75rem; color:var(--stage-ink); background:rgba(255,255,255,.05);
  border:1px solid var(--stage-grid); border-radius:2px; padding:.28rem .55rem; cursor:pointer;
}
.btn:hover{background:rgba(255,255,255,.11)}
.btn:focus-visible{outline:2px solid var(--accent); outline-offset:1px}
.btn[aria-pressed="true"]{background:var(--accent); border-color:var(--accent); color:#062224}
.btn:disabled{opacity:.45; cursor:default}
.key{display:flex; align-items:center; gap:.45rem}
.key .swatch{width:1.1rem; height:0; border-top:2px solid; flex:none}
.key .box{width:.7rem; height:.7rem; border-radius:1px; flex:none}
.key label{display:flex; align-items:center; gap:.45rem; cursor:pointer}
.key input{accent-color:var(--accent); margin:0}

/* ---- panneau d'import ---- */
.accueil{
  position:absolute; inset:0; display:grid; place-items:center; padding:1.5rem;
  background:var(--stage); overflow-y:auto;
}
.accueil[hidden]{display:none}
.depot{
  width:min(34rem, 100%); background:var(--panel); color:var(--ink);
  border:1px solid var(--line); border-radius:3px; padding:1.5rem;
}
.zone{
  display:block; position:relative; margin-top:.7rem;
  border:1.5px dashed var(--line-strong); border-radius:3px; padding:1.6rem 1rem;
  text-align:center; cursor:pointer; transition:border-color .12s, background .12s;
}
.zone:hover, .zone.actif{border-color:var(--accent); background:var(--accent-wash)}
.zone strong{display:block; font-size:1rem; margin-bottom:.2rem}
.zone span{color:var(--ink-soft); font-size:.82rem}
.zone input{position:absolute; width:1px; height:1px; opacity:0; pointer-events:none}
.champs{display:grid; grid-template-columns:repeat(2, 1fr); gap:.7rem; margin-top:1.1rem}
.champ{display:grid; gap:.2rem}
.champ label{font-size:.74rem; color:var(--ink-soft)}
.champ input, .champ select{
  font:inherit; font-size:.85rem; font-family:var(--pile-nombres);
  padding:.34rem .45rem; border:1px solid var(--line-strong); border-radius:2px;
  background:var(--paper); color:var(--ink); width:100%;
}
.champ input:focus-visible, .champ select:focus-visible{outline:2px solid var(--accent); outline-offset:0}
details{margin-top:.9rem}
details summary{cursor:pointer; font-size:.8rem; color:var(--accent)}
details summary:focus-visible{outline:2px solid var(--accent); outline-offset:2px}
.lancer{
  margin-top:1.1rem; width:100%; font:inherit; font-weight:600; font-size:.9rem;
  padding:.6rem; border:0; border-radius:2px; background:var(--accent); color:#042023; cursor:pointer;
}
.lancer:disabled{opacity:.5; cursor:default}
.hors-serveur{
  margin-top:1rem; padding:.85rem 1rem; border-left:2px solid var(--oxide);
  background:var(--paper); font-size:.82rem; color:var(--ink-soft);
}
.hors-serveur[hidden]{display:none}
.hors-serveur p{margin:0 0 .5rem}
.hors-serveur p:last-child{margin-bottom:0}
.hors-serveur strong{color:var(--ink)}
.hors-serveur code{
  display:inline-block; font-family:var(--pile-nombres); font-size:.8rem;
  background:var(--panel); border:1px solid var(--line); border-radius:2px; padding:.15rem .4rem;
  color:var(--ink); overflow-wrap:anywhere;
}
.etat{margin-top:.8rem; font-size:.82rem; color:var(--ink-soft); min-height:1.2em}
.etat.erreur{color:var(--erreur)}
.progres{
  height:2px; margin-top:.7rem; background:var(--line); overflow:hidden; display:none;
}
.progres.actif{display:block}
.progres i{display:block; height:100%; width:35%; background:var(--accent); animation:glisse 1.1s linear infinite}
@keyframes glisse{from{transform:translateX(-100%)} to{transform:translateX(390%)}}

/* ---- rail de lecture ---- */
.rail{border-left:1px solid var(--line); background:var(--panel); overflow-y:auto; min-height:0}
@media (max-width:960px){ .rail{overflow-y:visible; border-left:0; border-top:1px solid var(--line)} }
.rail section{padding:1.1rem 1.25rem; border-bottom:1px solid var(--line)}
.rail section:last-child{border-bottom:0}
.rail h2{margin:0; font-size:.95rem; font-weight:600; letter-spacing:-.005em}
.rail .step{display:flex; align-items:baseline; gap:.5rem; margin:0 0 .7rem}
.rail .step .n{
  font-family:var(--pile-nombres); font-size:.7rem; color:var(--accent);
  border:1px solid var(--accent); border-radius:2px; padding:0 .28rem; flex:none;
}
.facts{display:grid; gap:.75rem}
.fact{
  display:grid; grid-template-columns:1fr auto; gap:.1rem .6rem; align-items:baseline;
  padding-bottom:.75rem; border-bottom:1px solid var(--line);
}
.fact:last-child{padding-bottom:0; border-bottom:0}
.fact .k{color:var(--ink-soft); font-size:.79rem}
.fact .conf{font-size:.7rem; color:var(--ink-faint)}
.fact .val{
  grid-column:1 / -1; font-family:var(--pile-nombres);
  font-variant-numeric:tabular-nums; font-size:.92rem; color:var(--ink); overflow-wrap:anywhere;
}
table{width:100%; border-collapse:collapse; font-size:.82rem}
th,td{text-align:left; padding:.3rem 0; vertical-align:baseline}
th{font-weight:500; color:var(--ink-soft); font-size:.78rem}
.scroller{overflow-x:auto}
.perf th:first-child{width:45%}
.perf td{
  font-family:var(--pile-nombres); font-variant-numeric:tabular-nums;
  text-align:right; white-space:nowrap; padding-left:.7rem;
}
.perf thead th{border-bottom:1px solid var(--line-strong); padding-bottom:.35rem; text-align:right}
.perf thead th:first-child{text-align:left}
.perf tbody tr + tr th, .perf tbody tr + tr td{border-top:1px solid var(--line)}

.headline{display:flex; align-items:baseline; gap:.5rem; margin:.2rem 0 .5rem}
.headline .value{
  font-family:var(--pile-nombres); font-size:2rem; font-weight:500;
  letter-spacing:-.02em; line-height:1; color:var(--accent);
}
.headline .unit{color:var(--ink-soft); font-size:.85rem}
.note{margin:.5rem 0 0; color:var(--ink-soft); font-size:.78rem}
.warnings{margin:0; padding-left:1.05rem; font-size:.8rem; color:var(--ink-soft)}
.warnings li{margin:.35rem 0}
.warnings li::marker{color:var(--warn)}
.map{width:100%; height:auto; display:block; border:1px solid var(--line); border-radius:2px; background:#fff}
.caption{margin:.45rem 0 0; font-size:.75rem; color:var(--ink-faint)}
.uncert{border-left:2px solid var(--oxide); padding:.1rem 0 .1rem .7rem; font-size:.8rem; color:var(--ink-soft)}
.fichiers{display:grid; gap:.4rem; margin:0; padding:0; list-style:none; font-size:.82rem}
.fichiers a{color:var(--accent); text-decoration:none; font-family:var(--pile-nombres)}
.fichiers a:hover{text-decoration:underline}
.vide{color:var(--ink-soft); font-size:.85rem}

/* ---- instrument : la carte d'occupation et les deux curseurs ---- */
.instrument{
  border-bottom:1px solid var(--line); background:var(--paper);
  padding:.9rem 1.1rem 1rem; display:grid; gap:.75rem;
  grid-template-columns:minmax(0,1fr) 19rem;
}
.instrument h2{
  font-size:var(--t-section, 20px); font-weight:600; margin:0 0 .35rem; line-height:1.25;
}
.carte-boite{position:relative; border:1px solid var(--line); background:var(--fond)}
#carte-canevas{display:block; width:100%; height:auto; cursor:crosshair}
.instrument{grid-template-columns:minmax(0,1fr) 21rem}
.lecture{
  position:absolute; top:.5rem; right:.5rem; background:var(--fond);
  border:1px solid var(--line); padding:.4rem .55rem; min-width:8.5rem;
  font-family:var(--pile-nombres); font-variant-numeric:tabular-nums;
  font-size:var(--t-nombre, 13px); line-height:1.55; pointer-events:none;
}
.lecture dt{display:inline-block; width:1.4rem; color:var(--ink-soft); font-family:var(--pile-texte)}
.lecture dd{display:inline; margin:0}
.lecture div{white-space:nowrap}
.bascules{display:flex; flex-wrap:wrap; gap:.3rem; margin-top:.5rem}
.bascule{
  font:inherit; font-size:var(--t-controle, 13px); font-weight:500;
  color:var(--ink); background:transparent; border:1px solid var(--line-strong);
  padding:.25rem .55rem; cursor:pointer;
}
.bascule[aria-pressed="true"]{background:var(--mesure); border-color:var(--mesure); color:var(--fond)}
.bascule:focus-visible{outline:2px solid var(--mesure); outline-offset:1px}

.reglages{display:grid; gap:.9rem; align-content:start}
.reglage{display:grid; gap:.25rem}
.reglage .titre{
  display:flex; justify-content:space-between; align-items:baseline; gap:.5rem;
  font-size:var(--t-controle, 13px); font-weight:500;
}
.reglage output{
  font-family:var(--pile-nombres); font-variant-numeric:tabular-nums; font-weight:600;
}
.reglage input[type=range]{width:100%; accent-color:var(--mesure); margin:0}
.reglage .bornes{
  display:flex; justify-content:space-between;
  font-size:var(--t-note, 11px); color:var(--ink-soft);
  font-family:var(--pile-nombres); font-variant-numeric:tabular-nums;
}
.reglage.franchi output{color:var(--limite)}
.reglage.franchi input[type=range]{accent-color:var(--limite)}
.alerte-limite{
  font-size:var(--t-note, 11px); color:var(--limite); font-weight:500;
  border-left:2px solid var(--limite); padding-left:.45rem;
}
.alerte-limite[hidden]{display:none}

.vivant{width:100%; border-collapse:collapse; font-size:var(--t-nombre, 13px)}
.vivant th{
  text-align:left; font-family:var(--pile-texte); font-weight:500;
  color:var(--ink-soft); font-size:var(--t-controle, 13px); padding:.22rem 0;
}
.vivant td{
  text-align:right; font-family:var(--pile-nombres); font-variant-numeric:tabular-nums;
  padding:.22rem 0 .22rem .6rem; white-space:nowrap;
}
.vivant td.src{
  text-align:left; font-family:var(--pile-texte); font-size:var(--t-note, 11px);
  color:var(--ink-soft); padding-left:.5rem;
}
.vivant tr + tr th, .vivant tr + tr td{border-top:1px solid var(--line)}
/* Provenance : la couleur ne porte pas seule. */
.p-mesure{font-weight:600}
.p-declare{font-weight:400; border-left:2px solid var(--declare); padding-left:.35rem}
.p-defaut{font-weight:400; color:var(--ink-soft); font-style:italic}
.c-low{border-left:2px solid var(--limite); padding-left:.35rem}
.miniplot{width:100%; height:auto; display:block; border:1px solid var(--line); margin-top:.2rem}
@media (max-width:60rem){ .instrument{grid-template-columns:minmax(0,1fr)} }

@media (prefers-reduced-motion:reduce){ *{animation:none !important; transition:none !important} }
</style>

<div class="app">
  <header class="bar">
    <h1 id="titre">Inspecteur de roue</h1>
    <span class="file" id="fichier"></span>
    <div class="verdicts" id="verdicts"></div>
    <button type="button" class="bascule" id="bascule-theme"
            aria-pressed="false" title="Basculer clair / sombre">sombre</button>
  </header>

  <main>

    <section class="instrument" id="instrument" hidden>
      <div>
        <h2>Carte d'occupation f(r, z)</h2>
        <div class="carte-boite">
          <canvas id="carte-canevas" width="900" height="460"></canvas>
          <dl class="lecture" id="lecture">
            <div><dt>r</dt><dd id="lec-r">&mdash;</dd></div>
            <div><dt>z</dt><dd id="lec-z">&mdash;</dd></div>
            <div><dt>f</dt><dd id="lec-f">&mdash;</dd></div>
          </dl>
        </div>
        <div class="bascules" id="bascules"></div>
      </div>

      <div class="reglages">
        <div class="reglage" id="reglage-regime">
          <div class="titre"><span>Regime</span><output id="sortie-regime"></output></div>
          <input type="range" id="curseur-regime" min="200" max="4000" step="10">
          <div class="bornes"><span id="borne-min"></span><span id="borne-max"></span></div>
          <p class="alerte-limite" id="alerte-cavitation" hidden></p>
        </div>

        <div class="reglage" id="reglage-beta">
          <div class="titre"><span>beta2 &plusmn; 1&deg;</span><output id="sortie-beta"></output></div>
          <input type="range" id="curseur-beta" min="-1" max="1" step="0.05" value="0">
          <div class="bornes"><span>&minus;1&deg;</span><span>+1&deg;</span></div>
        </div>

        <table class="vivant"><tbody id="valeurs-vives"></tbody></table>
        <canvas class="miniplot" id="miniplot" width="360" height="180"></canvas>
      </div>
    </section>
    <section class="stage" id="stage">
      <canvas id="gl"></canvas>
      <div class="fallback">
        Ce navigateur n'expose pas WebGL&nbsp;2 : la vue 3D ne peut pas s'afficher.
        Les resultats de l'analyse restent lisibles dans le panneau de droite.
      </div>
      <div class="overlay tools" id="outils"></div>
      <div class="overlay legend" id="legende"></div>
      <div class="overlay hint">
        glisser&nbsp;: pivoter &nbsp;·&nbsp; molette&nbsp;: zoomer<br>
        maj + glisser&nbsp;: deplacer
      </div>

      <div class="accueil" id="accueil" hidden>
        <form class="depot" id="formulaire">
          <p class="eyebrow">Importer une roue</p>
          <label class="zone" id="zone">
            <strong id="zone-titre">Deposez un fichier 3D ici</strong>
            <span id="zone-detail">ou cliquez pour le choisir &mdash; .stl, .obj, .ply, .off, .dxf</span>
            <input type="file" id="fichier-entree" accept=".stl,.obj,.ply,.off,.dxf,.3ds,.step,.stp,.iges,.igs">
          </label>

          <div class="champs">
            <div class="champ">
              <label for="unite">Unite du fichier</label>
              <select id="unite">
                <option value="cm" selected>centimetres</option>
                <option value="mm">millimetres</option>
                <option value="m">metres</option>
                <option value="in">pouces</option>
              </select>
            </div>
            <div class="champ">
              <label for="rpm">Regimes (tr/min)</label>
              <input id="rpm" value="1000 2000 3000" inputmode="numeric">
            </div>
          </div>

          <details>
            <summary>Reglages avances</summary>
            <div class="champs">
              <div class="champ">
                <label for="r-asp">Rayon d'aspiration (cm)</label>
                <input id="r-asp" placeholder="detecte" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="pales">Nombre de pales</label>
                <input id="pales" placeholder="detecte" inputmode="numeric">
              </div>
              <div class="champ">
                <label for="beta1">beta1 (deg)</label>
                <input id="beta1" placeholder="detecte" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="beta2">beta2 (deg)</label>
                <input id="beta2" placeholder="detecte" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="rotation">Sens de rotation</label>
                <select id="rotation">
                  <option value="">a indiquer</option>
                  <option value="horaire">horaire (vu de +Z)</option>
                  <option value="antihoraire">anti-horaire (vu de +Z)</option>
                </select>
              </div>
              <div class="champ">
                <label for="temperature">Temperature (&deg;C)</label>
                <input id="temperature" value="20" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="altitude">Altitude (m)</label>
                <input id="altitude" value="0" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="hauteur">Hauteur d'aspiration (m)</label>
                <input id="hauteur" value="0" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="pertes">Pertes d'aspiration (m)</label>
                <input id="pertes" value="0.5" inputmode="decimal">
              </div>
              <div class="champ">
                <label for="grille">Grille (nr &times; nz)</label>
                <input id="grille" value="200" inputmode="numeric">
              </div>
            </div>
          </details>

          <div class="hors-serveur" id="hors-serveur" hidden>
            <p><strong>Cette page ne peut pas analyser toute seule.</strong></p>
            <p>
              Elle est la façade de l'application ; le calcul se fait en Python, par un
              serveur local. Ouverte depuis un fichier, elle n'a personne a qui parler.
            </p>
            <p>Dans un terminal, a la racine du projet&nbsp;:</p>
            <code>python3 -m impeller_analyzer.serve</code>
            <p>
              puis ouvrez <code>http://127.0.0.1:__PORT__/</code> — le navigateur s'y ouvre
              tout seul. Vous y retrouverez ce meme ecran, en etat de marche.
            </p>
          </div>

          <button class="lancer" id="lancer" type="submit" disabled>Analyser</button>
          <div class="progres" id="progres"><i></i></div>
          <p class="etat" id="etat"></p>
        </form>
      </div>
    </section>

    <aside class="rail" id="rail"></aside>
  </main>
</div>

<script>
"use strict";
const BOOT = __DONNEES__;
const SERVEUR = __SERVEUR__;

/* ---------- outils ---------- */
function el(tag, attrs = {}, ...children){
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (key === "class") node.className = value;
    else node.setAttribute(key, value);
  }
  for (const child of children) if (child !== null && child !== undefined) node.append(child);
  return node;
}
function bytesOf(b64){
  const bin = atob(b64), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

/* ---------- petites matrices 4x4 ---------- */
const M4 = {
  multiply(a, b){
    const out = new Float32Array(16);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
      let s = 0;
      for (let k = 0; k < 4; k++) s += a[k*4+j] * b[i*4+k];
      out[i*4+j] = s;
    }
    return out;
  },
  perspective(fovy, aspect, near, far){
    const f = 1 / Math.tan(fovy / 2), d = near - far;
    return new Float32Array([f/aspect,0,0,0, 0,f,0,0, 0,0,(far+near)/d,-1, 0,0,2*far*near/d,0]);
  },
  lookAt(eye, target, up){
    const z = norm3([eye[0]-target[0], eye[1]-target[1], eye[2]-target[2]]);
    const x = norm3(cross3(up, z));
    const y = cross3(z, x);
    return new Float32Array([
      x[0],y[0],z[0],0, x[1],y[1],z[1],0, x[2],y[2],z[2],0,
      -(x[0]*eye[0]+x[1]*eye[1]+x[2]*eye[2]),
      -(y[0]*eye[0]+y[1]*eye[1]+y[2]*eye[2]),
      -(z[0]*eye[0]+z[1]*eye[1]+z[2]*eye[2]), 1]);
  },
  rotationZ(a){
    const c = Math.cos(a), s = Math.sin(a);
    return new Float32Array([c,s,0,0, -s,c,0,0, 0,0,1,0, 0,0,0,1]);
  }
};
function cross3(a, b){ return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]; }
function norm3(a){ const l = Math.hypot(a[0],a[1],a[2]) || 1; return [a[0]/l, a[1]/l, a[2]/l]; }

/* ---------- WebGL ---------- */
const canvas = document.getElementById("gl");
const stage = document.getElementById("stage");
const gl = canvas.getContext("webgl2", {antialias:true, alpha:false});
if (!gl) stage.classList.add("no-webgl");

const MESH_VS = `#version 300 es
in vec3 a_pos; in vec3 a_nrm; in vec3 a_col;
uniform mat4 u_proj, u_view, u_model;
out vec3 v_nrm; out vec3 v_col; out vec3 v_world;
void main(){
  vec4 world = u_model * vec4(a_pos, 1.0);
  v_world = world.xyz;
  v_nrm = mat3(u_model) * a_nrm;
  v_col = a_col;
  gl_Position = u_proj * u_view * world;
}`;
const MESH_FS = `#version 300 es
precision highp float;
in vec3 v_nrm; in vec3 v_col; in vec3 v_world;
uniform float u_clip;
uniform vec3 u_eye;
out vec4 fragColor;
void main(){
  if (u_clip > 0.5 && v_world.y > 0.0) discard;
  vec3 n = normalize(v_nrm);
  vec3 toEye = normalize(u_eye - v_world);
  if (dot(n, toEye) < 0.0) n = -n;
  float key  = max(dot(n, normalize(vec3(0.45, 0.7, 0.9))), 0.0);
  float fill = max(dot(n, normalize(vec3(-0.6, -0.3, 0.35))), 0.0);
  float rim  = pow(1.0 - max(dot(n, toEye), 0.0), 3.0);
  vec3 c = v_col * (0.24 + 0.72 * key + 0.20 * fill) + rim * 0.16;
  fragColor = vec4(pow(c, vec3(0.4545)), 1.0);
}`;
const LINE_VS = `#version 300 es
in vec3 a_pos;
uniform mat4 u_proj, u_view;
void main(){ gl_Position = u_proj * u_view * vec4(a_pos, 1.0); }`;
const LINE_FS = `#version 300 es
precision highp float;
uniform vec3 u_color;
out vec4 fragColor;
void main(){ fragColor = vec4(u_color, 1.0); }`;

function compile(src, type){
  const sh = gl.createShader(type);
  gl.shaderSource(sh, src); gl.compileShader(sh);
  if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(sh));
  return sh;
}
function program(vs, fs){
  const p = gl.createProgram();
  gl.attachShader(p, compile(vs, gl.VERTEX_SHADER));
  gl.attachShader(p, compile(fs, gl.FRAGMENT_SHADER));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
let meshProgram = null, lineProgram = null;
if (gl) {
  meshProgram = program(MESH_VS, MESH_FS);
  lineProgram = program(LINE_VS, LINE_FS);
  gl.enable(gl.DEPTH_TEST);
  gl.clearColor(0.094, 0.106, 0.129, 1);
}
function attach(data, prog, name, size, store){
  const buf = gl.createBuffer();
  store.push(buf);
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
  const loc = gl.getAttribLocation(prog, name);
  gl.enableVertexAttribArray(loc);
  gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0);
}

/* ---------- etat ---------- */
const PALETTE = [[0.42,0.45,0.49], [0.16,0.62,0.67], [0.60,0.55,0.46]];
let D = null;
let scene = null;          /* {vao, count, groups, buffers, vaos} */
let view = null, home = null, radius = 1;
let spin = 0, spinning = false, clip = false, last = 0, running = false;

function disposeScene(){
  if (!gl || !scene) return;
  for (const buf of scene.buffers) gl.deleteBuffer(buf);
  for (const vao of scene.vaos) gl.deleteVertexArray(vao);
  scene = null;
}

function buildScene(payload){
  disposeScene();
  if (!gl) return;
  const positions = new Float32Array(bytesOf(payload.positions).buffer);
  const faceClass = bytesOf(payload.classes);
  const normals = new Float32Array(positions.length);
  const colors = new Float32Array(positions.length);
  for (let f = 0; f < positions.length / 9; f++) {
    const o = f * 9;
    const ux = positions[o+3]-positions[o],  uy = positions[o+4]-positions[o+1], uz = positions[o+5]-positions[o+2];
    const vx = positions[o+6]-positions[o],  vy = positions[o+7]-positions[o+1], vz = positions[o+8]-positions[o+2];
    let nx = uy*vz-uz*vy, ny = uz*vx-ux*vz, nz = ux*vy-uy*vx;
    const len = Math.hypot(nx, ny, nz) || 1;
    nx /= len; ny /= len; nz /= len;
    const c = PALETTE[faceClass[f]] || PALETTE[2];
    for (let k = 0; k < 3; k++) {
      normals[o+k*3] = nx; normals[o+k*3+1] = ny; normals[o+k*3+2] = nz;
      colors[o+k*3] = c[0]; colors[o+k*3+1] = c[1]; colors[o+k*3+2] = c[2];
    }
  }
  const buffers = [], vaos = [];
  const vao = gl.createVertexArray();
  vaos.push(vao);
  gl.bindVertexArray(vao);
  attach(positions, meshProgram, "a_pos", 3, buffers);
  attach(normals, meshProgram, "a_nrm", 3, buffers);
  attach(colors, meshProgram, "a_col", 3, buffers);
  gl.bindVertexArray(null);

  const groups = payload.overlays.map(group => {
    const gvao = gl.createVertexArray();
    vaos.push(gvao);
    gl.bindVertexArray(gvao);
    attach(new Float32Array(group.points), lineProgram, "a_pos", 3, buffers);
    gl.bindVertexArray(null);
    const hex = group.color.replace("#", "");
    return {
      id: group.id, label: group.label, color: group.color, vao: gvao,
      count: group.points.length / 3, visible: true,
      rgb: [0, 2, 4].map(i => parseInt(hex.slice(i, i + 2), 16) / 255)
    };
  });
  scene = {vao, count: positions.length / 3, groups, buffers, vaos};
}

function frameCamera(payload){
  const lo = payload.bbox.min, hi = payload.bbox.max;
  radius = Math.max(Math.hypot(hi[0]-lo[0], hi[1]-lo[1]) / 2, (hi[2]-lo[2]) / 2, 1e-3);
  view = {yaw:-0.9, pitch:0.5, dist:radius * 3.4, target:[0, 0, (lo[2] + hi[2]) / 2]};
  home = {yaw:view.yaw, pitch:view.pitch, dist:view.dist, target:view.target.slice()};
  spin = 0;
}

function eyePosition(){
  const cp = Math.cos(view.pitch);
  return [
    view.target[0] + view.dist * cp * Math.cos(view.yaw),
    view.target[1] + view.dist * cp * Math.sin(view.yaw),
    view.target[2] + view.dist * Math.sin(view.pitch)
  ];
}
function resize(){
  const ratio = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(1, Math.round(canvas.clientWidth * ratio));
  const h = Math.max(1, Math.round(canvas.clientHeight * ratio));
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
}

function draw(now){
  requestAnimationFrame(draw);
  if (!gl) return;
  const dt = last ? Math.min((now - last) / 1000, 0.1) : 0;
  last = now;
  resize();
  gl.viewport(0, 0, canvas.width, canvas.height);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  if (!scene || !view) return;
  if (spinning) spin += dt * 1.1 * (D.resume.rotation_signe || D.resume.rotation_suggeree || 1);

  const eye = eyePosition();
  const proj = M4.perspective(Math.PI / 4.2, canvas.width / canvas.height, radius * 0.02, radius * 60);
  const camera = M4.lookAt(eye, view.target, [0, 0, 1]);

  gl.useProgram(meshProgram);
  gl.uniformMatrix4fv(gl.getUniformLocation(meshProgram, "u_proj"), false, proj);
  gl.uniformMatrix4fv(gl.getUniformLocation(meshProgram, "u_view"), false, camera);
  gl.uniformMatrix4fv(gl.getUniformLocation(meshProgram, "u_model"), false, M4.rotationZ(spin));
  gl.uniform1f(gl.getUniformLocation(meshProgram, "u_clip"), clip ? 1 : 0);
  gl.uniform3fv(gl.getUniformLocation(meshProgram, "u_eye"), new Float32Array(eye));
  gl.bindVertexArray(scene.vao);
  gl.drawArrays(gl.TRIANGLES, 0, scene.count);

  gl.useProgram(lineProgram);
  gl.uniformMatrix4fv(gl.getUniformLocation(lineProgram, "u_proj"), false, proj);
  gl.uniformMatrix4fv(gl.getUniformLocation(lineProgram, "u_view"), false, camera);
  const colorLoc = gl.getUniformLocation(lineProgram, "u_color");
  for (const group of scene.groups) {
    if (!group.visible) continue;
    gl.uniform3fv(colorLoc, new Float32Array(group.rgb));
    gl.bindVertexArray(group.vao);
    gl.drawArrays(gl.LINES, 0, group.count);
  }
  gl.bindVertexArray(null);
}

/* ---------- manipulation ---------- */
let pointer = null;
canvas.addEventListener("pointerdown", e => {
  if (!view) return;
  pointer = {x:e.clientX, y:e.clientY, pan:e.shiftKey || e.button === 1};
  canvas.setPointerCapture(e.pointerId);
});
canvas.addEventListener("pointermove", e => {
  if (!pointer || !view) return;
  const dx = e.clientX - pointer.x, dy = e.clientY - pointer.y;
  pointer.x = e.clientX; pointer.y = e.clientY;
  if (pointer.pan) {
    const k = view.dist * 0.0022;
    view.target[0] += (Math.sin(view.yaw) * dx) * k;
    view.target[1] += (-Math.cos(view.yaw) * dx) * k;
    view.target[2] += dy * k;
  } else {
    view.yaw -= dx * 0.008;
    view.pitch = Math.max(-1.52, Math.min(1.52, view.pitch + dy * 0.008));
  }
});
const release = () => { pointer = null; };
canvas.addEventListener("pointerup", release);
canvas.addEventListener("pointercancel", release);
canvas.addEventListener("wheel", e => {
  if (!view) return;
  e.preventDefault();
  view.dist = Math.max(radius * 0.6, Math.min(radius * 22, view.dist * Math.exp(e.deltaY * 0.0012)));
}, {passive:false});
let pinch = 0;
canvas.addEventListener("touchmove", e => {
  if (e.touches.length !== 2 || !view) return;
  e.preventDefault();
  const d = Math.hypot(e.touches[0].clientX - e.touches[1].clientX,
                       e.touches[0].clientY - e.touches[1].clientY);
  if (pinch) view.dist = Math.max(radius * 0.6, Math.min(radius * 22, view.dist * pinch / d));
  pinch = d;
}, {passive:false});
canvas.addEventListener("touchend", () => { pinch = 0; });

/* ---------- barre d'outils ---------- */
const outils = document.getElementById("outils");
const accueil = document.getElementById("accueil");
let spinButton = null;

function toolButton(label, action, pressed){
  const b = el("button", {class:"btn", type:"button"}, label);
  if (pressed !== undefined) b.setAttribute("aria-pressed", String(pressed));
  b.addEventListener("click", () => action(b));
  outils.append(b);
  return b;
}
toolButton("Vue de dessus (+Z)", () => { if (view) { view.yaw = -Math.PI/2; view.pitch = 1.45; } });
toolButton("Vue de face", () => { if (view) { view.yaw = -Math.PI/2; view.pitch = 0.06; } });
toolButton("Isometrique", () => { if (view) { view.yaw = home.yaw; view.pitch = home.pitch; } });
toolButton("Recadrer", () => { if (view) { view.dist = home.dist; view.target = home.target.slice(); } });
toolButton("Demi-coupe", b => { clip = !clip; b.setAttribute("aria-pressed", String(clip)); }, false);
spinButton = toolButton("Faire tourner", b => {
  spinning = !spinning;
  b.setAttribute("aria-pressed", String(spinning));
  b.textContent = spinning ? "Arreter" : "Faire tourner";
}, false);
if (SERVEUR) {
  toolButton("Changer de roue", () => {
    spinning = false;
    spinButton.setAttribute("aria-pressed", "false");
    spinButton.textContent = "Faire tourner";
    accueil.hidden = false;
  });
}

/* ---------- rail de lecture ---------- */
const rail = document.getElementById("rail");
const legende = document.getElementById("legende");
const titre = document.getElementById("titre");
const fichierNom = document.getElementById("fichier");
const verdicts = document.getElementById("verdicts");


/* ---------------------------------------------------------------------------
   L'instrument : la carte d'occupation, et les deux manipulations qui comptent.

   La carte est la lecture dont tout le reste decoule -- elle merite d'etre
   interrogeable et pas seulement regardable. Les deux curseurs recalculent dans
   la page, mais aucun ne reimplemente le modele : le regime passe par les lois
   de similitude que l'outil verifie deja a 0.00 %, et beta2 interpole entre
   trois courbes reellement calculees en Python.
--------------------------------------------------------------------------- */
const INSTRUMENT = {grille:null, vif:null, f:null, bascules:{}, rpm:0, dbeta:0};

function lireOctets(b64){
  const brut = atob(b64), out = new Uint8Array(brut.length);
  for (let i = 0; i < brut.length; i++) out[i] = brut.charCodeAt(i);
  return out;
}

/* Meme rampe que les figures : teinte unique, luminance monotone. Une echelle
   multicolore fabriquerait des frontieres que les donnees n'ont pas.

   Les deux bornes viennent du theme et non de valeurs en dur : sinon la carte
   garde un fond clair quand la page passe en sombre, et l'incoherence est
   exactement celle qu'on cherchait a supprimer entre la page et les figures. */
function composantes(nom, repli){
  const v = getComputedStyle(document.body).getPropertyValue(nom).trim();
  const m = /^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(v);
  return m ? [1,2,3].map(k => parseInt(m[k], 16)) : repli;
}
function rampe(t, bornes){
  t = Math.max(0, Math.min(1, t));
  const w = Math.pow(t, 0.85);
  const [a, b] = bornes;
  return [0,1,2].map(k => Math.round(a[k] + (b[k] - a[k]) * w));
}
/* La meme compression par morceaux que la figure PNG (plot.occupancy_scale) :
   la bande de pales, souvent tres faible, garde l'essentiel de la dynamique. */
function echelle(f){
  const s = INSTRUMENT.grille.seuils;
  if (f <= s.vide) return s.vide > 0 ? 0.15 * f / s.vide : 0;
  if (f >= s.solide) return 0.85 + 0.15 * (f - s.solide) / (1 - s.solide);
  return 0.15 + 0.70 * (f - s.vide) / (s.solide - s.vide);
}

function tracerCarte(){
  const g = INSTRUMENT.grille, cv = document.getElementById("carte-canevas");
  if (!g || !cv) return;
  const ctx = cv.getContext("2d");
  const pad = {g:52, d:16, h:14, b:34};
  const W = cv.width, H = cv.height;
  const px0 = pad.g, px1 = W - pad.d, py0 = pad.h, py1 = H - pad.b;
  const style = getComputedStyle(document.body);
  const encre = style.getPropertyValue("--encre").trim() || "#16232B";
  const grille = style.getPropertyValue("--grille").trim() || "#DDE3E0";
  const declare = style.getPropertyValue("--declare").trim() || "#5B4B8A";
  const limite = style.getPropertyValue("--limite").trim() || "#9B1D20";
  ctx.fillStyle = style.getPropertyValue("--fond").trim() || "#FBFBFA";
  ctx.fillRect(0, 0, W, H);

  const rMax = g.r[g.r.length - 1], zMin = g.z[0], zMax = g.z[g.z.length - 1];
  const X = r => px0 + (px1 - px0) * r / rMax;
  const Y = z => py1 - (py1 - py0) * (z - zMin) / (zMax - zMin);

  const bornes = [composantes("--fond", [251, 251, 250]),
                  composantes("--mesure", [29, 111, 106])];
  const img = ctx.createImageData(g.nr, g.nz);
  for (let iz = 0; iz < g.nz; iz++) {
    for (let ir = 0; ir < g.nr; ir++) {
      const f = INSTRUMENT.f[iz * g.nr + ir] / 255;
      const c = rampe(echelle(f), bornes);
      /* l'image se dessine du haut vers le bas, la grille du bas vers le haut */
      const o = ((g.nz - 1 - iz) * g.nr + ir) * 4;
      img.data[o] = c[0]; img.data[o+1] = c[1]; img.data[o+2] = c[2]; img.data[o+3] = 255;
    }
  }
  const tampon = document.createElement("canvas");
  tampon.width = g.nr; tampon.height = g.nz;
  tampon.getContext("2d").putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(tampon, px0, py0, px1 - px0, py1 - py0);

  ctx.strokeStyle = grille; ctx.lineWidth = 1;
  ctx.strokeRect(px0 + .5, py0 + .5, px1 - px0 - 1, py1 - py0 - 1);
  ctx.fillStyle = encre;
  ctx.font = '11px ' + style.getPropertyValue("--pile-nombres");
  ctx.textAlign = "center";
  for (let k = 0; k <= 4; k++) {
    const r = rMax * k / 4;
    ctx.fillText(r.toFixed(0), X(r), py1 + 14);
  }
  ctx.fillText("rayon r (mm)", (px0 + px1) / 2, py1 + 28);
  ctx.textAlign = "right";
  for (let k = 0; k <= 4; k++) {
    const z = zMin + (zMax - zMin) * k / 4;
    ctx.fillText(z.toFixed(0), px0 - 6, Y(z) + 3);
  }

  const rep = g.reperes || {};
  const trait = (couleur, tirets, dessin) => {
    ctx.save(); ctx.strokeStyle = couleur; ctx.lineWidth = 1.5;
    ctx.setLineDash(tirets); dessin(); ctx.restore();
  };
  const marque = (texte, x, y, couleur) => {
    ctx.save(); ctx.fillStyle = couleur; ctx.textAlign = "left";
    ctx.font = '11px ' + style.getPropertyValue("--pile-texte");
    ctx.fillText(texte, x + 3, y - 3); ctx.restore();
  };

  if (INSTRUMENT.bascules.axe) {
    trait(encre, [6, 4], () => {
      ctx.beginPath(); ctx.moveTo(X(0), py0); ctx.lineTo(X(0), py1); ctx.stroke();
    });
    marque("axe", X(0), py0 + 12, encre);
  }
  if (INSTRUMENT.bascules.rayons) {
    for (const [cle, libelle] of [["r1h","r1h"], ["r1s","r1s"], ["r2","r2"]]) {
      const v = rep[cle];
      if (!v) continue;
      trait(declare, [], () => {
        ctx.beginPath(); ctx.moveTo(X(v), py0); ctx.lineTo(X(v), py1); ctx.stroke();
      });
      marque(libelle, X(v), py1 - 4, declare);
    }
  }
  if (INSTRUMENT.bascules.plans) {
    for (const [cle, libelle] of [["z1","entree"], ["z2","sortie"]]) {
      const v = rep[cle];
      if (v === undefined || v === null) continue;
      trait(encre, [2, 3], () => {
        ctx.beginPath(); ctx.moveTo(px0, Y(v)); ctx.lineTo(px1, Y(v)); ctx.stroke();
      });
      marque(libelle, px0, Y(v), encre);
    }
    /* Le sens debitant : une fleche, pas une convention ecrite. */
    const x = px1 - 26, y0 = py0 + 26, y1 = py0 + 66;
    const bas = g.sens_debitant < 0;
    trait(encre, [], () => {
      ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(x - 4, bas ? y1 - 7 : y0 + 7);
      ctx.lineTo(x, bas ? y1 : y0); ctx.lineTo(x + 4, bas ? y1 - 7 : y0 + 7);
      ctx.stroke();
    });
  }
  if (INSTRUMENT.bascules.boucle && rep.boucle) {
    const b = rep.boucle;
    trait(limite, [5, 3], () => {
      for (const v of [b.r_interieur, b.r_exterieur]) {
        if (!v) continue;
        ctx.beginPath(); ctx.moveTo(X(v), py0); ctx.lineTo(X(v), py1); ctx.stroke();
      }
    });
    if (b.r_fusion) {
      trait(limite, [], () => {
        ctx.beginPath(); ctx.moveTo(X(b.r_fusion), py0); ctx.lineTo(X(b.r_fusion), py1); ctx.stroke();
      });
      marque("fusion des brins", X(b.r_fusion), py0 + 26, limite);
    }
    marque("bande a deux brins", X(b.r_interieur), py0 + 12, limite);
  }
}

function survolCarte(evenement){
  const g = INSTRUMENT.grille, cv = document.getElementById("carte-canevas");
  if (!g) return;
  const boite = cv.getBoundingClientRect();
  const echX = cv.width / boite.width, echY = cv.height / boite.height;
  const x = (evenement.clientX - boite.left) * echX;
  const y = (evenement.clientY - boite.top) * echY;
  const pad = {g:52, d:16, h:14, b:34};
  const px0 = pad.g, px1 = cv.width - pad.d, py0 = pad.h, py1 = cv.height - pad.b;
  const rMax = g.r[g.r.length - 1], zMin = g.z[0], zMax = g.z[g.z.length - 1];
  const dedans = x >= px0 && x <= px1 && y >= py0 && y <= py1;
  const mettre = (id, v) => { document.getElementById(id).textContent = v; };
  if (!dedans) { mettre("lec-r", "—"); mettre("lec-z", "—"); mettre("lec-f", "—"); return; }
  const r = rMax * (x - px0) / (px1 - px0);
  const z = zMin + (zMax - zMin) * (py1 - y) / (py1 - py0);
  const ir = Math.max(0, Math.min(g.nr - 1, Math.round(r / rMax * (g.nr - 1))));
  const iz = Math.max(0, Math.min(g.nz - 1, Math.round((z - zMin) / (zMax - zMin) * (g.nz - 1))));
  const f = INSTRUMENT.f[iz * g.nr + ir] / 255;
  mettre("lec-r", r.toFixed(1) + " mm");
  mettre("lec-z", (z >= 0 ? "+" : "") + z.toFixed(1) + " mm");
  mettre("lec-f", f.toFixed(3));
}

/* Regime : lois de similitude. H comme n^2, Q comme n, P comme n^3, NPSHr
   comme n^2. C'est le modele lui-meme, pas une approximation -- le controle de
   similitude du rapport le verifie a 0.00 %. */
function auRegime(base, rpm, reference){
  const n = rpm / reference;
  return {Q: base.Q * n, H: base.H * n * n, P: base.P * n * n * n,
          couple: base.couple * n * n, npshr: base.npshr * n * n,
          rendement: base.rendement};
}

/* beta2 : interpolation entre trois courbes reellement calculees. */
function familleInterpolee(ecart){
  const fam = INSTRUMENT.vif.familles;
  if (fam.length === 1) return fam[0];
  const centre = fam.find(f => Math.abs(f.ecart) < 1e-9) || fam[0];
  if (Math.abs(ecart) < 1e-9) return centre;
  const cible = fam.find(f => Math.sign(f.ecart) === Math.sign(ecart) && f.ecart !== 0);
  if (!cible) return centre;
  const t = Math.min(1, Math.abs(ecart) / Math.abs(cible.ecart));
  const melange = (a, b) => a.map((v, i) => v + (b[i] - v) * t);
  return {beta2: centre.beta2 + (cible.beta2 - centre.beta2) * t,
          Q: melange(centre.Q, cible.Q), H: melange(centre.H, cible.H),
          rendement: melange(centre.rendement, cible.rendement),
          npshr: melange(centre.npshr, cible.npshr), interpolee: t > 0 && t < 1};
}

function tracerMiniplot(famille){
  const cv = document.getElementById("miniplot");
  if (!cv) return;
  const ctx = cv.getContext("2d");
  const st = getComputedStyle(document.body);
  ctx.fillStyle = st.getPropertyValue("--fond").trim(); ctx.fillRect(0, 0, cv.width, cv.height);
  const pad = 24, W = cv.width, H = cv.height;
  const toutes = INSTRUMENT.vif.familles;
  const qMax = Math.max(...toutes.flatMap(f => f.Q)) || 1;
  const hMax = Math.max(...toutes.flatMap(f => f.H)) || 1;
  const X = q => pad + (W - pad - 8) * q / qMax;
  const Y = h => H - pad - (H - pad - 8) * h / hMax;
  ctx.strokeStyle = st.getPropertyValue("--grille").trim(); ctx.lineWidth = 1;
  ctx.strokeRect(pad + .5, 8.5, W - pad - 8, H - pad - 8);
  /* La bande beta2 +/- 1 degre : ce que le rapport annonce en mots, montre. */
  ctx.save(); ctx.globalAlpha = .30;
  ctx.fillStyle = st.getPropertyValue("--declare").trim();
  ctx.beginPath();
  const bas = toutes[0], haut = toutes[toutes.length - 1];
  bas.Q.forEach((q, i) => i ? ctx.lineTo(X(q), Y(bas.H[i])) : ctx.moveTo(X(q), Y(bas.H[i])));
  for (let i = haut.Q.length - 1; i >= 0; i--) ctx.lineTo(X(haut.Q[i]), Y(haut.H[i]));
  ctx.closePath(); ctx.fill(); ctx.restore();
  ctx.strokeStyle = st.getPropertyValue("--mesure").trim(); ctx.lineWidth = 2;
  ctx.beginPath();
  famille.Q.forEach((q, i) => i ? ctx.lineTo(X(q), Y(famille.H[i])) : ctx.moveTo(X(q), Y(famille.H[i])));
  ctx.stroke();
  ctx.fillStyle = st.getPropertyValue("--encre").trim();
  ctx.font = '10px ' + st.getPropertyValue("--pile-texte");
  ctx.textAlign = "left"; ctx.fillText("H (m)", 2, 14);
  ctx.textAlign = "right"; ctx.fillText("Q (m3/h)", W - 4, H - 6);
}

function rafraichirVif(){
  const v = INSTRUMENT.vif;
  if (!v || !v.point) return;
  const rpm = INSTRUMENT.rpm;
  const famille = familleInterpolee(INSTRUMENT.dbeta);
  /* Le facteur beta2 est lu sur la courbe interpolee, puis transporte au
     regime demande par les lois de similitude. Les deux effets se composent. */
  const centre = v.familles.find(f => Math.abs(f.ecart) < 1e-9) || v.familles[0];
  const iBep = centre.H.indexOf(Math.max(...centre.H.filter((h, i) => centre.Q[i] > 0)));
  const facteurH = (famille.H[iBep] || 0) / (centre.H[iBep] || 1);
  const base = Object.assign({}, v.point, {H: v.point.H * facteurH});
  const p = auRegime(base, rpm, v.rpm_reference);

  const provenance = v.beta2_provenance;
  const classe = provenance === "declare" ? "p-declare" : "p-mesure";
  const lignes = [
    ["Debit", p.Q.toFixed(1), "m3/h", "p-mesure", "mesure"],
    ["Hauteur", p.H.toFixed(2), "m", classe, provenance],
    ["Puissance", p.P.toFixed(2), "kW", classe, provenance],
    ["Couple", p.couple.toFixed(1), "N.m", classe, provenance],
    ["NPSH requis", p.npshr.toFixed(2), "m", "p-mesure", "mesure"],
    ["beta2", famille.beta2.toFixed(2), "deg", classe, provenance],
  ];
  const corps = document.getElementById("valeurs-vives");
  corps.replaceChildren();
  for (const [nom, valeur, unite, classeLigne, src] of lignes) {
    const tr = document.createElement("tr");
    const th = document.createElement("th"); th.textContent = nom;
    const td = document.createElement("td");
    td.className = classeLigne; td.textContent = valeur + " " + unite;
    const ts = document.createElement("td"); ts.className = "src"; ts.textContent = src;
    tr.append(th, td, ts); corps.append(tr);
  }

  document.getElementById("sortie-regime").textContent = rpm.toFixed(0) + " tr/min";
  const ecart = INSTRUMENT.dbeta;
  document.getElementById("sortie-beta").textContent =
    (ecart >= 0 ? "+" : "") + ecart.toFixed(2) + " deg";

  const franchi = v.rpm_limite > 0 && rpm > v.rpm_limite;
  document.getElementById("reglage-regime").classList.toggle("franchi", franchi);
  const alerte = document.getElementById("alerte-cavitation");
  alerte.hidden = !franchi;
  if (franchi) {
    alerte.textContent = "au-dela de " + v.rpm_limite + " tr/min : limite "
      + (v.limite_active || "de cavitation") + " franchie. Le NPSH disponible ne couvre plus "
      + "le NPSH requis a ce regime.";
  }
  tracerMiniplot(famille);
}

/* Le theme sombre reste disponible, il cesse d'etre le defaut. Le choix est
   garde d'une ouverture a l'autre ; sur un fichier ouvert depuis le disque,
   localStorage peut etre indisponible, et l'echec ne doit pas casser la page. */
function monterTheme(){
  const bouton = document.getElementById("bascule-theme");
  if (!bouton) return;
  let actuel = "clair";
  try { actuel = localStorage.getItem("theme-inspecteur") || "clair"; } catch (e) {}
  const appliquer = () => {
    document.documentElement.setAttribute("data-theme", actuel);
    bouton.setAttribute("aria-pressed", String(actuel === "sombre"));
    bouton.textContent = actuel === "sombre" ? "clair" : "sombre";
    if (INSTRUMENT.grille) tracerCarte();
    if (INSTRUMENT.vif) rafraichirVif();
  };
  bouton.addEventListener("click", () => {
    actuel = actuel === "sombre" ? "clair" : "sombre";
    try { localStorage.setItem("theme-inspecteur", actuel); } catch (e) {}
    appliquer();
  });
  appliquer();
}

function monterInstrument(payload){
  const section = document.getElementById("instrument");
  if (!section || !payload.grille || !payload.grille.nr) return;
  INSTRUMENT.grille = payload.grille;
  INSTRUMENT.vif = payload.vif && payload.vif.point ? payload.vif : null;
  INSTRUMENT.f = lireOctets(payload.grille.f);
  INSTRUMENT.bascules = {axe:true, rayons:true, plans:false,
                         boucle: !!(payload.grille.reperes || {}).boucle};
  section.hidden = false;

  const boite = document.getElementById("bascules");
  boite.replaceChildren();
  const choix = [["axe", "axe"], ["rayons", "r1h r1s r2"], ["plans", "plans et sens debitant"]];
  if (INSTRUMENT.bascules.boucle) choix.push(["boucle", "brins et fusion"]);
  for (const [cle, libelle] of choix) {
    const b = document.createElement("button");
    b.type = "button"; b.className = "bascule"; b.textContent = libelle;
    b.setAttribute("aria-pressed", String(!!INSTRUMENT.bascules[cle]));
    b.addEventListener("click", () => {
      INSTRUMENT.bascules[cle] = !INSTRUMENT.bascules[cle];
      b.setAttribute("aria-pressed", String(INSTRUMENT.bascules[cle]));
      tracerCarte();
    });
    boite.append(b);
  }
  const cv = document.getElementById("carte-canevas");
  cv.addEventListener("mousemove", survolCarte);
  cv.addEventListener("mouseleave", survolCarte);
  tracerCarte();

  const reglageBeta = document.getElementById("reglage-beta");
  if (!INSTRUMENT.vif) {
    document.getElementById("reglage-regime").hidden = true;
    reglageBeta.hidden = true;
    document.getElementById("miniplot").hidden = true;
    return;
  }
  const v = INSTRUMENT.vif;
  const curseur = document.getElementById("curseur-regime");
  curseur.min = Math.round(v.rpm_min); curseur.max = Math.round(v.rpm_max_affiche);
  curseur.value = Math.round(v.rpm_reference);
  INSTRUMENT.rpm = v.rpm_reference;
  document.getElementById("borne-min").textContent = curseur.min;
  document.getElementById("borne-max").textContent = curseur.max;
  curseur.addEventListener("input", () => {
    INSTRUMENT.rpm = Number(curseur.value); rafraichirVif();
  });
  const cbeta = document.getElementById("curseur-beta");
  reglageBeta.hidden = v.familles.length < 2;
  cbeta.addEventListener("input", () => {
    INSTRUMENT.dbeta = Number(cbeta.value); rafraichirVif();
  });
  rafraichirVif();
}

function facts(rows){
  /* Quatre colonnes depuis C4 : libelle, valeur, provenance, confiance. La
     provenance se marque aussi par la graisse et un filet, jamais par la seule
     couleur ; une confiance faible ajoute son filet rouge. */
  const box = el("div", {class:"facts"});
  const classes = {mesure:"p-mesure", declare:"p-declare", defaut:"p-defaut"};
  for (const row of rows) {
    const [libelle, valeur, provenance, confiance] = row;
    let classe = "val " + (classes[provenance] || "p-mesure");
    if (confiance === "faible") classe += " c-low";
    box.append(el("div", {class:"fact"},
      el("span", {class:"k"}, libelle),
      el("span", {class:"conf"}, [provenance, confiance].filter(Boolean).join(" \u00b7 ")),
      el("span", {class:classe}, confiance === "faible" ? "(" + valeur + ")" : valeur)));
  }
  return box;
}
function step(index, title){
  return el("div", {class:"step"}, el("span", {class:"n"}, index), el("h2", {}, title));
}

function buildRail(payload){
  rail.replaceChildren();
  rail.append(el("section", {}, step("1", "Geometrie extraite"), facts(payload.tables.geometry)));

  const head = el("tr", {}, el("th", {}, ""));
  for (const speed of payload.tables.speeds) head.append(el("th", {}, speed));
  const body = el("tbody");
  for (const row of payload.tables.performance) {
    const tr = el("tr", {}, el("th", {}, row[0]));
    for (const value of row.slice(1)) tr.append(el("td", {}, value));
    body.append(tr);
  }
  rail.append(el("section", {}, step("2", "Performances"),
    el("div", {class:"scroller"}, el("table", {class:"perf"}, el("thead", {}, head), body))));

  rail.append(el("section", {}, step("3", "Cavitation"),
    el("p", {class:"eyebrow"}, "Vitesse maximale sans cavitation"),
    el("div", {class:"headline"},
       el("span", {class:"value"}, String(payload.resume.vitesse_max)),
       el("span", {class:"unit"}, "tr/min")),
    el("p", {class:"note"}, "Limite active : " + payload.resume.limite + "."),
    el("p", {class:"note"}, payload.resume.hypotheses + ".")));

  if (payload.telechargements && payload.telechargements.length) {
    const list = el("ul", {class:"fichiers"});
    for (const item of payload.telechargements) {
      list.append(el("li", {}, el("a", {href:item.url, download:item.nom}, item.nom)));
    }
    rail.append(el("section", {}, el("p", {class:"eyebrow"}, "Rapport complet"), list));
  }

  if (payload.carte) {
    rail.append(el("section", {},
      el("p", {class:"eyebrow"}, "Carte d'occupation f(r, z)"),
      el("img", {class:"map", src:payload.carte,
                 alt:"Coupe meridienne de la roue, colorisee selon la fraction angulaire occupee par la matiere"}),
      el("p", {class:"caption"},
        "Coupe meridienne : jaune = moyeu et flasque, bande sombre = pales, fond = veine fluide. "
        + "C'est la lecture dont tout le reste decoule.")));
  }

  if (payload.avertissements.length) {
    const list = el("ul", {class:"warnings"});
    for (const message of payload.avertissements) list.append(el("li", {}, message));
    rail.append(el("section", {}, el("p", {class:"eyebrow"}, "Avertissements"), list));
  }

  rail.append(el("section", {}, el("p", {class:"uncert"},
    "Modele 1D ligne moyenne : hauteur ±18 %, debit ±25 %, NPSHr ±30 %. "
    + "A verifier par essai sur banc avant toute decision d'achat ou de dimensionnement.")));
}

function buildLegend(payload){
  legende.replaceChildren(
    el("div", {class:"key"}, el("span", {class:"box", style:"background:#6B737D"}), "moyeu, flasque"),
    el("div", {class:"key"}, el("span", {class:"box", style:"background:#29A0AB"}), "pales")
  );
  if (!scene) return;
  for (const group of scene.groups) {
    const input = el("input", {type:"checkbox", checked:"checked"});
    input.addEventListener("change", () => { group.visible = input.checked; });
    legende.append(el("div", {class:"key"},
      el("label", {}, input, el("span", {class:"swatch", style:"border-color:" + group.color}), group.label)));
  }
}

function buildHeader(payload){
  titre.textContent = "Roue " + payload.resume.type + " a " + payload.resume.pales + " pales";
  fichierNom.textContent = payload.nom + " · " + payload.faces.toLocaleString("fr-FR") + " triangles";
  verdicts.replaceChildren(
    el("span", {class:"chip"},
       el("span", {class:"dot", style:"background:var(--amber)"}), "Rotation ",
       el("b", {}, payload.resume.rotation)),
    el("span", {class:"chip"}, "β1 / β2 ",
       el("b", {}, payload.resume.beta1.toFixed(1) + " / " + payload.resume.beta2.toFixed(1) + "°")),
    el("span", {class:"chip"}, "Confiance ", el("b", {}, payload.resume.confiance))
  );
}

function render(payload){
  D = payload;
  spinning = false; clip = false; spin = 0;
  buildScene(payload);
  frameCamera(payload);
  buildHeader(payload);
  buildLegend(payload);
  buildRail(payload);
  monterInstrument(payload);
  monterTheme();
  for (const b of outils.querySelectorAll("button")) b.setAttribute("aria-pressed", "false");
  spinButton.textContent = "Faire tourner";
  var sens = payload.resume.rotation_signe || payload.resume.rotation_suggeree;
  spinButton.disabled = !sens;
  spinButton.title = payload.resume.rotation_signe
    ? ""
    : (sens ? "sens seulement suggere par la geometrie : a confirmer" : "sens de rotation non renseigne");
  accueil.hidden = true;
}

/* ---------- import (mode serveur) ---------- */
if (SERVEUR) {
  const zone = document.getElementById("zone");
  const entree = document.getElementById("fichier-entree");
  const lancer = document.getElementById("lancer");
  const etat = document.getElementById("etat");
  const progres = document.getElementById("progres");
  const zoneTitre = document.getElementById("zone-titre");
  const zoneDetail = document.getElementById("zone-detail");
  let choisi = null;

  /* Ouverte depuis le disque, la page n'a aucun serveur a interroger : autant le
     dire tout de suite plutot que de laisser le navigateur repondre
     « NetworkError » apres le clic. */
  const horsServeur = location.protocol === "file:";
  if (horsServeur) document.getElementById("hors-serveur").hidden = false;

  function message(texte, erreur){
    etat.textContent = texte || "";
    etat.classList.toggle("erreur", Boolean(erreur));
  }
  function accepter(file){
    if (!file) return;
    choisi = file;
    zoneTitre.textContent = file.name;
    zoneDetail.textContent = (file.size / 1048576).toFixed(2) + " Mo — cliquez pour en choisir un autre";
    lancer.disabled = horsServeur;
    message("");
  }
  entree.addEventListener("change", () => accepter(entree.files[0]));
  for (const type of ["dragenter", "dragover"]) {
    zone.addEventListener(type, e => { e.preventDefault(); zone.classList.add("actif"); });
  }
  for (const type of ["dragleave", "drop"]) {
    zone.addEventListener(type, e => { e.preventDefault(); zone.classList.remove("actif"); });
  }
  zone.addEventListener("drop", e => accepter(e.dataTransfer.files[0]));
  accueil.addEventListener("dragover", e => e.preventDefault());
  accueil.addEventListener("drop", e => { e.preventDefault(); accepter(e.dataTransfer.files[0]); });

  function nombre(id){
    const value = document.getElementById(id).value.trim();
    return value === "" ? null : value;
  }

  document.getElementById("formulaire").addEventListener("submit", async e => {
    e.preventDefault();
    if (!choisi || running) return;
    running = true;
    lancer.disabled = true;
    progres.classList.add("actif");
    message("Analyse en cours — quelques secondes selon la finesse du maillage…");

    const params = new URLSearchParams();
    params.set("nom", choisi.name);
    params.set("unite", document.getElementById("unite").value);
    params.set("rpm", document.getElementById("rpm").value.trim());
    const sens = document.getElementById("rotation").value;
    if (sens) params.set("rotation", sens);
    for (const [id, cle] of [["r-asp","r_aspiration"], ["pales","pales"], ["beta1","beta1"],
                             ["beta2","beta2"], ["temperature","temperature"], ["altitude","altitude"],
                             ["hauteur","hauteur"], ["pertes","pertes"], ["grille","grille"]]) {
      const value = nombre(id);
      if (value !== null) params.set(cle, value);
    }
    try {
      const reponse = await fetch("/analyse?" + params.toString(),
                                  {method:"POST", body: await choisi.arrayBuffer()});
      const data = await reponse.json();
      if (!reponse.ok || data.erreur) throw new Error(data.erreur || "erreur inattendue du serveur");
      render(data);
      message("");
    } catch (error) {
      /* Un echec de fetch remonte un TypeError sans detail utile : on le
         traduit en la seule cause plausible ici, le serveur absent. */
      const injoignable = error instanceof TypeError;
      message(injoignable
        ? "le serveur local ne repond pas. Verifiez qu'il tourne encore "
          + "(python3 -m impeller_analyzer.serve) et rechargez la page."
        : String(error.message || error), true);
    } finally {
      running = false;
      lancer.disabled = false;
      progres.classList.remove("actif");
    }
  });
}

/* ---------- demarrage ---------- */
if (BOOT) {
  render(BOOT);
} else {
  accueil.hidden = false;
  rail.append(el("section", {},
    el("p", {class:"eyebrow"}, "Aucune roue chargee"),
    el("p", {class:"vide"},
      "Deposez un fichier 3D dans la zone de gauche. L'analyse s'execute sur votre machine, "
      + "par le meme code que la ligne de commande : rien n'est envoye ailleurs.")));
}
if (gl) requestAnimationFrame(draw);
</script>
"""
