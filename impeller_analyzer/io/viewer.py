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
from ..mesh import TriMesh, cross, normalize, sub
from . import plot

#: Classes de facette rendues dans la vue.
CLASS_SOLID = 0  # moyeu, flasque : f >= F_SOLIDE juste sous la surface
CLASS_BLADE = 1  # pale : F_VIDE < f < F_SOLIDE
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
        elif value > config.F_VIDE:
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
    if blades is not None and blades.rotation_sign:
        groups.append(
            {
                "id": "rotation",
                "label": blades.rotation_label,
                "color": "#F0A202",
                "points": _arc_arrow(
                    topology.r_tip * 1.18, z_high + 0.12 * height, blades.rotation_sign
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


def _fill(payload_json: str, server: bool) -> str:
    """Injecte la charge utile, le mode et le port par defaut dans le gabarit."""
    return (
        _TEMPLATE.replace("__DONNEES__", payload_json)
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
        '<style>html,body{margin:0}img{max-width:100%}</style>\n'
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
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
:root{
  --paper:#F5F3EF; --panel:#FFFFFF; --ink:#1B1D21; --ink-soft:#5F636B; --ink-faint:#8B9099;
  --line:#E1DCD3; --line-strong:#CFC8BC;
  --accent:#0E7C86; --accent-wash:#E4F0F1;
  --oxide:#A9491A; --amber:#8F5E06; --steel:#6E747E;
  --stage:#181B21; --stage-grid:#242932; --stage-ink:#C9CDD4;
  --ok:#2F6E4F; --warn:#8A5A12; --erreur:#A32C1E;
  --rail:22.5rem;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --paper:#0F1216; --panel:#171A20; --ink:#E8E5DF; --ink-soft:#A0A6AF; --ink-faint:#767D87;
    --line:#252A33; --line-strong:#333944;
    --accent:#45B8C2; --accent-wash:#10272B;
    --oxide:#DE7F43; --amber:#E2AA36; --steel:#8B929C;
    --ok:#5CB68A; --warn:#D9A441; --erreur:#E8705E;
  }
}
:root[data-theme="dark"]{
  --paper:#0F1216; --panel:#171A20; --ink:#E8E5DF; --ink-soft:#A0A6AF; --ink-faint:#767D87;
  --line:#252A33; --line-strong:#333944;
  --accent:#45B8C2; --accent-wash:#10272B;
  --oxide:#DE7F43; --amber:#E2AA36; --steel:#8B929C;
  --ok:#5CB68A; --warn:#D9A441; --erreur:#E8705E;
}

*{box-sizing:border-box}
body{margin:0}
.app{
  height:100vh; overflow:hidden; background:var(--paper); color:var(--ink);
  font-family:"IBM Plex Sans", ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
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
.bar .file{font-family:"IBM Plex Mono", ui-monospace, monospace; font-size:.78rem; color:var(--ink-soft)}
.verdicts{display:flex; flex-wrap:wrap; gap:.4rem .5rem; margin-left:auto}
.chip{
  display:inline-flex; align-items:center; gap:.4rem; padding:.2rem .55rem;
  border:1px solid var(--line-strong); border-radius:2px; font-size:.78rem; color:var(--ink-soft);
  background:var(--paper);
}
.chip b{color:var(--ink); font-weight:600}
.chip .dot{width:.5rem; height:.5rem; border-radius:50%; flex:none}

/* ---- corps ---- */
main{flex:1; min-height:0; display:grid; grid-template-columns:1fr var(--rail)}
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
  font:inherit; font-size:.85rem; font-family:"IBM Plex Mono", ui-monospace, monospace;
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
  display:inline-block; font-family:"IBM Plex Mono", ui-monospace, monospace; font-size:.8rem;
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
  font-family:"IBM Plex Mono", monospace; font-size:.7rem; color:var(--accent);
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
  grid-column:1 / -1; font-family:"IBM Plex Mono", ui-monospace, monospace;
  font-variant-numeric:tabular-nums; font-size:.92rem; color:var(--ink); overflow-wrap:anywhere;
}
table{width:100%; border-collapse:collapse; font-size:.82rem}
th,td{text-align:left; padding:.3rem 0; vertical-align:baseline}
th{font-weight:500; color:var(--ink-soft); font-size:.78rem}
.scroller{overflow-x:auto}
.perf th:first-child{width:45%}
.perf td{
  font-family:"IBM Plex Mono", ui-monospace, monospace; font-variant-numeric:tabular-nums;
  text-align:right; white-space:nowrap; padding-left:.7rem;
}
.perf thead th{border-bottom:1px solid var(--line-strong); padding-bottom:.35rem; text-align:right}
.perf thead th:first-child{text-align:left}
.perf tbody tr + tr th, .perf tbody tr + tr td{border-top:1px solid var(--line)}

.headline{display:flex; align-items:baseline; gap:.5rem; margin:.2rem 0 .5rem}
.headline .value{
  font-family:"IBM Plex Mono", monospace; font-size:2rem; font-weight:500;
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
.fichiers a{color:var(--accent); text-decoration:none; font-family:"IBM Plex Mono", ui-monospace, monospace}
.fichiers a:hover{text-decoration:underline}
.vide{color:var(--ink-soft); font-size:.85rem}
@media (prefers-reduced-motion:reduce){ *{animation:none !important; transition:none !important} }
</style>

<div class="app">
  <header class="bar">
    <h1 id="titre">Inspecteur de roue</h1>
    <span class="file" id="fichier"></span>
    <div class="verdicts" id="verdicts"></div>
  </header>

  <main>
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
  if (spinning) spin += dt * 1.1 * (D.resume.rotation_signe || 1);

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

function facts(rows){
  const box = el("div", {class:"facts"});
  for (const row of rows) {
    box.append(el("div", {class:"fact"},
      el("span", {class:"k"}, row[0]),
      el("span", {class:"conf"}, row[2] || ""),
      el("span", {class:"val"}, row[1])));
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
  for (const b of outils.querySelectorAll("button")) b.setAttribute("aria-pressed", "false");
  spinButton.textContent = "Faire tourner";
  spinButton.disabled = !payload.resume.rotation_signe;
  spinButton.title = payload.resume.rotation_signe ? "" : "sens de rotation indetermine";
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
