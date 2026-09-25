"""Angles de pale lus sur une pale isolee (SPEC v2, 5.2 ; cahier C3.3 et C3.4).

En import global, la cambrure se lit mal : la coupe traverse plusieurs pales,
qu'il faut segmenter, et une aube en boucle la traverse deux fois. Une pale
importee **seule** leve les deux difficultes : chaque surface de courant la
coupe selon un profil, dont la ligne moyenne donne l'angle de pale.

Surfaces de courant. Sur une roue a refoulement radial, l'eau traverse les
aubes en s'eloignant de l'axe, dans des plans quasi perpendiculaires a celui-ci :
la surface de courant est le plan `z = constante`, et le profil se lit dans le
plan `(r, theta)`. L'angle de pale, mesure depuis la direction tangentielle, vaut

    beta = atan2(1, -s . r . dtheta/dr)

ou `s` est le signe de la rotation (+1 anti-horaire vu de +Z). Une aube courbee
vers l'arriere -- son bout trainant derriere son pied -- a beta < 90 degres ;
une aube radiale, 90 ; une aube courbee vers l'avant, plus de 90.

Aubes en boucle. Les niveaux se regroupent en **brins** selon le sens de leur
recul : une aube conventionnelle en a un, une aube en boucle deux, qui reculent
en sens opposes et se rejoignent pres du plan ou ils changent de sens. Chaque
brin est traite comme une grille d'aubes a part entiere, avec son bord d'attaque,
son bord de fuite et sa deviation ; ils ne sont pas fusionnes. Les niveaux de la
jonction, ou le profil revient sur lui-meme, sont ecartes : on n'y lit pas de
ligne moyenne.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import LOW, MEDIUM
from ..mesh import TriMesh
from .meridian import MeridianWalls

Vec3 = tuple[float, float, float]


@dataclass
class LevelProfile:
    """Profil de la pale dans un plan `z` : ligne moyenne et angles de bord."""

    z: float = 0.0
    r_min: float = 0.0
    r_max: float = 0.0
    sweep: float = 0.0  # dtheta/dr moyen, rad/m
    beta_le_deg: float = 0.0  # au bord interieur (attaque, ecoulement centrifuge)
    beta_te_deg: float = 0.0  # au bord exterieur (fuite)
    usable: bool = True
    discharges: bool | None = None  # un chemin mene-t-il du bord de fuite a la fente de sortie
    free_le: bool | None = None  # le bout interieur est-il un bord d'attaque libre (et non un pied)
    embedded: bool = False  # niveau noye dans une paroi : pied soude, pas une surface de courant
    camber: list = field(default_factory=list, repr=False)  # ligne moyenne (r, theta), theta a 2 pi pres


@dataclass
class Branch:
    """Un brin : suite de niveaux qui reculent dans le meme sens."""

    name: str = ""
    levels: list[LevelProfile] = field(default_factory=list)
    beta_le_deg: float = 0.0
    beta_te_deg: float = 0.0
    backward: bool = True  # courbe vers l'arriere pour le sens de rotation retenu

    @property
    def du(self) -> str:
        """Le nom precede de sa preposition : « de l'aube », « du brin bas »."""
        return "de l'aube" if self.name == "aube" else f"du {self.name}"

    @property
    def z_mean(self) -> float:
        return sum(level.z for level in self.levels) / len(self.levels) if self.levels else 0.0

    @property
    def discharging(self) -> int:
        """Nombre de niveaux dont le bord de fuite debouche vers la sortie."""
        return sum(1 for level in self.levels if level.discharges)

    @property
    def deviation_deg(self) -> float:
        """Deviation de la grille : ecart entre angle de fuite et angle d'attaque."""
        return self.beta_te_deg - self.beta_le_deg

    def to_dict(self) -> dict:
        return {
            "brin": self.name,
            "niveaux": len(self.levels),
            "z_moyen_m": self.z_mean,
            "beta_attaque_deg": self.beta_le_deg,
            "beta_fuite_deg": self.beta_te_deg,
            "deviation_deg": self.deviation_deg,
            "courbe_vers_l_arriere": self.backward,
            "niveaux_debouchant_vers_la_sortie": self.discharging,
        }


@dataclass
class AxialPush:
    """Sens ou un brin pousse l'eau le long de l'axe, la ou elle le traverse.

    Une aube dont l'azimut varie avec la hauteur est une vis : en tournant, elle
    pousse l'eau le long de l'axe dans le sens `-s . signe(dtheta/dz)`. C'est ce
    qui decide si un brin place sur l'arrivee de l'eau la conduit vers le brin
    qui refoule, ou la renvoie vers l'oeillard.
    """

    branch: str = ""
    passage: tuple[float, float] = (0.0, 0.0)  # m, rayons ou l'eau traverse le brin le long de l'axe
    toward_working: float = 0.0  # part des rayons ou il pousse l'eau vers le brin qui refoule
    helix_deg: float = 0.0  # angle moyen de la vis, depuis la tangente : petit, vis serree
    samples: int = 0
    upstream: bool = False  # le brin est entre l'entree et le brin qui refoule : l'eau le rencontre d'abord
    r_rms: float = 0.0  # m, rayon quadratique moyen du passage
    helix_rms_deg: float = 0.0  # angle de la vis a ce rayon : l'angle d'entree de la roue

    @property
    def area(self) -> float:
        """Section du passage, couronne entre ses deux rayons, m2."""
        return math.pi * (self.passage[1] ** 2 - self.passage[0] ** 2)

    def to_dict(self) -> dict:
        return {
            "brin": self.branch,
            "passage_r_m": list(self.passage),
            "part_poussant_vers_le_brin_refoulant": self.toward_working,
            "angle_d_helice_deg": self.helix_deg,
            "rayon_quadratique_moyen_m": self.r_rms,
            "angle_d_helice_au_rayon_moyen_deg": self.helix_rms_deg,
            "rayons_sondes": self.samples,
            "en_amont_du_brin_refoulant": self.upstream,
        }

    @property
    def is_inlet(self) -> bool:
        """Vrai si ce brin est l'entree de la roue : en amont, et poussant l'eau vers l'aval."""
        return self.upstream and self.toward_working >= config.AXIAL_PUSH_AGREEMENT


def _theta_at(camber: list[tuple[float, float]], r: float) -> float | None:
    """Azimut de la ligne moyenne au rayon `r` (interpolation lineaire)."""
    for (r0, t0), (r1, t1) in zip(camber, camber[1:]):
        if r0 <= r <= r1 and r1 > r0:
            return t0 + (t1 - t0) * (r - r0) / (r1 - r0)
    return None


def axial_twist(levels: list[LevelProfile], r: float) -> float | None:
    """dtheta/dz de la ligne moyenne au rayon `r`, sur les niveaux qui l'atteignent, rad/m."""
    points: list[tuple[float, float]] = []
    for level in sorted(levels, key=lambda l: l.z):
        theta = _theta_at(level.camber, r)
        if theta is None:
            continue
        if points:
            previous = points[-1][1]
            theta = previous + (theta - previous + math.pi) % (2.0 * math.pi) - math.pi
        points.append((level.z, theta))
    if len(points) < config.BRANCH_MIN_LEVELS:
        return None
    return _slope(points)


@dataclass
class IsolatedBladeAngles:
    """Angles publies, et ce qui les porte."""

    beta1_deg: float = 0.0
    beta2_deg: float = 0.0
    branches: list[Branch] = field(default_factory=list)
    levels: list[LevelProfile] = field(default_factory=list)
    working: Branch | None = None  # le brin qui refoule dans la fente de sortie
    r_le: float = 0.0  # m, rayon moyen du bord d'attaque du brin refoulant
    r_te: float = 0.0  # m, rayon moyen du bord de fuite du brin refoulant
    r_te_range: tuple[float, float] = (0.0, 0.0)
    r_le_range: tuple[float, float] = (0.0, 0.0)
    level_step: float = 0.0  # m, ecart entre deux plans de coupe : l'epaisseur d'une surface de courant
    le_height: float = 0.0  # m, hauteur du bord d'attaque : niveaux qui le portent x level_step
    te_height: float = 0.0  # m, hauteur du bord de fuite
    beta1_spread: tuple[float, float] = (0.0, 0.0)
    beta2_spread: tuple[float, float] = (0.0, 0.0)
    layout: str = ""  # disposition des brins : parallele, un_seul_refoule, indeterminee
    layout_detail: str = ""
    axial_push: AxialPush | None = None  # brin non refoulant : sens ou il pousse l'eau le long de l'axe
    walls: bool = False  # les parois ont-elles servi a dire ou l'eau sort
    confidence: str = LOW
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "beta1_deg": self.beta1_deg,
            "beta1_min_max_deg": list(self.beta1_spread),
            "beta2_deg": self.beta2_deg,
            "brin_refoulant": self.working.name if self.working else None,
            "rayon_d_attaque_m": self.r_le or None,
            "rayon_de_fuite_m": self.r_te or None,
            "rayon_de_fuite_min_max_m": list(self.r_te_range),
            "rayon_d_attaque_min_max_m": list(self.r_le_range),
            "hauteur_du_bord_d_attaque_m": self.le_height,
            "hauteur_du_bord_de_fuite_m": self.te_height,
            "beta2_min_max_deg": list(self.beta2_spread),
            "disposition_des_brins": self.layout or None,
            "poussee_axiale_du_brin_non_refoulant": (
                self.axial_push.to_dict() if self.axial_push else None
            ),
            "parois_consultees": self.walls,
            "brins": [branch.to_dict() for branch in self.branches],
            "niveaux_lus": sum(1 for level in self.levels if level.usable),
            "niveaux_ecartes": sum(1 for level in self.levels if not level.usable),
            "confiance": self.confidence,
        }


def _slope(pairs: list[tuple[float, float]]) -> float:
    """Pente de la regression lineaire de y sur x."""
    n = float(len(pairs))
    if n < 2:
        return 0.0
    mx = sum(x for x, _ in pairs) / n
    my = sum(y for _, y in pairs) / n
    var = sum((x - mx) ** 2 for x, _ in pairs)
    if var <= 0.0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in pairs) / var


def blade_angle(r: float, dtheta_dr: float, rotation_sign: int) -> float:
    """Angle de pale depuis la tangente, en degres : < 90 courbee vers l'arriere."""
    return math.degrees(math.atan2(1.0, -rotation_sign * r * dtheta_dr))


def _loops(segments, tolerance: float) -> list[list[tuple[float, float]]]:
    """Chaine les segments d'une coupe en contours fermes."""
    def key(p):
        return (round(p[0] / tolerance), round(p[1] / tolerance))

    neighbours: dict = {}
    where: dict = {}
    for a, b in segments:
        ka, kb = key(a), key(b)
        if ka == kb:
            continue
        neighbours.setdefault(ka, []).append(kb)
        neighbours.setdefault(kb, []).append(ka)
        where[ka], where[kb] = a, b
    seen: set = set()
    loops = []
    for start in neighbours:
        if start in seen:
            continue
        loop, previous, current = [start], None, start
        seen.add(start)
        while True:
            options = [k for k in neighbours[current] if k != previous and k not in seen]
            if not options:
                break
            previous, current = current, options[0]
            seen.add(current)
            loop.append(current)
        if len(loop) >= 3:
            loops.append([where[k] for k in loop])
    return loops


def _side_theta(side: list[tuple[float, float]], r: float) -> float | None:
    """Azimut d'une face du profil au rayon `r` (interpolation lineaire)."""
    for (r0, t0), (r1, t1) in zip(side, side[1:]):
        if (r0 - r) * (r1 - r) <= 0.0 and r0 != r1:
            return t0 + (t1 - t0) * (r - r0) / (r1 - r0)
    return None


def level_profile(mesh: TriMesh, axis: Vec3, z: float, rotation_sign: int) -> LevelProfile | None:
    """Ligne moyenne du profil de la pale dans le plan `z`, et ses angles de bord.

    Les segments de coupe sont chaines en un contour ferme, coupe en deux faces
    entre son point le plus proche de l'axe (attaque) et le plus eloigne
    (fuite). A chaque rayon, la ligne moyenne est au milieu des deux faces. Une
    regression sur les points de coupe melangeait intrados, extrados et face
    d'extremite : sur une aube de 4 mm a faible angle, la meme pale et son image
    miroir donnaient 18 et 35 degres. Les angles de bord sont lus sur la ligne
    moyenne, entre `CAMBER_EDGE_START` et `CAMBER_EDGE_END` de l'etendue radiale
    depuis chaque bord : l'arrondi du bord lui-meme n'en dit rien. L'angle local
    y est ajuste en droite et **extrapole au bord** : lu au milieu de la zone,
    il glissait vers l'angle de l'autre bord d'un huitieme de l'ecart
    beta2 - beta1, soit 5 degres sur beta2 pour la roue d'essai hel1.
    """
    segments = []
    for a, b, c in mesh.triangles():
        cut = []
        for p, q in ((a, b), (b, c), (c, a)):
            dp, dq = p[2] - z, q[2] - z
            if (dp > 0.0) != (dq > 0.0):
                t = dp / (dp - dq)
                cut.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
        if len(cut) == 2:
            segments.append((cut[0], cut[1]))
    if len(segments) < config.ISOLATED_MIN_POINTS:
        return None
    lo, hi = mesh.bounds()
    loops = _loops(segments, 1e-7 * max(hi[0] - lo[0], hi[1] - lo[1], 1e-9))
    if not loops:
        return None
    loop = max(loops, key=len)
    polar = []
    for x, y in loop:
        r, theta = math.hypot(x - axis[0], y - axis[1]), math.atan2(y - axis[1], x - axis[0])
        if polar:
            previous = polar[-1][1]
            theta = previous + (theta - previous + math.pi) % (2.0 * math.pi) - math.pi
        polar.append((r, theta))
    i_le = min(range(len(polar)), key=lambda i: polar[i][0])
    i_te = max(range(len(polar)), key=lambda i: polar[i][0])
    r_min, r_max = polar[i_le][0], polar[i_te][0]
    span = r_max - r_min
    if span <= 0.0:
        return None
    n = len(polar)
    # Deux faces, chacune de l'attaque vers la fuite.
    face_a = [polar[(i_le + k) % n] for k in range((i_te - i_le) % n + 1)]
    face_b = [polar[(i_le - k) % n] for k in range((i_le - i_te) % n + 1)]
    # Le contour n'encercle pas l'axe : les deux faces partagent la meme
    # determination de l'azimut au bord d'attaque.
    shift = face_b[0][1] - face_a[0][1]
    face_b = [(r, t - shift) for r, t in face_b]

    camber = []
    steps = config.CAMBER_STATIONS
    for k in range(1, steps):
        r = r_min + span * k / steps
        ta, tb = _side_theta(face_a, r), _side_theta(face_b, r)
        if ta is None or tb is None:
            return None
        camber.append((r, 0.5 * (ta + tb)))
    profile = LevelProfile(z=z, r_min=r_min, r_max=r_max, sweep=_slope(camber), camber=camber)
    # Angle local de la ligne moyenne, par differences centrees.
    local = [
        (b[0], blade_angle(b[0], (c[1] - a[1]) / (c[0] - a[0]), rotation_sign))
        for a, b, c in zip(camber, camber[1:], camber[2:])
    ]
    near_le = [p for p in local if config.CAMBER_EDGE_START * span <= p[0] - r_min <= config.CAMBER_EDGE_END * span]
    near_te = [p for p in local if config.CAMBER_EDGE_START * span <= r_max - p[0] <= config.CAMBER_EDGE_END * span]
    if len(near_le) < 2 or len(near_te) < 2:
        return None
    profile.beta_le_deg = _extrapolate(near_le, r_min)
    profile.beta_te_deg = _extrapolate(near_te, r_max)
    return profile


def _extrapolate(pairs: list[tuple[float, float]], x: float) -> float:
    """Valeur en `x` de la droite des moindres carres de y sur x."""
    mx = sum(p[0] for p in pairs) / len(pairs)
    my = sum(p[1] for p in pairs) / len(pairs)
    return my + _slope(pairs) * (x - mx)


def read_isolated_blade(
    mesh: TriMesh,
    axis: Vec3,
    rotation_sign: int,
    outlet_z: tuple[float, float],
    inlet_z: float,
    toroidal: bool,
    walls: MeridianWalls | None = None,
    outlet_radius: float = 0.0,
) -> IsolatedBladeAngles:
    """Angles d'une pale isolee a refoulement radial.

    `outlet_z` : etendue axiale de la fente de sortie, `outlet_radius` son rayon.
    `inlet_z` : cote du plan d'entree. `walls` : les parois de la roue dans le
    plan meridien ; avec elles, on sait quels bords de fuite debouchent vers la
    fente et quels bouts de profil sont des bords libres plutot que des pieds.
    """
    result = IsolatedBladeAngles()
    lo, hi = mesh.bounds()
    count = config.ISOLATED_LEVELS
    result.level_step = (hi[2] - lo[2]) / count
    for k in range(count):
        z = lo[2] + (k + 0.5) * (hi[2] - lo[2]) / count
        profile = level_profile(mesh, axis, z, rotation_sign)
        if profile is not None:
            result.levels.append(profile)
    if not result.levels:
        result.warnings.append("aucune coupe exploitable de la pale isolee : angles non lus")
        return result

    # Niveaux sans ligne moyenne lisible :
    # - les pointes de la boucle, profils de quelques millimetres dont les
    #   bords ne portent pas d'angle ;
    # - la jonction des brins, ou le profil fait un crochet : les deux niveaux
    #   qui encadrent un changement de sens de recul, et tout niveau dont le
    #   recul moyen s'effondre.
    # Niveaux noyes dans une paroi : la pale soudee plonge dans ses flasques, et
    # la coupe y rend son pied, pas une surface de courant. Sur la roue classique
    # de 8 pouces essayee, deux niveaux sur vingt-quatre, qui ajoutaient 2,5 mm a
    # la hauteur du bord de fuite.
    embedded = 0
    if walls is not None:
        for level in result.levels:
            radii = [level.r_min + (level.r_max - level.r_min) * (k + 0.5) / config.EMBEDDED_SAMPLES
                     for k in range(config.EMBEDDED_SAMPLES)]
            if sum(walls.wall(r, level.z) for r in radii) > config.EMBEDDED_FRACTION * len(radii):
                level.usable = False
                level.embedded = True
                embedded += 1
    widest = max(level.r_max - level.r_min for level in result.levels)
    typical = sorted(abs(level.sweep) for level in result.levels)[len(result.levels) // 2]
    for level in result.levels:
        if (level.r_max - level.r_min) < config.TIP_SPAN_FRACTION * widest:
            level.usable = False
        if abs(level.sweep) < config.JUNCTION_SWEEP_FRACTION * typical:
            level.usable = False
    for first, second in zip(result.levels, result.levels[1:]):
        if (first.sweep > 0.0) != (second.sweep > 0.0):
            first.usable = second.usable = False
    usable = [level for level in result.levels if level.usable]
    if not usable:
        result.warnings.append("aucun niveau de la pale n'a de recul net : angles non lus")
        return result

    # Brins : suites contigues de niveaux de meme sens de recul.
    for level in usable:
        sign = 1 if level.sweep > 0.0 else -1
        if not result.branches or (1 if result.branches[-1].levels[-1].sweep > 0.0 else -1) != sign:
            result.branches.append(Branch())
        result.branches[-1].levels.append(level)
    result.branches = [b for b in result.branches if len(b.levels) >= config.BRANCH_MIN_LEVELS]
    ordered = sorted(result.branches, key=lambda b: b.z_mean)
    for index, branch in enumerate(ordered):
        if len(ordered) == 1:
            branch.name = "aube"
        elif len(ordered) == 2:
            branch.name = ("brin bas", "brin haut")[index]
        else:
            branch.name = f"brin {index + 1}"
        branch.beta_le_deg = sum(l.beta_le_deg for l in branch.levels) / len(branch.levels)
        branch.beta_te_deg = sum(l.beta_te_deg for l in branch.levels) / len(branch.levels)
        branch.backward = -rotation_sign * branch.levels[0].sweep > 0.0
    if not result.branches:
        result.warnings.append("aucun brin assez etendu pour y lire des angles")
        return result

    if toroidal and len(result.branches) < 2:
        result.warnings.append(
            f"pale declaree toroidale, mais ses coupes ne montrent qu'un seul sens de recul : "
            "les deux brins d'une boucle devraient reculer en sens opposes. Les angles sont lus "
            "sur l'aube comme sur une aube conventionnelle."
        )
    if not toroidal and len(result.branches) > 1:
        result.warnings.append(
            f"pale declaree conventionnelle, mais ses coupes montrent {len(result.branches)} sens "
            "de recul opposes : c'est la signature d'une aube en boucle."
        )

    # Ou l'eau sort. Avec les parois, on cherche pour chaque niveau un chemin
    # de son bord de fuite a la fente, qui peut longer une paroi mais ne revient
    # jamais vers l'axe ; sans elles, on s'en tient a la hauteur de la fente.
    def in_slot(level: LevelProfile) -> bool:
        return outlet_z[0] <= level.z <= outlet_z[1]

    meridian = walls is not None and outlet_radius > 0.0
    result.walls = meridian
    for branch in result.branches:
        for level in branch.levels:
            if meridian:
                level.discharges = walls.reaches(level.r_max, level.z, outlet_radius, outlet_z)
                level.free_le = walls.edge_is_free(level.r_min, level.z, inward=True)
            else:
                level.discharges = in_slot(level)

    # Le brin qui refoule : celui dont le plus de niveaux debouchent. C'est lui,
    # et lui seul, que le modele 1D peut decrire.
    working = max(result.branches, key=lambda b: (
        b.discharging, sum(1 for l in b.levels if in_slot(l)), -abs(b.z_mean - inlet_z)))
    result.working = working
    low_confidence = False

    # Bord de fuite et bord d'attaque. Un plan horizontal n'est une surface de
    # courant que la ou l'ecoulement est radial : pres de l'oeillard, le bout
    # d'une coupe peut etre l'arete de la pale contre un flasque, et non un
    # bord. Avec les parois, beta2 se lit sur les niveaux qui debouchent, beta1
    # sur ceux dont le bout interieur est libre -- un bord incline est alors lu
    # sur toute sa hauteur. Sans elles, sur les niveaux qui atteignent le rayon
    # de sortie (resp. d'entree) de la pale.
    r_te_max = max(l.r_max for l in working.levels)
    r_le_min = min(l.r_min for l in working.levels)
    exits = [l for l in working.levels if l.r_max >= (1.0 - config.EDGE_REACH) * r_te_max]
    entries = [l for l in working.levels if l.r_min <= (1.0 + config.EDGE_REACH) * r_le_min]
    if meridian:
        debouchant = [l for l in working.levels if l.discharges]
        libres = [l for l in working.levels if l.free_le]
        if debouchant:
            exits = debouchant
        else:
            low_confidence = True
            result.warnings.append(
                f"aucun niveau {working.du} ne trouve de chemin vers la fente de sortie : les "
                "parois barrent tout. beta2 est lu au rayon de sortie de la pale, sans garantie "
                "que l'eau y passe."
            )
        if libres:
            entries = libres
    result.r_te = sum(l.r_max for l in exits) / len(exits)
    result.r_te_range = (min(l.r_max for l in exits), max(l.r_max for l in exits))
    result.r_le = sum(l.r_min for l in entries) / len(entries)
    result.r_le_range = (min(l.r_min for l in entries), max(l.r_min for l in entries))
    result.le_height = len(entries) * result.level_step
    result.te_height = len(exits) * result.level_step
    result.beta2_deg = sum(l.beta_te_deg for l in exits) / len(exits)
    result.beta1_deg = sum(l.beta_le_deg for l in entries) / len(entries)
    working.beta_le_deg, working.beta_te_deg = result.beta1_deg, result.beta2_deg
    # Le vrillage se juge sur les seuls niveaux qui portent le bord d'attaque :
    # un plan qui coupe la pale plus loin de l'axe y lit un angle d'aval, et
    # une aube sans vrillage passait pour vrillee de 20 degres.
    le = [l.beta_le_deg for l in entries]
    te = [l.beta_te_deg for l in exits]
    result.beta1_spread = (min(le), max(le))
    result.beta2_spread = (min(te), max(te))
    mm = config.MM_PER_M
    if max(le) - min(le) > config.LE_TWIST_WARN_DEG:
        low_confidence = True
        result.warnings.append(
            f"le bord d'attaque {working.du} est tres vrille : beta d'attaque de "
            f"{min(le):.0f} a {max(le):.0f} deg selon la hauteur"
            + (f", de r = {min(l.r_min for l in entries) * mm:.0f} a "
               f"{max(l.r_min for l in entries) * mm:.0f} mm" if len(entries) > 1 else "")
            + f". Son angle d'attaque moyen, {result.beta1_deg:.1f} deg, ne represente pas "
            "toute la hauteur : aucun debit d'adaptation unique n'y convient."
        )
    if max(te) - min(te) > config.LE_TWIST_WARN_DEG:
        low_confidence = True
        result.warnings.append(
            f"le bord de fuite {working.du} est tres vrille : beta de fuite de "
            f"{min(te):.0f} a {max(te):.0f} deg selon la hauteur, de r = "
            f"{result.r_te_range[0] * mm:.0f} a {result.r_te_range[1] * mm:.0f} mm. beta2 = "
            f"{result.beta2_deg:.1f} deg en est la moyenne sur la hauteur de la veine."
        )
    if outlet_radius > 0.0 and result.r_te < (1.0 - config.EDGE_REACH) * outlet_radius:
        result.warnings.append(
            f"le bord de fuite {working.du} est en moyenne a r = {result.r_te * mm:.1f} mm, la "
            f"fente de sortie a {outlet_radius * mm:.1f} mm : entre les deux, une couronne sans "
            "aube, ou l'eau conserve son moment cinetique. Pour la ligne moyenne, r2 est le "
            "rayon du bord de fuite, pas celui de la fente."
        )

    # Disposition des brins : par ou l'eau les traverse.
    others = [b for b in result.branches if b is not working]
    if others:
        def sens(branch: Branch) -> str:
            return "l arriere" if branch.backward else "l avant"

        courbure = "".join(
            f" Pour le sens de rotation retenu, le {working.name} est courbe vers {sens(working)}"
            f", le {other.name} vers {sens(other)}."
            for other in others[:1] if other.backward != working.backward
        )
        if not meridian:
            result.layout = "indeterminee"
            result.layout_detail = (
                f"aube en boucle a {len(result.branches)} brins : sans piece de paroi (moyeu ou "
                "corps), l'outil ne peut dire si l'eau les traverse en serie ou en parallele. "
                f"Les angles publies sont ceux du {working.name}, qui occupe le plus la hauteur "
                f"de la fente de sortie.{courbure}"
            )
        elif all(2 * b.discharging >= len(b.levels) for b in others):
            result.layout = "parallele"
            result.layout_detail = (
                f"aube en boucle a {len(result.branches)} brins, qui debouchent tous vers la "
                f"fente de sortie : l'eau les traverse en parallele.{courbure} Le modele de "
                f"ligne moyenne ne decrit qu'une grille d'aubes ; il est applique au "
                f"{working.name}, dont le plus de niveaux debouchent, et l'apport des autres "
                "brins ne s'y calcule pas."
            )
        else:
            result.layout = "un_seul_refoule"
            muets = [b for b in others if 2 * b.discharging < len(b.levels)]
            detail = (
                f"aube en boucle a {len(result.branches)} brins, dont un seul refoule : le "
                f"{working.name} debouche vers la fente de sortie, le "
                f"{' et le '.join(b.name for b in muets)} non -- les parois ferment sa "
                "chambre sur l'exterieur. Il est sur le chemin de l'eau venue de l'oeillard, "
                f"qui le traverse pour gagner le {working.name} : en serie, en amont.{courbure}"
            )
            push = _axial_push(muets[0], working, walls, rotation_sign)
            if push is not None:
                push.upstream = abs(muets[0].z_mean - inlet_z) < abs(working.z_mean - inlet_z)
            result.axial_push = push
            if push is not None:
                a, b = (x * mm for x in push.passage)
                ou = f"la ou l'eau le traverse le long de l'axe (r = {a:.0f} a {b:.0f} mm)"
                if push.toward_working >= config.AXIAL_PUSH_AGREEMENT:
                    detail += (
                        f" {ou.capitalize()}, son pas helicoidal pousse l'eau vers le "
                        f"{working.name} pour le sens de rotation retenu : il s'oppose a son "
                        f"retour vers l'oeillard (vis a {push.helix_deg:.0f} deg de la tangente)."
                    )
                elif push.toward_working <= 1.0 - config.AXIAL_PUSH_AGREEMENT:
                    low_confidence = True
                    detail += (
                        f" {ou.capitalize()}, son pas helicoidal repousse l'eau vers l'oeillard "
                        "pour le sens de rotation retenu, a contre-sens du debit. Verifiez le "
                        f"sens de rotation : dans l'autre, il la pousserait vers le {working.name}."
                    )
                else:
                    detail += (
                        f" {ou.capitalize()}, son pas helicoidal change de sens : il pousse l'eau "
                        f"vers le {working.name} sur {push.toward_working:.0%} des rayons "
                        "seulement."
                    )
                if walls.rotating:
                    detail += (
                        f" Au-dela de r = {b:.0f} mm, sa chambre, fermee et entrainee avec la "
                        "roue, tourne en bloc avec l'eau qu'elle contient : le brin n'y fait pas "
                        "travailler l'eau."
                    )
            if push is not None and push.is_inlet:
                detail += (
                    " Le modele de ligne moyenne traite les deux brins comme une seule roue : "
                    f"l'eau y entre par la vis du {push.branch}, sur le passage, et en sort par le "
                    f"bord de fuite du {working.name}. Ses angles de profil, lus dans des plans "
                    "que l'eau n'y suit pas, sont indicatifs."
                )
            else:
                detail += (
                    f" Le modele de ligne moyenne ne decrit que le {working.name} : l'effet de "
                    "l'autre brin n'y est pas compte, et ses angles de profil, lus dans des plans "
                    "que l'eau n'y suit pas, sont indicatifs."
                )
            result.layout_detail = detail

    result.confidence = LOW if low_confidence else MEDIUM
    for branch in result.branches:
        result.notes.append(
            f"{branch.name} ({len(branch.levels)} niveaux, z moyen "
            f"{branch.z_mean * mm:.1f} mm) : beta d'attaque "
            f"{branch.beta_le_deg:.1f} deg, de fuite {branch.beta_te_deg:.1f} deg, deviation "
            f"{branch.deviation_deg:+.1f} deg ; courbe vers "
            f"{'l arriere' if branch.backward else 'l avant'} pour le sens de rotation retenu"
            + (f" ; {branch.discharging} niveau(x) sur {len(branch.levels)} debouchent vers la "
               "fente de sortie." if meridian else ".")
        )
    ecartes = sum(1 for level in result.levels if not level.usable and not level.embedded)
    if ecartes:
        result.notes.append(
            f"{ecartes} niveau(x) ecarte(s) -- pointes de la boucle ou jonction des brins : le "
            "profil y est trop court ou revient sur lui-meme, aucune ligne moyenne ne s'y lit."
        )
    if embedded:
        result.notes.append(
            f"{embedded} niveau(x) noye(s) dans une paroi, ecarte(s) : la pale y plonge dans son "
            "flasque, et la coupe rend son pied soude, pas une surface de courant."
        )
    return result


def _axial_push(
    branch: Branch, working: Branch, walls: MeridianWalls, rotation_sign: int
) -> AxialPush | None:
    """Sens ou `branch` pousse l'eau le long de l'axe, la ou elle le traverse.

    L'eau traverse le brin non refoulant le long de l'axe pour gagner le brin
    qui refoule. Le passage est lu dans le plan qui separe les deux brins : les
    rayons que n'y occupe aucune paroi. Sur hel1, l'ouverture du disque
    intermediaire, de r = 35 a 62 mm.
    """
    if branch.z_mean > working.z_mean:
        z_cut = 0.5 * (max(l.z for l in working.levels) + min(l.z for l in branch.levels))
    else:
        z_cut = 0.5 * (min(l.z for l in working.levels) + max(l.z for l in branch.levels))
    r_lo = min(l.r_min for l in branch.levels)
    r_hi = max(l.r_max for l in branch.levels)
    count = config.AXIAL_PUSH_SAMPLES
    radii = [r_lo + (k + 0.5) * (r_hi - r_lo) / count for k in range(count)]
    passage = [r for r in radii if not walls.wall(r, z_cut)]
    toward = 1 if working.z_mean > branch.z_mean else -1
    pushes, helix = [], []
    for r in passage:
        twist = axial_twist(branch.levels, r)
        if not twist:
            continue
        pushes.append(-rotation_sign * (1 if twist > 0.0 else -1) == toward)
        helix.append(math.degrees(math.atan2(1.0, r * abs(twist))))
    if not pushes:
        return None
    push = AxialPush(
        branch=branch.name,
        passage=(min(passage), max(passage)),
        toward_working=sum(pushes) / len(pushes),
        helix_deg=sum(helix) / len(helix),
        samples=len(pushes),
    )
    push.r_rms = math.sqrt(0.5 * (push.passage[0] ** 2 + push.passage[1] ** 2))
    twist = axial_twist(branch.levels, push.r_rms)
    push.helix_rms_deg = (
        math.degrees(math.atan2(1.0, push.r_rms * abs(twist))) if twist else push.helix_deg
    )
    return push
