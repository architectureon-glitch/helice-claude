"""Ligne de cambrure, angles de pale, sens de rotation et de refoulement (SPEC 4.2 a 4.4)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import HIGH, LOW, MEDIUM, ConfidenceMap, worst
from ..numeric import linear_interp, smoothing_spline
from .sections import MeridionalCurve, Profile, Section
from .topology import AXIAL, CENTRIFUGAL, MIXED, Topology

Point2 = tuple[float, float]

COUNTERCLOCKWISE = "anti-horaire"
NOT_SUPPLIED = "a indiquer (--rotation)"
CLOCKWISE = "horaire"


@dataclass
class Camber:
    """Ligne de cambrure d'un profil deroule, dans le repere de la corde."""

    leading_edge: Point2 = (0.0, 0.0)  # (tangentiel, meridien)
    trailing_edge: Point2 = (0.0, 0.0)
    chord: float = 0.0
    stations: list[float] = field(default_factory=list)  # abscisse le long de la corde
    offsets: list[float] = field(default_factory=list)  # ecart a la corde, lisse
    thickness: list[float] = field(default_factory=list)
    slopes: list[float] = field(default_factory=list)  # dm/dt le long de la cambrure
    meridional: list[float] = field(default_factory=list)  # m de chaque station
    max_thickness: float = 0.0
    valid: bool = False


@dataclass
class SectionAngles:
    """Grandeurs de pale extraites d'une coupe."""

    span: float = 0.0
    radius: float = 0.0
    beta1_deg: float = 0.0
    beta2_deg: float = 0.0
    beta_mean_deg: float = 0.0
    chord: float = 0.0
    max_thickness: float = 0.0
    pitch: float = 0.0
    slope_sign: int = 0
    n_profiles: int = 0
    chord_to_thickness: float = 0.0
    family: int = 0  # rang de la famille dans la coupe, 0 = la plus enroulee
    n_families: int = 1
    meridional_extent: float = 0.0  # etendue en m du profil le long de la veine
    radial_span: float = 0.0  # portee radiale en m, du plus petit au plus grand rayon
    camber_ratio: float = 0.0  # fleche maximale de la ligne de cambrure / corde

    def to_dict(self) -> dict:
        """Vue serialisable en JSON (SI, angles en degres)."""
        return {
            "envergure": self.span,
            "rayon_m": self.radius,
            "beta1_deg": self.beta1_deg,
            "beta2_deg": self.beta2_deg,
            "corde_m": self.chord,
            "epaisseur_max_m": self.max_thickness,
            "pas_helicoidal_m": self.pitch,
            "profils": self.n_profiles,
            "corde_sur_epaisseur": self.chord_to_thickness,
            "famille": self.family,
            "familles": self.n_families,
            "etendue_meridienne_m": self.meridional_extent,
            "portee_radiale_m": self.radial_span,
            "cambrure_relative": self.camber_ratio,
        }


@dataclass
class BladeGeometry:
    """Resultat complet de la phase 4."""

    sections: list[SectionAngles] = field(default_factory=list)
    beta1_deg: float = 0.0
    beta2_deg: float = 0.0
    reference_radius: float = 0.0
    chord: float = 0.0
    max_thickness: float = 0.0
    pitch: float = 0.0
    wrap_consistency: float = 0.0  # enroulement mesure / enroulement implique par beta
    n_families: int = 1  # familles de profils par coupe : 1 pour une aube simple
    n_effective_blades: int = 0  # surfaces de pale vues par l'ecoulement sur un tour
    rotation_sign: int = 0  # +1 anti-horaire vu de +Z, -1 horaire, 0 non renseigne
    rotation_label: str = ""
    observed_rotation_sign: int = 0  # ce que la geometrie suggere, a titre indicatif
    observed_rotation_label: str = ""
    forced_rotation: bool = False
    helix_coefficient: float = 0.0  # k = dz/dtheta au rayon de reference, en m/rad
    forced_beta: bool = False
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "beta1_deg": self.beta1_deg,
            "beta2_deg": self.beta2_deg,
            "rayon_de_reference_m": self.reference_radius,
            "corde_m": self.chord,
            "epaisseur_max_m": self.max_thickness,
            "pas_helicoidal_m": self.pitch,
            "coherence_enroulement": self.wrap_consistency,
            "familles_de_profils": self.n_families,
            "aubes_effectives": self.n_effective_blades,
            "sens_de_rotation": self.rotation_label,
            "signe_de_rotation": self.rotation_sign,
            "sens_de_rotation_impose": self.forced_rotation,
            "sens_suggere_par_la_geometrie": self.observed_rotation_label or None,
            "signe_suggere_par_la_geometrie": self.observed_rotation_sign,
            "coefficient_helicoidal_k_m_par_rad": self.helix_coefficient,
            "angles_imposes": self.forced_beta,
            "coupes": [section.to_dict() for section in self.sections],
            "confiance": dict(self.confidence),
            "avertissements": list(self.warnings),
            "remarques": list(self.notes),
        }


# ---------------------------------------------------------------------------
# 4.2 Ligne de cambrure
# ---------------------------------------------------------------------------
def _farthest_pair(points: list[Point2]) -> tuple[int, int]:
    """Indices des deux points les plus eloignes du profil (corde maximale)."""
    best = (0, 0)
    best_distance = -1.0
    for i in range(len(points)):
        xi, yi = points[i]
        for j in range(i + 1, len(points)):
            dx = points[j][0] - xi
            dy = points[j][1] - yi
            distance = dx * dx + dy * dy
            if distance > best_distance:
                best_distance = distance
                best = (i, j)
    return best


def _arc_between(points: list[Point2], start: int, end: int) -> list[Point2]:
    """Portion du contour ferme allant de l'indice `start` a l'indice `end`."""
    out = []
    index = start
    while True:
        out.append(points[index])
        if index == end:
            break
        index = (index + 1) % len(points)
    return out


def _resample_by_arclength(arc: list[Point2], count: int) -> list[Point2]:
    """Rechantillonne une polyligne a abscisse curviligne normalisee reguliere."""
    cumulative = [0.0]
    for index in range(1, len(arc)):
        cumulative.append(cumulative[-1] + math.dist(arc[index - 1], arc[index]))
    total = cumulative[-1]
    if total <= 0.0:
        return [arc[0]] * count
    out = []
    segment = 0
    for step in range(count):
        target = total * (step + 0.5) / count
        while segment < len(cumulative) - 2 and cumulative[segment + 1] < target:
            segment += 1
        span = cumulative[segment + 1] - cumulative[segment]
        ratio = (target - cumulative[segment]) / span if span > 0.0 else 0.0
        a, b = arc[segment], arc[segment + 1]
        out.append((a[0] + ratio * (b[0] - a[0]), a[1] + ratio * (b[1] - a[1])))
    return out


def camber_line(profile: Profile, leading_edge_at_max: bool) -> Camber:
    """Ligne de cambrure d'un profil deroule (SPEC 4.2).

    Le bord d'attaque est celui des deux bouts de la corde maximale qui se
    trouve du cote de l'aspiration : `m` maximal pour une coupe cylindrique
    (m = z, aspiration vers +Z), `m` minimal pour une surface de courant
    centrifuge (m croit avec le rayon, le fluide entre au petit rayon).

    **Ecart assume a la SPEC 4.2.** La SPEC prend, a chaque station, le milieu
    du segment intrados/extrados perpendiculaire a la corde.  Ce milieu n'est
    sur la cambrure que si les deux faces y sont paralleles.  C'est le cas d'une
    pale d'helice d'epaisseur constante en z, et l'extraction retrouve alors
    l'angle a 0.01 degre pres.  Ce n'est pas le cas d'une aube centrifuge
    d'epaisseur **tangentielle** constante -- la regle de dessin la plus
    courante : son epaisseur perpendiculaire croit avec le rayon, les deux faces
    divergent, et le milieu de la coupe derive de 3 a 4 degres pres des bords,
    la ou beta1 et beta2 sont justement mesures.

    Les deux faces sont donc appariees par **abscisse curviligne normalisee**
    entre bord d'attaque et bord de fuite, ce qui est exact pour les deux
    conventions d'epaisseur (la construction reste le milieu d'un segment
    intrados/extrados, seule change la facon de l'apparier).  Les stations sont
    ensuite reportees sur la corde, et la suite du traitement -- ecretage des
    bouts, lissage, angles sur les 10 premiers et 10 derniers pour cent -- est
    celle de la SPEC.  L'epaisseur reportee est la distance entre les deux
    points apparies.
    """
    camber = Camber()
    points = profile.points()
    if len(points) < config.MIN_PROFILE_POINTS:
        return camber

    i, j = _farthest_pair(points)
    a, b = points[i], points[j]
    if (a[1] >= b[1]) == leading_edge_at_max:
        leading, trailing = a, b
        leading_index, trailing_index = i, j
    else:
        leading, trailing = b, a
        leading_index, trailing_index = j, i
    camber.leading_edge = leading
    camber.trailing_edge = trailing

    dx = trailing[0] - leading[0]
    dy = trailing[1] - leading[1]
    chord = math.hypot(dx, dy)
    if chord <= 0.0:
        return camber
    camber.chord = chord
    e1 = (dx / chord, dy / chord)
    e2 = (-e1[1], e1[0])

    side_a = _resample_by_arclength(_arc_between(points, leading_index, trailing_index), config.N_STATIONS)
    side_b = _resample_by_arclength(
        list(reversed(_arc_between(points, trailing_index, leading_index))), config.N_STATIONS
    )

    local = [
        ((p[0] - leading[0]) * e1[0] + (p[1] - leading[1]) * e1[1],
         (p[0] - leading[0]) * e2[0] + (p[1] - leading[1]) * e2[1])
        for p in points
    ]
    raw: list[tuple[float, float, float]] = []
    for pa, pb in zip(side_a, side_b):
        middle = (0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]))
        u = (middle[0] - leading[0]) * e1[0] + (middle[1] - leading[1]) * e1[1]
        v = (middle[0] - leading[0]) * e2[0] + (middle[1] - leading[1]) * e2[1]
        # L'epaisseur reste celle de la SPEC : l'etendue de la coupe
        # perpendiculaire a la corde. L'appariement par abscisse curviligne est
        # decale d'une longueur de face de bout entre les deux faces, ce qui
        # laisse la cambrure exacte mais gonflerait l'epaisseur mesuree entre
        # points apparies.
        raw.append((u, v, _cut_extent(local, u)))
    raw.sort()
    stations = [raw[0][0]]
    offsets = [raw[0][1]]
    thickness = [raw[0][2]]
    for u, v, width in raw[1:]:
        if u > stations[-1]:
            stations.append(u)
            offsets.append(v)
            thickness.append(width)
    if len(stations) < 4:
        return camber
    camber.thickness = thickness
    camber.max_thickness = max(thickness)

    kept = _trim_ends(stations, thickness, chord, camber.max_thickness)
    stations = [stations[i] for i in kept]
    offsets = [offsets[i] for i in kept]
    # L'epaisseur reportee est celle de la partie exploitable : aux deux bouts,
    # les points apparies encadrent la face de bout et leur distance est la
    # diagonale du profil, pas son epaisseur.
    camber.thickness = [thickness[i] for i in kept]
    camber.max_thickness = max(camber.thickness)

    # Lissage : le residu vise est une fraction de l'amplitude de la cambrure --
    # l'echelle propre de la grandeur lissee -- et non de la corde. Sur une aube
    # centrifuge, dont la corde deroulee vaut vingt fois l'ecart a la corde, un
    # residu cale sur la corde depasse largement le bruit reel et rabat la
    # spline vers la droite : la pente aux deux stations extremes, la ou beta1
    # et beta2 sont mesures, s'en trouve faussee de plusieurs degres.
    # L'epaisseur maximale sert de plancher pour une cambrure quasi droite.
    amplitude = max(offsets) - min(offsets)
    scale = max(amplitude, camber.max_thickness)
    target = (config.SPLINE_SMOOTH * scale) ** 2 * len(stations)
    spline = smoothing_spline(stations, offsets, target)
    smoothed = [spline(u) for u in stations]
    # La pente est prise par ajustement parabolique local sur la cambrure lissee
    # plutot que par derivation directe de la spline : une spline naturelle
    # impose une derivee seconde nulle a ses noeuds extremes, ce qui fausserait
    # la pente aux deux stations ou beta1 et beta2 sont justement mesures.
    derivatives = _local_derivatives(stations, smoothed)

    camber.stations = stations
    camber.offsets = smoothed
    camber.slopes = []
    camber.meridional = []
    for u, v, dv in zip(stations, smoothed, derivatives):
        dt = e1[0] + dv * e2[0]
        dm = e1[1] + dv * e2[1]
        camber.slopes.append(dm / dt if dt != 0.0 else math.copysign(math.inf, dm))
        camber.meridional.append(leading[1] + u * e1[1] + v * e2[1])
    camber.valid = True
    return camber


def _cut_extent(local: list[Point2], station: float) -> float:
    """Etendue de la coupe du profil perpendiculaire a la corde, a l'abscisse donnee."""
    crossings: list[float] = []
    for index in range(len(local)):
        u0, v0 = local[index]
        u1, v1 = local[(index + 1) % len(local)]
        if (u0 - station > 0.0) != (u1 - station > 0.0):
            ratio = (station - u0) / (u1 - u0)
            crossings.append(v0 + ratio * (v1 - v0))
    if len(crossings) < 2:
        return 0.0
    return max(crossings) - min(crossings)


def _local_derivatives(x: list[float], y: list[float]) -> list[float]:
    """Derivee en chaque point par la parabole passant par ses deux voisins.

    Exacte pour une parabole, y compris aux extremites et a pas irregulier.
    """
    length = len(x)
    if length < 3:
        if length == 2 and x[1] != x[0]:
            slope = (y[1] - y[0]) / (x[1] - x[0])
            return [slope, slope]
        return [0.0] * length
    out = []
    for index in range(length):
        i = min(max(index, 1), length - 2)
        x0, x1, x2 = x[i - 1], x[i], x[i + 1]
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        d01, d12 = x1 - x0, x2 - x1
        if d01 == 0.0 or d12 == 0.0 or (x2 - x0) == 0.0:
            out.append(0.0)
            continue
        t = x[index]
        out.append(
            y0 * (2.0 * t - x1 - x2) / ((x0 - x1) * (x0 - x2))
            + y1 * (2.0 * t - x0 - x2) / (-d01 * d12)
            + y2 * (2.0 * t - x0 - x1) / ((x2 - x0) * d12)
        )
    return out


def _taper_rates(stations: list[float], thickness: list[float]) -> list[float]:
    """Vitesse d'ouverture de l'epaisseur le long de la corde, par differences centrees."""
    length = len(stations)
    rates = []
    for index in range(length):
        low = max(0, index - 1)
        high = min(length - 1, index + 1)
        span = stations[high] - stations[low]
        rates.append((thickness[high] - thickness[low]) / span if span > 0.0 else 0.0)
    return rates


def _trim_ends(
    stations: list[float], thickness: list[float], chord: float, max_thickness: float
) -> list[int]:
    """Indices des stations exploitables, une fois les deux bouts ecartes.

    Aux deux extremites, les deux points apparies tombent de part et d'autre de
    la meme face de bout : l'epaisseur mesuree y varie brutalement -- elle part
    de la diagonale du profil sur un bord carre, de zero sur un bord effile --
    a une vitesse d'un ordre de grandeur superieure a l'effilement normal.
    C'est ce contraste, quel qu'en soit le signe, et non l'epaisseur elle-meme,
    qui delimite la zone a ecarter : l'epaisseur d'une aube dont beta varie du
    moyeu au carter varie elle aussi, tout a fait legitimement.
    """
    length = len(stations)
    if chord <= 0.0 or max_thickness <= 0.0:
        return list(range(length))
    limit = config.CAMBER_TRIM_TAPER * max_thickness / chord
    rates = _taper_rates(stations, thickness)
    trim = config.CAMBER_END_TRIM * max_thickness
    first = 0
    while first < length and (stations[first] < trim or abs(rates[first]) > limit):
        first += 1
    last = length - 1
    while last > first and (stations[last] > chord - trim or abs(rates[last]) > limit):
        last -= 1
    kept = list(range(first, last + 1))
    if len(kept) < config.CAMBER_MIN_STATIONS:
        return list(range(length))
    return kept


def beta_windows(camber: Camber) -> tuple[list[int], list[int]]:
    """Stations retenues pour beta1 (10 premiers % de corde) et beta2 (10 derniers).

    Les stations de la cambrure ont deja ete ecretees aux deux bouts par
    `camber_line` ; les fenetres sont donc prises depuis les extremites de la
    partie conservee.
    """
    stations = camber.stations
    window = config.BETA_CHORD_FRACTION * camber.chord
    start, end = stations[0], stations[-1]
    head = [i for i, u in enumerate(stations) if u <= start + window] or [0]
    tail = [i for i, u in enumerate(stations) if u >= end - window] or [len(stations) - 1]
    return head, tail


# ---------------------------------------------------------------------------
# 4.3 Angles
# ---------------------------------------------------------------------------
def _clamp_beta(value_deg: float) -> float:
    """Borne un angle de pale aux limites numeriques du modele."""
    return max(config.BETA_MIN_DEG, min(config.BETA_MAX_DEG, value_deg))


def _beta_from_slope(slope: float, curve: MeridionalCurve, reference_radius: float, meridional: float) -> float:
    """Angle de pale, en degres, depuis la direction tangentielle.

    `tan(beta) = dm / (r dtheta)`.  Le profil est deroule avec
    `t = reference * theta`, donc `dtheta = dt / reference` et
    `tan(beta) = (dm/dt) * reference / r`.  `reference_radius` doit etre
    exactement le rayon ayant servi au deroulement du profil, sans quoi le
    rapport introduit un biais systematique ; pour une coupe cylindrique
    `r = reference` et l'on retrouve `dz / (r dtheta)`.
    """
    radius = curve.radius_at(meridional)
    if radius <= 0.0:
        radius = reference_radius
    return math.degrees(math.atan(abs(slope) * reference_radius / radius))


def section_angle_families(
    section: Section, leading_edge_at_max: bool, n_blades: int
) -> list[tuple[SectionAngles, list[Camber]]]:
    """Angles de pale d'une coupe, une entree par famille de profils.

    Une aube en boucle donne plusieurs familles sur une meme coupe ; les
    moyenner ensemble donnerait un angle qui ne decrit aucune des deux.
    """
    families = section.families(n_blades)
    out = []
    for rank, profiles in enumerate(families):
        angles, cambers = section_angles(section, leading_edge_at_max, profiles)
        angles.family = rank
        angles.n_families = len(families)
        out.append((angles, cambers))
    return out


def section_angles(
    section: Section, leading_edge_at_max: bool, profiles: list[Profile] | None = None
) -> tuple[SectionAngles, list[Camber]]:
    """Angles de pale d'une coupe, medianes sur les profils retenus."""
    result = SectionAngles(span=section.span, radius=section.reference_radius)
    usable = section.usable() if profiles is None else list(profiles)
    pairs = [(profile, camber_line(profile, leading_edge_at_max)) for profile in usable]
    pairs = [(profile, camber) for profile, camber in pairs if camber.valid]
    cambers = [camber for _, camber in pairs]
    result.n_profiles = len(cambers)
    if not pairs or section.curve is None:
        return result, cambers

    beta1_values: list[float] = []
    beta2_values: list[float] = []
    beta_mean_values: list[float] = []
    slope_signs: list[int] = []
    for profile, camber in pairs:
        # Le rayon de deroulement du profil est celui qui doit servir a la
        # conversion dt -> dtheta : tout autre rayon biaise beta.
        reference = profile.reference or section.reference_radius
        # La moyenne sur la fenetre rend l'angle au milieu de celle-ci, ce qui
        # rabat beta1 et beta2 vers la moyenne et aplatit le vrillage : le biais
        # va de 0.7 degre sur une aube 20/25 a 4.6 degres sur une aube 40/65.
        # L'evaluer au bord de la fenetre par une droite des moindres carres le
        # reduit d'un degre au-dela de 50 degres, mais en ajoute deux sous 20 --
        # ou sont les aubes de pompe. La moyenne reste donc, et le biais est
        # documente plutot que deplace.
        head_indices, tail_indices = beta_windows(camber)
        beta1_values.append(
            _mean(
                _beta_from_slope(camber.slopes[i], section.curve, reference, camber.meridional[i])
                for i in head_indices
            )
        )
        beta2_values.append(
            _mean(
                _beta_from_slope(camber.slopes[i], section.curve, reference, camber.meridional[i])
                for i in tail_indices
            )
        )
        beta_mean_values.append(
            _mean(
                _beta_from_slope(slope, section.curve, reference, m)
                for slope, m in zip(camber.slopes, camber.meridional)
            )
        )
        slope_signs.append(1 if _mean(camber.slopes) >= 0.0 else -1)

    result.beta1_deg = _clamp_beta(_median(beta1_values))
    result.beta2_deg = _clamp_beta(_median(beta2_values))
    result.beta_mean_deg = _clamp_beta(_median(beta_mean_values))
    result.chord = _median([c.chord for c in cambers])
    result.max_thickness = _median([c.max_thickness for c in cambers])
    result.chord_to_thickness = result.chord / result.max_thickness if result.max_thickness > 0.0 else math.inf
    result.pitch = 2.0 * math.pi * result.radius * math.tan(math.radians(result.beta_mean_deg))
    result.slope_sign = 1 if sum(slope_signs) >= 0 else -1
    result.meridional_extent = _median([
        max(profile.meridional) - min(profile.meridional) for profile, _ in pairs
    ])
    result.radial_span = _median([
        max(profile.radius) - min(profile.radius) for profile, _ in pairs
    ])
    # Fleche relative de la ligne de cambrure. C'est elle qui porte la portance a
    # incidence nulle : pour un profil mince en arc de cercle, la theorie des
    # profils minces donne alpha_0 = -2 h/c et donc Cl(0) = 4 pi h/c. Une lecture
    # geometrique, la ou un modele d'helice courant demande a l'utilisateur de
    # taper un Cl de dessin qu'il ne connait pas.
    result.camber_ratio = _median([
        max(abs(offset) for offset in camber.offsets) / camber.chord
        for camber in cambers
        if camber.offsets and camber.chord > 0.0
    ])
    return result, cambers


def _cut_profile(
    polygon: list[Point2], origin: Point2, normal: Point2
) -> tuple[Point2, float] | None:
    """Coupe le contour ferme par la droite passant par `origin` de normale `normal`.

    Renvoie le milieu des deux intersections encadrant `origin` et l'epaisseur
    mesuree, ou `None` si la droite ne traverse pas proprement le profil.
    """
    nx, ny = normal
    # Coordonnee le long de la droite, et cote par rapport a elle.
    hits: list[tuple[float, Point2]] = []
    for index in range(len(polygon)):
        ax, ay = polygon[index]
        bx, by = polygon[(index + 1) % len(polygon)]
        side_a = (ax - origin[0]) * nx + (ay - origin[1]) * ny
        side_b = (bx - origin[0]) * nx + (by - origin[1]) * ny
        if (side_a > 0.0) == (side_b > 0.0):
            continue
        ratio = side_a / (side_a - side_b) if side_a != side_b else 0.5
        point = (ax + ratio * (bx - ax), ay + ratio * (by - ay))
        along = (point[0] - origin[0]) * (-ny) + (point[1] - origin[1]) * nx
        hits.append((along, point))
    if len(hits) < 2:
        return None
    hits.sort()
    # On retient l'intersection immediatement de part et d'autre de l'origine :
    # sur un profil fortement cambre, la droite peut recouper le contour plus
    # loin, et prendre les extremes reviendrait a mesurer l'enveloppe.
    after = next((h for h in hits if h[0] >= 0.0), None)
    before = next((h for h in reversed(hits) if h[0] < 0.0), None)
    if after is None or before is None:
        low, high = hits[0], hits[-1]
    else:
        low, high = before, after
    centre = (0.5 * (low[1][0] + high[1][0]), 0.5 * (low[1][1] + high[1][1]))
    return centre, abs(high[0] - low[0])


def _mean(values) -> float:
    """Moyenne arithmetique d'un iterable non vide."""
    data = list(values)
    return sum(data) / len(data) if data else 0.0


def _median(values) -> float:
    """Mediane, robuste a une pale aberrante."""
    data = sorted(values)
    if not data:
        return 0.0
    middle = len(data) // 2
    if len(data) % 2:
        return data[middle]
    return 0.5 * (data[middle - 1] + data[middle])


# ---------------------------------------------------------------------------
# 4.4 Sens de rotation
# ---------------------------------------------------------------------------
def _add_twist_note(geometry: BladeGeometry) -> None:
    """Dit que la lecture aplatit le vrillage, et de combien.

    beta1 et beta2 sont pris comme la moyenne sur les dix premiers et dix
    derniers pour cent de corde ; une moyenne de fenetre rend la valeur au
    milieu de celle-ci, pas a son bord. La lecture rabat donc les deux angles
    vers la moyenne : beta1 ressort trop grand, beta2 trop petit, et l'ecart
    croit avec le vrillage -- 0.7 degre sur une aube 20/25, 4.6 sur une aube
    40/65. Le biais n'est pas corrige : une correction demanderait de caler une
    loi d'aube, et celle des roues de synthese n'est pas celle des roues
    reelles. Il est donc chiffre et dit, pour que la marge soit connue.
    """
    twist = geometry.beta2_deg - geometry.beta1_deg
    if twist < config.TWIST_NOTE_DEG or geometry.forced_beta:
        return
    geometry.notes.append(
        f"vrillage lu : {twist:.1f} degres du bord d'attaque au bord de fuite. La lecture prend "
        f"chaque angle en moyenne sur un dixieme de corde, ce qui rabat les deux extremites vers "
        f"la moyenne et restitue {config.TWIST_RECOVERY_MIN * 100.0:.0f} a "
        f"{config.TWIST_RECOVERY_MAX * 100.0:.0f} % du vrillage reel sur des roues d'angles "
        f"connus. Le vrillage reel est donc plutot de "
        f"{twist / config.TWIST_RECOVERY_MAX:.1f} a {twist / config.TWIST_RECOVERY_MIN:.1f} "
        f"degres : beta1 un peu plus petit et beta2 un peu plus grand que ceux du tableau. Biais "
        "mesure et non corrige."
    )


def wrap_consistency(
    sections: list[Section], topology: Topology, beta1_deg: float, beta2_deg: float
) -> float:
    """Enroulement mesure rapporte a celui qu'impliquent les angles lus.

    Une aube centrifuge dont l'angle vaut `beta` s'enroule de
    `theta = ln(r2/r1) / tan(beta)` : c'est l'integration de
    `dtheta = dr / (r tan(beta))`. Le rapport entre l'enroulement reellement
    mesure sur les profils et celui-la vaut donc un, a la variation de `beta` le
    long de l'aube pres.

    Loin de un, les deux lectures se contredisent, et c'est le signe que les
    coupes n'ont pas rendu de vrais profils -- typiquement des **fragments**,
    quand la surface de courant effleure l'aube au lieu de la traverser. Le
    controle ne coute rien et ne suppose aucun seuil arbitraire : il confronte
    la mesure a elle-meme.

    L'enroulement retenu est le **plus grand** des profils, non leur mediane :
    la question posee est « une vraie coupe d'aube existe-t-elle ? », et une
    seule suffit a repondre oui. La mediane, elle, se laisse noyer par les
    fragments des qu'ils sont nombreux -- sur un maillage decime a 20 % elle
    tombe a 7 degres pendant que les vrais profils en font 148, et le controle
    condamnerait une lecture pourtant juste a un demi-degre pres.

    Renvoie 0 quand le controle n'est pas applicable.
    """
    if topology.machine_type not in (CENTRIFUGAL, MIXED):
        return 0.0
    r_1s, r_2 = topology.r_1s, topology.r_2
    if r_1s <= 0.0 or r_2 <= r_1s:
        return 0.0
    beta = math.radians(_clamp_beta(0.5 * (beta1_deg + beta2_deg)))
    if math.tan(beta) <= 0.0:
        return 0.0
    attendu = math.log(r_2 / r_1s) / math.tan(beta)
    mesures = [
        profile.reference_wrap()
        for section in sections
        for profile in section.usable()
    ]
    if not mesures or attendu <= 0.0:
        return 0.0
    return max(mesures) / attendu


def forced_rotation_warning(geometry: "BladeGeometry") -> str:
    """Reserve a emettre quand le sens impose contredit celui que lit la geometrie.

    Chaine vide s'il n'y a rien a dire : sens non impose, aucune suggestion, ou
    les deux d'accord.  A appeler **apres** que la lecture des angles est
    arretee : sur une aube en boucle, la suggestion vient des normales et non de
    la cambrure, et la citer trop tot nommerait le sens d'une lecture ecartee.
    """
    if not geometry.forced_rotation or not geometry.observed_rotation_sign:
        return ""
    if geometry.observed_rotation_sign == geometry.rotation_sign:
        return ""
    return (
        "le sens impose est l'inverse de ce que suggere la geometrie "
        f"({geometry.observed_rotation_label}). C'est le sens impose qui est retenu ; "
        "verifiez qu'il correspond bien a la piece, la suggestion pouvant se tromper "
        "sur une aube quasi radiale ou une roue a la frontiere de deux familles."
    )


def rotation_sense(machine_type: str, slope_sign: int) -> int:
    """Signe de omega, +1 anti-horaire vu de +Z (cote aspiration), -1 horaire.

    Cas axial et mixte : la nappe moyenne est localement helicoidale,
    `z = z0 + k*theta` ; une telle surface tournant a omega translate a
    `v_z = -k*omega`.  Refouler vers -Z impose `v_z < 0`, donc
    `signe(omega) = signe(k)`, et `k` a le signe de la pente `dm/dt` de la
    cambrure.

    Cas centrifuge : les aubes de pompe sont incurvees vers l'arriere, l'aube
    fuit le sens de rotation quand r augmente, donc
    `signe(omega) = -signe(dtheta/dr)`.  Le long de la cambrure `dtheta/dr` a le
    signe de la pente `dm/dt` (m croit avec r), d'ou le signe oppose.
    """
    if slope_sign == 0:
        return 0
    return slope_sign if machine_type in (AXIAL, MIXED) else -slope_sign


def rotation_label(sign: int) -> str:
    """Libelle du sens de rotation, vu du cote aspiration (+Z)."""
    if sign > 0:
        return f"{COUNTERCLOCKWISE} (vu de +Z, cote aspiration)"
    if sign < 0:
        return f"{CLOCKWISE} (vu de +Z, cote aspiration)"
    return NOT_SUPPLIED


def rotation_sign_from_name(name: str | None) -> int | None:
    """Traduit `horaire` / `antihoraire` en signe, `None` si rien n'est demande."""
    if name is None:
        return None
    key = name.strip().lower().replace("-", "").replace("_", "")
    if key in ("horaire", "cw", "sensHoraire".lower()):
        return -1
    if key in ("antihoraire", "ccw", "trigonometrique", "trigo"):
        return 1
    raise ValueError(f"sens de rotation inconnu : {name!r}")


# ---------------------------------------------------------------------------
# Chaine complete
# ---------------------------------------------------------------------------
def analyse(
    sections: list[Section],
    topology: Topology,
    forced_beta1_deg: float | None = None,
    forced_beta2_deg: float | None = None,
    forced_rotation: int | None = None,
) -> BladeGeometry:
    """Extraction des angles de pale et du sens de rotation (SPEC 4.2 a 4.4)."""
    geometry = BladeGeometry()
    leading_edge_at_max = topology.machine_type != CENTRIFUGAL

    per_section: list[SectionAngles] = []
    family_counts: list[int] = []
    for section in sections:
        families = [
            angles
            for angles, _ in section_angle_families(
                section, leading_edge_at_max, topology.blades.n_blades
            )
            if angles.n_profiles
        ]
        if not families:
            continue
        family_counts.append(len(families))
        # La famille retenue pour le modele 1D est celle de plus grande portee
        # **radiale** : c'est la surface qui conduit l'ecoulement de l'ouie au
        # refoulement. Les autres familles d'une coupe sont des accidents
        # locaux -- sur la roue toroidale de reference, le bourrelet ou les deux
        # brins fusionnent, qui ne couvre que les seize derniers millimetres de
        # rayon et donnerait un angle sans rapport avec le guidage.
        per_section.append(max(families, key=lambda angles: angles.radial_span))
    geometry.sections = per_section
    geometry.n_families = int(_median(family_counts)) if family_counts else 1
    geometry.n_effective_blades = topology.blades.n_blades * geometry.n_families
    if geometry.n_families > 1:
        geometry.notes.append(
            f"{geometry.n_families} familles de profils par coupe : l'ecoulement voit "
            f"{geometry.n_effective_blades} surfaces de pale par tour pour "
            f"{topology.blades.n_blades} aubes. beta1 et beta2 sont lus sur la famille de plus "
            "grande portee **radiale**, celle qui conduit l'ecoulement de l'ouie au refoulement ; "
            "les autres sont des accidents locaux. Le glissement, lui, se calcule au refoulement, "
            f"ou les brins d'une meme aube ont fusionne : il garde {topology.blades.n_blades} "
            "passages."
        )

    if not per_section:
        geometry.confidence.set("angles_de_pale", LOW)
        geometry.confidence.set("sens_de_rotation", LOW)
        geometry.warnings.append(
            "aucune coupe exploitable : les angles de pale n'ont pas pu etre extraits. "
            "Fournissez --beta1, --beta2 et --blades pour poursuivre le calcul."
        )
    else:
        if topology.machine_type == CENTRIFUGAL:
            # Les coupes sont des surfaces de courant reperees par l'envergure :
            # le modele 1D retient la surface a mi-envergure.
            spans = [s.span for s in per_section]
            geometry.beta1_deg = linear_interp(spans, [s.beta1_deg for s in per_section], 0.5)
            geometry.beta2_deg = linear_interp(spans, [s.beta2_deg for s in per_section], 0.5)
            geometry.reference_radius = topology.r_1
            geometry.notes.append(
                "roue centrifuge : les coupes sont des surfaces de courant meridiennes ; "
                "beta1 et beta2 du modele 1D sont pris a mi-envergure"
            )
        else:
            radii = [s.radius for s in per_section]
            geometry.beta1_deg = linear_interp(radii, [s.beta1_deg for s in per_section], topology.r_1)
            geometry.beta2_deg = linear_interp(radii, [s.beta2_deg for s in per_section], topology.r_1)
            geometry.reference_radius = topology.r_1
        geometry.chord = _median([s.chord for s in per_section])
        geometry.max_thickness = _median([s.max_thickness for s in per_section])
        geometry.pitch = 2.0 * math.pi * geometry.reference_radius * math.tan(
            math.radians(_clamp_beta(0.5 * (geometry.beta1_deg + geometry.beta2_deg)))
        )

        slope_sign = 1 if sum(s.slope_sign for s in per_section) >= 0 else -1
        geometry.observed_rotation_sign = rotation_sense(topology.machine_type, slope_sign)
        geometry.observed_rotation_label = rotation_label(geometry.observed_rotation_sign)
        geometry.helix_coefficient = slope_sign * geometry.reference_radius * math.tan(
            math.radians(_clamp_beta(geometry.beta2_deg))
        )

        geometry.wrap_consistency = wrap_consistency(
            sections, topology, geometry.beta1_deg, geometry.beta2_deg
        )
        _add_twist_note(geometry)
        geometry.confidence.set("angles_de_pale", _angle_confidence(per_section, geometry))
        geometry.confidence.set(
            "sens_de_rotation", _rotation_confidence(topology, geometry, per_section)
        )
        _add_rotation_notes(topology, geometry)

    # Le sens de rotation est une **entree**, pas un resultat. La geometrie le
    # suggere, mais elle ne le tranche pas : sur une roue a rapport r2/r1s pose a
    # cheval sur la frontiere mixte / centrifuge, les deux familles donnent des
    # sens opposes, et sur une aube quasi radiale la lecture n'a aucune marge.
    # L'utilisateur, lui, a la piece sous les yeux.
    if forced_rotation:
        geometry.rotation_sign = 1 if forced_rotation > 0 else -1
        geometry.forced_rotation = True
        geometry.confidence.set("sens_de_rotation", HIGH)
        geometry.notes.append(
            f"sens de rotation impose par l'utilisateur : {rotation_label(geometry.rotation_sign)}"
        )
        # La comparaison entre le sens impose et celui que suggere la geometrie
        # n'est pas faite ici : sur une aube en boucle, la suggestion est reprise
        # sur les normales apres coup, et une comparaison faite maintenant
        # citerait celle de la cambrure -- ecartee. Elle est rendue par
        # `forced_rotation_warning`, que l'analyse appelle une fois la lecture
        # des angles arretee.
    else:
        geometry.rotation_sign = 0
        geometry.confidence.set("sens_de_rotation", LOW)
        if geometry.observed_rotation_sign:
            geometry.notes.append(
                f"la geometrie suggere {geometry.observed_rotation_label}, a titre indicatif "
                "seulement : le sens retenu doit etre donne par --rotation"
            )
    geometry.rotation_label = rotation_label(geometry.rotation_sign)

    if forced_beta1_deg is not None or forced_beta2_deg is not None:
        if forced_beta1_deg is not None:
            geometry.beta1_deg = _clamp_beta(float(forced_beta1_deg))
        if forced_beta2_deg is not None:
            geometry.beta2_deg = _clamp_beta(float(forced_beta2_deg))
        geometry.forced_beta = True
        geometry.confidence.set("angles_de_pale", HIGH)
        geometry.warnings.append(
            f"angles de pale imposes par l'utilisateur : beta1 = {geometry.beta1_deg:.1f} deg, "
            f"beta2 = {geometry.beta2_deg:.1f} deg"
        )
    return geometry


def _angle_confidence(per_section: list[SectionAngles], geometry: BladeGeometry) -> str:
    """Confiance des angles, d'apres le nombre et la qualite des coupes."""
    level = HIGH
    if len(per_section) < config.MIN_BLADE_SECTIONS:
        level = LOW
        geometry.warnings.append(
            f"seules {len(per_section)} coupes exploitables sur {config.N_SECTIONS} : "
            "angles de pale peu fiables"
        )
    elif len(per_section) < config.N_SECTIONS:
        level = MEDIUM
        geometry.warnings.append(
            f"{len(per_section)} coupes exploitables sur {config.N_SECTIONS} demandees"
        )
    worst_ratio = min(s.chord_to_thickness for s in per_section)
    if worst_ratio < config.CHORD_THICKNESS_MIN:
        level = worst(level, MEDIUM)
        geometry.warnings.append(
            f"profil mal conditionne (corde / epaisseur = {worst_ratio:.1f} < "
            f"{config.CHORD_THICKNESS_MIN}) : la ligne de cambrure est peu significative"
        )
    spread = max(s.beta2_deg for s in per_section) - min(s.beta2_deg for s in per_section)
    if spread > config.BETA_MAX_DEG / 2.0:
        level = worst(level, MEDIUM)
        geometry.warnings.append(
            f"beta2 varie de {spread:.0f} deg entre moyeu et carter : vrillage tres marque, "
            "le modele 1D au rayon quadratique moyen est une approximation grossiere"
        )
    return level


def _rotation_confidence(
    topology: Topology, geometry: BladeGeometry, per_section: list[SectionAngles]
) -> str:
    """Confiance du sens de rotation (SPEC 4.4)."""
    level = geometry.confidence.get_level("angles_de_pale", HIGH)
    if topology.machine_type == CENTRIFUGAL:
        if geometry.beta2_deg > config.BETA2_RADIAL_DEG:
            level = LOW
            geometry.warnings.append(
                f"beta2 = {geometry.beta2_deg:.0f} deg > {config.BETA2_RADIAL_DEG} deg : aubes quasi "
                "radiales, le sens de rotation est ambigu"
            )
        signs = {s.slope_sign for s in per_section}
        if len(signs) > 1:
            level = LOW
            geometry.warnings.append(
                "le sens de courbure des aubes change d'une coupe a l'autre : "
                "sens de rotation indetermine"
            )
    else:
        if geometry.pitch > config.HELIX_PITCH_RATIO_MAX * geometry.reference_radius:
            level = worst(level, MEDIUM)
            geometry.warnings.append(
                f"pas helicoidal {geometry.pitch * config.MM_PER_M:.0f} mm superieur a "
                f"{config.HELIX_PITCH_RATIO_MAX} fois le rayon : nappe presque plane, "
                "le sens de rotation deduit est peu marque"
            )
    return level


def _add_rotation_notes(topology: Topology, geometry: BladeGeometry) -> None:
    """Remarques a reporter avec le sens de rotation."""
    if topology.machine_type == CENTRIFUGAL:
        geometry.notes.append(
            "la suggestion geometrique suppose des aubes incurvees vers l'arriere, cas de la "
            "quasi-totalite des pompes. Pour des aubes incurvees vers l'avant (rare en pompe, "
            "courant sur un ventilateur a cage d'ecureuil), elle est a inverser."
        )
    else:
        geometry.notes.append(
            "sens deduit du signe de k = dz/dtheta : une nappe helicoidale tournant a omega "
            "translate a v_z = -k*omega, et le refoulement vers -Z impose signe(omega) = signe(k)."
        )


def discharge_direction(
    topology: Topology,
    geometry: BladeGeometry,
    cm2: float = 0.0,
    cu2: float = 0.0,
) -> dict:
    """Sens de sortie du liquide (SPEC 4.4, dernier point).

    `cm2` et `cu2` viennent de la phase 5 ; sans eux seule la composante
    meridienne est renseignee.
    """
    if topology.machine_type == CENTRIFUGAL:
        meridional = "radial (+r)"
        meridional_angle = 0.0
    elif topology.machine_type == AXIAL:
        meridional = "axial (-Z)"
        meridional_angle = 90.0
    else:
        # Roue mixte : inclinaison de la vitesse meridienne, deduite de la
        # geometrie de la section de sortie.
        axial_extent = max(topology.b_2, 0.0)
        radial_extent = max(topology.r_2s - topology.r_2h, 0.0)
        meridional_angle = math.degrees(math.atan2(axial_extent, radial_extent)) if (axial_extent or radial_extent) else config.MIXED_DISCHARGE_DEG
        meridional = f"mixte, {meridional_angle:.0f} deg entre -Z et +r"

    alpha2 = math.degrees(math.atan2(cm2, cu2)) if cu2 > 0.0 else None
    return {
        "composante_meridienne": meridional,
        "inclinaison_meridienne_deg": meridional_angle,
        "composante_tangentielle": f"dans le sens de rotation ({geometry.rotation_label})",
        "cu2_m_par_s": cu2,
        "cm2_m_par_s": cm2,
        "alpha2_deg": alpha2,
    }
