"""Coupes de pale et profils deroules (SPEC 4.1).

La SPEC decrit des coupes cylindriques, ce qui est la bonne surface pour une
helice axiale : les lignes de courant y sont des cylindres `r = constante`, et
le profil deroule `(r*theta, z)` a bien une corde.

Pour une roue centrifuge cette coupe degenere : a rayon constant, une aube
radiale se reduit a un rectangle sans corde, et `tan(beta) = dr / (r dtheta)`
-- la formule que la SPEC donne elle-meme pour ce cas (4.3) -- n'est pas
calculable puisque `r` y est constant.  Ce module travaille donc sur la surface
generale dont la coupe cylindrique est un cas particulier : la **surface de
courant meridienne**, engendree par une courbe `(r, z)` prise entre le bord
moyeu et le bord carter de la zone de pales.  La coordonnee meridienne est
l'abscisse curviligne le long de cette courbe, et

    tan(beta) = dm / (r dtheta)

redonne `dz / (r dtheta)` pour une courbe verticale (helice axiale) et
`dr / (r dtheta)` pour une courbe horizontale (roue centrifuge).

Une seule mecanique d'extraction sert aux deux familles : le champ scalaire
`phi` s'annule sur la surface visee, et l'iso-zero est suivie sur le maillage.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from .. import config
from ..mesh import TriMesh
from .occupancy import OccupancyMap
from .topology import Topology, CENTRIFUGAL, clean_blade_mask

Point2 = tuple[float, float]


@dataclass
class MeridionalCurve:
    """Courbe `(r, z)` du plan meridien engendrant une surface de revolution."""

    radii: list[float] = field(default_factory=list)
    heights: list[float] = field(default_factory=list)
    arclength: list[float] = field(default_factory=list)
    vertical: bool = False  # cas particulier r = constante (cylindre)

    @classmethod
    def cylinder(cls, radius: float, z_low: float, z_high: float) -> "MeridionalCurve":
        """Cylindre `r = radius` : la coordonnee meridienne est z."""
        return cls(radii=[radius, radius], heights=[z_low, z_high], arclength=[z_low, z_high], vertical=True)

    @classmethod
    def from_points(cls, points: Sequence[Point2]) -> "MeridionalCurve":
        """Courbe generale, echantillonnee en rayon croissant."""
        ordered = sorted(points)
        radii = [p[0] for p in ordered]
        heights = [p[1] for p in ordered]
        arclength = [0.0]
        for index in range(1, len(ordered)):
            arclength.append(
                arclength[-1]
                + math.hypot(radii[index] - radii[index - 1], heights[index] - heights[index - 1])
            )
        return cls(radii=radii, heights=heights, arclength=arclength, vertical=False)

    def _bracket(self, radius: float) -> tuple[int, float]:
        """Segment contenant `radius` et fraction d'avancement dedans."""
        radii = self.radii
        if radius <= radii[0]:
            return 0, 0.0
        if radius >= radii[-1]:
            return len(radii) - 2, 1.0
        low, high = 0, len(radii) - 1
        while high - low > 1:
            mid = (low + high) // 2
            if radii[mid] <= radius:
                low = mid
            else:
                high = mid
        span = radii[low + 1] - radii[low]
        return low, (radius - radii[low]) / span if span > 0.0 else 0.0

    def level(self, radius: float, height: float) -> float:
        """Champ scalaire nul sur la surface, positif au-dessus / a l'exterieur."""
        if self.vertical:
            return radius - self.radii[0]
        index, ratio = self._bracket(radius)
        z_curve = self.heights[index] + ratio * (self.heights[index + 1] - self.heights[index])
        return height - z_curve

    def meridional(self, radius: float, height: float) -> float:
        """Coordonnee meridienne d'un point de la surface."""
        if self.vertical:
            return height
        index, ratio = self._bracket(radius)
        return self.arclength[index] + ratio * (self.arclength[index + 1] - self.arclength[index])

    def reference_radius(self) -> float:
        """Rayon representatif de la surface (milieu de sa plage de rayons)."""
        return 0.5 * (self.radii[0] + self.radii[-1])

    def radius_at(self, meridional: float) -> float:
        """Rayon du point de la surface d'abscisse curviligne donnee."""
        if self.vertical:
            return self.radii[0]
        arclength = self.arclength
        if meridional <= arclength[0]:
            return self.radii[0]
        if meridional >= arclength[-1]:
            return self.radii[-1]
        low, high = 0, len(arclength) - 1
        while high - low > 1:
            mid = (low + high) // 2
            if arclength[mid] <= meridional:
                low = mid
            else:
                high = mid
        span = arclength[low + 1] - arclength[low]
        ratio = (meridional - arclength[low]) / span if span > 0.0 else 0.0
        return self.radii[low] + ratio * (self.radii[low + 1] - self.radii[low])


@dataclass
class Profile:
    """Un profil de pale deroule dans le plan (tangentiel, meridien)."""

    tangential: list[float] = field(default_factory=list)  # t = reference * theta, en m
    meridional: list[float] = field(default_factory=list)  # m, en m
    radius: list[float] = field(default_factory=list)  # rayon local de chaque point, en m
    reference: float = 0.0  # rayon ayant servi au deroulement, en m

    def __len__(self) -> int:
        return len(self.tangential)

    def wrap(self) -> float:
        """Etendue angulaire du profil, en radians (via t = r_ref * theta)."""
        if not self.tangential:
            return 0.0
        return max(self.tangential) - min(self.tangential)

    def points(self) -> list[Point2]:
        """Le profil comme polygone 2D ferme."""
        return list(zip(self.tangential, self.meridional))


@dataclass
class Section:
    """Ensemble des profils obtenus sur une meme surface de courant."""

    span: float = 0.0  # fraction d'envergure, 0 = moyeu, 1 = carter
    reference_radius: float = 0.0
    curve: MeridionalCurve | None = None
    profiles: list[Profile] = field(default_factory=list)

    def usable(self) -> list[Profile]:
        """Profils exploitables : assez de points et enroulement raisonnable."""
        return [
            profile for profile in self.profiles
            if len(profile) >= config.MIN_PROFILE_POINTS
            and profile.reference_wrap() < config.SECTION_MAX_WRAP
        ]


def _profile_reference_wrap(self: Profile) -> float:
    """Etendue angulaire en radians du profil."""
    if not self.tangential or self.reference <= 0.0:
        return 0.0
    return self.wrap() / self.reference


Profile.reference_wrap = _profile_reference_wrap  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Construction des surfaces de coupe
# ---------------------------------------------------------------------------
def blade_zone_columns(occupancy: OccupancyMap) -> dict[int, tuple[int, int]]:
    """Pour chaque colonne radiale, la plage de lignes occupee par la zone de pales."""
    blade = clean_blade_mask(occupancy.blade_mask())
    spans: dict[int, tuple[int, int]] = {}
    for ir in range(occupancy.nr):
        rows = [iz for iz in range(occupancy.nz) if blade[iz][ir]]
        if rows:
            spans[ir] = (min(rows), max(rows))
    return spans


def blade_zone_rows(occupancy: OccupancyMap) -> dict[int, tuple[int, int]]:
    """Pour chaque ligne, la plage de colonnes occupee par la zone de pales."""
    blade = clean_blade_mask(occupancy.blade_mask())
    spans: dict[int, tuple[int, int]] = {}
    for iz in range(occupancy.nz):
        columns = [ir for ir in range(occupancy.nr) if blade[iz][ir]]
        if columns:
            spans[iz] = (min(columns), max(columns))
    return spans


def section_surfaces(
    occupancy: OccupancyMap,
    topology: Topology,
    count: int = config.N_SECTIONS,
) -> list[tuple[float, MeridionalCurve]]:
    """Les `count` surfaces de courant, du moyeu au carter (SPEC 4.1).

    Roue axiale ou mixte : cylindres `r = constante`, repartis entre le moyeu et
    le tip avec le retrait `SECTION_MARGIN` aux deux bouts, exactement comme la
    SPEC le decrit.  Roue centrifuge : courbes `z(r)` interpolees entre le bord
    moyeu et le bord carter de la zone de pales.
    """
    if topology.machine_type != CENTRIFUGAL:
        r_low = topology.r_1h if topology.r_1h > 0.0 else occupancy.r_centres[0]
        r_high = topology.r_tip
        margin = config.SECTION_MARGIN * (r_high - r_low)
        z_low, z_high = occupancy.z_min, occupancy.z_max
        surfaces = []
        for index in range(count):
            span = index / (count - 1) if count > 1 else 0.5
            radius = (r_low + margin) + span * ((r_high - margin) - (r_low + margin))
            surfaces.append((span, MeridionalCurve.cylinder(radius, z_low, z_high)))
        return surfaces

    columns = blade_zone_columns(occupancy)
    if not columns:
        return []
    surfaces = []
    for index in range(count):
        span = index / (count - 1) if count > 1 else 0.5
        inner = config.SECTION_MARGIN + span * (1.0 - 2.0 * config.SECTION_MARGIN)
        points = []
        for ir in sorted(columns):
            low, high = columns[ir]
            z_low = occupancy.z_centres[low]
            z_high = occupancy.z_centres[high]
            points.append((occupancy.r_centres[ir], z_low + inner * (z_high - z_low)))
        if len(points) >= 2:
            surfaces.append((span, MeridionalCurve.from_points(points)))
    return surfaces


# ---------------------------------------------------------------------------
# Extraction de l'iso-zero sur le maillage
# ---------------------------------------------------------------------------
def _loops_from_segments(
    segments: list[tuple[tuple[int, int], tuple[int, int]]]
) -> list[list[tuple[int, int]]]:
    """Chaine les segments en boucles fermees, via l'arete du maillage porteuse.

    Deux triangles voisins partagent l'arete traversee, donc le meme point de
    coupe : chainer par identifiant d'arete evite toute comparaison de
    coordonnees flottantes.
    """
    adjacency: dict[tuple[int, int], list[int]] = {}
    for index, (a, b) in enumerate(segments):
        adjacency.setdefault(a, []).append(index)
        adjacency.setdefault(b, []).append(index)
    used = [False] * len(segments)
    loops: list[list[tuple[int, int]]] = []
    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        chain = [segments[start][0], segments[start][1]]
        for direction in (1, 0):
            while True:
                end = chain[-1] if direction else chain[0]
                nxt = None
                for candidate in adjacency.get(end, ()):
                    if used[candidate]:
                        continue
                    a, b = segments[candidate]
                    other = b if a == end else a
                    used[candidate] = True
                    nxt = other
                    break
                if nxt is None:
                    break
                if direction:
                    chain.append(nxt)
                else:
                    chain.insert(0, nxt)
                if nxt == (chain[0] if direction else chain[-1]):
                    break
        loops.append(chain)
    return loops


def extract_section(mesh: TriMesh, curve: MeridionalCurve) -> Section:
    """Coupe le maillage par la surface de revolution et deroule les profils."""
    vertices = mesh.vertices
    scale = max(mesh.scale_length(), config.OCCUPANCY_MIN_CELL)
    eps = config.SLICE_EPS_REL * scale

    def level_at(point: Sequence[float]) -> float:
        return curve.level(math.hypot(point[0], point[1]), point[2])

    phi = []
    for v in vertices:
        value = level_at(v)
        phi.append(eps if -eps < value < eps else value)

    points: dict[tuple[int, int], tuple[float, float, float]] = {}
    segments: list[tuple[tuple[int, int], tuple[int, int]]] = []

    def crossing(i: int, j: int) -> tuple[int, int]:
        """Point de coupe sur l'arete (i, j), recale par bissection."""
        key = (i, j) if i < j else (j, i)
        if key in points:
            return key
        a, b = vertices[key[0]], vertices[key[1]]
        fa, fb = phi[key[0]], phi[key[1]]
        t = fa / (fa - fb) if fa != fb else 0.5
        low, high = 0.0, 1.0
        for _ in range(config.SECTION_REFINE_STEPS):
            point = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
            value = level_at(point)
            if (value > 0.0) == (fa > 0.0):
                low = t
            else:
                high = t
            t = 0.5 * (low + high)
        points[key] = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
        return key

    for i, j, k in mesh.faces:
        fi, fj, fk = phi[i], phi[j], phi[k]
        if (fi > 0.0) == (fj > 0.0) == (fk > 0.0):
            continue
        keys = []
        for a, b in ((i, j), (j, k), (k, i)):
            if (phi[a] > 0.0) != (phi[b] > 0.0):
                keys.append(crossing(a, b))
        if len(keys) == 2 and keys[0] != keys[1]:
            segments.append((keys[0], keys[1]))

    section = Section(reference_radius=curve.reference_radius(), curve=curve)
    for chain in _loops_from_segments(segments):
        if len(chain) < config.MIN_PROFILE_POINTS:
            continue
        raw = [points[key] for key in chain]
        radii = [math.hypot(p[0], p[1]) for p in raw]
        reference = sum(radii) / len(radii)
        angles = _unwrap([math.atan2(p[1], p[0]) for p in raw])
        profile = Profile(
            tangential=[reference * angle for angle in angles],
            meridional=[curve.meridional(radii[index], raw[index][2]) for index in range(len(raw))],
            radius=radii,
            reference=reference,
        )
        section.profiles.append(profile)
    return section


def _unwrap(angles: list[float]) -> list[float]:
    """Deroule une suite d'angles pour la rendre continue le long du contour."""
    out = [angles[0]]
    for value in angles[1:]:
        delta = value - out[-1]
        while delta > math.pi:
            delta -= 2.0 * math.pi
        while delta < -math.pi:
            delta += 2.0 * math.pi
        out.append(out[-1] + delta)
    return out


def extract_sections(
    mesh: TriMesh,
    occupancy: OccupancyMap,
    topology: Topology,
    count: int = config.N_SECTIONS,
) -> list[Section]:
    """Les `count` coupes de la roue, du moyeu au carter."""
    sections = []
    for span, curve in section_surfaces(occupancy, topology, count):
        section = extract_section(mesh, curve)
        section.span = span
        sections.append(section)
    return sections
