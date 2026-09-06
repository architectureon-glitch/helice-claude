"""Generateurs de geometries synthetiques a reponse analytique connue.

Ces primitives servent aux tests unitaires de toutes les phases et a la campagne
de validation de la phase 8 : chaque geometrie est construite a partir de
grandeurs imposees (nombre de pales, angles de pale, rayons), qui sont donc la
verite de reference a laquelle l'extraction est comparee.

Toutes les dimensions sont en metres, l'axe de rotation est Z.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

from .io import repair as repair_module
from .mesh import TriMesh

Vec3 = tuple[float, float, float]


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------
def box(size_x: float, size_y: float, size_z: float, center: Sequence[float] = (0.0, 0.0, 0.0)) -> TriMesh:
    """Parallelepipede aligne sur les axes, centre sur `center`."""
    hx, hy, hz = size_x / 2.0, size_y / 2.0, size_z / 2.0
    cx, cy, cz = center
    vertices = [
        (cx - hx, cy - hy, cz - hz), (cx + hx, cy - hy, cz - hz),
        (cx + hx, cy + hy, cz - hz), (cx - hx, cy + hy, cz - hz),
        (cx - hx, cy - hy, cz + hz), (cx + hx, cy - hy, cz + hz),
        (cx + hx, cy + hy, cz + hz), (cx - hx, cy + hy, cz + hz),
    ]
    faces = [
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return TriMesh(vertices, faces)


def cube(edge: float, center: Sequence[float] = (0.0, 0.0, 0.0)) -> TriMesh:
    """Cube d'arete `edge`."""
    return box(edge, edge, edge, center)


def revolve(profile: Sequence[Sequence[float]], segments: int = 180) -> TriMesh:
    """Solide de revolution autour de Z, engendre par un contour ferme (r, z).

    `profile` est une polyligne fermee du demi-plan (r >= 0, z) ; le dernier
    point est relie au premier.  Les aretes situees sur l'axe (r = 0 aux deux
    extremites) ne produisent pas de facette.  L'orientation des normales est
    corrigee en fin de construction.
    """
    n = len(profile)
    if n < 3:
        raise ValueError("un profil de revolution demande au moins trois points")
    vertices: list[Vec3] = []
    rings: list[list[int]] = []
    for r, z in profile:
        if r <= 0.0:
            index = len(vertices)
            vertices.append((0.0, 0.0, float(z)))
            rings.append([index] * segments)
            continue
        ring: list[int] = []
        for s in range(segments):
            angle = 2.0 * math.pi * s / segments
            ring.append(len(vertices))
            vertices.append((r * math.cos(angle), r * math.sin(angle), float(z)))
        rings.append(ring)

    faces: list[tuple[int, int, int]] = []
    for i in range(n):
        j = (i + 1) % n
        ring_a, ring_b = rings[i], rings[j]
        degenerate_a = profile[i][0] <= 0.0
        degenerate_b = profile[j][0] <= 0.0
        if degenerate_a and degenerate_b:
            continue
        for s in range(segments):
            t = (s + 1) % segments
            a0, a1 = ring_a[s], ring_a[t]
            b0, b1 = ring_b[s], ring_b[t]
            if degenerate_a:
                faces.append((a0, b0, b1))
            elif degenerate_b:
                faces.append((a0, a1, b0))
            else:
                faces.append((a0, a1, b1))
                faces.append((a0, b1, b0))
    mesh = TriMesh(vertices, faces)
    mesh.merge_vertices()
    repair_module.fix_normals(mesh)
    return mesh


def cylinder(radius: float, height: float, segments: int = 180, z_center: float = 0.0) -> TriMesh:
    """Cylindre plein d'axe Z."""
    half = height / 2.0
    return revolve(
        [(0.0, z_center - half), (radius, z_center - half), (radius, z_center + half), (0.0, z_center + half)],
        segments,
    )


def tube(r_inner: float, r_outer: float, height: float, segments: int = 180, z_center: float = 0.0) -> TriMesh:
    """Couronne cylindrique (tube) d'axe Z."""
    half = height / 2.0
    return revolve(
        [
            (r_inner, z_center - half), (r_outer, z_center - half),
            (r_outer, z_center + half), (r_inner, z_center + half),
        ],
        segments,
    )


def _closed_box_from_grids(grid_a: list[list[Vec3]], grid_b: list[list[Vec3]]) -> TriMesh:
    """Surface fermee entre deux nappes de meme topologie (nu x nv).

    `grid_a` et `grid_b` sont les deux faces opposees ; les quatre bandes de bord
    referment le volume.  Le resultat est un solide topologiquement equivalent a
    une boite, donc etanche.
    """
    nu = len(grid_a)
    nv = len(grid_a[0])
    vertices: list[Vec3] = []
    index_a = [[0] * nv for _ in range(nu)]
    index_b = [[0] * nv for _ in range(nu)]
    for i in range(nu):
        for j in range(nv):
            index_a[i][j] = len(vertices)
            vertices.append(grid_a[i][j])
    for i in range(nu):
        for j in range(nv):
            index_b[i][j] = len(vertices)
            vertices.append(grid_b[i][j])

    faces: list[tuple[int, int, int]] = []

    def quad(p0: int, p1: int, p2: int, p3: int) -> None:
        faces.append((p0, p1, p2))
        faces.append((p0, p2, p3))

    for i in range(nu - 1):
        for j in range(nv - 1):
            quad(index_a[i][j], index_a[i + 1][j], index_a[i + 1][j + 1], index_a[i][j + 1])
            quad(index_b[i][j], index_b[i][j + 1], index_b[i + 1][j + 1], index_b[i + 1][j])
    for i in range(nu - 1):
        quad(index_a[i][0], index_b[i][0], index_b[i + 1][0], index_a[i + 1][0])
        quad(index_a[i][nv - 1], index_a[i + 1][nv - 1], index_b[i + 1][nv - 1], index_b[i][nv - 1])
    for j in range(nv - 1):
        quad(index_a[0][j], index_a[0][j + 1], index_b[0][j + 1], index_b[0][j])
        quad(index_a[nu - 1][j], index_b[nu - 1][j], index_b[nu - 1][j + 1], index_a[nu - 1][j + 1])

    mesh = TriMesh(vertices, faces)
    mesh.merge_vertices()
    repair_module.fix_normals(mesh)
    return mesh


def combine(meshes: Sequence[TriMesh]) -> TriMesh:
    """Assemble plusieurs solides en un maillage unique (union non booleenne).

    Les solides restent des composantes connexes distinctes, orientees vers
    l'exterieur : la carte d'occupation, qui travaille par nombre d'enroulement,
    les traite bien comme une union (voir `geometry.occupancy`).
    """
    vertices: list[Vec3] = []
    faces: list[tuple[int, int, int]] = []
    for mesh in meshes:
        base = len(vertices)
        vertices.extend(mesh.vertices)
        faces.extend((base + i, base + j, base + k) for i, j, k in mesh.faces)
    return TriMesh(vertices, faces)


# ---------------------------------------------------------------------------
# Roues
# ---------------------------------------------------------------------------
def axial_blade_surface(
    r: float,
    theta: float,
    beta_deg: float,
    sense: int,
    theta_0: float,
    z_center: float,
) -> float:
    """Cote z de la nappe moyenne d'une pale a angle de pale constant.

    `z = z_center + sense * r * tan(beta) * (theta - theta_0)` : par
    construction `tan(beta) = dz / (r dtheta)` a tout rayon, donc l'angle de pale
    vaut exactement `beta_deg` partout (SPEC 4.3).
    """
    return z_center + sense * r * math.tan(math.radians(beta_deg)) * (theta - theta_0)


def axial_impeller(
    n_blades: int = 4,
    beta_deg: float = 30.0,
    r_hub: float = 0.030,
    r_tip: float = 0.100,
    blade_wrap_deg: float = 60.0,
    thickness: float = 0.004,
    sense: int = 1,
    hub_margin: float = 0.010,
    n_radial: int = 14,
    n_chord: int = 18,
    hub_segments: int = 180,
) -> TriMesh:
    """Helice axiale a pales d'angle constant `beta_deg`, moyeu cylindrique.

    `sense = +1` donne `k = dz/dtheta > 0`, donc une rotation attendue
    anti-horaire vue de +Z ; `sense = -1` donne l'inverse (SPEC 4.4).
    """
    wrap = math.radians(blade_wrap_deg)
    half_t = thickness / 2.0
    z_extent = r_tip * math.tan(math.radians(beta_deg)) * wrap / 2.0
    hub_half = z_extent + half_t + hub_margin
    parts = [cylinder(r_hub, 2.0 * hub_half, hub_segments, 0.0)]

    root = r_hub * 0.9  # la pale penetre le moyeu pour garantir l'union
    for blade in range(n_blades):
        theta_0 = 2.0 * math.pi * blade / n_blades
        grid_low: list[list[Vec3]] = []
        grid_high: list[list[Vec3]] = []
        for i in range(n_radial):
            r = root + (r_tip - root) * i / (n_radial - 1)
            row_low: list[Vec3] = []
            row_high: list[Vec3] = []
            for j in range(n_chord):
                theta = theta_0 - wrap / 2.0 + wrap * j / (n_chord - 1)
                z_mid = axial_blade_surface(r, theta, beta_deg, sense, theta_0, 0.0)
                x, y = r * math.cos(theta), r * math.sin(theta)
                row_low.append((x, y, z_mid - half_t))
                row_high.append((x, y, z_mid + half_t))
            grid_low.append(row_low)
            grid_high.append(row_high)
        parts.append(_closed_box_from_grids(grid_low, grid_high))
    return combine(parts)


def _centrifugal_theta(
    r: float,
    r1: float,
    r2: float,
    beta1_deg: float,
    beta2_deg: float,
    meridional_slope: Callable[[float], float] | None = None,
    steps: int = 400,
) -> float:
    """Angle polaire de la ligne de cambrure d'une aube centrifuge.

    Integration de `dtheta = dm / (r tan(beta(r)))` avec `beta` variant
    lineairement de `beta1` a `beta2` et `dm = sqrt(1 + (dz/dr)^2) dr`
    l'abscisse curviligne meridienne.  L'angle de pale de la geometrie produite,
    au sens de la SPEC 4.3 mesure sur la surface de courant, vaut donc exactement
    `beta1` en `r1` et `beta2` en `r2` -- y compris quand la veine descend, ce
    qui est le cas des que la roue a un oeillard axial.
    """
    if r <= r1:
        return 0.0
    total = 0.0
    n = max(2, int(steps * (r - r1) / (r2 - r1)))
    for i in range(n):
        ra = r1 + (r - r1) * i / n
        rb = r1 + (r - r1) * (i + 1) / n
        rc = 0.5 * (ra + rb)
        s = (rc - r1) / (r2 - r1)
        beta = math.radians(beta1_deg + (beta2_deg - beta1_deg) * s)
        slope = meridional_slope(rc) if meridional_slope is not None else 0.0
        total += (rb - ra) * math.sqrt(1.0 + slope * slope) / (rc * math.tan(beta))
    return total


def centrifugal_impeller(
    n_blades: int = 6,
    beta1_deg: float = 22.0,
    beta2_deg: float = 25.0,
    r1: float = 0.035,
    r2: float = 0.090,
    b1: float = 0.020,
    b2: float = 0.010,
    eye_height: float = 0.035,
    thickness: float = 0.004,
    edge_taper: float = 0.0,
    sense: int = 1,
    shroud_thickness: float = 0.006,
    n_radial: int = 20,
    n_span: int = 5,
    hub_segments: int = 180,
) -> TriMesh:
    """Roue centrifuge a aubes en arc, avec oeillard axial et flasque arriere.

    La veine meridienne descend de `eye_height` (bord d'attaque, cote aspiration
    +Z) au plan de sortie z = 0, et sa largeur passe de `b1` a `b2` : le maximum
    de z de la zone de pales est donc bien atteint au rayon d'oeillard `r1`, et
    le rapport `r2 / r1s` classe la roue (SPEC 3.2 et 3.3).
    """

    def z_hi(r: float) -> float:
        s = (r - r1) / (r2 - r1)
        return eye_height * (1.0 - s)

    def width(r: float) -> float:
        s = (r - r1) / (r2 - r1)
        return b1 + (b2 - b1) * s

    def z_lo(r: float) -> float:
        return z_hi(r) - width(r)

    def mid_slope(r: float) -> float:
        """Pente dz/dr de la surface de courant a mi-envergure."""
        step = (r2 - r1) * 1e-4
        z_mid_a = z_hi(r - step) - 0.5 * width(r - step)
        z_mid_b = z_hi(r + step) - 0.5 * width(r + step)
        return (z_mid_b - z_mid_a) / (2.0 * step)

    # Moyeu / flasque arriere : solide de revolution sous la veine.
    bottom = z_lo(r2) - shroud_thickness
    profile: list[tuple[float, float]] = [(0.0, bottom)]
    for i in range(n_radial):
        r = r1 + (r2 - r1) * i / (n_radial - 1)
        profile.append((r, bottom))
    for i in range(n_radial - 1, -1, -1):
        r = r1 + (r2 - r1) * i / (n_radial - 1)
        profile.append((r, z_lo(r)))
    profile.append((0.0, z_lo(r1)))
    parts = [revolve(profile, hub_segments)]

    def half_angle_at(r: float) -> float:
        """Demi-epaisseur angulaire de l'aube au rayon r.

        `edge_taper` amincit l'aube pres du bord d'attaque et du bord de fuite
        par une loi elliptique, comme le veut le dessin usuel d'une aube :
        laisse a zero, l'aube est coupee carre aux deux bouts.
        """
        half = 0.5 * thickness / r
        if edge_taper <= 0.0:
            return half
        s = (r - r1) / (r2 - r1)
        edge = min(s, 1.0 - s) / edge_taper
        if edge >= 1.0:
            return half
        return half * math.sqrt(max(0.0, 1.0 - (1.0 - edge) ** 2))

    for blade in range(n_blades):
        base = 2.0 * math.pi * blade / n_blades
        grid_a: list[list[Vec3]] = []
        grid_b: list[list[Vec3]] = []
        for i in range(n_radial):
            r = r1 + (r2 - r1) * i / (n_radial - 1)
            theta_c = base + sense * _centrifugal_theta(r, r1, r2, beta1_deg, beta2_deg, mid_slope)
            half = half_angle_at(r)
            row_a: list[Vec3] = []
            row_b: list[Vec3] = []
            lo = z_lo(r) - shroud_thickness * 0.5  # penetration dans le flasque
            hi = z_hi(r)
            for j in range(n_span):
                z = lo + (hi - lo) * j / (n_span - 1)
                row_a.append((r * math.cos(theta_c - half), r * math.sin(theta_c - half), z))
                row_b.append((r * math.cos(theta_c + half), r * math.sin(theta_c + half), z))
            grid_a.append(row_a)
            grid_b.append(row_b)
        parts.append(_closed_box_from_grids(grid_a, grid_b))
    return combine(parts)
