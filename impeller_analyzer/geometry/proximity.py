"""Distance point-maillage et distance de Hausdorff (equivalent de `trimesh.proximity`).

Sert au controle croise de la phase 3 : apres rotation de `2*pi/N`, un maillage
reellement periodique d'ordre N doit se superposer a lui-meme.  Le hachage
spatial des triangles rend la requete quasi constante en nombre de triangles.
"""

from __future__ import annotations

import math
from typing import Sequence

from .. import config
from ..mesh import TriMesh


def closest_point_on_triangle(
    point: Sequence[float],
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
) -> tuple[float, float, float]:
    """Point du triangle (a, b, c) le plus proche de `point` (regions de Voronoi)."""
    abx, aby, abz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    acx, acy, acz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    apx, apy, apz = point[0] - a[0], point[1] - a[1], point[2] - a[2]

    d1 = abx * apx + aby * apy + abz * apz
    d2 = acx * apx + acy * apy + acz * apz
    if d1 <= 0.0 and d2 <= 0.0:
        return (a[0], a[1], a[2])

    bpx, bpy, bpz = point[0] - b[0], point[1] - b[1], point[2] - b[2]
    d3 = abx * bpx + aby * bpy + abz * bpz
    d4 = acx * bpx + acy * bpy + acz * bpz
    if d3 >= 0.0 and d4 <= d3:
        return (b[0], b[1], b[2])

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3) if d1 != d3 else 0.0
        return (a[0] + v * abx, a[1] + v * aby, a[2] + v * abz)

    cpx, cpy, cpz = point[0] - c[0], point[1] - c[1], point[2] - c[2]
    d5 = abx * cpx + aby * cpy + abz * cpz
    d6 = acx * cpx + acy * cpy + acz * cpz
    if d6 >= 0.0 and d5 <= d6:
        return (c[0], c[1], c[2])

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6) if d2 != d6 else 0.0
        return (a[0] + w * acx, a[1] + w * acy, a[2] + w * acz)

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        denominator = (d4 - d3) + (d5 - d6)
        w = (d4 - d3) / denominator if denominator != 0.0 else 0.0
        return (
            b[0] + w * (c[0] - b[0]),
            b[1] + w * (c[1] - b[1]),
            b[2] + w * (c[2] - b[2]),
        )

    denominator = va + vb + vc
    if denominator == 0.0:
        return (a[0], a[1], a[2])
    v = vb / denominator
    w = vc / denominator
    return (a[0] + abx * v + acx * w, a[1] + aby * v + acy * w, a[2] + abz * v + acz * w)


class TriangleGrid:
    """Hachage spatial uniforme des triangles d'un maillage."""

    def __init__(self, mesh: TriMesh, cells: int = config.PROXIMITY_CELLS):
        self.mesh = mesh
        lo, hi = mesh.bounds()
        self.origin = lo
        largest = max(hi[k] - lo[k] for k in range(3)) or 1.0
        self.cell = largest / max(1, cells)
        self.buckets: dict[tuple[int, int, int], list[int]] = {}
        for index in range(len(mesh.faces)):
            a, b, c = mesh.triangle(index)
            low = [min(a[k], b[k], c[k]) for k in range(3)]
            high = [max(a[k], b[k], c[k]) for k in range(3)]
            ranges = [
                range(self._axis_index(low[k], k), self._axis_index(high[k], k) + 1)
                for k in range(3)
            ]
            for ix in ranges[0]:
                for iy in ranges[1]:
                    for iz in ranges[2]:
                        self.buckets.setdefault((ix, iy, iz), []).append(index)

    def _axis_index(self, value: float, axis: int) -> int:
        """Index de cellule le long d'un axe."""
        return int(math.floor((value - self.origin[axis]) / self.cell))

    def distance(self, point: Sequence[float]) -> float:
        """Distance exacte du point a la surface du maillage."""
        base = tuple(self._axis_index(point[k], k) for k in range(3))
        best = math.inf
        ring = 0
        max_ring = config.PROXIMITY_CELLS * 2 + 2
        while ring <= max_ring:
            candidates: set[int] = set()
            for dx in range(-ring, ring + 1):
                for dy in range(-ring, ring + 1):
                    for dz in range(-ring, ring + 1):
                        # Seule la couche externe du cube est nouvelle a ce tour.
                        if ring and max(abs(dx), abs(dy), abs(dz)) != ring:
                            continue
                        bucket = self.buckets.get((base[0] + dx, base[1] + dy, base[2] + dz))
                        if bucket:
                            candidates.update(bucket)
            for index in candidates:
                a, b, c = self.mesh.triangle(index)
                q = closest_point_on_triangle(point, a, b, c)
                d = math.dist(point, q)
                if d < best:
                    best = d
            if best < ring * self.cell:
                break
            ring += 1
        return best


def hausdorff_distance(
    mesh_a: TriMesh,
    mesh_b: TriMesh,
    samples: int = config.HAUSDORFF_SAMPLES,
) -> float:
    """Distance de Hausdorff symetrique approchee entre deux maillages.

    Les points sont echantillonnes sur les surfaces (ponderation par l'aire,
    suite deterministe), puis leur distance a l'autre surface est calculee
    exactement.  L'erreur est donc celle de l'echantillonnage des points de
    depart, pas de la mesure de distance.
    """
    grid_a = TriangleGrid(mesh_a)
    grid_b = TriangleGrid(mesh_b)
    worst = 0.0
    for point in mesh_a.sample_surface(samples):
        worst = max(worst, grid_b.distance(point))
    for point in mesh_b.sample_surface(samples):
        worst = max(worst, grid_a.distance(point))
    return worst
