"""Point dans un solide : parite des traversees d'un rayon (equivalent de `trimesh.contains`).

Sert aux controles du mode composants qui portent sur des **volumes** et non
sur des boites : deux copies de pale se recoupent-elles, une pale entre-t-elle
dans le moyeu. Une boite englobante ne repond a aucune des deux questions --
celle d'une pale en boucle fait presque le tour de l'axe, celle d'un corps de
roue fermee contient toutes les pales par construction.

Le rayon part le long de +X. Les triangles sont ranges dans une grille du plan
(y, z) : seuls ceux dont l'ombre contient le point sont examines. L'origine du
rayon est decalee d'une fraction irrationnelle de cellule, pour qu'il ne passe
jamais exactement par une arete ou un sommet, ce qui compterait une traversee
deux fois ou pas du tout.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

from .. import config
from ..mesh import TriMesh


class SolidTester:
    """Teste l'appartenance de points a un solide ferme."""

    def __init__(self, mesh: TriMesh, cells: int = config.PROXIMITY_CELLS):
        self.mesh = mesh
        lo, hi = mesh.bounds()
        self.lo, self.hi = lo, hi
        span = max(hi[1] - lo[1], hi[2] - lo[2]) or 1.0
        self.cell = span / max(1, cells)
        # Decalage du rayon : petit devant la cellule, irrationnel pour eviter
        # les aretes et les sommets du maillage.
        self.jitter = (self.cell * 1e-4 * math.sqrt(2.0), self.cell * 1e-4 * math.sqrt(3.0))
        self.buckets: dict[tuple[int, int], list[int]] = {}
        for index in range(len(mesh.faces)):
            a, b, c = mesh.triangle(index)
            y0, y1 = min(a[1], b[1], c[1]), max(a[1], b[1], c[1])
            z0, z1 = min(a[2], b[2], c[2]), max(a[2], b[2], c[2])
            for iy in range(self._index(y0, 1), self._index(y1, 1) + 1):
                for iz in range(self._index(z0, 2), self._index(z1, 2) + 1):
                    self.buckets.setdefault((iy, iz), []).append(index)

    def _index(self, value: float, axis: int) -> int:
        return int(math.floor((value - self.lo[axis]) / self.cell))

    def contains(self, point: Sequence[float]) -> bool:
        """Vrai si le point est a l'interieur du solide."""
        x = point[0]
        y = point[1] + self.jitter[0]
        z = point[2] + self.jitter[1]
        if not (self.lo[0] <= x <= self.hi[0] and self.lo[1] <= y <= self.hi[1]
                and self.lo[2] <= z <= self.hi[2]):
            return False
        crossings = 0
        for index in self.buckets.get((self._index(y, 1), self._index(z, 2)), ()):
            a, b, c = self.mesh.triangle(index)
            # Intersection du rayon (x, y, z) + t (1, 0, 0), t > 0, par
            # coordonnees barycentriques dans le plan (y, z).
            d = (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])
            if d == 0.0:
                continue
            u = ((y - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (z - a[2])) / d
            v = ((b[1] - a[1]) * (z - a[2]) - (y - a[1]) * (b[2] - a[2])) / d
            if u < 0.0 or v < 0.0 or u + v > 1.0:
                continue
            hit = a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0])
            if hit > x:
                crossings += 1
        return crossings % 2 == 1

    def fraction_inside(self, points: Iterable[Sequence[float]]) -> float:
        """Part des points a l'interieur du solide."""
        total = inside = 0
        for point in points:
            total += 1
            if self.contains(point):
                inside += 1
        return inside / total if total else 0.0
