"""Le plan meridien de la roue : ou sont les parois, et par ou l'eau peut sortir.

Les pieces de corps -- moyeu, flasques, disque intermediaire, coque -- sont des
solides de revolution. Leur trace dans un demi-plan `(r, z)` dit ou l'eau peut
passer. Deux questions s'y posent, que ni une boite ni une demi-droite ne
tranchent :

- le bord de fuite d'un brin d'aube debouche-t-il vers la sortie, ou sur une
  paroi ? Une demi-droite radiale tiree du bord de fuite heurtait le cone du
  fond de la roue d'essai hel1, alors que l'eau longe ce cone jusqu'a la
  fente. On cherche donc un **chemin** : l'eau peut monter ou descendre le long
  d'une paroi, mais ne revient jamais vers l'axe ;
- le bout d'un profil d'aube est-il un bord libre, ou un pied soude a une paroi ?

Une case de la grille est paroi si la matiere l'occupe a la majorite de trois
azimuts : une empreinte d'aube ou un percage local n'ouvre ni ne ferme un
passage.
"""

from __future__ import annotations

import math
from typing import Sequence

from .. import config
from .inclusion import SolidTester

Vec3 = tuple[float, float, float]

#: Azimuts de sondage, quelconques a dessein : ni sur la couture y = 0 des
#: solides de revolution d'AutoCAD, ni a un multiple simple de 2 pi / N.
_AZIMUTHS = (0.4103, 2.4991, 4.5687)


class MeridianWalls:
    """Parois de la roue dans le plan meridien, sondees a la demande."""

    def __init__(self, testers: Sequence[SolidTester], axis: Vec3, r_max: float,
                 z_range: tuple[float, float]):
        self.testers = list(testers)
        self.axis = axis
        self.cell = max(r_max, 1e-9) / config.MERIDIAN_CELLS
        self.z_low, self.z_high = z_range
        self._memo: dict[tuple[int, int], bool] = {}

    def _r(self, i: int) -> float:
        return (i + 0.5) * self.cell

    def _z(self, j: int) -> float:
        return (j + 0.5) * self.cell

    def _i(self, r: float) -> int:
        return int(math.floor(r / self.cell))

    def _j(self, z: float) -> int:
        return int(math.floor(z / self.cell))

    def _wall(self, i: int, j: int) -> bool:
        key = (i, j)
        hit = self._memo.get(key)
        if hit is None:
            r, z = self._r(i), self._z(j)
            votes = 0
            for angle in _AZIMUTHS:
                point = (self.axis[0] + r * math.cos(angle), self.axis[1] + r * math.sin(angle), z)
                if any(t.contains(point) for t in self.testers):
                    votes += 1
            hit = votes >= 2
            self._memo[key] = hit
        return hit

    def wall(self, r: float, z: float) -> bool:
        """Vrai si la paroi occupe le point `(r, z)` du plan meridien."""
        return self._wall(self._i(r), self._j(z))

    def edge_is_free(self, r: float, z: float, inward: bool) -> bool:
        """Le bout d'un profil en `(r, z)` est-il libre, ou appuye sur une paroi ?

        On sonde juste au-dela du bout, vers l'axe (`inward`) ou vers
        l'exterieur : une aube soudee penetre sa paroi, et la case voisine y est.
        """
        step = config.MERIDIAN_EDGE_CELLS * self.cell
        return not self.wall(r - step if inward else r + step, z)

    def reaches(self, r: float, z: float, target_r: float, z_band: tuple[float, float]) -> bool:
        """Un chemin mene-t-il de `(r, z)` a la bande `r = target_r`, `z` dans `z_band` ?

        Le chemin ne revient jamais vers l'axe : a chaque rayon, l'eau parcourt
        la colonne d'eau ou elle se trouve, entre deux parois, puis passe au
        rayon suivant par les cases libres. La colonne est bornee a l'etendue
        axiale des parois, de la bande et du point de depart : au-dela, rien ne
        guide l'eau vers la sortie.
        """
        j_low = min(self._j(self.z_low), self._j(z_band[0]), self._j(z)) - 1
        j_high = max(self._j(self.z_high), self._j(z_band[1]), self._j(z)) + 1
        i, i_end = self._i(r), self._i(target_r)
        start = self._j(z)
        # Le bout du profil peut tomber dans une case de paroi, l'aube y etant
        # soudee : on part de la premiere case libre au-dessus ou au-dessous.
        frontier = {j for j in (start, start - 1, start + 1) if not self._wall(i, j)}
        while frontier:
            column: set[int] = set()
            for j in frontier:
                if j in column:
                    continue
                k = j
                while k >= j_low and not self._wall(i, k):
                    column.add(k)
                    k -= 1
                k = j + 1
                while k <= j_high and not self._wall(i, k):
                    column.add(k)
                    k += 1
            if i >= i_end:
                return any(z_band[0] <= self._z(k) <= z_band[1] for k in column)
            i += 1
            frontier = {k for k in column if not self._wall(i, k)}
        return False


def free_fraction(walls: Sequence[SolidTester], points: Sequence[Vec3]) -> float:
    """Part des points qu'aucune paroi n'occupe."""
    if not points:
        return 1.0
    free = sum(1 for p in points if not any(t.contains(p) for t in walls))
    return free / len(points)
