"""Carte meridienne d'occupation angulaire `f(r, z)` (SPEC 2.2).

C'est la structure centrale du systeme : toutes les phases suivantes lisent
cette carte plutot que le maillage.

Principe.  Pour chaque hauteur `z` de la grille, le maillage est coupe par le
plan correspondant, ce qui donne un jeu de segments 2D.  Pour chacune des
`N_THETA` directions angulaires, un rayon partant de l'axe traverse ces segments
et le **nombre d'enroulement** (winding number) est accumule depuis l'exterieur
(ou il vaut 0) vers l'axe : la matiere est la ou il vaut au moins 1.  Compter le
nombre de directions occupees donne directement `f(r, z)`.

Ce comptage signe -- plutot qu'une simple parite -- est ce qui permet de traiter
correctement une roue livree comme une **union de solides qui s'interpenetrent**
(moyeu + pales), cas le plus frequent des exports CAO.

Cout : `GRID_NZ` coupes de plan, puis un balayage par segment et par secteur
angulaire.  Le remplissage de la grille passe par un tableau de differences, ce
qui evite de marcher les `GRID_NR` cellules pour chaque rayon.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field

from .. import config
from ..mesh import TriMesh

Interval = tuple[int, int]  # (index de cellule radiale de debut, fin exclue)


@dataclass
class OccupancyMap:
    """Grille meridienne `f(r, z)` et les intervalles radiaux occupes par secteur."""

    r_edges: list[float] = field(default_factory=list)
    z_edges: list[float] = field(default_factory=list)
    r_centres: list[float] = field(default_factory=list)
    z_centres: list[float] = field(default_factory=list)
    f: list[list[float]] = field(default_factory=list)  # f[iz][ir]
    intervals: list[list[tuple[Interval, ...]]] = field(default_factory=list)  # [iz][itheta]
    n_theta: int = config.N_THETA
    r_max: float = 0.0
    z_min: float = 0.0
    z_max: float = 0.0

    # -- acces ------------------------------------------------------------
    @property
    def nr(self) -> int:
        """Nombre de cellules radiales."""
        return len(self.r_centres)

    @property
    def nz(self) -> int:
        """Nombre de cellules axiales."""
        return len(self.z_centres)

    @property
    def dr(self) -> float:
        """Pas radial de la grille."""
        return self.r_max / self.nr if self.nr else 0.0

    @property
    def dz(self) -> float:
        """Pas axial de la grille."""
        return (self.z_max - self.z_min) / self.nz if self.nz else 0.0

    def radial_index(self, r: float) -> int:
        """Index de la cellule radiale contenant `r` (borne aux extremites)."""
        if self.dr <= 0.0:
            return 0
        return max(0, min(self.nr - 1, int(r / self.dr)))

    def axial_index(self, z: float) -> int:
        """Index de la cellule axiale contenant `z` (borne aux extremites)."""
        if self.dz <= 0.0:
            return 0
        return max(0, min(self.nz - 1, int((z - self.z_min) / self.dz)))

    def value(self, r: float, z: float) -> float:
        """Fraction angulaire occupee au point (r, z)."""
        if r > self.r_max or z < self.z_min or z > self.z_max:
            return 0.0
        return self.f[self.axial_index(z)][self.radial_index(r)]

    # -- classification ---------------------------------------------------
    def solid_mask(self) -> list[list[bool]]:
        """Cellules pleines : `f >= F_SOLIDE` (moyeu ou flasque)."""
        return [[value >= config.F_SOLIDE for value in row] for row in self.f]

    def blade_mask(self) -> list[list[bool]]:
        """Cellules de la zone de pales : `F_VIDE < f < F_SOLIDE`."""
        return [
            [config.F_VIDE < value < config.F_SOLIDE for value in row]
            for row in self.f
        ]

    def void_mask(self) -> list[list[bool]]:
        """Cellules de veine fluide ou exterieures : `f <= F_VIDE`."""
        return [[value <= config.F_VIDE for value in row] for row in self.f]

    def theta_signal(self, mask: list[list[bool]] | None = None) -> list[float]:
        """Signal `g(theta)` : occupation integree sur les cellules de `mask`.

        La somme est ponderee par le volume `r dr dz` de chaque cellule, de sorte
        qu'un secteur soit compte a la mesure de la matiere qu'il porte.  Sans
        `mask`, la zone de pales est utilisee (SPEC 3.1).
        """
        if mask is None:
            mask = self.blade_mask()
        signal = [0.0] * self.n_theta
        dr, dz = self.dr, self.dz
        for iz in range(self.nz):
            row = mask[iz]
            # Somme prefixe des poids des cellules retenues a cette hauteur.
            prefix = [0.0] * (self.nr + 1)
            for ir in range(self.nr):
                weight = self.r_centres[ir] * dr * dz if row[ir] else 0.0
                prefix[ir + 1] = prefix[ir] + weight
            if prefix[self.nr] <= 0.0:
                continue
            sectors = self.intervals[iz]
            for itheta in range(self.n_theta):
                total = 0.0
                for start, stop in sectors[itheta]:
                    total += prefix[stop] - prefix[start]
                signal[itheta] += total
        return signal

    def theta_spectrum(
        self,
        harmonics: int,
        mask: list[list[bool]] | None = None,
    ) -> list[float]:
        """Spectre angulaire moyen, module somme **cellule par cellule**.

        Le spectre de `theta_signal` -- qui integre d'abord sur toute la zone de
        pales puis transforme -- s'effondre des que les pales se recouvrent en
        projection : une roue a N pales dont l'enroulement approche `2*pi/N`
        remplit tous les secteurs, le signal integre devient presque constant et
        l'harmonique N disparait devant celle en 2N.  Le probleme vient de la
        phase : a rayon et hauteur fixes la pale occupe une bande angulaire
        etroite, mais sa position tourne avec la hauteur, et les contributions
        s'annulent a la sommation.

        On transforme donc **avant** de sommer, et on additionne les modules.
        Le coefficient de Fourier de chaque cellule est obtenu sans boucle sur
        les cellules : les secteurs occupes sont deja connus sous forme
        d'intervalles radiaux, dont un tableau de differences complexe donne
        directement le coefficient a tout rayon.
        """
        if mask is None:
            mask = self.blade_mask()
        orders = list(range(config.BLADES_MIN, harmonics + 1))
        if not orders or self.nz == 0:
            return [0.0] * (harmonics + 1)
        tables = {
            k: [cmath.exp(-2.0j * math.pi * k * i / self.n_theta) for i in range(self.n_theta)]
            for k in orders
        }
        amplitudes = [0.0] * (harmonics + 1)
        weight_total = 0.0
        stride = max(1, self.nz // config.SPECTRUM_MAX_ROWS)
        for iz in range(0, self.nz, stride):
            row = mask[iz]
            if not any(row):
                continue
            sectors = self.intervals[iz]
            diffs = {k: [0j] * (self.nr + 1) for k in orders}
            for itheta in range(self.n_theta):
                spans = sectors[itheta]
                if not spans:
                    continue
                for k in orders:
                    table = tables[k][itheta]
                    diff = diffs[k]
                    for start, stop in spans:
                        diff[start] += table
                        diff[stop] -= table
            running = {k: 0j for k in orders}
            for ir in range(self.nr):
                for k in orders:
                    running[k] += diffs[k][ir]
                if not row[ir]:
                    continue
                weight = self.r_centres[ir]
                weight_total += weight
                for k in orders:
                    amplitudes[k] += weight * abs(running[k]) / self.n_theta
        if weight_total > 0.0:
            amplitudes = [value / weight_total for value in amplitudes]
        return amplitudes

    def summary(self) -> dict:
        """Resume serialisable de la carte (les tableaux complets restent internes)."""
        blade = self.blade_mask()
        solid = self.solid_mask()
        return {
            "grille": {"nr": self.nr, "nz": self.nz, "n_theta": self.n_theta},
            "r_max_m": self.r_max,
            "z_min_m": self.z_min,
            "z_max_m": self.z_max,
            "cellules_pleines": sum(sum(1 for v in row if v) for row in solid),
            "cellules_pales": sum(sum(1 for v in row if v) for row in blade),
            "f_max": max((max(row) for row in self.f), default=0.0),
        }


def _merge_intervals(intervals: list[Interval]) -> tuple[Interval, ...]:
    """Fusionne des intervalles de cellules radiales qui se touchent."""
    if not intervals:
        return ()
    intervals.sort()
    merged = [list(intervals[0])]
    for start, stop in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], stop)
        else:
            merged.append([start, stop])
    return tuple((a, b) for a, b in merged)


def build_occupancy(
    mesh: TriMesh,
    nr: int = config.GRID_NR,
    nz: int = config.GRID_NZ,
    n_theta: int = config.N_THETA,
) -> OccupancyMap:
    """Construit la carte `f(r, z)` du maillage (suppose deja aligne sur Z).

    Le maillage doit avoir ses normales orientees vers l'exterieur : le signe du
    produit scalaire normale/rayon est ce qui distingue une entree d'une sortie
    de matiere.
    """
    lo, hi = mesh.bounds()
    z_min, z_max = lo[2], hi[2]
    r_max = max((math.hypot(v[0], v[1]) for v in mesh.vertices), default=0.0)
    if r_max <= config.OCCUPANCY_MIN_CELL or (z_max - z_min) <= config.OCCUPANCY_MIN_CELL:
        raise ValueError("maillage degenere : extension radiale ou axiale nulle")

    dr = r_max / nr
    dz = (z_max - z_min) / nz
    r_centres = [(i + 0.5) * dr for i in range(nr)]
    z_centres = [z_min + (k + 0.5) * dz for k in range(nz)]
    r_edges = [i * dr for i in range(nr + 1)]
    z_edges = [z_min + k * dz for k in range(nz + 1)]

    cos_theta = [math.cos(2.0 * math.pi * k / n_theta) for k in range(n_theta)]
    sin_theta = [math.sin(2.0 * math.pi * k / n_theta) for k in range(n_theta)]
    bins_per_turn = n_theta / (2.0 * math.pi)

    # Repartition des triangles par tranche : chaque triangle n'est examine que
    # pour les plans de coupe qu'il traverse reellement.
    buckets: list[list[int]] = [[] for _ in range(nz)]
    vertices = mesh.vertices
    for index, (i, j, k) in enumerate(mesh.faces):
        za, zb, zc = vertices[i][2], vertices[j][2], vertices[k][2]
        lo_z = za if za < zb else zb
        if zc < lo_z:
            lo_z = zc
        hi_z = za if za > zb else zb
        if zc > hi_z:
            hi_z = zc
        # La plage est elargie d'une tranche de chaque cote : un triangle dont
        # une arete affleure exactement le plan de coupe est traite comme
        # traversant par la regle de l'epsilon ci-dessous, alors que le calcul
        # d'indice le rejetterait a 1e-16 pres. Le test de signe ecarte ensuite
        # les triangles reellement hors plan, pour un cout negligeable.
        first = int(math.ceil((lo_z - z_min) / dz - 0.5)) - 1
        last = int(math.floor((hi_z - z_min) / dz - 0.5)) + 1
        if first < 0:
            first = 0
        if last > nz - 1:
            last = nz - 1
        for level in range(first, last + 1):
            buckets[level].append(index)

    eps = config.SLICE_EPS_REL * max(mesh.scale_length(), config.OCCUPANCY_MIN_CELL)
    grid: list[list[float]] = []
    all_intervals: list[list[tuple[Interval, ...]]] = []

    for level in range(nz):
        z_c = z_centres[level]
        crossings: list[list[tuple[float, int]]] = [[] for _ in range(n_theta)]

        for face in buckets[level]:
            i, j, k = mesh.faces[face]
            pa, pb, pc = vertices[i], vertices[j], vertices[k]
            da = pa[2] - z_c
            db = pb[2] - z_c
            dc = pc[2] - z_c
            if -eps < da < eps:
                da = eps
            if -eps < db < eps:
                db = eps
            if -eps < dc < eps:
                dc = eps
            if (da > 0.0) == (db > 0.0) == (dc > 0.0):
                continue

            # Normale non normalisee du triangle : seule sa projection compte.
            ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
            vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            if nx == 0.0 and ny == 0.0:
                continue

            points: list[tuple[float, float]] = []
            for (p0, d0), (p1, d1) in (((pa, da), (pb, db)), ((pb, db), (pc, dc)), ((pc, dc), (pa, da))):
                if (d0 > 0.0) != (d1 > 0.0):
                    t = d0 / (d0 - d1)
                    points.append((p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])))
            if len(points) != 2:
                continue
            (px, py), (qx, qy) = points
            ex, ey = qx - px, qy - py
            if ex == 0.0 and ey == 0.0:
                continue

            angle_p = math.atan2(py, px)
            delta = math.atan2(qy, qx) - angle_p
            if delta > math.pi:
                delta -= 2.0 * math.pi
            elif delta < -math.pi:
                delta += 2.0 * math.pi
            low_angle, high_angle = (angle_p, angle_p + delta) if delta >= 0.0 else (angle_p + delta, angle_p)

            # L'intervalle angulaire ne sert que de presélection, elargie d'un
            # secteur de chaque cote ; l'appartenance reelle est tranchee par le
            # signe du produit vectoriel, ce qui reste exact au passage de la
            # coupure +/-pi de l'arc tangente.
            first_bin = int(math.floor(low_angle * bins_per_turn)) - 1
            last_bin = int(math.ceil(high_angle * bins_per_turn)) + 1
            cross_pe = px * ey - py * ex  # produit vectoriel 2D (p, e)
            for raw in range(first_bin, last_bin + 1):
                bin_index = raw % n_theta
                cx, sy = cos_theta[bin_index], sin_theta[bin_index]
                side_p = cx * py - sy * px
                side_q = cx * qy - sy * qx
                # Convention semi-ouverte : un sommet exactement sur le rayon
                # n'appartient qu'a un seul des deux segments qui le portent,
                # donc il n'est jamais compte deux fois ni oublie.
                if not ((side_p <= 0.0 < side_q) or (side_q <= 0.0 < side_p)):
                    continue
                denominator = cx * ey - sy * ex
                if denominator == 0.0:
                    continue
                radius = cross_pe / denominator
                if radius <= 0.0:
                    continue
                crossings[bin_index].append((radius, 1 if (nx * cx + ny * sy) > 0.0 else -1))

        counts = [0] * (nr + 1)
        sectors: list[tuple[Interval, ...]] = []
        for bin_index in range(n_theta):
            hits = crossings[bin_index]
            if not hits:
                sectors.append(())
                continue
            hits.sort()
            winding = 0
            outer = None
            spans: list[Interval] = []
            for position in range(len(hits) - 1, -1, -1):
                radius, sign = hits[position]
                if winding >= 1 and outer is not None:
                    span = _cells_between(radius, outer, dr, nr)
                    if span is not None:
                        spans.append(span)
                winding += sign
                outer = radius
            if winding >= 1 and outer is not None:
                span = _cells_between(0.0, outer, dr, nr)
                if span is not None:
                    spans.append(span)
            merged = _merge_intervals(spans)
            sectors.append(merged)
            for start, stop in merged:
                counts[start] += 1
                counts[stop] -= 1
        all_intervals.append(sectors)

        row = [0.0] * nr
        running = 0
        for ir in range(nr):
            running += counts[ir]
            row[ir] = running / n_theta
        grid.append(row)

    return OccupancyMap(
        r_edges=r_edges,
        z_edges=z_edges,
        r_centres=r_centres,
        z_centres=z_centres,
        f=grid,
        intervals=all_intervals,
        n_theta=n_theta,
        r_max=r_max,
        z_min=z_min,
        z_max=z_max,
    )


def _cells_between(r_low: float, r_high: float, dr: float, nr: int) -> Interval | None:
    """Cellules radiales dont le centre tombe dans `[r_low, r_high]`."""
    first = int(math.ceil(r_low / dr - 0.5))
    last = int(math.floor(r_high / dr - 0.5))
    if first < 0:
        first = 0
    if last > nr - 1:
        last = nr - 1
    if first > last:
        return None
    return (first, last + 1)


def is_inside(mesh: TriMesh, point: tuple[float, float, float]) -> bool:
    """Test d'inclusion d'un point par nombre d'enroulement (rayon +X).

    Utilise pour verifier la carte d'occupation ; ce n'est pas le chemin de
    calcul de production, qui traite un plan entier a la fois.
    """
    x, y, z = point
    winding = 0
    for a, b, c in mesh.triangles():
        # Le rayon part de `point` vers +X : on teste l'appartenance au triangle
        # projete sur le plan (y, z), puis le signe de l'abscisse d'intersection.
        d0 = a[1] - y, a[2] - z
        d1 = b[1] - y, b[2] - z
        d2 = c[1] - y, c[2] - z
        s0 = d0[0] * d1[1] - d0[1] * d1[0]
        s1 = d1[0] * d2[1] - d1[1] * d2[0]
        s2 = d2[0] * d0[1] - d2[1] * d0[0]
        if not ((s0 >= 0.0 and s1 >= 0.0 and s2 >= 0.0) or (s0 <= 0.0 and s1 <= 0.0 and s2 <= 0.0)):
            continue
        total = s0 + s1 + s2
        if total == 0.0:
            continue
        hit_x = (s1 * a[0] + s2 * b[0] + s0 * c[0]) / total
        if hit_x <= x:
            continue
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx = uy * vz - uz * vy
        # Normale sortante dirigee vers +X : on sort de la matiere, donc le point
        # etait dedans pour cette traversee.
        winding += 1 if nx > 0.0 else -1
    return winding >= 1
