"""Boite a outils numerique en Python pur (bibliotheque standard seule).

Regroupe les briques dont le paquet a besoin et qui viendraient normalement de
`numpy`/`scipy` : diagonalisation d'une matrice symetrique, transformee de
Fourier discrete, spline cubique de lissage.  Chaque fonction est ecrite pour
les tailles de problemes du projet (matrices 3x3, signaux de 720 points,
splines de 25 noeuds), pas pour le calcul intensif.
"""

from __future__ import annotations

import cmath
import math
from typing import Sequence

from . import config


# ---------------------------------------------------------------------------
# Algebre lineaire
# ---------------------------------------------------------------------------
def jacobi_eigen(matrix: Sequence[Sequence[float]]) -> tuple[list[float], list[list[float]]]:
    """Diagonalise une matrice symetrique 3x3 par rotations de Jacobi.

    Renvoie `(valeurs_propres, vecteurs_propres)` triees par valeur propre
    croissante ; `vecteurs_propres[i]` est le vecteur associe a
    `valeurs_propres[i]`, unitaire.
    """
    n = len(matrix)
    a = [[float(matrix[i][j]) for j in range(n)] for i in range(n)]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    scale = max((abs(a[i][j]) for i in range(n) for j in range(n)), default=0.0)
    if scale == 0.0:
        return [0.0] * n, [list(row) for row in v]

    for _ in range(config.JACOBI_MAX_SWEEPS):
        off = math.sqrt(sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j))
        if off <= config.JACOBI_TOL * scale:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) <= config.JACOBI_TOL * scale:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp = v[k][p]
                    vkq = v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq

    values = [a[i][i] for i in range(n)]
    vectors = [[v[i][j] for i in range(n)] for j in range(n)]  # colonnes -> lignes
    order = sorted(range(n), key=lambda i: values[i])
    return [values[i] for i in order], [vectors[i] for i in order]


def solve_banded_symmetric(diagonals: list[list[float]], rhs: list[float]) -> list[float]:
    """Resout un systeme symetrique defini positif a bande, par Cholesky bande.

    `diagonals[k][i]` est le terme `A[i][i + k]`, pour `k` de 0 a la largeur de
    bande.  Utilise par la spline de lissage (bande de largeur 2).
    """
    bandwidth = len(diagonals) - 1
    n = len(rhs)
    lower = [[0.0] * n for _ in range(bandwidth + 1)]
    for i in range(n):
        total = diagonals[0][i]
        for k in range(1, bandwidth + 1):
            j = i - k
            if j >= 0:
                total -= lower[k][j] ** 2
        if total <= 0.0:
            total = abs(total) + config.JACOBI_TOL
        lower[0][i] = math.sqrt(total)
        for k in range(1, bandwidth + 1):
            row = i + k
            if row >= n:
                break
            total = diagonals[k][i]
            for m in range(1, bandwidth + 1 - k):
                left = i - m
                if left >= 0:
                    total -= lower[m][left] * lower[m + k][left]
            lower[k][i] = total / lower[0][i]

    y = [0.0] * n
    for i in range(n):
        total = rhs[i]
        for k in range(1, bandwidth + 1):
            j = i - k
            if j >= 0:
                total -= lower[k][j] * y[j]
        y[i] = total / lower[0][i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        total = y[i]
        for k in range(1, bandwidth + 1):
            j = i + k
            if j < n:
                total -= lower[k][i] * x[j]
        x[i] = total / lower[0][i]
    return x


# ---------------------------------------------------------------------------
# Analyse spectrale
# ---------------------------------------------------------------------------
def dft_magnitudes(signal: Sequence[float], harmonics: int) -> list[float]:
    """Amplitudes des harmoniques 0..`harmonics` d'un signal periodique reel.

    Transformee de Fourier discrete evaluee uniquement sur les harmoniques
    utiles (au plus `BLADES_MAX`), ce qui evite d'implementer une FFT complete
    pour un gain nul a cette taille.  L'amplitude est normalisee par la longueur
    du signal.
    """
    n = len(signal)
    if n == 0:
        return [0.0] * (harmonics + 1)
    out: list[float] = []
    for k in range(harmonics + 1):
        acc = 0j
        factor = -2.0j * math.pi * k / n
        for index, value in enumerate(signal):
            acc += value * cmath.exp(factor * index)
        out.append(abs(acc) / n)
    return out


# ---------------------------------------------------------------------------
# Spline cubique de lissage
# ---------------------------------------------------------------------------
class CubicSpline:
    """Spline cubique naturelle definie par ses valeurs et derivees secondes."""

    __slots__ = ("x", "y", "sigma")

    def __init__(self, x: Sequence[float], y: Sequence[float], sigma: Sequence[float]):
        self.x = list(x)
        self.y = list(y)
        self.sigma = list(sigma)

    def _segment(self, t: float) -> int:
        """Index du segment contenant `t` (extrapolation par les segments extremes)."""
        x = self.x
        if t <= x[0]:
            return 0
        if t >= x[-1]:
            return len(x) - 2
        low, high = 0, len(x) - 1
        while high - low > 1:
            mid = (low + high) // 2
            if x[mid] <= t:
                low = mid
            else:
                high = mid
        return low

    def __call__(self, t: float) -> float:
        """Valeur de la spline en `t`."""
        i = self._segment(t)
        h = self.x[i + 1] - self.x[i]
        d = t - self.x[i]
        s0, s1 = self.sigma[i], self.sigma[i + 1]
        slope = (self.y[i + 1] - self.y[i]) / h - h * (2.0 * s0 + s1) / 6.0
        return self.y[i] + d * slope + 0.5 * s0 * d * d + (s1 - s0) * d ** 3 / (6.0 * h)

    def derivative(self, t: float) -> float:
        """Derivee premiere de la spline en `t`."""
        i = self._segment(t)
        h = self.x[i + 1] - self.x[i]
        d = t - self.x[i]
        s0, s1 = self.sigma[i], self.sigma[i + 1]
        slope = (self.y[i + 1] - self.y[i]) / h - h * (2.0 * s0 + s1) / 6.0
        return slope + s0 * d + (s1 - s0) * d * d / (2.0 * h)


def smoothing_spline(x: Sequence[float], y: Sequence[float], residual_target: float) -> CubicSpline:
    """Spline cubique naturelle de lissage (formulation de Green & Silverman).

    Minimise `somme((y - g)^2) + lambda * integrale(g''^2)`, `lambda` etant
    ajuste par bissection pour que la somme des carres des residus atteigne
    `residual_target`.  `residual_target = 0` redonne l'interpolation exacte.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x et y doivent avoir la meme longueur")
    if n < 4:
        return _interpolating_spline(x, y)

    h = [x[i + 1] - x[i] for i in range(n - 1)]
    if min(h) <= 0.0:
        raise ValueError("les abscisses doivent etre strictement croissantes")
    m = n - 2  # noeuds interieurs

    # R : tridiagonale (m x m) ; Q^T Q : pentadiagonale (m x m).
    r_diag = [(h[j] + h[j + 1]) / 3.0 for j in range(m)]
    r_off = [h[j + 1] / 6.0 for j in range(m - 1)]

    def q_column(j: int) -> tuple[float, float, float]:
        """Colonne j de Q : coefficients sur les lignes j, j+1, j+2."""
        return (1.0 / h[j], -1.0 / h[j] - 1.0 / h[j + 1], 1.0 / h[j + 1])

    columns = [q_column(j) for j in range(m)]
    qtq0 = [sum(c * c for c in columns[j]) for j in range(m)]
    qtq1 = [columns[j][1] * columns[j + 1][0] + columns[j][2] * columns[j + 1][1] for j in range(m - 1)]
    qtq2 = [columns[j][2] * columns[j + 2][0] for j in range(m - 2)]
    qty = [sum(columns[j][k] * y[j + k] for k in range(3)) for j in range(m)]

    def fit(lam: float) -> tuple[list[float], list[float], float]:
        """Ajuste pour un lambda donne ; renvoie (valeurs, derivees secondes, residu)."""
        d0 = [r_diag[j] + lam * qtq0[j] for j in range(m)]
        d1 = [r_off[j] + lam * qtq1[j] for j in range(m - 1)] + [0.0]
        d2 = [lam * qtq2[j] for j in range(m - 2)] + [0.0, 0.0]
        sigma_inner = solve_banded_symmetric([d0, d1, d2], qty)
        g = list(y)
        for j in range(m):
            c = columns[j]
            for k in range(3):
                g[j + k] -= lam * c[k] * sigma_inner[j]
        sigma = [0.0] + sigma_inner + [0.0]
        residual = sum((y[i] - g[i]) ** 2 for i in range(n))
        return g, sigma, residual

    if residual_target <= 0.0:
        return _interpolating_spline(x, y)

    low, high = config.SPLINE_LAMBDA_MIN, config.SPLINE_LAMBDA_MAX
    _, _, residual_high = fit(high)
    if residual_high <= residual_target:
        g, sigma, _ = fit(high)
        return CubicSpline(x, g, sigma)
    best = fit(low)
    if best[2] >= residual_target:
        return CubicSpline(x, best[0], best[1])
    for _ in range(config.SPLINE_BISECTION_STEPS):
        mid = math.sqrt(low * high)
        g, sigma, residual = fit(mid)
        if residual > residual_target:
            high = mid
        else:
            low = mid
            best = (g, sigma, residual)
    return CubicSpline(x, best[0], best[1])


def _interpolating_spline(x: Sequence[float], y: Sequence[float]) -> CubicSpline:
    """Spline cubique naturelle passant exactement par les points."""
    n = len(x)
    if n < 3:
        return CubicSpline(x, y, [0.0] * n)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    m = n - 2
    d0 = [(h[j] + h[j + 1]) / 3.0 for j in range(m)]
    d1 = [h[j + 1] / 6.0 for j in range(m - 1)] + [0.0]
    rhs = [(y[j + 2] - y[j + 1]) / h[j + 1] - (y[j + 1] - y[j]) / h[j] for j in range(m)]
    sigma_inner = solve_banded_symmetric([d0, d1], rhs)
    return CubicSpline(x, y, [0.0] + sigma_inner + [0.0])


def linear_interp(x: Sequence[float], y: Sequence[float], t: float) -> float:
    """Interpolation lineaire de la table (x croissant, y), avec maintien aux bords."""
    if t <= x[0]:
        return y[0]
    if t >= x[-1]:
        return y[-1]
    low, high = 0, len(x) - 1
    while high - low > 1:
        mid = (low + high) // 2
        if x[mid] <= t:
            low = mid
        else:
            high = mid
    span = x[low + 1] - x[low]
    if span == 0.0:
        return y[low]
    ratio = (t - x[low]) / span
    return y[low] + ratio * (y[low + 1] - y[low])
