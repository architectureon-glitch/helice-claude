"""Maillage triangulaire unifie et operations geometriques de base.

C'est la structure d'echange interne du paquet : tous les lecteurs de la phase 1
produisent un `TriMesh`, toutes les phases suivantes le consomment.  Les
coordonnees sont **toujours en metres** : la conversion d'unite est faite une
seule fois, a l'import (SPEC 1).

L'implementation est en Python pur (bibliotheque standard uniquement) de facon a
rester executable sans dependance externe ; `numpy`/`trimesh` restent utilises
quand ils sont installes, mais ne sont jamais requis.
"""

from __future__ import annotations

import math
from typing import Iterable, Iterator, Sequence

from . import config

Vec3 = tuple[float, float, float]
Face = tuple[int, int, int]
Mat3 = tuple[Vec3, Vec3, Vec3]


def cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    """Produit vectoriel a x b."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    """Produit scalaire a . b."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    """Difference a - b."""
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: Sequence[float], b: Sequence[float]) -> Vec3:
    """Somme a + b."""
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a: Sequence[float], k: float) -> Vec3:
    """Produit du vecteur a par le scalaire k."""
    return (a[0] * k, a[1] * k, a[2] * k)


def norm(a: Sequence[float]) -> float:
    """Norme euclidienne de a."""
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def normalize(a: Sequence[float]) -> Vec3:
    """Vecteur unitaire colineaire a a (renvoie a inchange si a est nul)."""
    n = norm(a)
    if n == 0.0:
        return (a[0], a[1], a[2])
    return (a[0] / n, a[1] / n, a[2] / n)


def matvec(m: Mat3, v: Sequence[float]) -> Vec3:
    """Produit matrice 3x3 par vecteur."""
    return (dot(m[0], v), dot(m[1], v), dot(m[2], v))


def matmul(a: Mat3, b: Mat3) -> Mat3:
    """Produit de deux matrices 3x3."""
    return tuple(  # type: ignore[return-value]
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3)
    )


def transpose(m: Mat3) -> Mat3:
    """Transposee d'une matrice 3x3."""
    return tuple(tuple(m[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def rotation_matrix(axis: Sequence[float], angle: float) -> Mat3:
    """Matrice de rotation d'angle `angle` (rad) autour de `axis` (Rodrigues)."""
    x, y, z = normalize(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    t = 1.0 - c
    return (
        (t * x * x + c, t * x * y - s * z, t * x * z + s * y),
        (t * x * y + s * z, t * y * y + c, t * y * z - s * x),
        (t * x * z - s * y, t * y * z + s * x, t * z * z + c),
    )


def rotation_between(a: Sequence[float], b: Sequence[float]) -> Mat3:
    """Rotation minimale amenant le vecteur unitaire a sur le vecteur unitaire b."""
    u = normalize(a)
    v = normalize(b)
    axis = cross(u, v)
    sin_a = norm(axis)
    cos_a = dot(u, v)
    if sin_a < config.JACOBI_TOL:
        if cos_a > 0.0:
            return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        # Vecteurs opposes : demi-tour autour de n'importe quelle perpendiculaire.
        perp = cross(u, (1.0, 0.0, 0.0))
        if norm(perp) < config.JACOBI_TOL:
            perp = cross(u, (0.0, 1.0, 0.0))
        return rotation_matrix(perp, math.pi)
    return rotation_matrix(axis, math.atan2(sin_a, cos_a))


class TriMesh:
    """Maillage triangulaire indexe, en metres.

    Attributs
    ---------
    vertices : liste de sommets (x, y, z), en metres.
    faces    : liste de triplets d'indices, orientation directe vers l'exterieur
               apres reparation.
    """

    __slots__ = ("vertices", "faces", "_cache")

    def __init__(self, vertices: Iterable[Sequence[float]], faces: Iterable[Sequence[int]]):
        self.vertices: list[Vec3] = [(float(v[0]), float(v[1]), float(v[2])) for v in vertices]
        self.faces: list[Face] = [(int(f[0]), int(f[1]), int(f[2])) for f in faces]
        self._cache: dict[str, object] = {}

    # -- construction ------------------------------------------------------
    def copy(self) -> "TriMesh":
        """Copie profonde du maillage (le cache n'est pas recopie)."""
        return TriMesh(list(self.vertices), list(self.faces))

    def __repr__(self) -> str:  # pragma: no cover - confort d'inspection
        return f"TriMesh({len(self.vertices)} sommets, {len(self.faces)} triangles)"

    def invalidate(self) -> None:
        """Vide le cache des grandeurs derivees, apres modification en place."""
        self._cache.clear()

    # -- geometrie elementaire --------------------------------------------
    def triangle(self, index: int) -> tuple[Vec3, Vec3, Vec3]:
        """Les trois sommets du triangle `index`."""
        i, j, k = self.faces[index]
        return self.vertices[i], self.vertices[j], self.vertices[k]

    def triangles(self) -> Iterator[tuple[Vec3, Vec3, Vec3]]:
        """Iterateur sur les triangles, sous forme de triplets de sommets."""
        vertices = self.vertices
        for i, j, k in self.faces:
            yield vertices[i], vertices[j], vertices[k]

    def bounds(self) -> tuple[Vec3, Vec3]:
        """Boite englobante alignee sur les axes : (min, max)."""
        cached = self._cache.get("bounds")
        if cached is not None:
            return cached  # type: ignore[return-value]
        if not self.vertices:
            result = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        else:
            xs = [v[0] for v in self.vertices]
            ys = [v[1] for v in self.vertices]
            zs = [v[2] for v in self.vertices]
            result = ((min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs)))
        self._cache["bounds"] = result
        return result

    def extents(self) -> Vec3:
        """Dimensions de la boite englobante."""
        lo, hi = self.bounds()
        return (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])

    def scale_length(self) -> float:
        """Longueur caracteristique du maillage (diagonale de la bbox)."""
        return norm(self.extents())

    def area(self) -> float:
        """Aire totale de la surface."""
        cached = self._cache.get("area")
        if cached is not None:
            return cached  # type: ignore[return-value]
        total = 0.0
        for a, b, c in self.triangles():
            total += 0.5 * norm(cross(sub(b, a), sub(c, a)))
        self._cache["area"] = total
        return total

    def volume(self) -> float:
        """Volume signe, positif si les normales pointent vers l'exterieur.

        Theoreme de la divergence applique triangle par triangle ; exact pour un
        maillage ferme, approche sinon.
        """
        cached = self._cache.get("volume")
        if cached is not None:
            return cached  # type: ignore[return-value]
        total = 0.0
        for a, b, c in self.triangles():
            total += dot(a, cross(b, c))
        total /= 6.0
        self._cache["volume"] = total
        return total

    def face_normal(self, index: int) -> Vec3:
        """Normale unitaire du triangle `index`."""
        a, b, c = self.triangle(index)
        return normalize(cross(sub(b, a), sub(c, a)))

    def face_areas(self) -> list[float]:
        """Aire de chaque triangle."""
        return [0.5 * norm(cross(sub(b, a), sub(c, a))) for a, b, c in self.triangles()]

    def centroid(self) -> Vec3:
        """Barycentre du volume ferme (barycentre des sommets si volume nul)."""
        cached = self._cache.get("centroid")
        if cached is not None:
            return cached  # type: ignore[return-value]
        vol = 0.0
        acc = [0.0, 0.0, 0.0]
        for a, b, c in self.triangles():
            d = dot(a, cross(b, c)) / 6.0
            vol += d
            for axis in range(3):
                acc[axis] += d * 0.25 * (a[axis] + b[axis] + c[axis])
        if abs(vol) < config.OCCUPANCY_MIN_CELL and self.vertices:
            n = float(len(self.vertices))
            result = tuple(sum(v[axis] for v in self.vertices) / n for axis in range(3))
        elif not self.vertices:
            result = (0.0, 0.0, 0.0)
        else:
            result = (acc[0] / vol, acc[1] / vol, acc[2] / vol)
        self._cache["centroid"] = result
        return result  # type: ignore[return-value]

    def covariance(self, about: Sequence[float] | None = None) -> Mat3:
        """Matrice des moments d'ordre deux `integrale(x x^T dV)`.

        `about` est l'origine de reference (barycentre du volume par defaut).
        Decomposition en tetraedres (origine, a, b, c) ; les integrales sur le
        simplexe de reference valent 1/60 (termes carres) et 1/120 (termes
        croises).
        """
        origin = tuple(self.centroid()) if about is None else tuple(float(x) for x in about)
        acc = [[0.0] * 3 for _ in range(3)]
        for a0, b0, c0 in self.triangles():
            p1 = sub(a0, origin)
            p2 = sub(b0, origin)
            p3 = sub(c0, origin)
            det = dot(p1, cross(p2, p3))
            if det == 0.0:
                continue
            pts = (p1, p2, p3)
            for i in range(3):
                for j in range(3):
                    square = sum(pts[k][i] * pts[k][j] for k in range(3)) / 60.0
                    crossed = 0.0
                    for k in range(3):
                        for m in range(k + 1, 3):
                            crossed += pts[k][i] * pts[m][j] + pts[m][i] * pts[k][j]
                    acc[i][j] += det * (square + crossed / 120.0)
        return tuple(tuple(row) for row in acc)  # type: ignore[return-value]

    def inertia_tensor(self, about: Sequence[float] | None = None, density: float = 1.0) -> Mat3:
        """Tenseur d'inertie d'un solide homogene : `trace(C) I - C`."""
        cov = self.covariance(about)
        trace = cov[0][0] + cov[1][1] + cov[2][2]
        return tuple(  # type: ignore[return-value]
            tuple(density * ((trace - cov[i][i]) if i == j else -cov[i][j]) for j in range(3))
            for i in range(3)
        )

    # -- topologie ---------------------------------------------------------
    def edge_map(self) -> dict[tuple[int, int], list[int]]:
        """Table arete non orientee -> indices des faces qui la portent."""
        cached = self._cache.get("edge_map")
        if cached is not None:
            return cached  # type: ignore[return-value]
        edges: dict[tuple[int, int], list[int]] = {}
        for index, (i, j, k) in enumerate(self.faces):
            for a, b in ((i, j), (j, k), (k, i)):
                key = (a, b) if a < b else (b, a)
                edges.setdefault(key, []).append(index)
        self._cache["edge_map"] = edges
        return edges

    def boundary_edges(self) -> list[tuple[int, int]]:
        """Aretes portees par une seule face (bord libre du maillage)."""
        return [edge for edge, faces in self.edge_map().items() if len(faces) == 1]

    def is_watertight(self) -> bool:
        """Vrai si chaque arete est partagee par exactement deux faces."""
        if not self.faces:
            return False
        return all(len(faces) == 2 for faces in self.edge_map().values())

    def is_winding_consistent(self) -> bool:
        """Vrai si toute arete interne est parcourue en sens oppose par ses deux faces."""
        seen: dict[tuple[int, int], int] = {}
        for i, j, k in self.faces:
            for a, b in ((i, j), (j, k), (k, i)):
                seen[(a, b)] = seen.get((a, b), 0) + 1
        for (a, b), count in seen.items():
            if count > 1:
                return False
            if seen.get((b, a), 0) > 1:
                return False
        return True

    # -- transformations ---------------------------------------------------
    def apply_scale(self, factor: float) -> None:
        """Homothetie en place (conversion d'unite)."""
        self.vertices = [(v[0] * factor, v[1] * factor, v[2] * factor) for v in self.vertices]
        self.invalidate()

    def apply_translation(self, offset: Sequence[float]) -> None:
        """Translation en place."""
        dx, dy, dz = float(offset[0]), float(offset[1]), float(offset[2])
        self.vertices = [(v[0] + dx, v[1] + dy, v[2] + dz) for v in self.vertices]
        self.invalidate()

    def apply_rotation(self, matrix: Mat3) -> None:
        """Rotation en place par une matrice 3x3."""
        self.vertices = [matvec(matrix, v) for v in self.vertices]
        self.invalidate()

    def transformed(
        self,
        matrix: Mat3 | None = None,
        translation: Sequence[float] | None = None,
        factor: float | None = None,
    ) -> "TriMesh":
        """Copie transformee : homothetie, puis rotation, puis translation."""
        out = self.copy()
        if factor is not None:
            out.apply_scale(factor)
        if matrix is not None:
            out.apply_rotation(matrix)
        if translation is not None:
            out.apply_translation(translation)
        return out

    def flip(self) -> None:
        """Inverse l'orientation de toutes les faces."""
        self.faces = [(f[0], f[2], f[1]) for f in self.faces]
        self.invalidate()

    # -- nettoyage ---------------------------------------------------------
    def remove_degenerate_faces(self) -> int:
        """Supprime les triangles d'aire nulle ; renvoie le nombre supprime."""
        kept: list[Face] = []
        removed = 0
        for i, j, k in self.faces:
            if i == j or j == k or i == k:
                removed += 1
                continue
            a, b, c = self.vertices[i], self.vertices[j], self.vertices[k]
            if norm(cross(sub(b, a), sub(c, a))) <= 0.0:
                removed += 1
                continue
            kept.append((i, j, k))
        if removed:
            self.faces = kept
            self.invalidate()
        return removed

    def merge_vertices(self, tol: float = config.MERGE_TOL) -> int:
        """Fusionne les sommets distants de moins de `tol` ; renvoie le nombre fusionne.

        Hachage sur une grille de pas `tol` : deux sommets sont fusionnes s'ils
        tombent dans la meme cellule ou dans une cellule voisine et que leur
        distance reelle est inferieure a `tol`.
        """
        if tol <= 0.0 or not self.vertices:
            return 0
        buckets: dict[tuple[int, int, int], list[int]] = {}
        new_vertices: list[Vec3] = []
        remap = [0] * len(self.vertices)
        inv = 1.0 / tol
        neighbourhood = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
        tol2 = tol * tol
        for index, v in enumerate(self.vertices):
            key = (int(math.floor(v[0] * inv)), int(math.floor(v[1] * inv)), int(math.floor(v[2] * inv)))
            found = -1
            for dx, dy, dz in neighbourhood:
                bucket = buckets.get((key[0] + dx, key[1] + dy, key[2] + dz))
                if not bucket:
                    continue
                for candidate in bucket:
                    w = new_vertices[candidate]
                    dxv = w[0] - v[0]
                    dyv = w[1] - v[1]
                    dzv = w[2] - v[2]
                    if dxv * dxv + dyv * dyv + dzv * dzv <= tol2:
                        found = candidate
                        break
                if found >= 0:
                    break
            if found < 0:
                found = len(new_vertices)
                new_vertices.append(v)
                buckets.setdefault(key, []).append(found)
            remap[index] = found
        merged = len(self.vertices) - len(new_vertices)
        self.vertices = new_vertices
        self.faces = [(remap[i], remap[j], remap[k]) for i, j, k in self.faces]
        self.invalidate()
        self.remove_degenerate_faces()
        self.remove_duplicate_faces()
        return merged

    def remove_duplicate_faces(self) -> int:
        """Supprime les triangles repetes (meme triplet de sommets)."""
        seen: set[frozenset[int]] = set()
        kept: list[Face] = []
        removed = 0
        for face in self.faces:
            key = frozenset(face)
            if key in seen:
                removed += 1
                continue
            seen.add(key)
            kept.append(face)
        if removed:
            self.faces = kept
            self.invalidate()
        return removed

    def remove_unreferenced_vertices(self) -> int:
        """Supprime les sommets qu'aucune face n'utilise."""
        used = sorted({index for face in self.faces for index in face})
        if len(used) == len(self.vertices):
            return 0
        removed = len(self.vertices) - len(used)
        remap = {old: new for new, old in enumerate(used)}
        self.vertices = [self.vertices[old] for old in used]
        self.faces = [(remap[i], remap[j], remap[k]) for i, j, k in self.faces]
        self.invalidate()
        return removed

    # -- echantillonnage ---------------------------------------------------
    def sample_surface(self, count: int) -> list[Vec3]:
        """Echantillonne `count` points sur la surface, ponderes par l'aire.

        L'echantillonnage est deterministe (suite de van der Corput en base 2 et
        3) pour que deux appels sur le meme maillage donnent le meme resultat :
        la comparaison de Hausdorff de la phase 3 en depend.
        """
        areas = self.face_areas()
        total = sum(areas)
        if count <= 0 or total <= 0.0:
            return []
        cumulative: list[float] = []
        running = 0.0
        for a in areas:
            running += a
            cumulative.append(running)
        points: list[Vec3] = []
        for index in range(count):
            target = _van_der_corput(index + 1, 2) * total
            face = _bisect_left(cumulative, target)
            if face >= len(self.faces):
                face = len(self.faces) - 1
            u = _van_der_corput(index + 1, 3)
            v = _van_der_corput(index + 1, 5)
            if u + v > 1.0:
                u, v = 1.0 - u, 1.0 - v
            a, b, c = self.triangle(face)
            points.append(
                (
                    a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]),
                    a[1] + u * (b[1] - a[1]) + v * (c[1] - a[1]),
                    a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]),
                )
            )
        return points


def _van_der_corput(index: int, base: int) -> float:
    """Terme `index` de la suite de van der Corput en base `base`, dans [0, 1)."""
    result = 0.0
    denom = 1.0
    while index > 0:
        denom *= base
        index, rest = divmod(index, base)
        result += rest / denom
    return result


def _bisect_left(values: Sequence[float], target: float) -> int:
    """Index d'insertion de `target` dans la liste croissante `values`."""
    low, high = 0, len(values)
    while low < high:
        mid = (low + high) // 2
        if values[mid] < target:
            low = mid + 1
        else:
            high = mid
    return low
