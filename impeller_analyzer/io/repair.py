"""Reparation de maillage : orientation des normales et rebouchage des trous.

Equivalent en Python pur de `trimesh.repair.fix_normals` et `fill_holes` (SPEC
phase 1, point 3).  Quand `trimesh` est disponible, `loader` lui delegue ; sinon
ces fonctions prennent le relais et donnent le meme resultat sur les maillages
usuels.
"""

from __future__ import annotations

from .. import config
from ..mesh import TriMesh, cross, dot


def _face_components(mesh: TriMesh) -> list[list[int]]:
    """Partitionne les faces en composantes connexes par les aretes."""
    adjacency: dict[tuple[int, int], list[int]] = {}
    for index, (i, j, k) in enumerate(mesh.faces):
        for a, b in ((i, j), (j, k), (k, i)):
            key = (a, b) if a < b else (b, a)
            adjacency.setdefault(key, []).append(index)
    seen = [False] * len(mesh.faces)
    components: list[list[int]] = []
    for start in range(len(mesh.faces)):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component = []
        while stack:
            face = stack.pop()
            component.append(face)
            i, j, k = mesh.faces[face]
            for a, b in ((i, j), (j, k), (k, i)):
                key = (a, b) if a < b else (b, a)
                for neighbour in adjacency.get(key, ()):
                    if not seen[neighbour]:
                        seen[neighbour] = True
                        stack.append(neighbour)
        components.append(component)
    return components


def _component_volume(mesh: TriMesh, component: list[int]) -> float:
    """Volume signe d'une composante connexe de faces."""
    total = 0.0
    for face in component:
        a, b, c = mesh.triangle(face)
        total += dot(a, cross(b, c))
    return total / 6.0


def fix_normals(mesh: TriMesh) -> int:
    """Rend les orientations coherentes puis les tourne vers l'exterieur.

    Renvoie le nombre de faces retournees.  Chaque composante connexe est
    propagee par parcours en profondeur : une face voisine dont l'arete commune
    est parcourue dans le meme sens est retournee.  La composante entiere est
    ensuite retournee si son volume signe est negatif.
    """
    if not mesh.faces:
        return 0
    adjacency: dict[tuple[int, int], list[int]] = {}
    for index, (i, j, k) in enumerate(mesh.faces):
        for a, b in ((i, j), (j, k), (k, i)):
            key = (a, b) if a < b else (b, a)
            adjacency.setdefault(key, []).append(index)

    flipped = 0
    seen = [False] * len(mesh.faces)
    components: list[list[int]] = []
    for start in range(len(mesh.faces)):
        if seen[start]:
            continue
        seen[start] = True
        stack = [start]
        component = [start]
        while stack:
            face = stack.pop()
            i, j, k = mesh.faces[face]
            directed = ((i, j), (j, k), (k, i))
            for a, b in directed:
                key = (a, b) if a < b else (b, a)
                for neighbour in adjacency.get(key, ()):
                    if neighbour == face or seen[neighbour]:
                        continue
                    ni, nj, nk = mesh.faces[neighbour]
                    # Coherent si le voisin parcourt l'arete en sens oppose.
                    if (a, b) in ((ni, nj), (nj, nk), (nk, ni)):
                        mesh.faces[neighbour] = (ni, nk, nj)
                        flipped += 1
                    seen[neighbour] = True
                    component.append(neighbour)
                    stack.append(neighbour)
        components.append(component)

    mesh.invalidate()
    for component in components:
        if _component_volume(mesh, component) < 0.0:
            for face in component:
                i, j, k = mesh.faces[face]
                mesh.faces[face] = (i, k, j)
                flipped += 1
    mesh.invalidate()
    return flipped


def boundary_loops(mesh: TriMesh) -> list[list[int]]:
    """Boucles de sommets bordant les trous.

    Chaque boucle est orientee de sorte que les triangles de rebouchage
    `(v_i, v_{i+1}, centre)` soient coherents avec les faces existantes.
    """
    successor: dict[int, list[int]] = {}
    edge_faces = mesh.edge_map()
    for index, (i, j, k) in enumerate(mesh.faces):
        for a, b in ((i, j), (j, k), (k, i)):
            key = (a, b) if a < b else (b, a)
            if len(edge_faces.get(key, ())) == 1:
                # La face parcourt a -> b ; le bouchon doit parcourir b -> a.
                successor.setdefault(b, []).append(a)
    loops: list[list[int]] = []
    while successor:
        start = next(iter(successor))
        loop = [start]
        current = start
        while True:
            candidates = successor.get(current)
            if not candidates:
                loop = []
                break
            nxt = candidates.pop()
            if not candidates:
                successor.pop(current, None)
            if nxt == start:
                break
            if nxt in loop:
                loop = []
                break
            loop.append(nxt)
            current = nxt
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def fill_holes(mesh: TriMesh, max_edges: int = config.FILL_HOLES_MAX_EDGES) -> int:
    """Rebouche les trous d'au plus `max_edges` aretes ; renvoie le nombre rebouche.

    Un trou triangulaire est ferme par un seul triangle, les autres par un
    eventail centre sur un nouveau sommet place au barycentre de la boucle : le
    procede reste valide pour une boucle non plane.
    """
    filled = 0
    for loop in boundary_loops(mesh):
        if len(loop) > max_edges:
            continue
        if len(loop) == 3:
            mesh.faces.append((loop[0], loop[1], loop[2]))
            filled += 1
            continue
        n = float(len(loop))
        centre = (
            sum(mesh.vertices[v][0] for v in loop) / n,
            sum(mesh.vertices[v][1] for v in loop) / n,
            sum(mesh.vertices[v][2] for v in loop) / n,
        )
        centre_index = len(mesh.vertices)
        mesh.vertices.append(centre)
        for position in range(len(loop)):
            a = loop[position]
            b = loop[(position + 1) % len(loop)]
            mesh.faces.append((a, b, centre_index))
        filled += 1
    if filled:
        mesh.invalidate()
        mesh.remove_degenerate_faces()
    return filled


def repair(mesh: TriMesh, max_passes: int = config.REPAIR_MAX_PASSES) -> dict[str, int]:
    """Sequence complete : fusion, orientation, rebouchage, re-orientation.

    Renvoie un petit journal des operations, repris dans le rapport d'import.
    """
    log = {"merged_vertices": 0, "removed_faces": 0, "flipped_faces": 0, "filled_holes": 0}
    log["merged_vertices"] += mesh.merge_vertices(config.MERGE_TOL)
    log["removed_faces"] += mesh.remove_degenerate_faces()
    for _ in range(max_passes):
        log["flipped_faces"] += fix_normals(mesh)
        if mesh.is_watertight():
            break
        added = fill_holes(mesh)
        log["filled_holes"] += added
        if added == 0:
            break
    return log
