"""Ecriture de maillages, utilisee par les generateurs synthetiques et les tests.

Les fichiers sont ecrits dans l'unite d'import du projet : les coordonnees d'un
`TriMesh` etant en metres, `unit_factor` (m par unite de fichier) permet
d'exporter en centimetres, la convention d'entree de l'outil.
"""

from __future__ import annotations

import struct

from .. import config
from ..mesh import TriMesh, cross, normalize, sub


def _scaled(mesh: TriMesh, unit_factor: float) -> list[tuple[float, float, float]]:
    """Sommets convertis des metres vers l'unite de fichier demandee."""
    inv = 1.0 / unit_factor
    return [(v[0] * inv, v[1] * inv, v[2] * inv) for v in mesh.vertices]


def write_stl(mesh: TriMesh, path: str, unit_factor: float = config.UNIT_FACTOR, binary: bool = True) -> str:
    """Ecrit un STL binaire (par defaut) ou ASCII."""
    vertices = _scaled(mesh, unit_factor)
    if binary:
        with open(path, "wb") as handle:
            handle.write(b"impeller-analyzer".ljust(80, b"\0"))
            handle.write(struct.pack("<I", len(mesh.faces)))
            for i, j, k in mesh.faces:
                a, b, c = vertices[i], vertices[j], vertices[k]
                n = normalize(cross(sub(b, a), sub(c, a)))
                handle.write(struct.pack("<12fH", n[0], n[1], n[2], *a, *b, *c, 0))
        return path
    with open(path, "w", encoding="ascii") as handle:
        handle.write("solid impeller\n")
        for i, j, k in mesh.faces:
            a, b, c = vertices[i], vertices[j], vertices[k]
            n = normalize(cross(sub(b, a), sub(c, a)))
            handle.write(f"facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop\n")
            for p in (a, b, c):
                handle.write(f"    vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}\n")
            handle.write("  endloop\nendfacet\n")
        handle.write("endsolid impeller\n")
    return path


def write_obj(mesh: TriMesh, path: str, unit_factor: float = config.UNIT_FACTOR) -> str:
    """Ecrit un OBJ (indices en base 1)."""
    vertices = _scaled(mesh, unit_factor)
    with open(path, "w", encoding="ascii") as handle:
        handle.write("# impeller-analyzer\n")
        for v in vertices:
            handle.write(f"v {v[0]:.9g} {v[1]:.9g} {v[2]:.9g}\n")
        for i, j, k in mesh.faces:
            handle.write(f"f {i + 1} {j + 1} {k + 1}\n")
    return path


def write_off(mesh: TriMesh, path: str, unit_factor: float = config.UNIT_FACTOR) -> str:
    """Ecrit un OFF ASCII."""
    vertices = _scaled(mesh, unit_factor)
    with open(path, "w", encoding="ascii") as handle:
        handle.write(f"OFF\n{len(vertices)} {len(mesh.faces)} 0\n")
        for v in vertices:
            handle.write(f"{v[0]:.9g} {v[1]:.9g} {v[2]:.9g}\n")
        for face in mesh.faces:
            handle.write(f"3 {face[0]} {face[1]} {face[2]}\n")
    return path


def write_ply(mesh: TriMesh, path: str, unit_factor: float = config.UNIT_FACTOR, binary: bool = False) -> str:
    """Ecrit un PLY ASCII ou binaire little-endian."""
    vertices = _scaled(mesh, unit_factor)
    header = (
        "ply\n"
        f"format {'binary_little_endian' if binary else 'ascii'} 1.0\n"
        f"element vertex {len(vertices)}\n"
        "property float x\nproperty float y\nproperty float z\n"
        f"element face {len(mesh.faces)}\n"
        "property list uchar int vertex_indices\n"
        "end_header\n"
    )
    if binary:
        with open(path, "wb") as handle:
            handle.write(header.encode("ascii"))
            for v in vertices:
                handle.write(struct.pack("<3f", *v))
            for face in mesh.faces:
                handle.write(struct.pack("<B3i", 3, *face))
        return path
    with open(path, "w", encoding="ascii") as handle:
        handle.write(header)
        for v in vertices:
            handle.write(f"{v[0]:.9g} {v[1]:.9g} {v[2]:.9g}\n")
        for face in mesh.faces:
            handle.write(f"3 {face[0]} {face[1]} {face[2]}\n")
    return path


def write_dxf(mesh: TriMesh, path: str, unit_factor: float = config.UNIT_FACTOR) -> str:
    """Ecrit un DXF ASCII contenant un 3DFACE par triangle.

    Forme la plus portable d'un maillage en DXF : elle se relit sans ezdxf.
    """
    vertices = _scaled(mesh, unit_factor)
    out: list[str] = ["0", "SECTION", "2", "ENTITIES"]
    for i, j, k in mesh.faces:
        a, b, c = vertices[i], vertices[j], vertices[k]
        out.extend(["0", "3DFACE", "8", "0"])
        for corner, point in enumerate((a, b, c, c)):
            out.extend([str(10 + corner), f"{point[0]:.9g}"])
            out.extend([str(20 + corner), f"{point[1]:.9g}"])
            out.extend([str(30 + corner), f"{point[2]:.9g}"])
    out.extend(["0", "ENDSEC", "0", "EOF"])
    with open(path, "w", encoding="ascii") as handle:
        handle.write("\n".join(out) + "\n")
    return path


WRITERS = {
    ".stl": write_stl,
    ".obj": write_obj,
    ".off": write_off,
    ".ply": write_ply,
    ".dxf": write_dxf,
}
