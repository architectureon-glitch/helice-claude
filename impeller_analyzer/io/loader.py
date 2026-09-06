"""Lecture multi-format d'un maillage 3D vers un `TriMesh` unifie (SPEC phase 1).

Ordre de priorite des formats :

1. `.stl`, `.obj`, `.ply`, `.off` (et `.3ds` via `trimesh`) ;
2. `.step`, `.stp`, `.iges`, `.igs` par tessellation `cadquery`/OCC ;
3. `.dxf` (POLYFACE / 3DSOLID tesselle) par `ezdxf` ;
4. `.dwg` converti en DXF par l'ODA File Converter.

Les formats du groupe 1 sont lus par un parseur interne (bibliotheque standard
seule) ; `trimesh` est utilise en priorite quand il est installe, ce qui couvre
aussi `.3ds`.  Les groupes 2 a 4 dependent de leur bibliotheque respective et
levent une erreur explicite si elle manque.

Traitements systematiques a l'import : conversion d'unite, fusion des sommets
dupliques, reparation, rapport d'import.
"""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import tempfile
from dataclasses import dataclass, field, asdict

from .. import config
from ..confidence import HIGH, MEDIUM, ConfidenceMap
from ..mesh import TriMesh
from . import repair as repair_module


class ImportError_(Exception):
    """Erreur d'import lisible par l'utilisateur final."""


@dataclass
class ImportReport:
    """Rapport d'import (SPEC phase 1, point 4)."""

    path: str = ""
    extension: str = ""
    backend: str = ""  # "interne", "trimesh", "cadquery", "ezdxf", "oda+ezdxf"
    unit: str = ""
    unit_factor: float = config.UNIT_FACTOR
    n_vertices_raw: int = 0
    n_faces_raw: int = 0
    n_vertices: int = 0
    n_faces: int = 0
    merged_vertices: int = 0
    removed_faces: int = 0
    flipped_faces: int = 0
    filled_holes: int = 0
    watertight: bool = False
    winding_consistent: bool = False
    volume_m3: float = 0.0
    area_m2: float = 0.0
    bbox_min_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bbox_max_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    extents_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    warnings: list[str] = field(default_factory=list)
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        data = asdict(self)
        data["confidence"] = dict(self.confidence)
        return data


# ---------------------------------------------------------------------------
# Unites
# ---------------------------------------------------------------------------
def unit_factor(unit: str | float | None) -> tuple[str, float]:
    """Facteur de conversion vers le metre.

    Accepte un nom d'unite (`cm`, `mm`, `m`, `in`, ...) ou un facteur numerique
    explicite.  Par defaut : centimetre (SPEC 1.1).
    """
    if unit is None:
        return "cm", config.UNIT_FACTOR
    if isinstance(unit, (int, float)):
        value = float(unit)
        if value <= 0.0:
            raise ImportError_("le facteur d'unite doit etre strictement positif")
        return f"x{value:g}", value
    name = str(unit).strip().lower()
    if name in config.UNIT_FACTORS:
        return name, config.UNIT_FACTORS[name]
    try:
        value = float(name)
    except ValueError:
        known = ", ".join(sorted(config.UNIT_FACTORS))
        raise ImportError_(f"unite inconnue : {unit!r} (attendu : {known}, ou un facteur numerique)") from None
    if value <= 0.0:
        raise ImportError_("le facteur d'unite doit etre strictement positif")
    return f"x{value:g}", value


# ---------------------------------------------------------------------------
# Parseurs internes
# ---------------------------------------------------------------------------
def _polygon_to_triangles(indices: list[int]) -> list[tuple[int, int, int]]:
    """Triangulation en eventail d'un polygone convexe par convention."""
    return [(indices[0], indices[i], indices[i + 1]) for i in range(1, len(indices) - 1)]


def read_stl(path: str) -> TriMesh:
    """Lit un STL binaire ou ASCII."""
    with open(path, "rb") as handle:
        head = handle.read(84)
        handle.seek(0)
        payload = handle.read()
    is_binary = True
    if head[:5].lower().lstrip() .startswith(b"solid"):
        # Un STL ASCII commence par "solid" ; on confirme par la taille attendue.
        if len(payload) >= 84:
            count = struct.unpack("<I", payload[80:84])[0]
            is_binary = len(payload) == 84 + 50 * count
        else:
            is_binary = False
    if is_binary:
        count = struct.unpack("<I", payload[80:84])[0]
        expected = 84 + 50 * count
        if len(payload) < expected:
            raise ImportError_(f"STL binaire tronque : {len(payload)} octets pour {count} triangles annonces")
        vertices: list[tuple[float, float, float]] = []
        faces: list[tuple[int, int, int]] = []
        offset = 84
        for _ in range(count):
            values = struct.unpack_from("<12fH", payload, offset)
            offset += 50
            base = len(vertices)
            vertices.append((values[3], values[4], values[5]))
            vertices.append((values[6], values[7], values[8]))
            vertices.append((values[9], values[10], values[11]))
            faces.append((base, base + 1, base + 2))
        return TriMesh(vertices, faces)

    text = payload.decode("utf-8", errors="replace")
    vertices = []
    faces = []
    current: list[tuple[float, float, float]] = []
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0].lower() == "vertex" and len(parts) >= 4:
            current.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif parts[0].lower() == "endloop":
            if len(current) >= 3:
                base = len(vertices)
                vertices.extend(current)
                faces.extend(_polygon_to_triangles(list(range(base, base + len(current)))))
            current = []
    if not faces:
        raise ImportError_("STL ASCII sans triangle exploitable")
    return TriMesh(vertices, faces)


def read_obj(path: str) -> TriMesh:
    """Lit un OBJ (sommets `v` et faces `f`, les autres enregistrements sont ignores)."""
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "v" and len(parts) >= 4:
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif tag == "f" and len(parts) >= 4:
                indices = []
                for token in parts[1:]:
                    raw = token.split("/")[0]
                    if not raw:
                        continue
                    index = int(raw)
                    indices.append(index - 1 if index > 0 else len(vertices) + index)
                if len(indices) >= 3:
                    faces.extend(_polygon_to_triangles(indices))
    if not faces:
        raise ImportError_("OBJ sans face exploitable")
    return TriMesh(vertices, faces)


def read_off(path: str) -> TriMesh:
    """Lit un OFF ASCII."""
    tokens: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.split("#")[0]
            tokens.extend(line.split())
    if not tokens or not tokens[0].upper().endswith("OFF"):
        raise ImportError_("en-tete OFF absent")
    cursor = 1
    n_vertices = int(tokens[cursor])
    n_faces = int(tokens[cursor + 1])
    cursor += 3  # on saute le nombre d'aretes, non utilise
    vertices = []
    for _ in range(n_vertices):
        vertices.append((float(tokens[cursor]), float(tokens[cursor + 1]), float(tokens[cursor + 2])))
        cursor += 3
    faces: list[tuple[int, int, int]] = []
    for _ in range(n_faces):
        count = int(tokens[cursor])
        cursor += 1
        indices = [int(tokens[cursor + i]) for i in range(count)]
        cursor += count
        if count >= 3:
            faces.extend(_polygon_to_triangles(indices))
    return TriMesh(vertices, faces)


_PLY_TYPES = {
    "char": ("b", 1), "int8": ("b", 1),
    "uchar": ("B", 1), "uint8": ("B", 1),
    "short": ("h", 2), "int16": ("h", 2),
    "ushort": ("H", 2), "uint16": ("H", 2),
    "int": ("i", 4), "int32": ("i", 4),
    "uint": ("I", 4), "uint32": ("I", 4),
    "float": ("f", 4), "float32": ("f", 4),
    "double": ("d", 8), "float64": ("d", 8),
}


def read_ply(path: str) -> TriMesh:
    """Lit un PLY ASCII ou binaire (little/big endian), elements `vertex` et `face`."""
    with open(path, "rb") as handle:
        payload = handle.read()
    marker = b"end_header"
    position = payload.find(marker)
    if position < 0:
        raise ImportError_("en-tete PLY sans 'end_header'")
    header_end = payload.find(b"\n", position) + 1
    header = payload[:position].decode("ascii", errors="replace")

    fmt = "ascii"
    elements: list[tuple[str, int, list[tuple]]] = []
    for line in header.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "format":
            fmt = parts[1]
        elif parts[0] == "element":
            elements.append((parts[1], int(parts[2]), []))
        elif parts[0] == "property" and elements:
            if parts[1] == "list":
                elements[-1][2].append(("list", parts[2], parts[3], parts[4]))
            else:
                elements[-1][2].append(("scalar", parts[1], parts[2]))

    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []

    if fmt == "ascii":
        tokens = payload[header_end:].split()
        cursor = 0
        for name, count, properties in elements:
            for _ in range(count):
                values: dict[str, float] = {}
                indices: list[int] = []
                for prop in properties:
                    if prop[0] == "scalar":
                        values[prop[2]] = float(tokens[cursor])
                        cursor += 1
                    else:
                        length = int(tokens[cursor])
                        cursor += 1
                        indices = [int(tokens[cursor + i]) for i in range(length)]
                        cursor += length
                if name == "vertex":
                    vertices.append((values.get("x", 0.0), values.get("y", 0.0), values.get("z", 0.0)))
                elif name == "face" and len(indices) >= 3:
                    faces.extend(_polygon_to_triangles(indices))
    else:
        endian = "<" if "little" in fmt else ">"
        offset = header_end
        for name, count, properties in elements:
            for _ in range(count):
                values = {}
                indices = []
                for prop in properties:
                    if prop[0] == "scalar":
                        code, size = _PLY_TYPES[prop[1]]
                        values[prop[2]] = struct.unpack_from(endian + code, payload, offset)[0]
                        offset += size
                    else:
                        code, size = _PLY_TYPES[prop[1]]
                        length = struct.unpack_from(endian + code, payload, offset)[0]
                        offset += size
                        code, size = _PLY_TYPES[prop[2]]
                        indices = list(struct.unpack_from(f"{endian}{length}{code}", payload, offset))
                        offset += size * length
                if name == "vertex":
                    vertices.append((float(values.get("x", 0.0)), float(values.get("y", 0.0)), float(values.get("z", 0.0))))
                elif name == "face" and len(indices) >= 3:
                    faces.extend(_polygon_to_triangles(indices))
    if not faces:
        raise ImportError_("PLY sans face exploitable")
    return TriMesh(vertices, faces)


_NATIVE_READERS = {".stl": read_stl, ".obj": read_obj, ".ply": read_ply, ".off": read_off}


# ---------------------------------------------------------------------------
# Passerelles vers les bibliotheques externes
# ---------------------------------------------------------------------------
def _from_trimesh(obj) -> TriMesh:
    """Convertit une scene ou un maillage `trimesh` en `TriMesh`."""
    import trimesh  # import differe : dependance optionnelle

    if isinstance(obj, trimesh.Scene):
        obj = trimesh.util.concatenate([g for g in obj.geometry.values()])
    return TriMesh([tuple(v) for v in obj.vertices], [tuple(f) for f in obj.faces])


def read_with_trimesh(path: str) -> TriMesh:
    """Lit n'importe quel format supporte par `trimesh`."""
    try:
        import trimesh
    except Exception as exc:  # pragma: no cover - depend de l'environnement
        raise ImportError_(
            f"lecture de {os.path.basename(path)} : la bibliotheque 'trimesh' est requise "
            "pour ce format (pip install trimesh)"
        ) from exc
    return _from_trimesh(trimesh.load(path, force="mesh", process=False))


def read_brep(path: str) -> TriMesh:
    """Tesselle un STEP/IGES avec `cadquery` (OCC), tolerance `STEP_TESSELLATION`."""
    try:
        import cadquery  # type: ignore
    except Exception as exc:  # pragma: no cover - depend de l'environnement
        raise ImportError_(
            f"lecture de {os.path.basename(path)} : la bibliotheque 'cadquery' est requise pour "
            "les formats STEP/IGES (pip install cadquery). A defaut, exportez le modele en STL."
        ) from exc
    extension = os.path.splitext(path)[1].lower()
    shape = cadquery.importers.importStep(path) if extension in (".step", ".stp") else cadquery.importers.importShape("IGES", path)
    vertices_out: list[tuple[float, float, float]] = []
    faces_out: list[tuple[int, int, int]] = []
    for solid in shape.vals():
        tess = solid.tessellate(config.STEP_TESSELLATION)
        base = len(vertices_out)
        vertices_out.extend((float(p.x), float(p.y), float(p.z)) for p in tess[0])
        faces_out.extend((base + t[0], base + t[1], base + t[2]) for t in tess[1])
    if not faces_out:
        raise ImportError_("tessellation STEP/IGES vide")
    return TriMesh(vertices_out, faces_out)


def _dxf_pairs(path: str):
    """Couples (code de groupe, valeur) d'un DXF ASCII."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    for index in range(0, len(lines) - 1, 2):
        raw = lines[index].strip()
        if not raw.lstrip("-").isdigit():
            continue
        yield int(raw), lines[index + 1].strip()


def read_dxf_native(path: str) -> TriMesh:
    """Lit les 3DFACE et les maillages POLYFACE d'un DXF ASCII, sans dependance.

    Ce sont les deux formes sous lesquelles un maillage tesselle sort le plus
    souvent d'AutoCAD ; les MESH et les 3DSOLID non tesselles demandent `ezdxf`.
    """
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []

    entity = ""
    values: dict[int, str] = {}
    polyface_vertices: list[tuple[float, float, float]] = []
    polyface_faces: list[list[int]] = []
    in_polyface = False

    def flush() -> None:
        """Traite l'entite qui vient de se terminer."""
        nonlocal in_polyface
        if entity == "3DFACE":
            corners = []
            for corner in range(4):
                if 10 + corner not in values:
                    break
                corners.append(
                    (
                        float(values.get(10 + corner, 0.0)),
                        float(values.get(20 + corner, 0.0)),
                        float(values.get(30 + corner, 0.0)),
                    )
                )
            unique = [corners[0]] if corners else []
            for point in corners[1:]:
                if point != unique[-1]:
                    unique.append(point)
            if len(unique) >= 3:
                base = len(vertices)
                vertices.extend(unique)
                faces.extend(_polygon_to_triangles(list(range(base, base + len(unique)))))
        elif entity == "POLYLINE":
            in_polyface = int(values.get(70, "0") or 0) & 64 != 0
            polyface_vertices.clear()
            polyface_faces.clear()
        elif entity == "VERTEX" and in_polyface:
            flags = int(values.get(70, "0") or 0)
            if flags & 128 and not flags & 64:
                indices = [
                    abs(int(float(values[code]))) for code in (71, 72, 73, 74) if code in values
                ]
                polyface_faces.append([i for i in indices if i > 0])
            else:
                polyface_vertices.append(
                    (
                        float(values.get(10, 0.0)),
                        float(values.get(20, 0.0)),
                        float(values.get(30, 0.0)),
                    )
                )
        elif entity == "SEQEND" and in_polyface:
            base = len(vertices)
            vertices.extend(polyface_vertices)
            for polygon in polyface_faces:
                if len(polygon) >= 3:
                    faces.extend(
                        _polygon_to_triangles([base + index - 1 for index in polygon])
                    )
            in_polyface = False

    for code, value in _dxf_pairs(path):
        if code == 0:
            flush()
            entity = value.upper()
            values = {}
        else:
            values[code] = value
    flush()

    if not faces:
        raise ImportError_(
            "aucun 3DFACE ni maillage POLYFACE dans ce DXF : installez 'ezdxf' pour lire les "
            "entites MESH, ou tessellez les 3DSOLID avant export, ou exportez en STL"
        )
    return TriMesh(vertices, faces)


def read_dxf(path: str) -> TriMesh:
    """Lit les maillages POLYFACE / MESH / 3DSOLID tesselle d'un DXF avec `ezdxf`."""
    try:
        import ezdxf  # type: ignore
    except Exception:
        # Sans ezdxf, le lecteur interne couvre les 3DFACE et les POLYFACE, qui
        # sont les formes les plus courantes d'un maillage tesselle en DXF.
        return read_dxf_native(path)
    document = ezdxf.readfile(path)
    vertices_out: list[tuple[float, float, float]] = []
    faces_out: list[tuple[int, int, int]] = []

    def _append(points, polygons) -> None:
        base = len(vertices_out)
        vertices_out.extend((float(p[0]), float(p[1]), float(p[2])) for p in points)
        for polygon in polygons:
            indices = [base + int(i) for i in polygon]
            if len(indices) >= 3:
                faces_out.extend(_polygon_to_triangles(indices))

    for entity in document.modelspace():
        kind = entity.dxftype()
        if kind == "MESH":
            data = entity.get_data()
            _append(data.vertices, data.faces)
        elif kind == "POLYLINE" and entity.get_mode() in ("AcDbPolyFaceMesh", "AcDbPolygonMesh"):
            builder = entity.virtual_entities()
            for face in builder:
                points = [tuple(float(c) for c in p) for p in face.wcs_vertices()]
                if len(points) >= 3:
                    base = len(vertices_out)
                    vertices_out.extend(points)
                    faces_out.extend(_polygon_to_triangles(list(range(base, base + len(points)))))
        elif kind == "3DFACE":
            points = [tuple(float(c) for c in entity.dxf.get(name)) for name in ("vtx0", "vtx1", "vtx2", "vtx3")]
            unique = [points[0], points[1], points[2]]
            if points[3] != points[2]:
                unique.append(points[3])
            base = len(vertices_out)
            vertices_out.extend(unique)
            faces_out.extend(_polygon_to_triangles(list(range(base, base + len(unique)))))
    if not faces_out:
        return read_dxf_native(path)
    return TriMesh(vertices_out, faces_out)


def read_dwg(path: str) -> TriMesh:
    """Convertit un DWG en DXF via l'ODA File Converter puis delegue a `read_dxf`."""
    converter = shutil.which("ODAFileConverter") or shutil.which("ODAFileConverter.exe")
    if converter is None:
        raise ImportError_(
            "le format DWG demande l'ODA File Converter (executable 'ODAFileConverter' absent du PATH). "
            "Convertissez le fichier en DXF ou exportez-le en STL."
        )
    with tempfile.TemporaryDirectory() as workdir:
        source = os.path.join(workdir, "in")
        target = os.path.join(workdir, "out")
        os.makedirs(source)
        os.makedirs(target)
        shutil.copy(path, source)
        subprocess.run(
            [converter, source, target, "ACAD2018", "DXF", "0", "1"],
            check=True,
            capture_output=True,
        )
        produced = [f for f in os.listdir(target) if f.lower().endswith(".dxf")]
        if not produced:
            raise ImportError_("la conversion DWG -> DXF n'a produit aucun fichier")
        return read_dxf(os.path.join(target, produced[0]))


# ---------------------------------------------------------------------------
# Point d'entree
# ---------------------------------------------------------------------------
def read_raw(path: str, prefer_trimesh: bool = True) -> tuple[TriMesh, str]:
    """Lit le fichier sans traitement ; renvoie le maillage brut et le nom du lecteur."""
    if not os.path.isfile(path):
        raise ImportError_(f"fichier introuvable : {path}")
    extension = os.path.splitext(path)[1].lower()

    if extension in config.EXT_REFUSED:
        raise ImportError_(
            f"{os.path.basename(path)} : un fichier {extension} est un programme AutoLISP, "
            "pas un format geometrique. Exportez la roue en STL depuis votre logiciel de CAO "
            "(commande d'export maillage / STL), puis relancez l'analyse sur le .stl."
        )
    if extension in config.EXT_MESH_NATIVE:
        if prefer_trimesh:
            try:
                return read_with_trimesh(path), "trimesh"
            except ImportError_:
                pass
            except Exception:
                pass
        return _NATIVE_READERS[extension](path), "interne"
    if extension in config.EXT_MESH_TRIMESH_ONLY:
        return read_with_trimesh(path), "trimesh"
    if extension in config.EXT_CAD_BREP:
        return read_brep(path), "cadquery"
    if extension in config.EXT_DXF:
        return read_dxf(path), "dxf"
    if extension in config.EXT_DWG:
        return read_dwg(path), "oda+ezdxf"

    known = ", ".join(
        config.EXT_MESH_NATIVE + config.EXT_MESH_TRIMESH_ONLY + config.EXT_CAD_BREP + config.EXT_DXF + config.EXT_DWG
    )
    raise ImportError_(f"extension non geree : {extension!r} (formats lus : {known})")


def load_mesh(
    path: str,
    unit: str | float | None = None,
    do_repair: bool = True,
    prefer_trimesh: bool = True,
) -> tuple[TriMesh, ImportReport]:
    """Charge un fichier 3D et renvoie `(maillage en metres, rapport d'import)`.

    Sequence : lecture brute, conversion d'unite unique, fusion des sommets a
    `MERGE_TOL`, reparation (orientation + rebouchage), rapport.  Si le maillage
    reste non etanche, l'import continue mais `watertight` est faux et la
    confiance de toutes les grandeurs volumiques est plafonnee a `medium`
    (SPEC phase 1, point 3).
    """
    mesh, backend = read_raw(path, prefer_trimesh=prefer_trimesh)
    name, factor = unit_factor(unit)

    report = ImportReport(
        path=os.path.abspath(path),
        extension=os.path.splitext(path)[1].lower(),
        backend=backend,
        unit=name,
        unit_factor=factor,
        n_vertices_raw=len(mesh.vertices),
        n_faces_raw=len(mesh.faces),
    )
    if not mesh.faces:
        raise ImportError_(f"{os.path.basename(path)} : maillage vide")

    mesh.apply_scale(factor)  # conversion d'unite, une seule fois (SPEC 1)

    if do_repair:
        log = repair_module.repair(mesh)
        report.merged_vertices = log["merged_vertices"]
        report.removed_faces = log["removed_faces"]
        report.flipped_faces = log["flipped_faces"]
        report.filled_holes = log["filled_holes"]
    else:
        report.merged_vertices = mesh.merge_vertices(config.MERGE_TOL)

    report.n_vertices = len(mesh.vertices)
    report.n_faces = len(mesh.faces)
    report.watertight = mesh.is_watertight()
    report.winding_consistent = mesh.is_winding_consistent()
    report.volume_m3 = mesh.volume()
    report.area_m2 = mesh.area()

    lo, hi = mesh.bounds()
    report.bbox_min_mm = tuple(round(c * config.MM_PER_M, 6) for c in lo)  # type: ignore[assignment]
    report.bbox_max_mm = tuple(round(c * config.MM_PER_M, 6) for c in hi)  # type: ignore[assignment]
    report.extents_mm = tuple(
        round((hi[i] - lo[i]) * config.MM_PER_M, 6) for i in range(3)
    )  # type: ignore[assignment]

    report.confidence.set("volume", HIGH if report.watertight else MEDIUM)
    report.confidence.set("maillage", HIGH if report.watertight else MEDIUM)
    if not report.watertight:
        report.warnings.append(
            "maillage non etanche apres reparation : les grandeurs volumiques sont plafonnees "
            f"a la confiance 'medium' ({len(mesh.boundary_edges())} aretes de bord restantes)"
        )
    if report.volume_m3 <= 0.0:
        report.warnings.append(
            "volume signe negatif ou nul : orientation des normales incertaine, "
            "les tests d'inclusion peuvent etre errones"
        )
    return mesh, report
