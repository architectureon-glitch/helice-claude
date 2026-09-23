"""Stress de l'import : fichiers corrompus, variantes de format, coherence croisee.

Chaque cas est un fichier que l'on rencontre en vrai -- un STL binaire de
SolidWorks qui commence par `solid`, un OFF colore, un OBJ indexe a partir de
zero -- ou une corruption qui produisait une trace de pile au lieu d'une phrase.
La regle : un fichier invalide donne une `ImportError_` qui nomme la cause, un
fichier valide se lit, et tous les formats rendent la meme roue.
"""

from __future__ import annotations

import struct

from helpers import BaseTestCase

from impeller_analyzer import synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.io import loader, writer

TRIANGLE = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))


def binary_stl(triangles, header: bytes = b"", count: int | None = None, padding: bytes = b"") -> bytes:
    """STL binaire construit a la main, pour maitriser chaque octet."""
    body = b"".join(struct.pack("<12fH", 0, 0, 1, *t[0], *t[1], *t[2], 0) for t in triangles)
    announced = len(triangles) if count is None else count
    return header.ljust(80, b"\0")[:80] + struct.pack("<I", announced) + body + padding


class TestStlDetection(BaseTestCase):
    """Binaire ou ASCII : le mot `solid` en tete ne tranche pas a lui seul."""

    def write(self, name: str, payload: bytes) -> str:
        path = self.path(name)
        with open(path, "wb") as handle:
            handle.write(payload)
        return path

    def test_ascii_precede_d_espaces(self):
        """Un STL ASCII indente se lit en ASCII, et non comme un binaire tronque."""
        text = (
            "  solid t\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\n"
            "vertex 0 1 0\nendloop\nendfacet\nendsolid\n"
        )
        mesh, _ = loader.read_raw(self.write("a.stl", text.encode()), prefer_trimesh=False)
        self.assertEqual(len(mesh.faces), 1)

    def test_binaire_solidworks_avec_bourrage(self):
        """En-tete `solid...` et octets de bourrage : c'est un binaire, il se lit."""
        payload = binary_stl([TRIANGLE] * 4, header=b"solid SolidWorks", padding=b"\0\0\0\0")
        mesh, _ = loader.read_raw(self.write("b.stl", payload), prefer_trimesh=False)
        self.assertEqual(len(mesh.faces), 4)

    def test_texte_quelconque(self):
        """Du texte sans `solid` n'est pas annonce comme un binaire tronque."""
        path = self.write("t.stl", b"ceci n'est pas un stl\n" * 10)
        with self.assertRaises(loader.ImportError_) as caught:
            loader.read_raw(path, prefer_trimesh=False)
        self.assertIn("texte", str(caught.exception))

    def test_binaire_tronque(self):
        path = self.write("tr.stl", binary_stl([TRIANGLE] * 10, count=20))
        with self.assertRaises(loader.ImportError_) as caught:
            loader.read_raw(path, prefer_trimesh=False)
        self.assertIn("tronque", str(caught.exception))


class TestIntegrite(BaseTestCase):
    """Indices hors bornes et coordonnees non finies sont refuses a l'import."""

    def write_text(self, name: str, text: str) -> str:
        path = self.path(name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        return path

    def assertRefused(self, path: str, fragment: str) -> None:
        with self.assertRaises(loader.ImportError_) as caught:
            loader.load_mesh(path, prefer_trimesh=False)
        self.assertIn(fragment, str(caught.exception))

    def test_obj_indice_hors_bornes(self):
        self.assertRefused(self.write_text("h.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 99\n"), "n'existe pas")

    def test_obj_indexe_a_partir_de_zero(self):
        self.assertRefused(self.write_text("z.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n"), "n'existe pas")

    def test_off_indice_hors_bornes(self):
        self.assertRefused(self.write_text("h.off", "OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 7\n"), "n'existe pas")

    def test_coordonnee_infinie(self):
        path = self.path("inf.stl")
        with open(path, "wb") as handle:
            handle.write(binary_stl([((float("inf"), 0, 0), (1, 0, 0), (0, 1, 0))] * 4))
        self.assertRefused(path, "non finie")

    def test_coordonnee_nan_en_ascii(self):
        text = (
            "solid t\nfacet normal 0 0 1\nouter loop\nvertex nan 0 0\nvertex 1 0 0\n"
            "vertex 0 1 0\nendloop\nendfacet\nendsolid\n"
        )
        self.assertRefused(self.write_text("nan.stl", text), "non finie")

    def test_dxf_binaire(self):
        path = self.path("b.dxf")
        with open(path, "wb") as handle:
            handle.write(b"AutoCAD Binary DXF\r\n\x1a\x00" + b"\x00" * 64)
        self.assertRefused(path, "DXF binaire")


class TestVariantesDeFormat(BaseTestCase):
    """Variantes legales que les parseurs lisaient de travers."""

    def write_text(self, name: str, text: str) -> str:
        path = self.path(name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        return path

    def test_off_avec_couleurs_par_face(self):
        """La couleur apres les indices ne decale plus la lecture des faces suivantes."""
        path = self.write_text(
            "c.off", "OFF\n4 2 0\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n3 0 1 2 255 0 0\n3 0 2 3 0 255 0\n"
        )
        mesh, _ = loader.read_raw(path, prefer_trimesh=False)
        self.assertEqual(mesh.faces, [(0, 1, 2), (0, 2, 3)])

    def test_off_compte_sur_la_ligne_d_en_tete(self):
        path = self.write_text("e.off", "OFF 3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n")
        mesh, _ = loader.read_raw(path, prefer_trimesh=False)
        self.assertEqual(mesh.faces, [(0, 1, 2)])

    def test_off_tronque(self):
        path = self.write_text("t.off", "OFF\n4 2 0\n0 0 0\n1 0 0\n")
        with self.assertRaises(loader.ImportError_) as caught:
            loader.read_raw(path, prefer_trimesh=False)
        self.assertIn("tronque", str(caught.exception))

    def test_ply_avec_coordonnees_de_texture(self):
        """Seule la liste `vertex_indices` decrit la face, pas `texcoord`."""
        path = self.write_text(
            "t.ply",
            "ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\n"
            "property float z\nelement face 1\nproperty list uchar int vertex_indices\n"
            "property list uchar float texcoord\nend_header\n0 0 0\n1 0 0\n0 1 0\n"
            "3 0 1 2 6 0 0 1 0 0 1\n",
        )
        mesh, _ = loader.read_raw(path, prefer_trimesh=False)
        self.assertEqual(mesh.faces, [(0, 1, 2)])


class TestCoherenceEntreFormats(BaseTestCase):
    """La meme roue ecrite dans sept formats rend exactement la meme analyse."""

    def test_sept_formats_une_seule_roue(self):
        mesh = synthetic.centrifugal_impeller()
        paths = [
            writer.write_stl(mesh, self.path("r.stl")),
            writer.write_stl(mesh, self.path("ra.stl"), binary=False),
            writer.write_obj(mesh, self.path("r.obj")),
            writer.write_off(mesh, self.path("r.off")),
            writer.write_ply(mesh, self.path("r.ply")),
            writer.write_ply(mesh, self.path("rb.ply"), binary=True),
            writer.write_dxf(mesh, self.path("r.dxf")),
        ]
        signatures = set()
        for path in paths:
            data = run(path, Options(speeds=(1450.0,), grid_nr=80, grid_nz=80, n_theta=360)).to_dict()
            signatures.add(
                (
                    data["topologie"]["pales"]["nombre_de_pales"],
                    round(data["topologie"]["r_2_m"], 6),
                    round(data["pales"]["beta2_deg"], 3),
                    round(data["import"]["volume_m3"] * 1e6, 3),
                )
            )
        self.assertEqual(len(signatures), 1, signatures)


class TestFacteurDUnite(BaseTestCase):
    """Un facteur d'unite non fini rendait toutes les coordonnees NaN, sans erreur."""

    def test_facteurs_refuses(self):
        for unit in ("nan", "inf", "-1", "0", float("nan"), float("inf")):
            with self.subTest(unit=unit), self.assertRaises(loader.ImportError_):
                loader.unit_factor(unit)

    def test_facteurs_acceptes(self):
        self.assertEqual(loader.unit_factor("mm"), ("mm", 0.001))
        self.assertEqual(loader.unit_factor("0.001")[1], 0.001)
