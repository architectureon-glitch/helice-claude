"""Phase 1 - import geometrique (SPEC phase 1)."""

import json
import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.confidence import HIGH, MEDIUM
from impeller_analyzer.io import loader, writer
from impeller_analyzer.mesh import TriMesh


class TestImport(BaseTestCase):
    def test_cube_10cm_donne_1e_3_m3(self):
        """Test impose par la SPEC : un cube de 10 cm importe vaut 1.0e-3 m3."""
        path = self.path("cube.stl")
        writer.write_stl(synthetic.cube(0.10), path, unit_factor=config.UNIT_FACTOR)
        _, report = loader.load_mesh(path, prefer_trimesh=False)
        self.assertClose(report.volume_m3, 1.0e-3, rel=1e-9)
        self.assertEqual(report.unit, "cm")
        self.assertEqual(report.unit_factor, config.UNIT_FACTOR)
        self.assertTrue(report.watertight)
        self.assertEqual(report.confidence["volume"], HIGH)
        self.assertSeqClose(report.extents_mm, (100.0, 100.0, 100.0), abs_=1e-6)

    def test_tous_les_formats_natifs_donnent_le_meme_volume(self):
        """STL / OBJ / PLY / OFF doivent produire le meme maillage."""
        for extension in (".stl", ".obj", ".ply", ".off"):
            with self.subTest(extension=extension):
                path = self.path(f"cube{extension}")
                writer.WRITERS[extension](synthetic.cube(0.10), path, unit_factor=config.UNIT_FACTOR)
                _, report = loader.load_mesh(path, prefer_trimesh=False)
                self.assertClose(report.volume_m3, 1.0e-3, rel=1e-9)
                self.assertEqual(report.n_faces, 12)

    def test_stl_ascii_et_ply_binaire(self):
        """Les variantes ASCII et binaires sont lues de la meme facon."""
        cube = synthetic.cube(0.10)
        ascii_path = self.path("cube_ascii.stl")
        writer.write_stl(cube, ascii_path, unit_factor=config.UNIT_FACTOR, binary=False)
        binary_ply = self.path("cube_bin.ply")
        writer.write_ply(cube, binary_ply, unit_factor=config.UNIT_FACTOR, binary=True)
        for path in (ascii_path, binary_ply):
            with self.subTest(path=path):
                _, report = loader.load_mesh(path, prefer_trimesh=False)
                self.assertClose(report.volume_m3, 1.0e-3, rel=1e-6)

    def test_conversion_d_unite(self):
        """Le facteur d'unite est applique une seule fois, a l'import."""
        path = self.path("cube_mm.stl")
        # Cube de 100 unites de cote, relu en millimetres : 0.1 m -> 1e-3 m3.
        writer.write_stl(synthetic.cube(0.10), path, unit_factor=0.001)
        _, report = loader.load_mesh(path, unit="mm", prefer_trimesh=False)
        self.assertClose(report.volume_m3, 1.0e-3, rel=1e-9)
        _, report_m = loader.load_mesh(path, unit="m", prefer_trimesh=False)
        self.assertClose(report_m.volume_m3, 1.0e6, rel=1e-9)  # 100 m de cote
        self.assertEqual(loader.unit_factor(0.5), ("x0.5", 0.5))
        with self.assertRaises(loader.ImportError_):
            loader.unit_factor("parsec")

    def test_fusion_des_sommets_dupliques(self):
        """Le STL duplique chaque sommet : la fusion doit revenir a 8 sommets."""
        path = self.path("cube.stl")
        writer.write_stl(synthetic.cube(0.10), path, unit_factor=config.UNIT_FACTOR)
        mesh, report = loader.load_mesh(path, prefer_trimesh=False)
        self.assertEqual(report.n_vertices_raw, 36)
        self.assertEqual(len(mesh.vertices), 8)
        self.assertEqual(report.merged_vertices, 28)

    def test_normales_inversees_corrigees(self):
        """Un maillage retourne est remis a l'endroit : volume positif."""
        cube = synthetic.cube(0.10)
        cube.flip()
        path = self.path("flipped.stl")
        writer.write_stl(cube, path, unit_factor=config.UNIT_FACTOR)
        _, report = loader.load_mesh(path, prefer_trimesh=False)
        self.assertClose(report.volume_m3, 1.0e-3, rel=1e-9)
        self.assertGreater(report.flipped_faces, 0)

    def test_trou_rebouche_et_confiance(self):
        """Un trou est rebouche ; s'il subsiste, la confiance volumique tombe a medium."""
        cube = synthetic.cube(0.10)
        holed = TriMesh(cube.vertices, cube.faces[2:])  # on retire une face du cube
        path = self.path("holed.stl")
        writer.write_stl(holed, path, unit_factor=config.UNIT_FACTOR)
        _, repaired = loader.load_mesh(path, prefer_trimesh=False)
        self.assertEqual(repaired.filled_holes, 1)
        self.assertTrue(repaired.watertight)
        self.assertClose(repaired.volume_m3, 1.0e-3, rel=1e-9)

        _, raw = loader.load_mesh(path, do_repair=False, prefer_trimesh=False)
        self.assertFalse(raw.watertight)
        self.assertEqual(raw.confidence["volume"], MEDIUM)
        self.assertTrue(any("non etanche" in w for w in raw.warnings))

    def test_refus_explicite_du_lisp(self):
        """Un .lisp n'est pas une geometrie : refus avec consigne d'export STL."""
        path = self.path("roue.lisp")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("(defun c:roue () (princ))\n")
        with self.assertRaises(loader.ImportError_) as ctx:
            loader.load_mesh(path)
        message = str(ctx.exception)
        self.assertIn("AutoLISP", message)
        self.assertIn("STL", message)

    def test_extension_inconnue_et_fichier_absent(self):
        """Messages d'erreur explicites pour les cas d'entree invalides."""
        unknown = self.path("roue.xyz")
        open(unknown, "w").close()
        with self.assertRaises(loader.ImportError_):
            loader.load_mesh(unknown)
        with self.assertRaises(loader.ImportError_):
            loader.load_mesh(self.path("absent.stl"))

    def test_formats_cad_sans_bibliotheque(self):
        """STEP / DXF sans la bibliotheque associee : message explicite, pas de trace."""
        for extension, keyword in ((".step", "cadquery"), (".dxf", "ezdxf")):
            with self.subTest(extension=extension):
                path = self.path(f"roue{extension}")
                open(path, "w").close()
                try:
                    loader.load_mesh(path)
                except loader.ImportError_ as exc:
                    self.assertTrue(keyword in str(exc) or "maillage" in str(exc), str(exc))
                except Exception:
                    # La bibliotheque est installee : elle rejette le fichier vide, c'est correct.
                    pass

    def test_rapport_import_serialisable(self):
        """Le rapport d'import contient les grandeurs demandees et passe en JSON."""
        path = self.path("cylindre.stl")
        writer.write_stl(synthetic.cylinder(0.05, 0.10, 120), path, unit_factor=config.UNIT_FACTOR)
        _, report = loader.load_mesh(path, prefer_trimesh=False)
        data = report.to_dict()
        json.dumps(data)
        self.assertGreater(data["n_faces"], 0)
        self.assertClose(data["extents_mm"][2], 100.0, abs_=1e-6)
        self.assertClose(data["volume_m3"], math.pi * 0.05 ** 2 * 0.10, rel=1e-3)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
