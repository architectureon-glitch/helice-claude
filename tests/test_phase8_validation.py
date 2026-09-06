"""Phase 8 - campagne de validation et elargissement des formats (SPEC phase 8)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic, validation
from impeller_analyzer.geometry import topology as topo
from impeller_analyzer.io import loader, writer

#: Grille des controles : celle de production, la precision visee etant de 2 %.
GRID = (150, 150)


class TestCheckArithmetic(BaseTestCase):
    def test_ecart_relatif_et_absolu(self):
        """Un controle exprime son ecart en relatif ou en absolu, selon la grandeur."""
        relative = validation.Check("r", 1.02, 1.00, 0.05)
        self.assertClose(relative.deviation, 0.02, rel=1e-12)
        self.assertTrue(relative.passed)

        absolute = validation.Check("beta", 27.0, 25.0, 1.5, "deg", relative=False)
        self.assertClose(absolute.deviation, 2.0, rel=1e-12)
        self.assertFalse(absolute.passed)
        self.assertIn("ECHEC", absolute.line())
        self.assertIn("deg", absolute.to_dict()["unite"])

    def test_valeur_attendue_nulle(self):
        """Une valeur attendue nulle n'explose pas."""
        self.assertClose(validation.Check("x", 0.0, 0.0, 0.0).deviation, 0.0, abs_=0.0)
        self.assertEqual(validation.Check("x", 1.0, 0.0, 0.5).deviation, math.inf)

    def test_rapport_agrege_les_controles(self):
        """Un rapport est conforme si tous ses controles le sont."""
        report = validation.Report("essai")
        report.checks.append(validation.Check("a", 1.0, 1.0, 0.0))
        self.assertTrue(report.passed)
        report.checks.append(validation.Check("b", 2.0, 1.0, 0.1))
        self.assertFalse(report.passed)
        self.assertIn("ECHEC", report.render())
        self.assertFalse(report.to_dict()["conforme"])


class TestSyntheticGeometries(BaseTestCase):
    """SPEC 8.1 : geometries a reponse analytique connue."""

    def test_helice_a_pales_planes(self):
        """N, r1s, beta et le sens de rotation sont retrouves dans les tolerances."""
        report = validation.check_synthetic(
            "helice",
            synthetic.axial_impeller(n_blades=4, beta_deg=30.0, r_hub=0.030, r_tip=0.100),
            {
                "n_blades": 4,
                "r_1s": 0.100,
                "beta1_deg": 30.0,
                "beta2_deg": 30.0,
                "machine_type": topo.AXIAL,
                "rotation_sign": 1,
            },
            *GRID,
        )
        self.assertTrue(report.passed, report.render())

    def test_roue_centrifuge_a_aubes_en_arc(self):
        """Meme controle sur une roue centrifuge."""
        report = validation.check_synthetic(
            "roue centrifuge",
            synthetic.centrifugal_impeller(n_blades=6, beta1_deg=22.0, beta2_deg=25.0),
            {
                "n_blades": 6,
                "r_1s": 0.035,
                "r_2": 0.090,
                "beta1_deg": 22.0,
                "beta2_deg": 25.0,
                "machine_type": topo.CENTRIFUGAL,
                "rotation_sign": -1,
            },
            *GRID,
        )
        self.assertTrue(report.passed, report.render())

    def test_ecart_detecte(self):
        """Une attente fausse est bien signalee comme non conforme."""
        report = validation.check_synthetic(
            "helice",
            synthetic.axial_impeller(n_blades=4, beta_deg=30.0),
            {"n_blades": 5, "beta1_deg": 30.0, "beta2_deg": 30.0},
            80,
            80,
        )
        self.assertFalse(report.passed)


class TestInvarianceAndRobustness(BaseTestCase):
    def test_invariance_par_rotation_et_translation(self):
        """SPEC 8.3 : rotation de 37 deg et translation, resultats identiques a 0.5 % pres."""
        report = validation.check_invariance(
            synthetic.centrifugal_impeller(), angle_deg=37.0, nr=GRID[0], nz=GRID[1]
        )
        self.assertTrue(report.passed, report.render())
        for check in report.checks:
            self.assertLessEqual(check.deviation, config.VALID_INVARIANCE_TOL)

    def test_robustesse_a_la_decimation(self):
        """SPEC 8.4 : maillage decime a 20 %, beta2 conserve a moins de 3 degres."""
        report = validation.check_robustness(
            synthetic.centrifugal_impeller(), ratio=config.VALID_DECIMATION_RATIO, nr=GRID[0], nz=GRID[1]
        )
        self.assertTrue(report.passed, report.render())

    def test_decimation_conserve_le_solide(self):
        """La decimation atteint la cible, reste etanche et conserve le volume."""
        mesh = synthetic.centrifugal_impeller()
        coarse = mesh.decimate(config.VALID_DECIMATION_RATIO)
        ratio = len(coarse.faces) / len(mesh.faces)
        self.assertClose(ratio, config.VALID_DECIMATION_RATIO, abs_=0.02)
        self.assertTrue(coarse.is_watertight())
        self.assertClose(coarse.volume(), mesh.volume(), rel=0.20)


class TestReferenceCase(BaseTestCase):
    """SPEC 8.2 : cas de reference exterieur."""

    GEOMETRY = {
        "r_1s_m": 0.0275, "r_1h_m": 0.0120, "r_1_m": 0.02122, "r_2_m": 0.1000,
        "b_2_m": 0.0080, "beta1_deg": 20.0, "beta2_deg": 25.0, "n_blades": 6,
    }

    def _case(self, head: float) -> dict:
        return {
            "designation": "essai",
            "geometrie": dict(self.GEOMETRY),
            "point_publie": {"regime_tr_min": 2900.0, "Q_m3_h": 12.5, "H_m": head},
        }

    def test_mecanique_du_controle(self):
        """Le controle compare la hauteur au BEP a la valeur publiee."""
        measured = validation.check_reference_case(self._case(50.0)).checks[0].measured
        self.assertGreater(measured, 0.0)
        # Une valeur publiee egale au calcul passe ; un ecart de 40 % ne passe pas.
        self.assertTrue(validation.check_reference_case(self._case(measured)).passed)
        self.assertFalse(validation.check_reference_case(self._case(measured * 1.40)).passed)
        self.assertTrue(
            validation.check_reference_case(
                self._case(measured * (1.0 + 0.5 * config.VALID_REFERENCE_TOL))
            ).passed
        )

    def test_chargement_du_fichier_livre(self):
        """Le cas de reference livre se charge et se traite."""
        import os

        path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "examples",
            "cas_de_reference.json",
        )
        case = validation.load_reference_case(path)
        self.assertIn("A REMPLACER", case["source"])
        report = validation.check_reference_case(case)
        self.assertEqual(len(report.checks), 1)
        self.assertGreater(report.checks[0].measured, 0.0)

    def test_geometrie_inexploitable(self):
        """Une geometrie qui ne pompe pas est signalee, sans exception."""
        case = self._case(50.0)
        case["geometrie"]["beta2_deg"] = config.BETA_MIN_DEG
        report = validation.check_reference_case(case)
        self.assertFalse(report.passed)


class TestExtendedFormats(BaseTestCase):
    def test_dxf_lu_sans_bibliotheque_externe(self):
        """Un DXF de 3DFACE est lu par le parseur interne, sans ezdxf."""
        path = self.path("cube.dxf")
        writer.write_dxf(synthetic.cube(0.10), path, unit_factor=config.UNIT_FACTOR)
        mesh = loader.read_dxf_native(path)
        self.assertEqual(len(mesh.faces), 12)
        _, report = loader.load_mesh(path)
        self.assertClose(report.volume_m3, 1.0e-3, rel=1e-9)
        self.assertEqual(report.backend, "dxf")

    def test_dxf_sans_maillage(self):
        """Un DXF sans entite exploitable donne un message explicite."""
        path = self.path("vide.dxf")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("0\nSECTION\n2\nENTITIES\n0\nENDSEC\n0\nEOF\n")
        with self.assertRaises(loader.ImportError_) as ctx:
            loader.read_dxf_native(path)
        self.assertIn("3DFACE", str(ctx.exception))

    def test_extensions_declarees(self):
        """Les quatre familles de formats de la SPEC sont declarees."""
        for extension in (".stl", ".obj", ".ply", ".off"):
            self.assertIn(extension, config.EXT_MESH_NATIVE)
        self.assertIn(".3ds", config.EXT_MESH_TRIMESH_ONLY)
        for extension in (".step", ".stp", ".iges"):
            self.assertIn(extension, config.EXT_CAD_BREP)
        self.assertIn(".dxf", config.EXT_DXF)
        self.assertIn(".dwg", config.EXT_DWG)


class TestCampaign(BaseTestCase):
    def test_campagne_par_defaut(self):
        """La campagne livree est conforme dans son ensemble."""
        reports = validation.default_cases(*GRID)
        self.assertEqual(len(reports), 5)
        for report in reports:
            with self.subTest(controle=report.title):
                self.assertTrue(report.passed, report.render())
        self.assertTrue(all(isinstance(r.to_dict(), dict) for r in reports))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
