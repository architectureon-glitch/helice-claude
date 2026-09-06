"""Phase 3 - topologie de la roue (SPEC phase 3)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.confidence import HIGH, LOW
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.geometry import topology as topo
from impeller_analyzer.geometry.proximity import TriangleGrid, hausdorff_distance
from impeller_analyzer.mesh import rotation_matrix

#: Les cartes d'occupation sont couteuses : on les construit une fois pour toutes.
_CACHE: dict[tuple, tuple] = {}


def prepared(key, builder, nr=120, nz=120):
    """Maillage aligne et carte d'occupation, mis en cache pour la suite de tests."""
    if key not in _CACHE:
        mesh, _ = axis_mod.align_to_z(builder())
        _CACHE[key] = (mesh, occ_mod.build_occupancy(mesh, nr=nr, nz=nz))
    return _CACHE[key]


class TestBladeCount(BaseTestCase):
    def test_nombre_de_pales_par_fft(self):
        """L'harmonique dominante de g(theta) donne le nombre de pales."""
        for n in (3, 4, 5, 7):
            with self.subTest(n_blades=n):
                mesh, occupancy = prepared(
                    ("axial", n), lambda n=n: synthetic.axial_impeller(n_blades=n)
                )
                count = topo.count_blades(occupancy, mesh=None)
                self.assertEqual(count.n_blades, n)
                self.assertEqual(count.confidence, HIGH)
                self.assertGreater(count.ratio_to_competing, config.FFT_RATIO_MIN)

    def test_nombre_de_pales_d_une_roue_centrifuge(self):
        """Meme detection sur une roue a aubes en arc."""
        _, occupancy = prepared(("centri", 6), synthetic.centrifugal_impeller)
        count = topo.count_blades(occupancy, mesh=None)
        self.assertEqual(count.n_blades, 6)
        self.assertEqual(count.confidence, HIGH)

    def test_controle_croise_par_hausdorff(self):
        """La rotation de 2*pi/N doit superposer le maillage a lui-meme."""
        mesh, occupancy = prepared(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        count = topo.count_blades(occupancy, mesh=mesh)
        self.assertIsNotNone(count.hausdorff_relative)
        self.assertLess(count.hausdorff_relative, config.SYM_TOL)
        self.assertEqual(count.confidence, HIGH)

    def test_hausdorff_detecte_une_mauvaise_periodicite(self):
        """Une rotation d'ordre errone donne une distance tres superieure au seuil."""
        mesh = synthetic.axial_impeller(n_blades=4)
        rotated = mesh.transformed(rotation_matrix((0.0, 0.0, 1.0), 2.0 * math.pi / 3.0))
        distance = hausdorff_distance(mesh, rotated, samples=800)
        r_tip = max(math.hypot(v[0], v[1]) for v in mesh.vertices)
        self.assertGreater(distance / r_tip, config.SYM_TOL)

    def test_distance_point_maillage_exacte(self):
        """Le hachage spatial rend la distance exacte a la surface."""
        grid = TriangleGrid(synthetic.cylinder(0.050, 0.100, 180))
        self.assertClose(grid.distance((0.080, 0.0, 0.0)), 0.030, rel=2e-3)
        self.assertClose(grid.distance((0.0, 0.0, 0.0)), 0.050, rel=2e-3)

    def test_nombre_de_pales_impose_par_l_utilisateur(self):
        """--blades prime toujours et remonte la confiance."""
        _, occupancy = prepared(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        count = topo.count_blades(occupancy, mesh=None, forced=9)
        self.assertEqual(count.n_blades, 9)
        self.assertTrue(count.forced)
        self.assertEqual(count.confidence, HIGH)
        self.assertTrue(any("impose" in w for w in count.warnings))

    def test_solide_de_revolution_sans_pale(self):
        """Sans zone de pales, la detection le dit au lieu d'inventer un nombre."""
        mesh = synthetic.cylinder(0.050, 0.100, 360)
        occupancy = occ_mod.build_occupancy(mesh, nr=60, nz=60)
        count = topo.count_blades(occupancy, mesh=None)
        self.assertEqual(count.n_blades, 0)
        self.assertEqual(count.confidence, LOW)


class TestRadii(BaseTestCase):
    def test_rayons_d_une_helice_axiale(self):
        """r_tip, r_1h et r_1s sont retrouves a mieux que 2 %."""
        _, occupancy = prepared(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        result = topo.characteristic_radii(occupancy)
        self.assertClose(result.r_tip, 0.100, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_1s, 0.100, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_1h, 0.030, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_1, math.sqrt((result.r_1s ** 2 + result.r_1h ** 2) / 2.0), rel=1e-12)
        self.assertGreater(result.z_1, result.z_2)
        self.assertFalse(result.closed_impeller)

    def test_classification_axiale(self):
        """r2/r1s proche de 1 : roue axiale."""
        _, occupancy = prepared(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        result = topo.characteristic_radii(occupancy)
        self.assertEqual(result.machine_type, topo.AXIAL)
        self.assertLess(result.ratio_r2_r1s, config.R_RATIO_AXIAL_MAX)
        # Section de sortie annulaire pour une roue axiale.
        expected = math.pi * (result.r_2s ** 2 - result.r_2h ** 2) * config.TAU_2
        self.assertClose(result.area_2, expected, rel=1e-12)

    def test_classification_centrifuge(self):
        """r2/r1s au-dela de 1.80 : roue centrifuge, section de sortie cylindrique."""
        _, occupancy = prepared(("centri", 6), synthetic.centrifugal_impeller)
        result = topo.characteristic_radii(occupancy)
        self.assertEqual(result.machine_type, topo.CENTRIFUGAL)
        self.assertGreaterEqual(result.ratio_r2_r1s, config.R_RATIO_MIXED_MAX)
        self.assertClose(result.r_1s, 0.035, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_2, 0.090, rel=config.VALID_GEOM_TOL)
        expected = 2.0 * math.pi * result.r_2 * result.b_2 * config.TAU_2
        self.assertClose(result.area_2, expected, rel=1e-12)

    def test_classification_mixte(self):
        """Un rapport intermediaire donne bien une roue mixte."""
        _, occupancy = prepared(
            ("mixte",),
            lambda: synthetic.centrifugal_impeller(r1=0.060, r2=0.090, b1=0.018, b2=0.012, eye_height=0.030),
        )
        result = topo.characteristic_radii(occupancy)
        self.assertEqual(result.machine_type, topo.MIXED)
        self.assertGreaterEqual(result.ratio_r2_r1s, config.R_RATIO_AXIAL_MAX)
        self.assertLess(result.ratio_r2_r1s, config.R_RATIO_MIXED_MAX)

    def test_rayon_d_aspiration_utilisateur_prime(self):
        """--r-aspiration (en cm) remplace la valeur detectee et remonte la confiance."""
        _, occupancy = prepared(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        result = topo.apply_user_suction_radius(topo.characteristic_radii(occupancy), 4.5)
        self.assertClose(result.r_aspiration, 0.045, rel=1e-12)
        self.assertClose(result.r_1s, 0.045, rel=1e-12)
        self.assertEqual(result.r_aspiration_source, "utilisateur")
        self.assertEqual(result.confidence["rayons"], HIGH)
        self.assertTrue(any("impose" in w for w in result.warnings))
        with self.assertRaises(ValueError):
            topo.apply_user_suction_radius(topo.characteristic_radii(occupancy), -1.0)

    def test_ilots_parasites_ecartes(self):
        """Les ilots minuscules de la zone de pales sont retires du masque."""
        mask = [[False] * 40 for _ in range(40)]
        for iz in range(10, 30):
            for ir in range(10, 30):
                mask[iz][ir] = True  # gros ilot : conserve
        mask[2][2] = True  # cellule isolee : ecartee
        cleaned = topo.clean_blade_mask(mask)
        self.assertTrue(cleaned[20][20])
        self.assertFalse(cleaned[2][2])


class TestCrossCheck(BaseTestCase):
    def test_vitesse_specifique_et_famille(self):
        """n_q classe la machine selon les bornes 35 et 80 de la SPEC."""
        self.assertClose(topo.specific_speed(1450.0, 0.025, 20.0), 1450.0 * math.sqrt(0.025) / 20.0 ** 0.75, rel=1e-12)
        self.assertEqual(topo.type_from_specific_speed(20.0), topo.CENTRIFUGAL)
        self.assertEqual(topo.type_from_specific_speed(50.0), topo.MIXED)
        self.assertEqual(topo.type_from_specific_speed(120.0), topo.AXIAL)

    def test_incoherence_signalee(self):
        """Une classification geometrique en desaccord avec n_q produit un avertissement."""
        result = topo.Topology()
        result.machine_type = topo.CENTRIFUGAL
        result.ratio_r2_r1s = 2.5
        self.assertIsNone(topo.cross_check_type(result, 20.0))
        message = topo.cross_check_type(result, 120.0)
        self.assertIsNotNone(message)
        self.assertIn("incoherence", message)
        self.assertIn("axiale", message)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
