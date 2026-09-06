"""Phase 2 - detection d'axe, recentrage et carte d'occupation (SPEC phase 2)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.confidence import HIGH, LOW
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.io import plot
from impeller_analyzer.mesh import rotation_matrix


class TestAxis(BaseTestCase):
    def test_axe_d_un_cylindre_est_z(self):
        """Un solide de revolution donne deux valeurs propres egales et l'axe Z."""
        result = axis_mod.detect_axis(synthetic.cylinder(0.05, 0.20, 180))
        self.assertSeqClose(result.axis, (0.0, 0.0, 1.0), abs_=1e-9)
        self.assertLess(result.pair_ratio, config.AXIS_TOL)
        self.assertEqual(result.confidence, HIGH)
        self.assertClose(result.angle_to_z_deg, 0.0, abs_=1e-9)

    def test_axe_d_une_roue_a_pales(self):
        """Une roue a N pales garde deux valeurs propres egales : axe detecte a coup sur."""
        for mesh in (synthetic.axial_impeller(), synthetic.centrifugal_impeller()):
            with self.subTest(mesh=repr(mesh)):
                result = axis_mod.detect_axis(mesh)
                self.assertEqual(result.confidence, HIGH)
                self.assertClose(result.angle_to_z_deg, 0.0, abs_=1e-6)

    def test_maillage_incline_est_realigne_avec_avertissement(self):
        """Ecart > AXIS_WARN_DEG : avertissement, realignement, axe ramene sur Z."""
        mesh = synthetic.axial_impeller()
        mesh.apply_rotation(rotation_matrix((1.0, 0.0, 0.0), math.radians(20.0)))
        mesh.apply_translation((0.30, -0.20, 0.50))
        detected = axis_mod.detect_axis(mesh)
        self.assertClose(detected.angle_to_z_deg, 20.0, abs_=1e-3)
        self.assertTrue(any("realigne" in w for w in detected.warnings))

        aligned, result = axis_mod.align_to_z(mesh, detected)
        self.assertTrue(result.realigned)
        self.assertSeqClose(aligned.centroid(), (0.0, 0.0, 0.0), abs_=1e-9)
        self.assertClose(axis_mod.detect_axis(aligned).angle_to_z_deg, 0.0, abs_=1e-6)

    def test_faible_ecart_ne_declenche_pas_de_rotation(self):
        """En deca du seuil, la convention Z est appliquee sans rotation parasite."""
        mesh = synthetic.axial_impeller()
        mesh.apply_rotation(rotation_matrix((1.0, 0.0, 0.0), math.radians(2.0)))
        aligned, result = axis_mod.align_to_z(mesh)
        self.assertLess(result.angle_to_z_deg, config.AXIS_WARN_DEG)
        self.assertFalse(result.realigned)

    def test_corps_isotrope_ne_donne_pas_d_axe(self):
        """Un cube a trois valeurs propres egales : confiance low, axe Z conserve."""
        result = axis_mod.detect_axis(synthetic.cube(0.10))
        self.assertEqual(result.confidence, LOW)
        self.assertSeqClose(result.axis, (0.0, 0.0, 1.0), abs_=1e-12)
        self.assertTrue(any("isotrope" in w for w in result.warnings))

    def test_recentrage_place_l_origine_au_barycentre(self):
        """Apres recentrage l'origine est sur l'axe, a la hauteur du barycentre."""
        mesh = synthetic.cylinder(0.05, 0.20, 120, z_center=0.37)
        mesh.apply_translation((0.10, 0.20, 0.0))
        aligned, _ = axis_mod.align_to_z(mesh)
        self.assertSeqClose(aligned.centroid(), (0.0, 0.0, 0.0), abs_=1e-9)


class TestOccupancy(BaseTestCase):
    def test_cylindre_plein_donne_f_egal_un(self):
        """Test impose par la SPEC : f = 1 partout dans l'emprise du cylindre."""
        mesh = synthetic.cylinder(0.05, 0.10, 720)
        occupancy = occ_mod.build_occupancy(mesh, nr=50, nz=50)
        for row in occupancy.f:
            for value in row:
                self.assertClose(value, 1.0, abs_=1e-12)
        # Hors emprise : au-dela du rayon maximal et de la hauteur, f est nul.
        self.assertEqual(occupancy.value(0.06, 0.0), 0.0)
        self.assertEqual(occupancy.value(0.01, 0.20), 0.0)

    def test_tube_creux_separe_plein_et_vide(self):
        """Une couronne donne f = 1 entre les deux rayons et 0 a l'interieur."""
        mesh = synthetic.tube(0.020, 0.050, 0.10, 720)
        occupancy = occ_mod.build_occupancy(mesh, nr=100, nz=40)
        self.assertClose(occupancy.value(0.010, 0.0), 0.0, abs_=1e-12)
        self.assertClose(occupancy.value(0.035, 0.0), 1.0, abs_=1e-12)
        self.assertClose(occupancy.value(0.019, 0.0), 0.0, abs_=1e-12)

    def test_union_de_solides_interpenetres(self):
        """Le comptage signe traite bien une union : pas d'inversion dans le recouvrement."""
        mesh = synthetic.combine([
            synthetic.cylinder(0.030, 0.10, 360),
            synthetic.box(0.020, 0.020, 0.20),  # traverse le cylindre de part en part
        ])
        occupancy = occ_mod.build_occupancy(mesh, nr=60, nz=60)
        # Sur l'axe, dans la zone de recouvrement, la matiere est bien pleine.
        self.assertClose(occupancy.value(0.002, 0.0), 1.0, abs_=1e-12)
        # Au-dessus du cylindre il ne reste que la barre carree de 20 mm : a
        # 12 mm de l'axe le cercle sort du carre, l'occupation est partielle.
        partial = occupancy.value(0.012, 0.080)
        self.assertGreater(partial, config.F_VIDE)
        self.assertLess(partial, config.F_SOLIDE)

    def test_coherence_avec_le_test_d_inclusion_direct(self):
        """La carte doit concorder avec un test d'inclusion point par point."""
        mesh = synthetic.tube(0.020, 0.050, 0.10, 360)
        occupancy = occ_mod.build_occupancy(mesh, nr=40, nz=20)
        for r, z, expected in ((0.010, 0.0, False), (0.035, 0.0, True), (0.045, 0.03, True)):
            with self.subTest(r=r, z=z):
                point = (r * math.cos(0.7), r * math.sin(0.7), z)
                self.assertEqual(occ_mod.is_inside(mesh, point), expected)
                self.assertEqual(occupancy.value(r, z) > 0.5, expected)

    def test_masques_et_signal_angulaire(self):
        """Les trois masques partitionnent la grille ; g(theta) porte la periodicite N."""
        from impeller_analyzer.numeric import dft_magnitudes

        mesh, _ = axis_mod.align_to_z(synthetic.axial_impeller(n_blades=5))
        occupancy = occ_mod.build_occupancy(mesh, nr=100, nz=100)
        solid = occupancy.solid_mask()
        blade = occupancy.blade_mask()
        void = occupancy.void_mask()
        for iz in range(occupancy.nz):
            for ir in range(occupancy.nr):
                self.assertEqual(sum((solid[iz][ir], blade[iz][ir], void[iz][ir])), 1)

        signal = occupancy.theta_signal()
        self.assertEqual(len(signal), config.N_THETA)
        magnitudes = dft_magnitudes(signal, config.BLADES_MAX)
        dominant = max(range(config.BLADES_MIN, config.BLADES_MAX + 1), key=lambda k: magnitudes[k])
        self.assertEqual(dominant, 5)

    def test_carte_png_produite(self):
        """La visualisation PNG de f(r, z) est produite et bien formee."""
        mesh, _ = axis_mod.align_to_z(synthetic.axial_impeller())
        occupancy = occ_mod.build_occupancy(mesh, nr=60, nz=60)
        path = self.path("carte.png")
        plot.occupancy_png(occupancy, path)
        with open(path, "rb") as handle:
            payload = handle.read()
        self.assertEqual(payload[:8], b"\x89PNG\r\n\x1a\n")
        self.assertGreater(len(payload), 1000)

    def test_echelle_de_couleur_calee_sur_les_seuils(self):
        """L'echelle d'affichage separe nettement vide, pales et plein."""
        self.assertClose(plot.occupancy_scale(0.0), 0.0, abs_=1e-12)
        self.assertClose(plot.occupancy_scale(config.F_VIDE), config.PLOT_BAND_VIDE, abs_=1e-12)
        self.assertClose(plot.occupancy_scale(config.F_SOLIDE), config.PLOT_BAND_PLEIN, abs_=1e-12)
        self.assertClose(plot.occupancy_scale(1.0), 1.0, abs_=1e-12)
        self.assertGreater(plot.occupancy_scale(0.05), plot.occupancy_scale(config.F_VIDE))

    def test_maillage_degenere_rejete(self):
        """Un maillage sans extension radiale ne peut pas produire de carte."""
        flat = synthetic.box(0.0, 0.0, 0.0)
        with self.assertRaises(ValueError):
            occ_mod.build_occupancy(flat, nr=10, nz=10)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
