"""Cote aspiration et roues a flasque : ce que l'outil doit mesurer, ce qu'il doit ignorer.

Une roue de pompe aspire par le centre et refoule lateralement. Deux pieges la
guettent a l'import :

* le fichier peut etre a l'envers -- rien dans un STL n'impose l'aspiration vers
  +Z -- et la roue est alors lue depuis son cote refoulement ;
* le disque arriere, ou le flasque avant, deborde souvent les aubes, et les
  confondre avec elles fausse les rayons caracteristiques.
"""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.confidence import HIGH, LOW
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.geometry import topology as topo
from impeller_analyzer.io import writer
from impeller_analyzer.mesh import rotation_matrix

GRID = 120
FLIP = rotation_matrix((1.0, 0.0, 0.0), math.pi)


def occupancy_of(mesh, grid=GRID):
    """Carte d'occupation d'un maillage aligne sur Z."""
    aligned, _ = axis_mod.align_to_z(mesh)
    return aligned, occ_mod.build_occupancy(aligned, nr=grid, nz=grid)


class TestSuctionSide(BaseTestCase):
    def test_asymetrie_signe_la_veine(self):
        """La veine d'une roue radiale part du petit rayon ; celle d'une helice, non."""
        radiale = synthetic.centrifugal_impeller()
        _, droite = occupancy_of(radiale)
        _, envers = occupancy_of(radiale.transformed(matrix=FLIP))
        _, axiale = occupancy_of(
            synthetic.axial_impeller(n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0)
        )
        self.assertLess(axis_mod.meridional_asymmetry(droite), -config.SUCTION_ASYMMETRY_MIN)
        self.assertGreater(axis_mod.meridional_asymmetry(envers), config.SUCTION_ASYMMETRY_MIN)
        self.assertLess(abs(axis_mod.meridional_asymmetry(axiale)), config.SUCTION_ASYMMETRY_MIN)

    def test_detection_du_cote_aspiration(self):
        """Une roue radiale a l'envers est reconnue comme telle et signalee."""
        roue = synthetic.centrifugal_impeller(front_shroud=True)
        _, droite = occupancy_of(roue)
        _, envers = occupancy_of(roue.transformed(matrix=FLIP))

        conforme = axis_mod.detect_suction_side(droite)
        self.assertEqual(conforme.sign, 1)
        self.assertFalse(conforme.flipped)
        self.assertEqual(conforme.confidence, HIGH)
        self.assertEqual(conforme.warnings, [])

        retournee = axis_mod.detect_suction_side(envers)
        self.assertEqual(retournee.sign, -1)
        self.assertTrue(retournee.flipped)
        self.assertEqual(retournee.confidence, HIGH)
        self.assertTrue(any("a l'envers" in w for w in retournee.warnings))

    def test_helice_axiale_indeterminee(self):
        """Sur une helice, la geometrie ne tranche pas : la convention est conservee."""
        _, occupancy = occupancy_of(
            synthetic.axial_impeller(n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0)
        )
        result = axis_mod.detect_suction_side(occupancy)
        self.assertEqual(result.sign, 0)
        self.assertFalse(result.flipped)
        self.assertEqual(result.confidence, LOW)
        self.assertTrue(any("--aspiration -z" in note for note in result.notes))

    def test_cote_impose_par_l_utilisateur(self):
        """--aspiration prime sur la detection."""
        _, occupancy = occupancy_of(synthetic.centrifugal_impeller())
        forced = axis_mod.detect_suction_side(occupancy, axis_mod.SUCTION_MINUS_Z)
        self.assertEqual(forced.sign, -1)
        self.assertTrue(forced.flipped)
        self.assertTrue(forced.forced)
        kept = axis_mod.detect_suction_side(occupancy, axis_mod.SUCTION_PLUS_Z)
        self.assertEqual(kept.sign, 1)
        self.assertFalse(kept.flipped)

    def test_retournement_est_un_deplacement_rigide(self):
        """Retourner le maillage ne change ni son volume ni son etancheite."""
        mesh = synthetic.centrifugal_impeller()
        flipped = axis_mod.flip_axis(mesh)
        self.assertClose(flipped.volume(), mesh.volume(), rel=1e-12)
        self.assertEqual(flipped.is_watertight(), mesh.is_watertight())
        self.assertClose(flipped.bounds()[1][2], -mesh.bounds()[0][2], rel=1e-12)


class TestShroudedGeometry(BaseTestCase):
    def test_roue_fermee_reconnue(self):
        """Le flasque avant est detecte, et l'oeillard n'est pas confondu avec lui."""
        _, occupancy = occupancy_of(synthetic.centrifugal_impeller(front_shroud=True))
        result = topo.characteristic_radii(occupancy)
        self.assertTrue(result.closed_impeller)
        self.assertEqual(result.machine_type, topo.CENTRIFUGAL)
        # r_1s est le rayon des pales au plan d'aspiration, pas celui du flasque.
        self.assertClose(result.r_1s, 0.035, rel=config.VALID_GEOM_TOL)

    def test_flasque_plat_l_oeillard_fixe_le_rayon_d_aspiration(self):
        """Quand les pales courent sous le flasque, r_1s est le percement, pas leur bout."""
        roue = synthetic.centrifugal_impeller(
            front_shroud=True, flat_shroud=True, r1=0.035, r2=0.090
        )
        _, occupancy = occupancy_of(roue)
        result = topo.characteristic_radii(occupancy)
        self.assertTrue(result.closed_impeller)

        # Au plan d'entree les aubes vont d'un bord a l'autre : les lire donnerait
        # le rayon exterieur de la roue, soit deux fois et demie le bon rayon.
        blade = topo.clean_blade_mask(occupancy.blade_mask())
        rows = [
            iz for iz in range(occupancy.nz)
            if sum(1 for cell in blade[iz] if cell) >= config.BLADE_ROW_MIN_CELLS
        ]
        sur_les_pales = occupancy.r_centres[topo._blade_outer_index(blade[rows[-1]], occupancy.f[rows[-1]])]
        self.assertClose(sur_les_pales, 0.090, rel=0.02)

        self.assertClose(result.r_1s, 0.035, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_2, 0.090, rel=config.VALID_GEOM_TOL)
        self.assertEqual(result.machine_type, topo.CENTRIFUGAL)
        self.assertEqual(result.warnings, [])

    def test_flasque_conique_l_oeillard_reste_au_bord_d_attaque(self):
        """Sur un flasque incline, les deux regles coincident : pas de regression."""
        _, occupancy = occupancy_of(synthetic.centrifugal_impeller(front_shroud=True))
        result = topo.characteristic_radii(occupancy)
        self.assertClose(result.r_1s, 0.035, rel=config.VALID_GEOM_TOL)
        self.assertEqual(result.warnings, [])

    def test_le_plateau_arriere_n_est_pas_un_moyeu_de_sortie(self):
        """Sur une roue fermee classee mixte, b2 ne doit pas se refermer sur rien.

        Au plan de fuite, la recherche de moyeu attrape le plateau arriere, qui
        occupe tout le rayon : `r_2h` vaut alors presque `r_2s` et la formule
        annulaire `b_2 = r_2s - r_2h` s'effondre -- 1.2 mm, soit le pas de la
        grille, pour une sortie qui en fait 25. La section de refoulement se lit
        la ou elle est ouverte : en hauteur, au bout des aubes.
        """
        roue = synthetic.centrifugal_impeller(
            n_blades=5, beta1_deg=10.0, beta2_deg=20.0,
            r1=0.0927, r2=0.1670, b1=0.030, b2=0.025, eye_height=0.050,
            thickness=0.012, front_shroud=True, shroud_thickness=0.005,
            n_radial=36, n_span=8, hub_segments=200,
        )
        _, occupancy = occupancy_of(roue, grid=140)
        result = topo.characteristic_radii(occupancy)

        self.assertEqual(result.machine_type, topo.MIXED)  # la branche visee
        self.assertGreater(result.r_2h, 0.9 * result.r_2s)  # le plateau pris pour un moyeu
        self.assertClose(result.b_2, 0.025, rel=0.15)
        self.assertGreater(result.area_2, 0.015)  # et non les 11 cm2 du cas degenere

    def test_disque_arriere_debordant(self):
        """Un disque arriere plus large que les aubes ne doit pas fixer r_2."""
        roue = synthetic.combine([
            synthetic.centrifugal_impeller(r1=0.035, r2=0.090),
            synthetic.revolve([(0.0, -0.017), (0.105, -0.017), (0.105, -0.011), (0.0, -0.011)], 180),
        ])
        _, occupancy = occupancy_of(roue)
        result = topo.characteristic_radii(occupancy)
        self.assertClose(result.r_tip, 0.105, rel=0.02)  # la matiere va bien jusque-la
        self.assertClose(result.r_blade_tip, 0.090, rel=config.VALID_GEOM_TOL)
        self.assertClose(result.r_2, result.r_blade_tip, rel=1e-12)
        self.assertEqual(result.machine_type, topo.CENTRIFUGAL)


class TestEndToEnd(BaseTestCase):
    def _run(self, mesh, **kwargs):
        path = self.path("roue.stl")
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(grid_nr=GRID, grid_nz=GRID, n_theta=360,
                                 symmetry_check=False, speeds=(1450.0,), **kwargs))

    def test_roue_a_l_envers_donne_le_meme_resultat(self):
        """Une rotation rigide du fichier ne doit rien changer au diagnostic."""
        roue = synthetic.centrifugal_impeller(front_shroud=True)
        droite = self._run(roue)
        envers = self._run(roue.transformed(matrix=FLIP))

        self.assertFalse(droite.suction.flipped)
        self.assertTrue(envers.suction.flipped)
        for name, first, second in (
            ("type", droite.topology.machine_type, envers.topology.machine_type),
            ("sortie", droite.discharge["composante_meridienne"],
             envers.discharge["composante_meridienne"]),
            ("rotation", droite.blades.rotation_sign, envers.blades.rotation_sign),
        ):
            with self.subTest(grandeur=name):
                self.assertEqual(first, second)
        self.assertClose(envers.topology.r_1s, droite.topology.r_1s, rel=0.02)
        self.assertClose(envers.topology.r_2, droite.topology.r_2, rel=0.02)

    def test_sortie_laterale_annoncee(self):
        """Une roue qui aspire au centre et refoule en peripherie sort radialement."""
        for label, mesh in (
            ("semi-ouverte", synthetic.centrifugal_impeller()),
            ("fermee", synthetic.centrifugal_impeller(front_shroud=True)),
            ("fermee a l'envers", synthetic.centrifugal_impeller(front_shroud=True).transformed(matrix=FLIP)),
        ):
            with self.subTest(roue=label):
                result = self._run(mesh)
                self.assertEqual(result.topology.machine_type, topo.CENTRIFUGAL)
                self.assertIn("radial", result.discharge["composante_meridienne"])
                self.assertIsNotNone(result.discharge["alpha2_deg"])

    def test_avertissement_dans_le_rapport(self):
        """Le retournement est signale a l'utilisateur, pas fait en douce."""
        result = self._run(synthetic.centrifugal_impeller().transformed(matrix=FLIP))
        self.assertTrue(any("a l'envers" in message for message in result.warnings))
        self.assertEqual(result.confidence["cote_aspiration"], HIGH)
        self.assertIn("cote_aspiration", result.to_dict())

    def test_option_aspiration_imposee(self):
        """--aspiration +z desactive la correction automatique."""
        roue = synthetic.centrifugal_impeller().transformed(matrix=FLIP)
        auto = self._run(roue)
        impose = self._run(roue, suction=axis_mod.SUCTION_PLUS_Z)
        self.assertTrue(auto.suction.flipped)
        self.assertFalse(impose.suction.flipped)
        self.assertEqual(auto.topology.machine_type, topo.CENTRIFUGAL)
        self.assertEqual(impose.topology.machine_type, topo.AXIAL)  # lue a l'envers, comme demande


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
