"""Aubes qui se referment sur elles-memes : ce que l'outil doit refuser de conclure.

Une helice toroidale n'a pas d'aube au sens de la SPEC : sa pale est une boucle,
sans bord d'attaque ni bord de fuite distincts. Toute la phase 4 -- cambrure,
beta1, beta2, sens de rotation -- et donc toute l'hydraulique qui en decoule
reposent sur l'hypothese inverse. L'outil doit reconnaitre la forme et le dire,
plutot que de sortir des chiffres qui ne decrivent rien.
"""

import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.confidence import LOW
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import blade_loops as loops
from impeller_analyzer.geometry import blade_normals as normals
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.geometry import topology as topo
from impeller_analyzer.io import writer

GRID = 110


def occupancy_of(mesh, grid=GRID, n_theta=240):
    aligned, _ = axis_mod.align_to_z(mesh)
    return occ_mod.build_occupancy(aligned, nr=grid, nz=grid, n_theta=n_theta)


class TestDetection(BaseTestCase):
    def test_helice_toroidale_reconnue(self):
        """Deux troncons en hauteur a azimut fixe : c'est une boucle."""
        occupancy = occupancy_of(synthetic.toroidal_propeller(n_blades=3))
        mask = topo.clean_blade_mask(occupancy.blade_mask())
        result = loops.detect_looped_blades(occupancy, 3, mask)

        self.assertTrue(result.looped)
        self.assertGreaterEqual(result.peak_fraction, config.LOOP_DOUBLE_FRACTION)
        # La bande dedoublee couvre la portee de l'aube, pas seulement le moyeu.
        self.assertGreater(result.r_outer, 2.0 * result.r_inner)
        self.assertTrue(result.warnings)

    def test_roues_a_aubes_simples_ne_sont_pas_signalees(self):
        """Aucune fausse alerte sur les trois familles ordinaires."""
        for label, mesh in (
            ("axiale", synthetic.axial_impeller(
                n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0)),
            ("centrifuge ouverte", synthetic.centrifugal_impeller()),
            ("centrifuge fermee", synthetic.centrifugal_impeller(front_shroud=True)),
        ):
            with self.subTest(roue=label):
                occupancy = occupancy_of(mesh)
                mask = topo.clean_blade_mask(occupancy.blade_mask())
                n_blades = topo.count_blades(occupancy).n_blades
                result = loops.detect_looped_blades(occupancy, n_blades, mask)
                self.assertFalse(result.looped)
                self.assertLess(result.peak_fraction, config.LOOP_DOUBLE_FRACTION)
                self.assertEqual(result.warnings, [])

    def test_le_nombre_de_boucles_est_bien_compte(self):
        """Le comptage d'aubes n'est pas trompe par les deux brins d'une boucle."""
        for n_blades in (3, 5):
            with self.subTest(pales=n_blades):
                occupancy = occupancy_of(synthetic.toroidal_propeller(n_blades=n_blades))
                self.assertEqual(topo.count_blades(occupancy).n_blades, n_blades)


class TestAnglesParNormales(BaseTestCase):
    """La methode de repli doit etre etalonnee sur des roues dont on sait les angles."""

    def _angles(self, mesh, grid=GRID):
        aligned, _ = axis_mod.align_to_z(mesh)
        occupancy = occ_mod.build_occupancy(aligned, nr=grid, nz=grid, n_theta=360)
        return normals.analyse(aligned, occupancy, topo.analyse(occupancy, None))

    def test_etalonnage_sur_angles_connus(self):
        """Sur des roues d'angles imposes, la lecture tient dans une fourchette etroite."""
        for beta1, beta2 in ((22.0, 25.0), (35.0, 50.0), (40.0, 65.0), (15.0, 20.0)):
            with self.subTest(beta=(beta1, beta2)):
                got = self._angles(
                    synthetic.centrifugal_impeller(beta1_deg=beta1, beta2_deg=beta2)
                )
                for measured, expected in ((got.beta1_deg, beta1), (got.beta2_deg, beta2)):
                    # Le biais est du cote bas et n'a jamais depasse 5 degres.
                    self.assertGreaterEqual(measured, expected - 5.0)
                    self.assertLessEqual(measured, expected + 1.0)

    def test_le_sens_d_enroulement_suit_la_roue(self):
        """Retourner le sens des aubes doit retourner le signe lu."""
        direct = self._angles(synthetic.centrifugal_impeller(sense=1))
        inverse = self._angles(synthetic.centrifugal_impeller(sense=-1))
        self.assertEqual(direct.slope_sign, -inverse.slope_sign)
        # La derive azimutale doit etre franche, pas un signe arrache au bruit.
        self.assertGreater(abs(direct.wrap_drift_deg), config.WRAP_SENSE_MIN_DEG)
        self.assertGreater(abs(inverse.wrap_drift_deg), config.WRAP_SENSE_MIN_DEG)

    def test_le_sens_ne_depend_pas_du_pas_de_grille(self):
        """Le sens d'enroulement doit etre le meme d'une grille a l'autre.

        Lu sur la composante radiale des normales il basculait : sur une aube
        helicoidale cette composante est minuscule et son signe n'est que du
        bruit. Lu sur la derive azimutale, il tient.
        """
        roue = synthetic.toroidal_propeller(n_blades=3)
        signes = {self._angles(roue, grid=g).slope_sign for g in (100, 130, 160)}
        self.assertEqual(len(signes), 1, f"sens instable selon la grille : {signes}")

    def test_insensible_au_pas_de_grille(self):
        """Deux grilles differentes doivent donner le meme angle."""
        roue = synthetic.centrifugal_impeller(beta1_deg=22.0, beta2_deg=25.0)
        fine = self._angles(roue, grid=150)
        grossiere = self._angles(roue, grid=100)
        self.assertClose(fine.beta2_deg, grossiere.beta2_deg, abs_=1.0)

    def test_roue_axiale_renvoyee_a_la_cambrure(self):
        """Sur une helice axiale la methode se declare hors sujet plutot que de repondre."""
        got = self._angles(
            synthetic.axial_impeller(n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0)
        )
        self.assertEqual(got.beta2_deg, 0.0)
        self.assertTrue(any("axiale" in note for note in got.notes))


class TestRoueToroidale(BaseTestCase):
    """Le generateur de roue de pompe a aubes en boucle, distinct de l'helice nue."""

    def test_topologie_de_boucle_reconnue(self):
        """Deux brins au milieu, un seul au bout : c'est ce qu'il doit produire.

        Le maillage passe par la reparation, comme dans le pipeline : les rubans
        vrilles du generateur sortent avec des normales incoherentes, et lue
        brute, la carte placait les brins sous le plateau. La boucle n'etait
        alors reconnue que grace aux franges du plateau, comptees a tort comme
        des brins.
        """
        from impeller_analyzer.io import repair

        roue = synthetic.toroidal_impeller(n_blades=3)
        repair.repair(roue)
        occupancy = occupancy_of(roue, grid=150, n_theta=300)
        mask = topo.clean_blade_mask(occupancy.blade_mask())
        result = loops.detect_looped_blades(occupancy, 3, mask)
        self.assertTrue(result.looped)
        self.assertGreaterEqual(result.peak_fraction, config.LOOP_DOUBLE_FRACTION)
        self.assertGreater(result.merge_radius, 0.0)

    def test_enveloppe_respectee(self):
        """Les rubans sont verticaux : ils ne debordent pas en rayon."""
        import math

        roue = synthetic.toroidal_impeller(r2=0.1670)
        rayon = max(math.hypot(v[0], v[1]) for v in roue.vertices)
        self.assertLessEqual(rayon, 0.1670 * 1.005)

    def test_redresser_les_aubes_raccourcit_le_canal(self):
        """Le gain d'une aube en boucle redressee passe par le frottement.

        A beta2 faible le canal s'allonge -- `L = (r2-r1)/sin(beta)` -- et la
        boucle presente deux fois ses faces sur toute cette longueur. C'est la
        que se joue son rendement, non dans la diffusion, qui reste ici proche
        de zero d'un bout a l'autre.
        """
        from impeller_analyzer.analysis import Options as O

        mesures = {}
        for beta2 in (8.0, 32.0):
            roue = synthetic.toroidal_impeller(beta1_deg=9.0, beta2_deg=beta2,
                                               b1=0.044, b2=0.042, thickness=0.014)
            path = self.path(f"m{beta2:.0f}.stl")
            writer.write_stl(roue, path, unit_factor=config.UNIT_FACTOR)
            result = run(path, O(unit="cm", speeds=(1450.0,), grid_nr=140, grid_nz=140,
                                 n_theta=300, symmetry_check=False, rotation=1,
                                 blades=3, beta1_deg=9.0, beta2_deg=beta2))
            mesures[beta2] = result.channel_losses
        couche, droite = mesures[8.0], mesures[32.0]
        self.assertGreater(couche.slenderness, 2.0 * droite.slenderness)
        self.assertGreater(couche.wetted_area, 2.0 * droite.wetted_area)
        self.assertGreater(droite.efficiency, couche.efficiency)


class TestConsequences(BaseTestCase):
    def _run(self, mesh, **kwargs):
        path = self.path("helice.stl")
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(grid_nr=GRID, grid_nz=GRID, n_theta=240,
                                 symmetry_check=False, speeds=(1450.0,), **kwargs))

    def test_les_angles_sont_repris_sur_les_normales(self):
        """La boucle bascule la lecture des angles, elle ne l'abandonne pas."""
        result = self._run(synthetic.toroidal_propeller(n_blades=3))
        self.assertIsNotNone(result.blade_loops)
        self.assertTrue(result.blade_loops.looped)
        self.assertTrue(any("boucle" in message for message in result.warnings))

        if result.topology.machine_type == topo.AXIAL:
            self.skipTest("roue axiale : la cambrure s'applique, les normales ne sont pas requises")
        self.assertIsNotNone(result.blade_normals)
        angles = result.blade_normals
        self.assertGreater(angles.beta1_deg, 0.0)
        self.assertGreater(angles.beta2_deg, 0.0)
        self.assertIn(angles.slope_sign, (-1, 1))
        # Les angles retenus sont bien ceux des normales, pas ceux de la cambrure.
        self.assertAlmostEqual(result.blades.beta1_deg, angles.beta1_deg, places=6)
        self.assertAlmostEqual(result.blades.beta2_deg, angles.beta2_deg, places=6)
        self.assertEqual(result.confidence.get_level("angles_de_pale"), angles.confidence)

    def test_les_angles_imposes_priment_sur_les_normales(self):
        """--beta1/--beta2 doit court-circuiter la lecture par les normales."""
        result = self._run(synthetic.toroidal_propeller(n_blades=3), beta1_deg=18.0, beta2_deg=27.0)
        self.assertAlmostEqual(result.blades.beta1_deg, 18.0, places=6)
        self.assertAlmostEqual(result.blades.beta2_deg, 27.0, places=6)

    def test_un_seul_angle_impose_laisse_l_autre_aux_normales(self):
        """Imposer beta2 ne doit pas faire retomber beta1 sur la cambrure.

        La cambrure ne s'applique pas a une aube en boucle : sur la roue reelle
        elle donnait 87 degres pour un beta1 que les normales lisent a 10, et le
        point de fonctionnement disparaissait purement et simplement.
        """
        roue = synthetic.toroidal_propeller(n_blades=3)
        libre = self._run(roue)
        if libre.blade_normals is None or libre.blade_normals.beta1_deg <= 0.0:
            self.skipTest("lecture par les normales indisponible sur cette roue")

        partiel = self._run(roue, beta2_deg=27.0)
        self.assertAlmostEqual(partiel.blades.beta2_deg, 27.0, places=6)
        self.assertAlmostEqual(partiel.blades.beta1_deg, libre.blades.beta1_deg, places=6)

    def test_sans_moyeu_le_rayon_de_sortie_est_celui_des_aubes(self):
        """La moyenne quadratique moyeu-carter n'a pas de sens sans moyeu.

        Elle degenere alors en r_2s / racine(2). Sur la roue toroidale de
        reference elle donnait 118 mm pour des aubes allant a 167, et la valeur
        basculait selon que le rapport r2/r1s tombait d'un cote ou de l'autre de
        la frontiere mixte / centrifuge, a un millieme pres -- soit un tiers sur
        u2 et deux tiers sur la hauteur, au gre du pas de grille.
        """
        import math

        roue = synthetic.centrifugal_impeller(
            front_shroud=True, flat_shroud=True, r1=0.055, r2=0.090, b1=0.014, b2=0.010
        )
        aligned, _ = axis_mod.align_to_z(roue)
        occupancy = occ_mod.build_occupancy(aligned, nr=120, nz=120, n_theta=240)
        result = topo.analyse(occupancy, None)

        self.assertEqual(result.machine_type, topo.MIXED)  # la branche visee
        self.assertEqual(result.r_2h, 0.0)  # et sans moyeu au plan de sortie
        self.assertClose(result.r_2, result.r_blade_tip, rel=1e-12)
        self.assertGreater(result.r_2, 2.0 * result.r_2s / math.sqrt(2.0))
        self.assertTrue(any("moyeu au plan de sortie" in w for w in result.warnings))

    def test_la_geometrie_reste_exploitable(self):
        """Axe, nombre d'aubes et rayons ne dependent pas de la forme des aubes."""
        result = self._run(synthetic.toroidal_propeller(n_blades=3, r_tip=0.100))
        self.assertEqual(result.topology.blades.n_blades, 3)
        self.assertClose(result.topology.r_tip, 0.100, rel=0.05)
        self.assertLess(result.axis.angle_to_z_deg, 1.0)

    def test_la_sensibilite_a_beta2_est_chiffree(self):
        """Quand un degre sur beta2 fait bouger la hauteur, il faut le dire, pas le taire."""
        from impeller_analyzer.hydraulics import meanline

        result = self._run(synthetic.toroidal_propeller(n_blades=3))
        self.assertGreaterEqual(result.head_sensitivity, 0.0)
        if result.head_sensitivity > config.BETA_SENSITIVITY_ALERT:
            self.assertTrue(any("hypersensible" in m for m in result.warnings))
        # Une roue ordinaire, aux aubes franchement inclinees, ne doit pas la declencher.
        sage = self._run(synthetic.centrifugal_impeller(beta1_deg=22.0, beta2_deg=25.0))
        self.assertLess(sage.head_sensitivity, config.BETA_SENSITIVITY_ALERT)
        self.assertFalse(any("hypersensible" in m for m in sage.warnings))
        del meanline

    def test_la_table_de_confiance_suit_le_texte(self):
        """Annoncer « ordres de grandeur » et coter « moyenne » serait se contredire."""
        result = self._run(synthetic.toroidal_propeller(n_blades=3))
        if result.head_sensitivity <= config.BETA_SENSITIVITY_ALERT:
            self.skipTest("hauteur peu sensible a beta2 sur cette roue")
        for quantity in ("hauteur", "puissance", "couple", "rendement"):
            with self.subTest(grandeur=quantity):
                self.assertEqual(result.confidence.get_level(quantity), LOW)
        # Debit et NPSHr ne passent pas par cu2 : ils ne sont pas degrades.
        for quantity in ("debit", "npshr"):
            with self.subTest(grandeur=quantity):
                self.assertNotEqual(result.confidence.get_level(quantity), LOW)

    def test_le_rapport_porte_la_mention_non_applicable(self):
        """Le tableau des performances doit etre desamorce dans le rapport."""
        from impeller_analyzer.io import report

        result = self._run(synthetic.toroidal_propeller(n_blades=3))
        chemin = report.write_markdown(result, self.workdir, source="helice.stl")
        with open(chemin, encoding="utf-8") as flux:
            texte = flux.read()
        self.assertIn("reserves", texte)
        self.assertIn("toroidal", texte)


if __name__ == "__main__":
    unittest.main()
