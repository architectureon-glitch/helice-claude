"""Phase 4 - coupes, ligne de cambrure, angles de pale et sens de rotation (SPEC phase 4)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.confidence import HIGH, LOW
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import blade_angles as ba
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.geometry import sections as sec
from impeller_analyzer.geometry import topology as topo

_CACHE: dict = {}


def analysed(key, builder, nr=150, nz=150):
    """Chaine complete phases 2 a 4, mise en cache pour la suite de tests."""
    if key not in _CACHE:
        mesh, _ = axis_mod.align_to_z(builder())
        occupancy = occ_mod.build_occupancy(mesh, nr=nr, nz=nz)
        topology = topo.analyse(occupancy, None)
        sections = sec.extract_sections(mesh, occupancy, topology)
        geometry = ba.analyse(sections, topology)
        _CACHE[key] = (mesh, occupancy, topology, sections, geometry)
    return _CACHE[key]


class TestSections(BaseTestCase):
    def test_une_coupe_par_pale(self):
        """Chaque surface de coupe rend autant de profils fermes que de pales."""
        _, _, topology, sections, _ = analysed(
            ("axial", 4), lambda: synthetic.axial_impeller(n_blades=4)
        )
        self.assertEqual(len(sections), config.N_SECTIONS)
        for section in sections:
            with self.subTest(span=section.span):
                self.assertEqual(len(section.usable()), 4)

    def test_coupes_cylindriques_pour_une_roue_axiale(self):
        """Roue axiale : les surfaces de courant sont des cylindres r = constante."""
        _, occupancy, topology, sections, _ = analysed(
            ("axial", 4), lambda: synthetic.axial_impeller(n_blades=4)
        )
        radii = [section.reference_radius for section in sections]
        self.assertTrue(all(section.curve.vertical for section in sections))
        self.assertLess(radii[0], radii[-1])
        self.assertGreater(radii[0], topology.r_1h)
        self.assertLess(radii[-1], topology.r_tip)

    def test_surfaces_de_courant_pour_une_roue_centrifuge(self):
        """Roue centrifuge : la coupe suit la veine, du moyeu au carter."""
        _, _, topology, sections, _ = analysed(("centri", 6), synthetic.centrifugal_impeller)
        self.assertEqual(topology.machine_type, topo.CENTRIFUGAL)
        for section in sections:
            self.assertFalse(section.curve.vertical)
        # La surface au moyeu est plus basse que celle au carter : a rayon egal,
        # le champ z - z_courbe y est donc plus grand.
        low = sections[0].curve
        high = sections[-1].curve
        middle = 0.5 * (low.radii[0] + low.radii[-1])
        self.assertGreater(low.level(middle, 0.0), high.level(middle, 0.0))

    def test_profil_deroule_est_ferme_et_coherent(self):
        """Les points d'un profil sont sur la surface de coupe et forment une boucle."""
        _, _, _, sections, _ = analysed(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        section = sections[5]
        profile = section.usable()[0]
        for radius in profile.radius:
            self.assertClose(radius, section.reference_radius, rel=1e-3)
        self.assertGreater(len(profile), config.MIN_PROFILE_POINTS)
        self.assertClose(profile.reference_wrap(), math.radians(60.0), abs_=math.radians(1.0))


class TestCamber(BaseTestCase):
    def test_cambrure_d_un_profil_droit(self):
        """Sur un profil a cambrure nulle, la ligne de cambrure est la corde."""
        # Bande rectangulaire de 100 x 4 mm dans le plan (tangentiel, meridien).
        points = []
        for x in range(0, 101, 2):
            points.append((x * 1e-3, 0.002))
        for x in range(100, -1, -2):
            points.append((x * 1e-3, -0.002))
        profile = sec.Profile(
            tangential=[p[0] for p in points],
            meridional=[p[1] for p in points],
            radius=[0.05] * len(points),
            reference=0.05,
        )
        camber = ba.camber_line(profile, leading_edge_at_max=True)
        self.assertTrue(camber.valid)
        self.assertClose(camber.max_thickness, 0.004, rel=0.05)
        # La corde joint deux coins opposes de la bande, elle est donc legerement
        # inclinee ; la cambrure, elle, doit rester la mediane horizontale, donc
        # de pente nulle dans le plan (tangentiel, meridien).
        for slope in camber.slopes:
            self.assertClose(slope, 0.0, abs_=1e-3)
        # Elle s'ecarte donc lineairement de la corde, a la pente de l'inclinaison
        # de celle-ci : 4 mm d'epaisseur pour 100 mm de corde.
        drift = (camber.offsets[-1] - camber.offsets[0]) / (camber.stations[-1] - camber.stations[0])
        self.assertClose(abs(drift), 0.04, rel=0.05)

    def test_derivees_locales_exactes_sur_une_parabole(self):
        """L'ajustement parabolique local est exact, y compris aux extremites."""
        x = [0.0, 1.0, 2.0, 3.5, 5.0]
        y = [value * value for value in x]
        got = ba._local_derivatives(x, y)
        self.assertSeqClose(got, [2.0 * value for value in x], abs_=1e-9)

    def test_ecretage_des_bouts(self):
        """Les stations dont l'epaisseur s'ouvre brutalement sont ecartees."""
        stations = [float(i) for i in range(20)]
        thickness = [0.1, 0.6] + [1.0] * 16 + [0.6, 0.1]
        kept = ba._trim_ends(stations, thickness, chord=20.0, max_thickness=1.0)
        self.assertNotIn(0, kept)
        self.assertNotIn(19, kept)
        self.assertIn(10, kept)


class TestBladeAngles(BaseTestCase):
    def test_helice_a_pales_planes_calees_a_30(self):
        """Test impose par la SPEC : 4 pales, beta = 30 deg a tous les rayons."""
        _, _, topology, _, geometry = analysed(
            ("axial", 4), lambda: synthetic.axial_impeller(n_blades=4, beta_deg=30.0)
        )
        self.assertEqual(topology.blades.n_blades, 4)
        self.assertClose(topology.r_1h, 0.030, rel=config.VALID_GEOM_TOL)
        self.assertClose(topology.r_tip, 0.100, rel=config.VALID_GEOM_TOL)
        self.assertClose(geometry.beta1_deg, 30.0, abs_=1.0)
        self.assertClose(geometry.beta2_deg, 30.0, abs_=1.0)
        for section in geometry.sections:
            with self.subTest(radius=section.radius):
                self.assertClose(section.beta1_deg, 30.0, abs_=1.0)
                self.assertClose(section.beta2_deg, 30.0, abs_=1.0)

    def test_autres_calages_axiaux(self):
        """L'extraction reste exacte de 15 a 45 degres."""
        for beta in (15.0, 45.0):
            with self.subTest(beta=beta):
                _, _, _, _, geometry = analysed(
                    ("axial-beta", beta),
                    lambda beta=beta: synthetic.axial_impeller(n_blades=5, beta_deg=beta),
                )
                self.assertClose(geometry.beta1_deg, beta, abs_=config.VALID_BETA_DEG)
                self.assertClose(geometry.beta2_deg, beta, abs_=config.VALID_BETA_DEG)

    def test_aubes_centrifuges_en_arc(self):
        """Roue centrifuge : beta1 et beta2 retrouves a 2 degres pres."""
        for beta1, beta2 in ((22.0, 25.0), (25.0, 25.0), (30.0, 20.0)):
            with self.subTest(beta1=beta1, beta2=beta2):
                _, _, _, _, geometry = analysed(
                    ("centri-beta", beta1, beta2),
                    lambda a=beta1, b=beta2: synthetic.centrifugal_impeller(beta1_deg=a, beta2_deg=b),
                )
                self.assertClose(geometry.beta1_deg, beta1, abs_=config.VALID_BETA_DEG)
                self.assertClose(geometry.beta2_deg, beta2, abs_=config.VALID_BETA_DEG)

    def test_corde_epaisseur_et_pas(self):
        """Corde, epaisseur maximale et pas helicoidal sont coherents."""
        _, _, topology, _, geometry = analysed(
            ("axial", 4), lambda: synthetic.axial_impeller(n_blades=4)
        )
        self.assertGreater(geometry.chord, 0.0)
        self.assertClose(geometry.max_thickness, 0.004, rel=0.25)
        expected = 2.0 * math.pi * geometry.reference_radius * math.tan(
            math.radians(0.5 * (geometry.beta1_deg + geometry.beta2_deg))
        )
        self.assertClose(geometry.pitch, expected, rel=1e-9)
        for section in geometry.sections:
            self.assertGreater(section.chord_to_thickness, config.CHORD_THICKNESS_MIN)


class TestRotationSense(BaseTestCase):
    def test_sens_axial_suit_le_signe_de_k(self):
        """Helice axiale : signe(omega) = signe(k), et s'inverse avec la geometrie."""
        _, _, _, _, direct = analysed(
            ("axial-sens", 1), lambda: synthetic.axial_impeller(n_blades=4, sense=1)
        )
        _, _, _, _, reverse = analysed(
            ("axial-sens", -1), lambda: synthetic.axial_impeller(n_blades=4, sense=-1)
        )
        self.assertEqual(direct.rotation_sign, 1)
        self.assertIn(ba.COUNTERCLOCKWISE, direct.rotation_label)
        self.assertEqual(reverse.rotation_sign, -1)
        self.assertIn(ba.CLOCKWISE, reverse.rotation_label)
        self.assertEqual(direct.confidence["sens_de_rotation"], HIGH)

    def test_sens_centrifuge_est_oppose_a_la_courbure(self):
        """Aube incurvee vers l'arriere : signe(omega) = -signe(dtheta/dr)."""
        _, _, _, _, direct = analysed(("centri-sens", 1), lambda: synthetic.centrifugal_impeller(sense=1))
        _, _, _, _, reverse = analysed(("centri-sens", -1), lambda: synthetic.centrifugal_impeller(sense=-1))
        self.assertEqual(direct.rotation_sign, -1)
        self.assertEqual(reverse.rotation_sign, 1)
        self.assertTrue(any("incurvees vers l'avant" in note for note in direct.notes))

    def test_regle_de_sens_isolee(self):
        """La regle de signe est directement verifiable, hors geometrie."""
        self.assertEqual(ba.rotation_sense(topo.AXIAL, 1), 1)
        self.assertEqual(ba.rotation_sense(topo.MIXED, -1), -1)
        self.assertEqual(ba.rotation_sense(topo.CENTRIFUGAL, 1), -1)
        self.assertEqual(ba.rotation_sense(topo.CENTRIFUGAL, -1), 1)
        self.assertEqual(ba.rotation_sense(topo.AXIAL, 0), 0)

    def test_sens_de_sortie_du_liquide(self):
        """Le sens de sortie reporte la composante meridienne et l'angle absolu."""
        _, _, topology, _, geometry = analysed(("centri", 6), synthetic.centrifugal_impeller)
        discharge = ba.discharge_direction(topology, geometry, cm2=3.0, cu2=12.0)
        self.assertIn("radial", discharge["composante_meridienne"])
        self.assertClose(discharge["alpha2_deg"], math.degrees(math.atan2(3.0, 12.0)), rel=1e-12)
        self.assertIn(geometry.rotation_label, discharge["composante_tangentielle"])

        _, _, axial_topology, _, axial_geometry = analysed(
            ("axial", 4), lambda: synthetic.axial_impeller(n_blades=4)
        )
        axial = ba.discharge_direction(axial_topology, axial_geometry)
        self.assertIn("-Z", axial["composante_meridienne"])
        self.assertIsNone(axial["alpha2_deg"])


class TestManualFallback(BaseTestCase):
    def test_angles_imposes_par_l_utilisateur(self):
        """--beta1 / --beta2 remplacent l'extraction et remontent la confiance."""
        _, _, topology, sections, _ = analysed(("axial", 4), lambda: synthetic.axial_impeller(n_blades=4))
        forced = ba.analyse(sections, topology, forced_beta1_deg=18.0, forced_beta2_deg=27.0)
        self.assertClose(forced.beta1_deg, 18.0, rel=1e-12)
        self.assertClose(forced.beta2_deg, 27.0, rel=1e-12)
        self.assertTrue(forced.forced_beta)
        self.assertEqual(forced.confidence["angles_de_pale"], HIGH)

    def test_absence_de_coupe_exploitable(self):
        """Sans coupe exploitable, l'outil le dit et invite au repli manuel."""
        geometry = ba.analyse([], topo.Topology())
        self.assertEqual(geometry.confidence["angles_de_pale"], LOW)
        self.assertTrue(any("--beta1" in w for w in geometry.warnings))
        forced = ba.analyse([], topo.Topology(), forced_beta1_deg=20.0, forced_beta2_deg=24.0)
        self.assertClose(forced.beta1_deg, 20.0, rel=1e-12)
        self.assertEqual(forced.confidence["angles_de_pale"], HIGH)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
