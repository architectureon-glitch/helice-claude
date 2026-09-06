"""Phase 5 - modele hydraulique de ligne moyenne (SPEC phase 5)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.confidence import HIGH, LOW, MEDIUM, ConfidenceMap
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import blade_angles as ba
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.geometry import sections as sec
from impeller_analyzer.geometry import topology as topo
from impeller_analyzer.hydraulics import meanline, similarity

#: Roue de reference des tests hydrauliques, valeurs geometriques imposees.
REFERENCE = meanline.MeanlineInput(
    r_1=0.0250,
    r_1s=0.0350,
    r_1h=0.0000,
    r_2=0.0900,
    area_1=math.pi * 0.0350 ** 2 * config.TAU_1,
    area_2=2.0 * math.pi * 0.0900 * 0.0100 * config.TAU_2,
    beta1_deg=22.0,
    beta2_deg=25.0,
    n_blades=6,
)

_CACHE: dict = {}


def analysed(key, builder):
    """Chaine phases 2 a 5 sur une geometrie synthetique, mise en cache."""
    if key not in _CACHE:
        mesh, _ = axis_mod.align_to_z(builder())
        occupancy = occ_mod.build_occupancy(mesh, nr=150, nz=150)
        topology = topo.analyse(occupancy, None)
        geometry = ba.analyse(sec.extract_sections(mesh, occupancy, topology), topology)
        data = meanline.MeanlineInput.from_geometry(topology, geometry)
        _CACHE[key] = (topology, geometry, data, meanline.curves_for_speeds(data))
    return _CACHE[key]


class TestBuildingBlocks(BaseTestCase):
    def test_glissement_de_wiesner(self):
        """sigma = 1 - sqrt(sin b2) / N^0.7."""
        for beta, blades in ((25.0, 6), (35.0, 4), (12.0, 9)):
            with self.subTest(beta=beta, blades=blades):
                expected = 1.0 - math.sqrt(math.sin(math.radians(beta))) / blades ** 0.7
                self.assertClose(meanline.slip_factor(beta, blades), expected, rel=1e-12)
        self.assertClose(meanline.slip_factor(25.0, 0), 1.0, rel=1e-12)

    def test_vitesse_de_rotation(self):
        """omega = 2 pi n / 60."""
        self.assertClose(meanline.angular_velocity(3000.0), 2.0 * math.pi * 3000.0 / 60.0, rel=1e-12)

    def test_point_nominal_a_incidence_nulle(self):
        """u1 = omega r1, cm1 = u1 tan(b1), Q_n = cm1 A1 (SPEC 5.1)."""
        omega = meanline.angular_velocity(1450.0)
        flow, u1, cm1 = meanline.nominal_flow(REFERENCE, omega)
        self.assertClose(u1, omega * REFERENCE.r_1, rel=1e-12)
        self.assertClose(cm1, u1 * math.tan(math.radians(REFERENCE.beta1_deg)), rel=1e-12)
        self.assertClose(flow, cm1 * REFERENCE.area_1, rel=1e-12)

    def test_hauteur_d_euler_avec_glissement(self):
        """cu2 = sigma u2 - cm2 / tan(b2), H_th = u2 cu2 / g (SPEC 5.2)."""
        omega = meanline.angular_velocity(1450.0)
        flow, _, _ = meanline.nominal_flow(REFERENCE, omega)
        slip = meanline.slip_factor(REFERENCE.beta2_deg, REFERENCE.n_blades)
        head, cm2, cu2 = meanline.euler_head(REFERENCE, omega, flow, slip)
        u2 = omega * REFERENCE.r_2
        self.assertClose(cm2, flow / REFERENCE.area_2, rel=1e-12)
        self.assertClose(cu2, slip * u2 - cm2 / math.tan(math.radians(REFERENCE.beta2_deg)), rel=1e-12)
        self.assertClose(head, u2 * cu2 / config.G, rel=1e-12)


class TestCurve(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.curve = meanline.build_curve(REFERENCE, 1450.0)

    def test_balayage_conforme(self):
        """29 points de Q = 0 a 1.40 Q_n, le nominal tombant sur l'un d'eux."""
        self.assertEqual(len(self.curve.points), config.Q_SWEEP_POINTS)
        flows = [point.flow for point in self.curve.points]
        self.assertClose(flows[0], 0.0, abs_=1e-6 * self.curve.flow_nominal)
        self.assertClose(flows[-1], config.Q_SWEEP_MAX * self.curve.flow_nominal, rel=1e-12)
        self.assertClose(
            self.curve.points[self.curve.nominal_index].flow, self.curve.flow_nominal, rel=1e-9
        )

    def test_perte_de_frottement_calee_au_nominal(self):
        """k_f donne 6 % de H_th au point nominal (SPEC 5.3)."""
        nominal = self.curve.nominal_point()
        self.assertClose(
            nominal.loss_friction,
            config.K_FROTTEMENT_REL * nominal.head_theoretical,
            rel=1e-9,
        )
        self.assertClose(
            self.curve.friction_coefficient,
            config.K_FROTTEMENT_REL * nominal.head_theoretical / nominal.flow ** 2,
            rel=1e-9,
        )

    def test_perte_d_incidence_nulle_au_nominal(self):
        """L'incidence est nulle au point nominal, positive ailleurs."""
        nominal = self.curve.nominal_point()
        self.assertClose(nominal.loss_incidence, 0.0, abs_=1e-12)
        for point in self.curve.points:
            self.assertGreaterEqual(point.loss_incidence, 0.0)
        self.assertGreater(self.curve.points[0].loss_incidence, 0.0)
        self.assertGreater(self.curve.points[-1].loss_incidence, 0.0)

    def test_hauteur_est_bien_h_th_moins_les_pertes(self):
        """H(Q) = H_th(Q) - h_frottement - h_incidence."""
        for point in self.curve.points:
            self.assertClose(
                point.head,
                point.head_theoretical - point.loss_friction - point.loss_incidence,
                rel=1e-12,
            )

    def test_courbe_decroissante_pour_des_aubes_arrieres(self):
        """Aubes incurvees vers l'arriere : H_th decroit avec le debit."""
        heads = [point.head_theoretical for point in self.curve.points]
        for previous, following in zip(heads, heads[1:]):
            self.assertLessEqual(following, previous + 1e-12)

    def test_puissances_et_couple(self):
        """P_hydraulique = rho g Q H, couple = P_arbre / omega (SPEC 5.4)."""
        for point in self.curve.points:
            self.assertClose(
                point.hydraulic_power,
                config.RHO * config.G * point.flow * max(point.head, 0.0),
                rel=1e-12,
            )
            self.assertClose(point.torque, point.shaft_power / self.curve.omega, rel=1e-12)
            self.assertGreaterEqual(point.shaft_power, point.hydraulic_power - 1e-9)

    def test_bep_et_rendement_de_reference(self):
        """Le BEP est le maximum de rendement, egal a eta_h eta_vol eta_mec."""
        best = self.curve.best_efficiency_point()
        target = config.ETA_H * config.ETA_VOL * config.ETA_MEC
        self.assertClose(best.efficiency, target, rel=1e-9)
        for point in self.curve.points:
            self.assertLessEqual(point.efficiency, best.efficiency + 1e-12)
        self.assertGreater(best.flow, 0.0)
        self.assertGreater(best.head, 0.0)

    def test_vitesse_specifique_au_bep(self):
        """n_q est calculee au point de meilleur rendement."""
        best = self.curve.best_efficiency_point()
        self.assertClose(
            self.curve.specific_speed,
            1450.0 * math.sqrt(best.flow) / best.head ** 0.75,
            rel=1e-12,
        )

    def test_geometrie_inexploitable(self):
        """Une geometrie incomplete produit un avertissement, pas une exception."""
        empty = meanline.build_curve(meanline.MeanlineInput(), 1450.0)
        self.assertEqual(empty.points, [])
        self.assertTrue(any("geometrie insuffisante" in w for w in empty.warnings))

        # Angles tels que la roue ne peut pas fonctionner en pompe.
        inert = meanline.MeanlineInput(
            r_1=0.05, r_1s=0.05, r_1h=0.0, r_2=0.05,
            area_1=0.01, area_2=0.01, beta1_deg=30.0, beta2_deg=30.0, n_blades=4,
        )
        curve = meanline.build_curve(inert, 1450.0)
        self.assertIsNone(curve.nominal_point())
        self.assertTrue(any("Euler negative" in w for w in curve.warnings))


class TestSimilarity(BaseTestCase):
    def test_lois_de_similitude(self):
        """Q ~ n, H ~ n^2, P ~ n^3 : ecart inferieur a 2 % (SPEC 5.5)."""
        curves = meanline.curves_for_speeds(REFERENCE, (1000.0, 2000.0, 3000.0))
        check = similarity.check(curves[0], curves[2])
        self.assertTrue(check.passed)
        self.assertLess(check.worst_deviation, config.SIMILARITY_TOL)
        self.assertClose(check.flow_scaled, check.flow_direct, rel=1e-9)
        self.assertClose(check.head_scaled, check.head_direct, rel=1e-9)
        self.assertClose(check.power_scaled, check.power_direct, rel=1e-9)

    def test_extrapolations_elementaires(self):
        """Les trois lois d'echelle sont appliquees telles quelles."""
        self.assertClose(similarity.scale_flow(2.0, 1000.0, 3000.0), 6.0, rel=1e-12)
        self.assertClose(similarity.scale_head(2.0, 1000.0, 3000.0), 18.0, rel=1e-12)
        self.assertClose(similarity.scale_power(2.0, 1000.0, 3000.0), 54.0, rel=1e-12)

    def test_ecart_detecte(self):
        """Un ecart au-dela du seuil est signale."""
        curves = meanline.curves_for_speeds(REFERENCE, (1000.0, 3000.0))
        curves[1].points[curves[1].nominal_index].head *= 1.10
        check = similarity.check(curves[0], curves[1])
        self.assertFalse(check.passed)
        self.assertTrue(any("erreur de calcul" in w for w in check.warnings))


class TestOnSyntheticImpellers(BaseTestCase):
    def test_roue_centrifuge_complete(self):
        """Chaine complete sur une roue centrifuge : grandeurs plausibles et coherentes."""
        topology, _, data, curves = analysed(("centri",), synthetic.centrifugal_impeller)
        self.assertTrue(data.valid())
        for curve in curves:
            with self.subTest(rpm=curve.rpm):
                best = curve.best_efficiency_point()
                self.assertGreater(best.head, 0.0)
                self.assertGreater(best.flow, 0.0)
                self.assertClose(
                    best.efficiency, config.ETA_H * config.ETA_VOL * config.ETA_MEC, rel=1e-9
                )
        # Controle croise du type de roue par la vitesse specifique (SPEC 3.3).
        self.assertEqual(
            topo.type_from_specific_speed(curves[0].specific_speed), topology.machine_type
        )
        self.assertIsNone(topo.cross_check_type(topology, curves[0].specific_speed))

    def test_helice_axiale_complete(self):
        """Une helice dont l'angle varie du bord d'attaque au bord de fuite fournit une hauteur."""
        topology, _, _, curves = analysed(
            ("axial",),
            lambda: synthetic.axial_impeller(
                n_blades=4, beta_deg=20.0, beta2_deg=35.0, blade_wrap_deg=70.0
            ),
        )
        self.assertEqual(topology.machine_type, topo.AXIAL)
        best = curves[0].best_efficiency_point()
        self.assertGreater(best.head, 0.0)
        self.assertEqual(
            topo.type_from_specific_speed(curves[0].specific_speed), topo.AXIAL
        )

    def test_similitude_sur_geometrie_reelle(self):
        """La similitude tient aussi sur la geometrie extraite d'un maillage."""
        _, _, _, curves = analysed(("centri",), synthetic.centrifugal_impeller)
        self.assertTrue(similarity.check(curves[0], curves[2]).passed)

    def test_propagation_de_la_confiance(self):
        """La confiance hydraulique est la plus faible des confiances geometriques."""
        topology, geometry, _, _ = analysed(("centri",), synthetic.centrifugal_impeller)
        levels = meanline.confidence_of(topology, geometry)
        self.assertEqual(set(levels), {"debit", "hauteur", "puissance", "couple", "rendement"})
        self.assertEqual(levels.overall(), levels["hauteur"])

        degraded = topo.Topology()
        degraded.confidence = ConfidenceMap({"rayons": LOW, "sections": HIGH, "nombre_de_pales": HIGH})
        blade = ba.BladeGeometry()
        blade.confidence = ConfidenceMap({"angles_de_pale": MEDIUM})
        self.assertEqual(meanline.confidence_of(degraded, blade).overall(), LOW)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
