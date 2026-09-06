"""Phase 6 - cavitation : NPSHr, NPSHa et vitesse maximale (SPEC phase 6)."""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config
from impeller_analyzer.hydraulics import cavitation, meanline

#: Meme roue de reference que la phase 5.
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


def reference_curve(rpm=1450.0):
    """Courbe caracteristique de la roue de reference, NPSHr renseigne."""
    return cavitation.apply_to_curve(meanline.build_curve(REFERENCE, rpm))


class TestFluidState(BaseTestCase):
    def test_pression_atmospherique(self):
        """p_atm = 101325 (1 - 2.25577e-5 h)^5.2559."""
        self.assertClose(cavitation.atmospheric_pressure(0.0), config.P_ATM_SEA_LEVEL, rel=1e-12)
        expected = config.P_ATM_SEA_LEVEL * (1.0 - config.ATM_LAPSE_COEF * 1500.0) ** config.ATM_EXPONENT
        self.assertClose(cavitation.atmospheric_pressure(1500.0), expected, rel=1e-12)
        self.assertLess(cavitation.atmospheric_pressure(2000.0), cavitation.atmospheric_pressure(0.0))

    def test_pression_de_vapeur_par_antoine(self):
        """La correlation d'Antoine redonne la valeur de reference a 20 degres C."""
        got = cavitation.vapour_pressure(20.0)
        self.assertClose(got, config.P_VAP, rel=0.01)
        expected = 10.0 ** (config.ANTOINE_A - config.ANTOINE_B / (config.ANTOINE_C + 20.0))
        self.assertClose(got, expected * config.ANTOINE_MMHG_TO_PA, rel=1e-12)
        # A 100 degres C la vapeur atteint la pression atmospherique.
        self.assertClose(cavitation.vapour_pressure(100.0), config.P_ATM_SEA_LEVEL, rel=0.01)

    def test_masse_volumique_par_table(self):
        """La table de masse volumique est interpolee lineairement."""
        self.assertClose(cavitation.water_density(20.0), config.RHO, rel=1e-12)
        self.assertClose(cavitation.water_density(22.5), 0.5 * (998.20 + 997.05), rel=1e-12)
        self.assertLess(cavitation.water_density(80.0), cavitation.water_density(20.0))

    def test_npsha_par_defaut(self):
        """Valeur imposee par la SPEC : NPSHa = 9.61 m dans les conditions par defaut."""
        site = cavitation.installation()
        self.assertClose(site.npsha, 9.61, abs_=0.01)
        expected = (site.p_atm - site.p_vap) / (site.rho * config.G) + config.HAUTEUR_ASPIRATION - config.PERTES_ASPIRATION
        self.assertClose(site.npsha, expected, rel=1e-12)

    def test_npsha_depend_de_l_installation(self):
        """Altitude, temperature, hauteur et pertes agissent dans le bon sens."""
        base = cavitation.installation().npsha
        self.assertLess(cavitation.installation(altitude=1500.0).npsha, base)
        self.assertLess(cavitation.installation(temperature_c=70.0).npsha, base)
        self.assertGreater(cavitation.installation(suction_height=3.0).npsha, base)
        self.assertLess(cavitation.installation(suction_losses=2.0).npsha, base)
        self.assertIn("en aspiration", cavitation.installation(suction_height=-2.0).describe())


class TestNPSHr(BaseTestCase):
    def test_methode_cinematique(self):
        """NPSHr = lc cm1^2/(2g) + lw w1s^2/(2g)."""
        got = cavitation.npshr_kinematic(3.0, 12.0)
        expected = (config.LAMBDA_C * 9.0 + config.LAMBDA_W * 144.0) / (2.0 * config.G)
        self.assertClose(got, expected, rel=1e-12)

    def test_methode_vitesse_specifique_d_aspiration(self):
        """NPSHr = (n sqrt(Q) / n_ss)^(4/3)."""
        got = cavitation.npshr_suction_speed(1450.0, 0.020)
        expected = (1450.0 * math.sqrt(0.020) / config.N_SS) ** config.NSS_EXPONENT
        self.assertClose(got, expected, rel=1e-12)
        self.assertClose(cavitation.npshr_suction_speed(1450.0, 0.0), 0.0, abs_=1e-15)

    def test_valeur_retenue_est_la_plus_penalisante(self):
        """Les deux methodes sont reportees, le maximum est retenu."""
        curve = reference_curve()
        for point in curve.points:
            self.assertClose(
                point.npshr, max(point.npshr_kinematic, point.npshr_suction_speed), rel=1e-12
            )
            self.assertGreaterEqual(point.npshr, point.npshr_kinematic - 1e-12)
            self.assertGreaterEqual(point.npshr, point.npshr_suction_speed - 1e-12)

    def test_npshr_croit_avec_la_vitesse_comme_n_carre(self):
        """A coefficient de debit constant, NPSHr varie comme n^2."""
        slow = reference_curve(1000.0).nominal_point()
        fast = reference_curve(2000.0).nominal_point()
        self.assertClose(fast.npshr, 4.0 * slow.npshr, rel=1e-9)


class TestMaximumSpeed(BaseTestCase):
    def test_limite_npsh(self):
        """n_max = n_ref sqrt(NPSHa / (marge NPSHr)), arrondie a la dizaine inferieure."""
        curve = reference_curve(1450.0)
        site = cavitation.installation()
        limit = cavitation.maximum_speed(curve, site)
        reference = curve.nominal_point()
        expected = 1450.0 * math.sqrt(site.npsha / (config.MARGE_NPSH * reference.npshr))
        self.assertClose(limit.rpm_max_npsh, expected, rel=1e-12)
        self.assertEqual(limit.active_limit, cavitation.NPSH_LIMIT)
        self.assertEqual(limit.rpm_max, int(math.floor(expected / config.RPM_ROUNDING) * config.RPM_ROUNDING))
        self.assertEqual(limit.rpm_max % 10, 0)

    def test_limite_vitesse_relative(self):
        """Quand w1s est la contrainte active, c'est elle qui fixe la vitesse."""
        # Un grand rayon d'oeillard fait exploser w1s bien avant le NPSH.
        wide = meanline.MeanlineInput(
            r_1=math.sqrt((0.30 ** 2 + 0.10 ** 2) / 2.0), r_1s=0.30, r_1h=0.10, r_2=0.32,
            area_1=math.pi * (0.30 ** 2 - 0.10 ** 2) * config.TAU_1,
            area_2=2.0 * math.pi * 0.32 * 0.08 * config.TAU_2,
            beta1_deg=30.0, beta2_deg=50.0, n_blades=5,
        )
        curve = cavitation.apply_to_curve(meanline.build_curve(wide, 1000.0))
        limit = cavitation.maximum_speed(curve, cavitation.installation(suction_height=40.0))
        expected = 1000.0 * config.W1S_MAX / curve.nominal_point().w1s
        self.assertClose(limit.rpm_max_w1s, expected, rel=1e-12)
        self.assertEqual(limit.active_limit, cavitation.SPEED_LIMIT)
        self.assertLessEqual(limit.rpm_max, limit.rpm_max_w1s)

    def test_resultat_independant_du_regime_de_reference(self):
        """La vitesse maximale ne depend pas du regime depuis lequel on l'extrapole."""
        site = cavitation.installation()
        values = {
            cavitation.maximum_speed(reference_curve(rpm), site).rpm_max
            for rpm in (1000.0, 1450.0, 3000.0)
        }
        self.assertEqual(len(values), 1)

    def test_installation_impossible(self):
        """Un NPSH disponible negatif est signale, sans exception."""
        site = cavitation.installation(suction_height=-12.0)
        self.assertLess(site.npsha, 0.0)
        limit = cavitation.maximum_speed(reference_curve(), site)
        self.assertEqual(limit.rpm_max, 0)
        self.assertTrue(any("NPSH disponible negatif" in w for w in limit.warnings))

    def test_avertissement_si_le_regime_analyse_cavite(self):
        """Un regime superieur a la limite produit un avertissement explicite."""
        site = cavitation.installation()
        limit = cavitation.maximum_speed(reference_curve(6000.0), site)
        self.assertLess(limit.rpm_max, 6000.0)
        self.assertTrue(any("cavite" in w for w in limit.warnings))

    def test_courbe_sans_point_nominal(self):
        """Sans point nominal la vitesse maximale n'est pas calculable."""
        empty = meanline.build_curve(meanline.MeanlineInput(), 1450.0)
        limit = cavitation.maximum_speed(empty, cavitation.installation())
        self.assertEqual(limit.rpm_max, 0)
        self.assertEqual(limit.confidence, "low")

    def test_hypotheses_reportees(self):
        """Les hypotheses d'installation accompagnent la vitesse maximale."""
        site = cavitation.installation(altitude=800.0, temperature_c=35.0, suction_height=1.5)
        limit = cavitation.maximum_speed(reference_curve(), site)
        for fragment in ("800", "35", "+1.5", "en charge"):
            self.assertIn(fragment, limit.assumptions)
        self.assertIn("NPSHa_m", site.to_dict())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
