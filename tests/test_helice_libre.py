"""Helice libre : poussee, rendement propulsif, et le plafond qui les juge.

Le modele est un bilan par element de pale et quantite de mouvement, dont la
geometrie -- calage, corde, cambrure, epaisseur -- est lue sur la piece et non
demandee a l'utilisateur.  Ce qui se verifie ici n'est donc pas un chiffre
absolu, qui demanderait une courbe d'essai, mais les **invariants** que tout
modele d'helice doit respecter : la forme de la courbe de poussee, l'invariance
d'echelle des coefficients, et surtout le plafond du disque actif ideal, qu'un
rendement calcule ne peut pas franchir.
"""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.hydraulics import propulsion
from impeller_analyzer.io import writer

GRID = 120
N_THETA = 240
RPM = 1450.0


class PropulsionTestCase(BaseTestCase):
    """Lit une helice de synthese et l'analyse en helice libre."""

    def lire(self, mesh, name="helice.stl", **options):
        """Analyse geometrique complete d'un maillage de synthese."""
        path = self.path(name)
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(
            unit="cm", speeds=(RPM,), grid_nr=GRID, grid_nz=GRID, n_theta=N_THETA,
            symmetry_check=False, rotation=1, **options
        ))

    def helice(self, **kwargs):
        """Une helice axiale fine, du domaine ou le modele a un sens."""
        parametres = dict(n_blades=4, beta_deg=22.0, blade_wrap_deg=15.0,
                          r_hub=0.030, r_tip=0.100, thickness=0.003)
        parametres.update(kwargs)
        result = self.lire(synthetic.axial_impeller(**parametres))
        self.assertEqual(result.topology.machine_type, "axiale")
        return propulsion.analyse(
            result.topology, result.blades, rpm=RPM, speed=3.0, fluid="eau"
        )


class TestPlafondDeFroude(PropulsionTestCase):
    """Le disque actif ideal est un plafond, pas une reference indicative."""

    def test_le_rendement_ne_depasse_jamais_le_plafond(self):
        """Sur toute la courbe, eta <= eta_Froude.

        C'est l'invariant qui separe un modele d'helice d'un ajustement
        empirique : la quantite de mouvement borne le rendement pour une poussee
        donnee, quelles que soient les pales.  Un depassement serait une erreur
        de programme, et l'analyse le declare comme telle.
        """
        result = self.helice()
        self.assertTrue(result.curve, "la courbe de poussee est vide")
        for point in result.curve:
            if point.efficiency <= 0.0 or point.froude_efficiency <= 0.0:
                continue
            with self.subTest(J=point.advance_ratio):
                self.assertLessEqual(
                    point.efficiency, point.froude_efficiency,
                    f"a J = {point.advance_ratio:.2f}, le rendement calcule "
                    f"({point.efficiency:.3f}) passe au-dessus du disque actif ideal "
                    f"({point.froude_efficiency:.3f})",
                )

    def test_sans_trainee_le_rendement_se_rapproche_du_plafond(self):
        """Trainee de profil mise a zero, l'ecart au plafond doit se refermer.

        Validation interne : le plafond de Froude ne compte ni la trainee de
        profil, ni la giration, ni la perte de bout de pale.  En retirant la
        premiere, l'ecart restant doit diminuer nettement -- sans s'annuler,
        puisque les deux autres subsistent.  Un modele dont l'ecart ne bougerait
        pas ne ferait pas passer la trainee par ou il faut.
        """
        def marge_maximale():
            result = self.helice()
            best = max((p for p in result.curve if p.efficiency > 0.0),
                       key=lambda p: p.efficiency)
            return best.ceiling_gap()

        complet = marge_maximale()
        base, induced, thickness = config.CD_BASE, config.CD_INDUCED_K, config.CD_THICKNESS_K
        config.CD_BASE = config.CD_INDUCED_K = config.CD_THICKNESS_K = 0.0
        try:
            sans_trainee = marge_maximale()
        finally:
            config.CD_BASE, config.CD_INDUCED_K, config.CD_THICKNESS_K = (
                base, induced, thickness
            )

        self.assertLess(sans_trainee, complet - 0.05,
                        "retirer la trainee de profil n'a pas rapproche du plafond")
        self.assertGreater(sans_trainee, 0.0,
                           "l'ecart au plafond ne peut pas s'annuler : il reste la "
                           "giration et le bout de pale")

    def test_le_plafond_n_est_pas_defini_a_l_arret(self):
        """A vitesse d'avance nulle, le rendement propulsif n'a pas de sens."""
        self.assertEqual(propulsion.froude_efficiency(100.0, 998.0, 0.0, 0.03), 0.0)
        self.assertEqual(propulsion.froude_efficiency(100.0, 998.0, 5.0, 0.0), 0.0)
        # Poussee nulle : rien a comparer.
        self.assertEqual(propulsion.froude_efficiency(0.0, 998.0, 5.0, 0.03), 0.0)

    def test_le_plafond_baisse_quand_la_charge_monte(self):
        """Plus la poussee est forte a vitesse donnee, plus le plafond est bas."""
        precedent = 1.0
        for thrust in (50.0, 200.0, 800.0, 3200.0):
            value = propulsion.froude_efficiency(thrust, config.RHO, 4.0, 0.03)
            self.assertLess(value, precedent)
            self.assertGreater(value, 0.0)
            precedent = value


class TestCourbeDePoussee(PropulsionTestCase):
    """La forme de la courbe : ce qu'une helice fait quand elle avance."""

    def test_poussee_decroissante_puis_negative(self):
        """La poussee part d'un maximum a l'arret, decroit, puis change de signe.

        C'est la signature d'une helice a pas fixe : au-dela de son parametre
        d'avance de poussee nulle, elle freine.  Un modele d'helice qui ne rend
        pas cette forme-la ne decrit pas une helice.
        """
        result = self.helice()
        self.assertGreater(result.static_thrust, 0.0, "poussee statique nulle ou negative")
        poussees = [point.thrust for point in result.curve]
        self.assertAlmostEqual(poussees[0], result.static_thrust, delta=1e-6)
        for avant, apres in zip(poussees, poussees[1:]):
            self.assertLessEqual(apres, avant + 1e-9, "la poussee remonte avec l'avance")
        self.assertLess(poussees[-1], 0.0, "la poussee ne change jamais de signe")
        self.assertGreater(result.zero_thrust_advance, 0.0)

    def test_le_rendement_passe_par_un_maximum(self):
        """Nul a l'arret, nul a poussee nulle : un maximum entre les deux."""
        result = self.helice()
        rendements = [point.efficiency for point in result.curve]
        self.assertEqual(rendements[0], 0.0, "rendement propulsif non nul a l'arret")
        best = max(range(len(rendements)), key=lambda i: rendements[i])
        self.assertGreater(rendements[best], 0.0)
        self.assertNotIn(best, (0, len(rendements) - 1),
                         "le maximum de rendement tombe sur un bord du balayage")

    def test_coefficients_invariants_par_le_regime(self):
        """A J egal, CT et CP ne dependent pas de la vitesse de rotation.

        C'est la similitude : les coefficients sont adimensionnels par
        construction, et une helice tournant deux fois plus vite a la meme
        avance relative doit les retrouver identiques.  Un ecart signalerait une
        vitesse oubliee quelque part dans le bilan.
        """
        result = self.lire(synthetic.axial_impeller(
            n_blades=4, beta_deg=22.0, blade_wrap_deg=15.0,
            r_hub=0.030, r_tip=0.100, thickness=0.003,
        ))
        advance = 0.5
        mesures = []
        for rpm in (900.0, 1800.0):
            analyse = propulsion.analyse(result.topology, result.blades, rpm=rpm,
                                         speed=0.0, fluid="eau")
            diameter = analyse.diameter
            speed = advance * (rpm / config.SECONDS_PER_MINUTE) * diameter
            point = propulsion.operating_point(
                analyse.stations, speed, rpm, analyse.rho, diameter,
                analyse.hub_diameter, analyse.n_blades,
            )
            mesures.append((point.thrust_coefficient, point.power_coefficient))
        self.assertClose(mesures[0][0], mesures[1][0], rel=1e-6, msg="CT")
        self.assertClose(mesures[0][1], mesures[1][1], rel=1e-6, msg="CP")


class TestProfilLuSurLaPiece(PropulsionTestCase):
    """La portance vient de la cambrure mesuree, pas d'un coefficient saisi."""

    def test_angle_de_portance_nulle(self):
        """`alpha_0 = -2 h/c` : un profil cambre porte deja a incidence nulle."""
        self.assertEqual(propulsion.zero_lift_angle(0.0), 0.0)
        self.assertLess(propulsion.zero_lift_angle(0.04), 0.0)
        self.assertAlmostEqual(propulsion.zero_lift_angle(0.05), -0.10, places=12)
        cl, _ = propulsion.lift_coefficient(0.0, 0.05)
        self.assertGreater(cl, 0.0, "un profil cambre doit porter a incidence nulle")

    def test_decrochage_borne_la_portance(self):
        """Au-dela du decrochage, Cl est plafonne et l'etat est signale."""
        cl, stalled = propulsion.lift_coefficient(math.radians(45.0), 0.0)
        self.assertTrue(stalled)
        self.assertAlmostEqual(cl, config.CL_STALL)
        cl, stalled = propulsion.lift_coefficient(math.radians(3.0), 0.0)
        self.assertFalse(stalled)

    def test_la_trainee_croit_avec_l_epaisseur(self):
        """L'epaisseur relative mesuree alourdit la trainee de base."""
        mince = propulsion.drag_coefficient(config.CD_MIN_DRAG_CL, 0.05)
        epais = propulsion.drag_coefficient(config.CD_MIN_DRAG_CL, 0.15)
        self.assertGreater(epais, mince)

    def test_la_cambrure_est_lue_sur_la_geometrie(self):
        """Les coupes portent une fleche relative plausible, non nulle."""
        result = self.lire(synthetic.axial_impeller(
            n_blades=4, beta_deg=22.0, blade_wrap_deg=15.0,
            r_hub=0.030, r_tip=0.100, thickness=0.003,
        ))
        cambrures = [s.camber_ratio for s in result.blades.sections if s.camber_ratio > 0.0]
        self.assertTrue(cambrures, "aucune cambrure relative mesuree")
        for value in cambrures:
            self.assertLess(value, 0.30, "fleche relative invraisemblable")


class TestRefusDeRepondre(PropulsionTestCase):
    """Hors de son domaine, le modele dit non plutot que de sortir un chiffre."""

    def test_une_roue_centrifuge_n_est_pas_une_helice(self):
        """Le modele d'helice libre suppose un disque non carene et axial."""
        result = self.lire(synthetic.centrifugal_impeller(), name="centrifuge.stl")
        self.assertEqual(result.topology.machine_type, "centrifuge")
        analyse = propulsion.analyse(result.topology, result.blades, rpm=RPM,
                                     speed=3.0, fluid="eau")
        self.assertIsNone(analyse.point)
        self.assertEqual(analyse.confidence, "low")
        self.assertTrue(any("ne s'applique pas" in w for w in analyse.warnings))

    def test_l_analyse_propulsive_est_facultative(self):
        """Sans vitesse d'avance demandee, rien n'est calcule ni publie."""
        result = self.lire(synthetic.axial_impeller(
            n_blades=4, beta_deg=22.0, blade_wrap_deg=15.0,
            r_hub=0.030, r_tip=0.100, thickness=0.003,
        ))
        self.assertIsNone(result.propulsion)
        self.assertNotIn("propulsion", result.to_dict()["confiance"])

    def test_domaine_des_entrees(self):
        """Vitesse d'avance et fluide sont bornes comme le reste."""
        with self.assertRaises(ValueError) as capture:
            Options(propulsion_speed=-1.0).check()
        self.assertIn("vitesse d'avance", str(capture.exception))
        with self.assertRaises(ValueError) as capture:
            Options(fluid="mercure").check()
        self.assertIn("fluide inconnu", str(capture.exception))
        Options(propulsion_speed=10.0, fluid="air").check()


class TestChaineComplete(PropulsionTestCase):
    """De la CLI au rapport, l'analyse propulsive traverse tout l'outil."""

    def test_resultat_serialisable_et_rapport(self):
        """Le JSON porte l'analyse, et le rapport en tire une section."""
        from impeller_analyzer.io import report as report_module

        result = self.lire(
            synthetic.axial_impeller(n_blades=4, beta_deg=22.0, blade_wrap_deg=15.0,
                                     r_hub=0.030, r_tip=0.100, thickness=0.003),
            propulsion_speed=3.0,
        )
        self.assertIsNotNone(result.propulsion)
        data = result.to_dict()["propulsion"]
        self.assertIsNotNone(data["point"])
        self.assertGreater(data["poussee_statique_N"], 0.0)

        path = report_module.write_markdown(result, self.workdir, source="helice.stl")
        with open(path, encoding="utf-8") as handle:
            markdown = handle.read()
        self.assertIn("En helice libre", markdown)
        self.assertIn("Plafond ideal (Froude)", markdown)


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
