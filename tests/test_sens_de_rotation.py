"""Le sens de rotation est une entree, pas un resultat.

La geometrie le suggere ; elle ne le tranche pas. Deux cas suffisent a le
montrer : une roue dont le rapport r2/r1s tombe a cheval sur la frontiere
mixte / centrifuge, ou les deux familles appliquent des regles **opposees** ; et
une aube quasi radiale, ou la lecture n'a plus aucune marge. Dans les deux cas
l'outil annoncait un sens avec l'aplomb d'un resultat. Il le demande desormais,
et se contente de suggerer.
"""

import unittest
import urllib.parse

from helpers import BaseTestCase

from impeller_analyzer import config, serve, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.confidence import HIGH, LOW
from impeller_analyzer.geometry import blade_angles as ba
from impeller_analyzer.io import report, writer

GRID = 110


class TestLectureDuNom(BaseTestCase):
    def test_noms_acceptes(self):
        """Les deux sens, plus quelques ecritures courantes."""
        self.assertEqual(ba.rotation_sign_from_name("horaire"), -1)
        self.assertEqual(ba.rotation_sign_from_name("antihoraire"), 1)
        self.assertEqual(ba.rotation_sign_from_name("anti-horaire"), 1)
        self.assertEqual(ba.rotation_sign_from_name("CCW"), 1)
        self.assertEqual(ba.rotation_sign_from_name("cw"), -1)
        self.assertIsNone(ba.rotation_sign_from_name(None))

    def test_nom_inconnu_refuse(self):
        """Mieux vaut un refus net qu'un sens choisi au hasard."""
        with self.assertRaises(ValueError):
            ba.rotation_sign_from_name("nord")


class TestBoutEnBout(BaseTestCase):
    def _run(self, mesh, **kwargs):
        path = self.path("roue.stl")
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(grid_nr=GRID, grid_nz=GRID, n_theta=240,
                                 symmetry_check=False, speeds=(1450.0,), **kwargs))

    def test_sans_indication_le_sens_reste_a_fournir(self):
        """Rien n'est affirme, mais la suggestion geometrique est conservee."""
        result = self._run(synthetic.centrifugal_impeller())
        blades = result.blades
        self.assertEqual(blades.rotation_sign, 0)
        self.assertEqual(blades.rotation_label, ba.NOT_SUPPLIED)
        self.assertFalse(blades.forced_rotation)
        self.assertEqual(result.confidence.get_level("sens_de_rotation"), LOW)
        # La lecture geometrique reste faite, et reste juste.
        self.assertEqual(blades.observed_rotation_sign, -1)
        self.assertIn(ba.CLOCKWISE, blades.observed_rotation_label)

    def test_le_sens_indique_est_retenu(self):
        """Ce que l'utilisateur donne prime, et passe en confiance haute."""
        result = self._run(synthetic.centrifugal_impeller(), rotation=1)
        self.assertEqual(result.blades.rotation_sign, 1)
        self.assertTrue(result.blades.forced_rotation)
        self.assertIn(ba.COUNTERCLOCKWISE, result.blades.rotation_label)
        self.assertEqual(result.confidence.get_level("sens_de_rotation"), HIGH)

    def test_un_sens_contraire_a_la_suggestion_est_signale(self):
        """L'utilisateur decide, mais il doit savoir qu'il contredit la geometrie."""
        result = self._run(synthetic.centrifugal_impeller(), rotation=1)
        self.assertEqual(result.blades.observed_rotation_sign, -1)  # le contraire
        self.assertTrue(any("inverse de ce que suggere" in m for m in result.warnings))

    def test_le_sens_ne_change_rien_a_l_hydraulique(self):
        """Le modele ne depend que de |omega| et des angles : les courbes sont identiques."""
        roue = synthetic.centrifugal_impeller()
        libre = self._run(roue)
        impose = self._run(roue, rotation=-1)
        for a, b in zip(libre.curves, impose.curves):
            point_a, point_b = a.nominal_point(), b.nominal_point()
            self.assertIsNotNone(point_a)
            self.assertClose(point_b.head, point_a.head, rel=1e-12)
            self.assertClose(point_b.flow, point_a.flow, rel=1e-12)
            self.assertClose(point_b.shaft_power, point_a.shaft_power, rel=1e-12)

    def test_la_reserve_sur_la_frontiere_se_tait_quand_le_sens_est_donne(self):
        """Inutile de demander de verifier un sens que l'utilisateur vient d'affirmer.

        Sur une roue posee a cheval sur la frontiere mixte / centrifuge, la
        suggestion geometrique n'est pas fiable et l'outil le dit. Mais une fois
        le sens fourni, la reserve ne porte plus sur rien.
        """
        roue = synthetic.centrifugal_impeller(r1=0.050, r2=0.090)
        libre = self._run(roue)
        if not libre.topology.rotation_ambiguity:
            self.skipTest("cette roue n'est pas a la frontiere des familles")
        self.assertTrue(any("pas fiable ici" in m for m in libre.warnings))
        impose = self._run(roue, rotation=1)
        self.assertFalse(any("pas fiable ici" in m for m in impose.warnings))

    def test_le_rapport_distingue_le_retenu_du_suggere(self):
        """Le tableau ne doit pas presenter une suggestion comme une mesure."""
        result = self._run(synthetic.centrifugal_impeller())
        chemin = report.write_markdown(result, self.workdir, source="roue.stl")
        with open(chemin, encoding="utf-8") as flux:
            texte = flux.read()
        self.assertIn("--rotation", texte)
        self.assertIn("suggere", texte)


class TestFormulaireWeb(BaseTestCase):
    def test_le_formulaire_transmet_le_sens(self):
        """L'application locale offre le meme choix que la ligne de commande."""
        for requete, attendu in (("", None), ("rotation=horaire", -1), ("rotation=antihoraire", 1)):
            with self.subTest(requete=requete or "(vide)"):
                query = urllib.parse.parse_qs(requete, keep_blank_values=True)
                self.assertEqual(serve.options_from_query(query).rotation, attendu)

    def test_sens_illisible_refuse_proprement(self):
        """Une valeur inconnue doit donner un message, pas une trace."""
        query = urllib.parse.parse_qs("rotation=nord")
        with self.assertRaises(serve.BadRequest):
            serve.options_from_query(query)


if __name__ == "__main__":
    unittest.main()
