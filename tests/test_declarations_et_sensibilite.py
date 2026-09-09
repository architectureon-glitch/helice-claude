"""Ce que l'outil a le droit d'affirmer, et ce qu'il doit laisser declarer.

Quatre corrections tenues ensemble par ce module, toutes de la meme famille :
l'outil affirmait quelque part plus qu'il ne savait.

* La **classification geometrique** repondait « centrifuge » en confiance haute
  sur une helice a aubes en boucle -- et cette erreur fermait le mode dont
  l'utilisateur avait besoin.  Elle se declare desormais, et une aube en boucle
  la fait passer en confiance faible d'elle-meme.
* Un **alesage traversant** etait rendu comme une absence de moyeu, `r1h = 0`,
  ce qui faisait compter le trou central comme section de passage.
* Une valeur **imposee en ligne de commande** se lisait comme une valeur
  extraite du maillage.
* L'incertitude annoncee par le modele suppose la geometrie juste ; rien ne
  disait ce qu'une erreur de lecture ferait vraiment.
"""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import (
    MACHINE_PROPELLER,
    MACHINE_PUMP,
    Options,
    run,
)
from impeller_analyzer.geometry import topology as topology_module
from impeller_analyzer.hydraulics import meanline, sensitivity
from impeller_analyzer.io import report as report_module
from impeller_analyzer.io import writer
from impeller_analyzer.provenance import DECLARE, MESURE

GRID = 130
N_THETA = 280


class DeclarationTestCase(BaseTestCase):
    """Lit une helice a aubes en boucle : la forme qui met la lecture en defaut."""

    def lire(self, mesh, name="piece.stl", **options):
        """Analyse complete d'un maillage de synthese."""
        path = self.path(name)
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(
            unit="cm", speeds=(1450.0,), grid_nr=GRID, grid_nz=GRID, n_theta=N_THETA,
            symmetry_check=False, rotation=1, **options
        ))

    def helice_toroidale(self, **options):
        """Helice a cinq pales en boucle, alesage traversant, comme la piece reelle."""
        return self.lire(synthetic.toroidal_propeller(
            n_blades=5, r_hub=0.045, r_tip=0.1336, hub_height=0.070, bore_radius=0.0267,
        ), **options)


class TestClassificationDeclarable(DeclarationTestCase):
    """C1 : le type de machine se declare, et une boucle interdit de l'affirmer."""

    def test_une_aube_en_boucle_abaisse_la_confiance_du_type(self):
        """Les criteres de classification supposent un canal meridien.

        Rapport r2/r1s, solidite, largeur de sortie : tous supposent un canal
        borde par le moyeu et le carter, que le fluide traverse une fois.  Une
        aube en boucle n'en a pas.  La classification garde une valeur -- il
        faut bien en publier une -- mais cesse de pouvoir etre affirmee.
        """
        result = self.helice_toroidale()
        self.assertTrue(result.blade_loops.looped, "la piece devait etre vue en boucle")
        self.assertEqual(result.confidence.get_level("type_de_roue"), "low")
        self.assertTrue(
            any("classification geometrique" in w for w in result.warnings),
            "l'abaissement de confiance doit etre explique",
        )

    def test_le_type_declare_prime_sur_la_lecture(self):
        """`--type-de-roue` remplace la classification et retablit la confiance."""
        lu = self.helice_toroidale()
        declare = self.helice_toroidale(wheel_type=topology_module.AXIAL)
        self.assertEqual(declare.topology.machine_type, topology_module.AXIAL)
        self.assertEqual(declare.confidence.get_level("type_de_roue"), "high")
        self.assertEqual(declare.provenance.get_source("type_de_roue"), DECLARE)
        self.assertEqual(lu.provenance.get_source("type_de_roue"), MESURE)

    def test_le_mode_helice_libre_ne_depend_plus_de_la_classification(self):
        """C'etait le blocage : une lecture fausse fermait le mode voulu."""
        auto = self.helice_toroidale(propulsion_speed=3.0)
        self.assertNotEqual(auto.topology.machine_type, topology_module.AXIAL)
        self.assertIsNone(auto.propulsion.point, "sans declaration, le modele doit refuser")

        declare = self.helice_toroidale(
            machine=MACHINE_PROPELLER, propulsion_speed=3.0
        )
        self.assertIsNotNone(declare.propulsion.point, "la declaration doit debloquer")
        self.assertGreater(declare.propulsion.point.thrust, 0.0)

    def test_la_pompe_carenee_refuse_l_analyse_propulsive(self):
        """Declarer une pompe carenee ferme le mode helice, meme avec une vitesse."""
        result = self.helice_toroidale(machine=MACHINE_PUMP, propulsion_speed=3.0)
        self.assertIsNone(result.propulsion)
        self.assertTrue(any("pompe carenee" in w for w in result.warnings))

    def test_declarations_hors_domaine(self):
        """Un modele ou un type inconnu est refuse avec la liste des valeurs."""
        with self.assertRaises(ValueError) as capture:
            Options(machine="turbine").check()
        self.assertIn("modele hydraulique inconnu", str(capture.exception))
        with self.assertRaises(ValueError) as capture:
            Options(wheel_type="helicoidale").check()
        self.assertIn("type de roue inconnu", str(capture.exception))


class TestAlesageEtMoyeu(DeclarationTestCase):
    """C2 : un trou traversant n'est ni un moyeu, ni une section de passage."""

    def test_alesage_reconnu_et_nomme(self):
        """La matiere commence a r > 0 sans qu'aucun moyeu plein n'existe."""
        result = self.helice_toroidale()
        topology = result.topology
        self.assertEqual(topology.hub_kind, topology_module.HUB_BORE)
        self.assertClose(topology.r_1h, 0.0267, rel=0.05)
        self.assertClose(topology.bore_radius, topology.r_1h, rel=1e-9)

    def test_la_section_d_entree_exclut_le_trou(self):
        """Compter l'alesage comme passage gonflerait le debit d'autant."""
        result = self.helice_toroidale()
        topology = result.topology
        attendu = math.pi * (topology.r_1s ** 2 - topology.r_1h ** 2) * config.TAU_1
        self.assertClose(topology.area_1, attendu, rel=1e-9)
        avec_le_trou = math.pi * topology.r_1s ** 2 * config.TAU_1
        self.assertGreater(avec_le_trou, topology.area_1 * 1.05,
                           "le trou doit peser sur la section")

    def test_confiance_moyenne_sans_moyeu_plein(self):
        """Un rayon interieur de matiere n'est pas un rayon de moyeu."""
        result = self.helice_toroidale()
        self.assertEqual(result.confidence.get_level("rayon_de_moyeu"), "medium")

    def test_moyeu_plein_reconnu_et_confiance_haute(self):
        """La meme helice sans alesage porte un vrai moyeu."""
        result = self.lire(synthetic.toroidal_propeller(
            n_blades=5, r_hub=0.045, r_tip=0.1336, hub_height=0.070,
        ), name="moyeu.stl")
        self.assertEqual(result.topology.hub_kind, topology_module.HUB_SOLID)
        self.assertGreater(result.topology.r_1h, 0.0)
        self.assertEqual(result.confidence.get_level("rayon_de_moyeu"), "high")

    def test_le_plateau_arriere_n_est_pas_un_moyeu_d_aspiration(self):
        """Une roue centrifuge ouverte : l'oeillard est libre malgre le flasque.

        Le flasque arriere est bien du plein depuis l'axe, mais a l'autre bout
        de la piece.  Le prendre pour le moyeu d'aspiration donnerait un rayon
        interieur superieur au rayon d'oeillard, et une section d'entree
        **negative**.
        """
        result = self.lire(synthetic.centrifugal_impeller(), name="centrifuge.stl")
        self.assertEqual(result.topology.hub_kind, topology_module.HUB_NONE)
        self.assertEqual(result.topology.r_1h, 0.0)
        self.assertGreater(result.topology.area_1, 0.0)


class TestProvenance(DeclarationTestCase):
    """C4 : une valeur imposee ne se lit pas comme une valeur extraite."""

    def test_le_tableau_porte_provenance_et_confiance(self):
        """Quatre colonnes : libelle, valeur, provenance, confiance."""
        result = self.helice_toroidale(blades=5, wheel_type=topology_module.AXIAL)
        table = report_module.geometry_table(result)
        self.assertTrue(table)
        for row in table:
            self.assertEqual(len(row), 4)
            self.assertIn(row[2], ("mesure", "declare", "defaut"))

    def test_declare_contre_mesure(self):
        """Ce que l'utilisateur impose est marque, le reste est marque mesure."""
        result = self.helice_toroidale(blades=5, beta1_deg=12.0, beta2_deg=20.0)
        source = result.provenance
        self.assertEqual(source.get_source("nombre_de_pales"), DECLARE)
        self.assertEqual(source.get_source("angles_de_pale"), DECLARE)
        self.assertEqual(source.get_source("sens_de_rotation"), DECLARE)
        self.assertEqual(source.get_source("rayons"), MESURE)

    def test_un_seul_angle_impose_n_est_ni_l_un_ni_l_autre(self):
        """beta1 impose et beta2 lu : la case ne peut pas dire « declare »."""
        result = self.helice_toroidale(beta1_deg=12.0)
        self.assertEqual(result.provenance.get_source("angles_de_pale"), "defaut")

    def test_la_provenance_est_serialisee(self):
        """Le JSON porte la table, a cote de celle des confiances."""
        result = self.helice_toroidale(blades=5)
        data = result.to_dict()
        self.assertIn("provenance", data)
        self.assertEqual(data["provenance"]["nombre_de_pales"], DECLARE)


class TestSensibilitePubliee(BaseTestCase):
    """C5 : ce qu'une erreur d'un degre fait vraiment aux resultats."""

    def entree(self, **remplacements):
        """Geometrie de reference du modele 1D, eventuellement retouchee."""
        import dataclasses

        path = self.path("roue.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), path,
                         unit_factor=config.UNIT_FACTOR)
        result = run(path, Options(
            unit="cm", speeds=(1450.0,), grid_nr=120, grid_nz=120, n_theta=240,
            symmetry_check=False, rotation=1,
        ))
        data = meanline.MeanlineInput.from_geometry(result.topology, result.blades)
        return dataclasses.replace(data, **remplacements) if remplacements else data

    def test_trois_grandeurs_trois_entrees(self):
        """Hauteur, debit et NPSHr, contre beta1, beta2 et le diametre."""
        report = sensitivity.analyse(self.entree(), 1450.0)
        self.assertEqual([row.quantity for row in report.rows],
                         ["hauteur", "debit", "npshr"])

    def test_l_elasticite_au_diametre_retrouve_les_lois_de_similitude(self):
        """H comme D2, Q comme D3, NPSHr comme D2 : la colonne se controle.

        C'est ce qui distingue cette colonne d'une simple curiosite : sa valeur
        est connue d'avance, et un ecart signalerait une erreur de programme
        plutot qu'une propriete de la roue.
        """
        report = sensitivity.analyse(self.entree(), 1450.0)
        attendu = {"hauteur": 2.0, "debit": 3.0, "npshr": 2.0}
        for row in report.rows:
            with self.subTest(grandeur=row.quantity):
                self.assertClose(row.diameter, attendu[row.quantity], abs_=0.02)

    def test_une_aube_couchee_declasse_la_hauteur(self):
        """Au-dela du seuil par degre, la confiance tombe d'elle-meme."""
        from impeller_analyzer.confidence import ConfidenceMap

        report = sensitivity.analyse(self.entree(beta2_deg=6.0), 1450.0)
        ligne = report.row("hauteur")
        if ligne is None:  # pas de point de fonctionnement : rien a declasser
            self.skipTest("la geometrie retouchee n'a pas de point de fonctionnement")
        self.assertTrue(ligne.alarming())
        confidence = ConfidenceMap()
        for quantity in ("hauteur", "puissance", "couple", "rendement"):
            confidence.set(quantity, "high")
        sensitivity.apply(report, confidence)
        self.assertEqual(confidence.get_level("hauteur"), "low")
        # La puissance et le couple derivent de la hauteur : ils suivent.
        self.assertEqual(confidence.get_level("puissance"), "low")
        self.assertEqual(confidence.get_level("couple"), "low")

    def test_une_roue_saine_n_est_pas_declassee(self):
        """Le seuil ne doit pas se declencher sur une geometrie ordinaire."""
        report = sensitivity.analyse(self.entree(), 1450.0)
        self.assertEqual(report.downgraded, [])
        for row in report.rows:
            self.assertFalse(row.alarming())

    def test_la_sensibilite_est_par_degre(self):
        """La mesure historique divisait par deux degres en en annoncant un.

        `head_sensitivity` rend desormais la variation **par degre**, ce que son
        avertissement disait deja.  Elle doit donc coincider avec la colonne
        beta2 du tableau.
        """
        data = self.entree()
        ancienne = meanline.head_sensitivity(data, 1450.0)
        ligne = sensitivity.analyse(data, 1450.0).row("hauteur")
        self.assertClose(ancienne, abs(ligne.beta2), rel=1e-9)

    def test_serialisation_sans_infini(self):
        """L'infini n'est pas du JSON : il se serialise en `None`."""
        report = sensitivity.analyse(self.entree(), 1450.0)
        data = report.to_dict()
        self.assertEqual(len(data["lignes"]), 3)
        for ligne in data["lignes"]:
            for key in ("par_degre_de_beta1", "par_degre_de_beta2",
                        "elasticite_au_diametre"):
                value = ligne[key]
                self.assertTrue(value is None or math.isfinite(value))


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
