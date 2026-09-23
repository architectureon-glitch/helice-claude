"""L'outil n'etudie que les helices toroidales, et compare chacune a sa jumelle normale.

Deux engagements, testes ici :

- une piece qui n'est pas toroidale est **refusee**, avec la raison -- la ligne
  de commande et la page levent le drapeau, l'appel direct au moteur le laisse
  baisse, le calcul d'une helice normale servant de reference ;
- chaque resultat vient avec le comparatif de la **meme geometrie a aubes
  normales** : un brin par aube au lieu d'une boucle. Ce que le modele compte de
  l'ecart, et ce qu'il n'en compte pas, est publie avec le tableau.
"""

import math
import unittest

from helpers import BaseTestCase
from test_composants_declares import ComposantsTestCase, decale
import test_pale_isolee as pale

from impeller_analyzer import components, config, synthetic
from impeller_analyzer.analysis import MACHINE_PROPELLER, Options, run, run_components
from impeller_analyzer.hydraulics import propulsion
from impeller_analyzer.io import report, writer

GRILLE = dict(grid_nr=80, grid_nz=80, n_theta=360, symmetry_check=False)


class ImportGlobal(BaseTestCase):
    def ecrire(self, mesh, nom):
        path = self.path(nom)
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return path


class TestRefusImportGlobal(ImportGlobal):
    """Import d'un seul fichier : la lecture tranche, une declaration la complete."""

    def test_roue_normale_refusee(self):
        result = run(self.ecrire(synthetic.centrifugal_impeller(), "normale.stl"),
                     Options(unit="cm", toroidal_only=True, **GRILLE))
        self.assertIn("pas toroidale", result.rejection)
        self.assertIn("analyse interrompue", result.warnings[0])
        self.assertEqual(result.curves, [])
        self.assertIsNone(result.comparison)

    def test_le_moteur_reste_general(self):
        """Sans le drapeau, le moteur calcule une roue normale : c'est la reference du comparatif."""
        result = run(self.ecrire(synthetic.centrifugal_impeller(), "normale.stl"),
                     Options(unit="cm", **GRILLE))
        self.assertIsNone(result.rejection)
        self.assertTrue(result.curves)
        self.assertIsNone(result.comparison)

    def test_declaree_toroidale_poursuivie_et_comparee(self):
        """La lecture peut se tromper : la declaration poursuit l'analyse, en le disant."""
        result = run(self.ecrire(synthetic.centrifugal_impeller(), "declaree.stl"),
                     Options(unit="cm", toroidal_only=True, topology_declared=True,
                             blade_topology=components.BLADE_TOROIDAL, **GRILLE))
        self.assertIsNone(result.rejection)
        self.assertIn("ne voit aucune aube en boucle", result.warnings[0])
        comparaison = result.comparison
        self.assertEqual(comparaison.machine, "pompe")
        lignes = {row.quantity: row for row in comparaison.rows}
        # La ligne moyenne ne voit pas la forme de l'aube : meme hauteur.
        hauteur = lignes["Hauteur au debit du cas"]
        self.assertAlmostEqual(hauteur.toroidal, hauteur.normal, places=9)
        # Le canal, si : deux brins frottent plus qu'un.
        surface = lignes["Surface mouillee des aubes et flasques"]
        self.assertGreater(surface.toroidal, surface.normal)
        rendement = lignes["Rendement de comparaison, au meilleur rendement"]
        self.assertLess(rendement.toroidal, rendement.normal)
        self.assertTrue(rendement.absolute)

    def test_declaree_conventionnelle_refusee(self):
        result = run(self.ecrire(synthetic.toroidal_propeller(), "tore.stl"),
                     Options(unit="cm", toroidal_only=True, topology_declared=True,
                             blade_topology=components.BLADE_CONVENTIONAL, **GRILLE))
        self.assertIn("declaree conventionnelle", result.rejection)


class TestRefusComposants(ComposantsTestCase):
    """Mode composants : la declaration, et le genre d'une pale etanche d'un seul tenant."""

    def pieces_(self, blade):
        return {
            components.SLOT_INLET: self.ecrire(decale(synthetic.cylinder(0.038, 0.001, z_center=0.050)), "e.stl"),
            components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.095, 0.096, 0.010, z_center=0.005)), "s.stl"),
            components.SLOT_HUB: self.ecrire(pale.parois(False), "c.stl"),
            components.SLOT_BLADE: self.ecrire(blade, "p.stl"),
        }

    def analyser(self, blade, **options):
        valeurs = dict(unit="cm", machine="pompe_carenee", blades=5, speeds=(1450.0,),
                       rotation=-1, component_paths=self.pieces_(blade), toroidal_only=True)
        valeurs.update(options)
        return run_components(Options(**valeurs))

    def test_declaree_conventionnelle(self):
        result = self.analyser(pale.aube(20.0, 30.0), blade_topology=components.BLADE_CONVENTIONAL)
        self.assertIn("declaree conventionnelle", result.rejection)
        self.assertIsNone(result.topology)

    def test_pale_de_genre_zero(self):
        """Une pale etanche d'un seul tenant, sans anse, n'est pas une boucle."""
        result = self.analyser(pale.aube(20.0, 30.0), blade_topology=components.BLADE_TOROIDAL)
        self.assertIn("n'est pas une boucle", result.rejection)

    def test_sans_drapeau_rien_n_est_refuse(self):
        result = self.analyser(pale.aube(20.0, 30.0), blade_topology=components.BLADE_CONVENTIONAL,
                               toroidal_only=False)
        self.assertIsNone(result.rejection)
        self.assertIsNone(result.comparison)


class TestComparatifRoueEnSerie(ComposantsTestCase):
    """La roue de type hel1 : entree par la vis du brin haut, sortie par le brin bas."""

    def resultat(self):
        paths = {
            components.SLOT_INLET: self.ecrire(decale(synthetic.cylinder(0.038, 0.001, z_center=0.050)), "e.stl"),
            components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.095, 0.096, 0.010, z_center=0.005)), "s.stl"),
            components.SLOT_HUB: self.ecrire(pale.parois(True), "c.stl"),
            components.SLOT_BLADE: self.ecrire(pale.boucle(pale.TestAntiRetour.PAS), "p.stl"),
        }
        # Deux rubans disjoints : genre 0, mais en deux morceaux -- le genre ne
        # tranche pas, et la declaration est retenue.
        return run_components(Options(
            unit="cm", machine="pompe_carenee", blades=5, speeds=(1450.0,), rotation=-1,
            blade_topology=components.BLADE_TOROIDAL, component_paths=paths, toroidal_only=True,
        ))

    def test_entree_de_la_normale_au_bord_d_attaque(self):
        result = self.resultat()
        self.assertIsNone(result.rejection)
        lignes = {row.quantity: row for row in result.comparison.rows}
        incidence = lignes["Incidence au bord d'attaque, au debit du cas"]
        self.assertLess(abs(incidence.toroidal), 1e-6, "le cas calcule est l'incidence nulle")
        self.assertGreater(abs(incidence.normal), 1.0)
        debit = lignes["Debit d'incidence nulle"]
        self.assertNotAlmostEqual(debit.toroidal, debit.normal, places=3)
        self.assertTrue(any("vis de son brin amont" in m for m in result.comparison.modelled))

    def test_variante_a_diametre_egal(self):
        """Bord de fuite a 90 mm, fente a 95,5 : la variante prolonge l'aube jusqu'a la fente."""
        comparaison = self.resultat().comparison
        self.assertIn("a diametre egal", comparaison.variant_title)
        a_vide = next(r for r in comparaison.variant_rows if r.quantity == "Hauteur a debit nul")
        attendu = (0.0955 / 0.090) ** 2
        self.assertLess(abs(a_vide.normal / a_vide.toroidal / attendu - 1.0), 0.05)

    def test_rapport(self):
        result = self.resultat()
        path = report.write_markdown(result, self.workdir, source="serie")
        with open(path, encoding="utf-8") as handle:
            texte = handle.read()
        self.assertIn("## Si l'helice etait normale", texte)
        self.assertIn("Ce que le modele compte", texte)
        self.assertIn("Ce qu'il ne compte pas", texte)
        self.assertIn("deg |", texte)


class TestComparatifHeliceLibre(ImportGlobal):
    """Helice libre : la boucle n'a pas de bout libre."""

    def test_perte_de_bout_supprimee_par_la_boucle(self):
        ouverte = propulsion.prandtl_loss(0.095, 0.02, 0.10, 0.3, 3)
        fermee = propulsion.prandtl_loss(0.095, 0.02, 0.10, 0.3, 3, closed_tip=True)
        pied_seul = propulsion.prandtl_loss(0.095, 0.02, 10.0, 0.3, 3)
        self.assertGreater(fermee, ouverte)
        self.assertAlmostEqual(fermee, pied_seul, places=9)

    def test_trois_colonnes(self):
        result = run(self.ecrire(synthetic.toroidal_propeller(), "tore.stl"),
                     Options(unit="cm", toroidal_only=True, machine=MACHINE_PROPELLER,
                             propulsion_speed=2.0, speeds=(1000.0,), **GRILLE))
        comparaison = result.comparison
        self.assertEqual(comparaison.machine, "helice_libre")
        self.assertEqual(len(comparaison.columns), 3)
        poussee = next(r for r in comparaison.rows if r.quantity == "Poussee")
        self.assertGreater(poussee.toroidal_alt, poussee.toroidal)
        self.assertGreater(poussee.toroidal, 0.0)
        # Le rendement propulsif publie est celui du bilan de la toroidale, bout conserve.
        rendement = next(r for r in comparaison.rows if r.quantity == "Rendement propulsif")
        self.assertAlmostEqual(rendement.toroidal, result.propulsion.point.efficiency * 100.0, places=6)


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
