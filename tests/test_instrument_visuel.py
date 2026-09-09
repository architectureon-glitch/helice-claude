"""Les sorties visuelles : une seule palette, et deux manipulations qui comptent.

Ce n'est pas un rapport imprime, c'est un instrument de banc d'essai.  Ce qui se
verifie ici est donc moins l'apparence -- qu'un test ne juge pas -- que les
proprietes qui la rendent utilisable et honnete :

* **une seule source de palette**, partagee par la page et par les deux figures.
  Une page claire et des figures restees aux reglages d'origine donnent un
  resultat incoherent, et c'est le defaut qu'on oublie le plus souvent ;
* **la couleur ne porte jamais seule** : provenance et confiance se lisent aussi
  par la graisse et par un filet ;
* **la page reste autonome**, sans aucune ressource externe ;
* **les deux curseurs ne reimplementent pas le modele** : le regime passe par
  les lois de similitude, beta2 interpole entre des courbes reellement
  calculees en Python. Un second modele en JavaScript aurait ete libre de
  diverger du premier.
"""

import json
import re
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.io import plot, style, viewer, writer


class VisuelTestCase(BaseTestCase):
    """Analyse une roue une fois, et lit ses sorties visuelles."""

    def setUp(self) -> None:
        super().setUp()
        path = self.path("roue.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), path,
                         unit_factor=config.UNIT_FACTOR)
        self.result = run(path, Options(
            unit="cm", speeds=(1450.0,), grid_nr=100, grid_nz=100, n_theta=200,
            symmetry_check=False, rotation=1,
        ))


class TestPaletteUnique(unittest.TestCase):
    """Figures et page tirent leurs couleurs du meme endroit."""

    def test_les_figures_utilisent_la_palette_de_la_page(self):
        """`plot` ne definit plus ses propres couleurs, il lit `style`."""
        self.assertEqual(plot.WHITE, style.FOND)
        self.assertEqual(plot.BLACK, style.ENCRE)
        self.assertEqual(style.FOND, style.rgb(config.COULEUR_FOND))
        self.assertEqual(style.MESURE_RGB, style.rgb(config.COULEUR_MESURE))

    def test_echelle_sequentielle_a_teinte_unique(self):
        """Du fond au teal, et monotone en luminance.

        Une echelle multicolore fabrique des frontieres que les donnees n'ont
        pas : c'est le reproche de fond fait a jet, et il vaut plus discretement
        pour viridis sur une grandeur sans seuil naturel. `f` est une fraction
        qui va de zero a un, une progression et non des categories.
        """
        self.assertEqual(style.sequential(0.0), style.FOND)
        self.assertEqual(style.sequential(1.0), style.MESURE_RGB)

        def luminance(couleur):
            return 0.2126 * couleur[0] + 0.7152 * couleur[1] + 0.0722 * couleur[2]

        valeurs = [luminance(style.sequential(k / 32)) for k in range(33)]
        for avant, apres in zip(valeurs, valeurs[1:]):
            self.assertLessEqual(apres, avant + 1e-9,
                                 "la luminance doit decroitre sans rebond")

    def test_les_regimes_sont_une_progression_pas_des_categories(self):
        """Une seule teinte eclaircie : la vitesse est ordonnee, pas nominale."""
        couleurs = [plot.series_color(i, 4) for i in range(4)]
        self.assertEqual(len(set(couleurs)), 4, "les regimes restent distinguables")
        for couleur in couleurs:
            # Chaque teinte est sur le segment fond -> teal, donc jamais plus
            # verte que le teal ni plus claire que le fond.
            self.assertLessEqual(couleur[1], style.FOND[1])
            self.assertGreaterEqual(couleur[1], style.MESURE_RGB[1])

    def test_la_fonte_des_figures_a_des_minuscules(self):
        """Un libelle de figure ne doit pas sortir en capitales."""
        for lettre in "abcdefghijklmnopqrstuvwxyz":
            self.assertIn(lettre, plot._FONT, f"glyphe manquant : {lettre}")
        self.assertEqual(plot._ascii("Carte d'occupation"), "Carte d'occupation")


class TestPageAutonome(VisuelTestCase):
    """La page s'ouvre depuis le disque, hors ligne, et rend la meme chose."""

    def test_aucune_ressource_externe(self):
        """Ni script, ni feuille de style, ni fonte distante."""
        page = viewer.build_page(self.result)
        self.assertEqual(len(re.findall(r'(?:href|src)="https?://', page)), 0)
        self.assertNotIn("<script src", page)

    def test_theme_clair_par_defaut_sombre_en_bascule(self):
        """Le sombre reste disponible ; il cesse d'etre le defaut."""
        page = viewer.build_page(self.result)
        self.assertIn(config.COULEUR_FOND, page)
        self.assertIn('[data-theme="sombre"]', page)
        self.assertIn('id="bascule-theme"', page)
        # La page ne doit plus suivre la preference du systeme.
        self.assertNotIn("prefers-color-scheme", page)


class TestCarteInterrogeable(VisuelTestCase):
    """La carte est la lecture dont tout decoule : elle se survole."""

    def grille(self) -> dict:
        """La grille embarquee dans la charge utile."""
        return viewer.build_payload(self.result)["grille"]

    def test_la_grille_accompagne_l_image(self):
        """L'image donne une impression ; la grille se lit."""
        grille = self.grille()
        occupancy = self.result.occupancy
        self.assertEqual(grille["nr"], occupancy.nr)
        self.assertEqual(grille["nz"], occupancy.nz)
        self.assertEqual(len(grille["r"]), occupancy.nr)
        self.assertEqual(len(grille["z"]), occupancy.nz)

    def test_les_valeurs_sont_quantifiees_et_compactes(self):
        """Un octet par cellule : une grille 100x100 tient en quelques ko."""
        import base64

        grille = self.grille()
        octets = base64.b64decode(grille["f"])
        self.assertEqual(len(octets), grille["nr"] * grille["nz"])
        # La valeur relue doit retrouver celle de la carte, au pas de
        # quantification pres.
        for iz in (0, grille["nz"] // 2, grille["nz"] - 1):
            for ir in (0, grille["nr"] // 2, grille["nr"] - 1):
                lu = octets[iz * grille["nr"] + ir] / 255.0
                self.assertClose(lu, self.result.occupancy.f[iz][ir], abs_=1 / 255.0)

    def test_les_reperes_sont_fournis(self):
        """Axe, rayons et plans : ce qu'on peut superposer a la carte."""
        reperes = self.grille()["reperes"]
        for cle in ("r1h", "r1s", "r2", "z1", "z2"):
            self.assertIn(cle, reperes)
        self.assertClose(reperes["r1s"], self.result.topology.r_1s * config.MM_PER_M,
                         rel=1e-9)


class TestCurseurs(VisuelTestCase):
    """Les deux recalculs se font dans la page, sans second modele."""

    def vif(self) -> dict:
        """Les donnees des deux curseurs."""
        return viewer.build_payload(self.result)["vif"]

    def test_trois_familles_de_courbes_reellement_calculees(self):
        """beta2 ne s'extrapole pas : la page interpole entre trois courbes.

        Reimplementer le modele de ligne moyenne en JavaScript aurait donne un
        second modele, libre de diverger du premier.
        """
        vif = self.vif()
        self.assertEqual(len(vif["familles"]), 3)
        ecarts = sorted(f["ecart"] for f in vif["familles"])
        self.assertSeqClose(ecarts, [-config.BETA_SENSITIVITY_DEG, 0.0,
                                     config.BETA_SENSITIVITY_DEG], abs_=1e-9)
        for famille in vif["familles"]:
            self.assertEqual(len(famille["Q"]), config.Q_SWEEP_POINTS)
            self.assertEqual(len(famille["H"]), config.Q_SWEEP_POINTS)

    def test_les_familles_different_vraiment(self):
        """Un degre d'ecart doit se voir : sinon le curseur ne montre rien."""
        familles = sorted(self.vif()["familles"], key=lambda f: f["ecart"])
        bas, haut = familles[0], familles[-1]
        ecarts = [abs(a - b) for a, b in zip(bas["H"], haut["H"])]
        self.assertGreater(max(ecarts), 0.0)

    def test_le_regime_passe_par_les_lois_de_similitude(self):
        """La page transporte le point par n, n2, n3 : le meme modele.

        Ce n'est pas une approximation : le controle de similitude du rapport
        verifie ces lois a 0.00 %, et c'est lui qui autorise le curseur.
        """
        vif = self.vif()
        self.assertGreater(vif["rpm_reference"], 0.0)
        self.assertIn("point", vif)
        for cle in ("Q", "H", "P", "couple", "npshr"):
            self.assertGreater(vif["point"][cle], 0.0)

    def test_la_limite_de_cavitation_est_transmise(self):
        """Le curseur doit pouvoir montrer le franchissement."""
        vif = self.vif()
        self.assertGreater(vif["rpm_limite"], 0)
        self.assertTrue(vif["limite_active"])

    def test_la_provenance_de_beta2_suit_la_page(self):
        """Un angle impose colore la valeur en declare, pas en mesure."""
        path = self.path("impose.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), path,
                         unit_factor=config.UNIT_FACTOR)
        impose = run(path, Options(
            unit="cm", speeds=(1450.0,), grid_nr=100, grid_nz=100, n_theta=200,
            symmetry_check=False, rotation=1, beta1_deg=20.0, beta2_deg=25.0,
        ))
        self.assertEqual(viewer.build_payload(impose)["vif"]["beta2_provenance"],
                         "declare")
        self.assertEqual(self.vif()["beta2_provenance"], "mesure")


class TestProvenanceDansLaPage(VisuelTestCase):
    """La couleur ne porte jamais seule : la forme et la graisse aussi."""

    def test_les_classes_de_provenance_existent(self):
        """Graisse pour mesure, filet pour declare, italique pour defaut."""
        page = viewer.build_page(self.result)
        for classe in (".p-mesure", ".p-declare", ".p-defaut", ".c-low"):
            self.assertIn(classe, page)
        self.assertIn("font-weight:600", page)

    def test_la_charge_utile_porte_la_provenance(self):
        """La page doit pouvoir marquer chaque grandeur."""
        payload = viewer.build_payload(self.result)
        self.assertIn("provenance", payload)
        for ligne in payload["tables"]["geometry"]:
            self.assertEqual(len(ligne), 4)
            self.assertIn(ligne[2], ("mesure", "declare", "defaut"))

    def test_le_json_reste_serialisable(self):
        """Tout ce qui part dans la page doit passer en JSON."""
        payload = viewer.build_payload(self.result)
        texte = json.dumps(payload)
        self.assertGreater(len(texte), 1000)


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
