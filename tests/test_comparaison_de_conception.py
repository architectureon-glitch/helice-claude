"""Pertes de canal : ce qui permet enfin de comparer deux conceptions.

Le modele de la SPEC pose chaque perte comme une fraction du point nominal --
le frottement vaut `K_FROTTEMENT_REL` fois la hauteur theorique, l'incidence est
nulle au nominal par definition, fuite et frottement de disque sont des
fractions fixes. Le rendement qui en sort est donc le meme a sept centiemes de
point pres pour des roues radicalement differentes : il dimensionne, il ne juge
pas. Les pertes de canal, calculees sur la geometrie, comblent ce manque.
"""

import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.hydraulics import losses
from impeller_analyzer.io import report, writer

GRID = 110


class TestGrandeursDeCanal(BaseTestCase):
    def test_longueur_croit_quand_l_aube_se_couche(self):
        """Une aube tres couchee fait un canal tres long : ds = dm / sin(beta)."""
        droite = losses.channel_length(0.04, 0.09, 45.0, 45.0)
        couchee = losses.channel_length(0.04, 0.09, 5.0, 5.0)
        self.assertGreater(couchee, 5.0 * droite)

    def test_diametre_hydraulique(self):
        """4A/P d'un canal rectangulaire, et zero si une dimension manque."""
        self.assertClose(losses.hydraulic_diameter(0.02, 0.02), 0.02, rel=1e-12)
        self.assertEqual(losses.hydraulic_diameter(0.0, 0.02), 0.0)

    def test_geometrie_incomplete_ne_produit_rien(self):
        """Mieux vaut ne rien dire que sortir un rendement invente."""
        result = losses.analyse(0.0, 0.09, 0.02, 0.01, 22.0, 25.0, 6, 10.0, 8.0, 12.0)
        self.assertEqual(result.efficiency, 0.0)
        self.assertTrue(result.warnings)


class TestPouvoirDiscriminant(BaseTestCase):
    def _losses(self, mesh, rpm=1450.0):
        path = self.path("roue.stl")
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(unit="cm", speeds=(rpm,), grid_nr=GRID, grid_nz=GRID,
                                 n_theta=240, symmetry_check=False, rotation=1))

    def test_la_roue_de_calage_retombe_sur_la_valeur_de_reference(self):
        """Par construction : c'est elle qui fixe ETA_COMPARAISON."""
        got = self._losses(synthetic.centrifugal_impeller()).channel_losses
        cible = config.ETA_H * config.ETA_VOL * config.ETA_MEC
        self.assertClose(got.efficiency, cible, abs_=0.005)

    def test_le_rendement_de_comparaison_distingue_les_formes(self):
        """Ce que le rendement de la SPEC ne sait pas faire."""
        ordinaire = self._losses(synthetic.centrifugal_impeller()).channel_losses
        toroidale = self._losses(synthetic.toroidal_propeller(n_blades=3)).channel_losses
        # Les deux doivent differer franchement, alors que le rendement SPEC ne bouge pas.
        self.assertGreater(abs(ordinaire.efficiency - toroidale.efficiency), 0.005)

    def test_independant_du_diametre_et_du_regime(self):
        """Un critere de conception ne doit juger que la forme."""
        roue = synthetic.centrifugal_impeller()
        petite = synthetic.centrifugal_impeller(
            r1=0.035 * 0.6, r2=0.090 * 0.6, b1=0.020 * 0.6, b2=0.010 * 0.6,
            eye_height=0.035 * 0.6, thickness=0.004 * 0.6,
        )
        grande = self._losses(roue).channel_losses
        reduite = self._losses(petite).channel_losses
        vite = self._losses(roue, rpm=2900.0).channel_losses
        self.assertClose(reduite.efficiency, grande.efficiency, abs_=0.01)
        self.assertClose(vite.efficiency, grande.efficiency, rel=1e-9)

    def test_le_rapport_porte_la_section(self):
        """Et il doit dire pourquoi ce chiffre n'est pas celui du tableau 2."""
        result = self._losses(synthetic.centrifugal_impeller())
        chemin = report.write_markdown(result, self.workdir, source="roue.stl")
        with open(chemin, encoding="utf-8") as flux:
            texte = flux.read()
        self.assertIn("Comparaison de conception", texte)
        self.assertIn("Rendement de comparaison", texte)
        self.assertIn("ne sert pas a", texte)


if __name__ == "__main__":
    unittest.main()
