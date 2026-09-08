"""Bilan d'Euler : d'ou vient la hauteur, et ce qu'aucune forme d'aube ne peut changer.

`H = u2 cu2 / g` n'est pas un modele mais un theoreme, tire de la conservation
du moment cinetique. Sa decomposition en trois termes -- centrifuge, diffusion,
energie cinetique -- dit par quel mecanisme l'energie passe au fluide. Elle ne
cree rien : la somme retombe sur Euler, et c'est ce que ces tests verifient.
"""

import math
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.hydraulics import energy
from impeller_analyzer.io import report, writer

GRID = 110


class TestIdentite(BaseTestCase):
    def test_la_somme_retombe_sur_euler(self):
        """La decomposition repartit, elle n'ajoute pas."""
        for u1, u2, cm1, cm2, cu2 in (
            (14.1, 25.4, 1.60, 1.02, 5.76),
            (5.0, 20.0, 3.00, 2.50, 14.0),
            (8.0, 12.0, 4.00, 3.00, 6.0),
        ):
            with self.subTest(u2=u2, cu2=cu2):
                got = energy.budget(u1, u2, cm1, cm2, cu2)
                self.assertClose(got.total, got.euler, rel=1e-12)
                self.assertLess(got.residual, 1e-12)

    def test_le_terme_centrifuge_ignore_les_aubes(self):
        """Il ne depend que des rayons et du regime."""
        a = energy.budget(u1=14.0, u2=25.0, cm1=2.0, cm2=1.0, cu2=6.0)
        b = energy.budget(u1=14.0, u2=25.0, cm1=5.0, cm2=4.0, cu2=18.0)
        self.assertClose(a.centrifugal, b.centrifugal, rel=1e-12)
        self.assertNotAlmostEqual(a.total, b.total, places=3)

    def test_diffusion_negative_quand_l_ecoulement_relatif_accelere(self):
        """w2 > w1 : le canal ne convertit rien en pression, il consomme."""
        lent = energy.budget(u1=14.0, u2=25.0, cm1=2.0, cm2=1.5, cu2=18.0)   # w2 petit
        vif = energy.budget(u1=14.0, u2=25.0, cm1=2.0, cm2=1.5, cu2=3.0)     # w2 grand
        self.assertGreater(lent.w1, lent.w2)
        self.assertGreater(lent.diffusion, 0.0)
        self.assertLess(vif.w2, 26.0)
        self.assertLess(vif.diffusion, 0.0)

    def test_budget_nul_sans_rotation(self):
        """A l'arret, aucun terme n'a de sens : la somme est nulle."""
        got = energy.budget(0.0, 0.0, 0.0, 0.0, 0.0)
        self.assertEqual(got.total, 0.0)
        self.assertEqual(got.shares(), (0.0, 0.0, 0.0))


class TestBoutEnBout(BaseTestCase):
    def _run(self, mesh, **kwargs):
        path = self.path("roue.stl")
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(grid_nr=GRID, grid_nz=GRID, n_theta=240,
                                 symmetry_check=False, speeds=(1450.0,), rotation=1, **kwargs))

    def test_le_bilan_accompagne_l_analyse(self):
        """Et il reste exact sur une geometrie reelle."""
        result = self._run(synthetic.centrifugal_impeller())
        self.assertIsNotNone(result.energy_budget)
        self.assertLess(result.energy_budget.residual, 1e-9)
        self.assertGreater(result.energy_budget.centrifugal, 0.0)

    def test_le_rapport_porte_la_section(self):
        """Le lecteur doit voir d'ou vient sa hauteur, pas seulement combien."""
        result = self._run(synthetic.centrifugal_impeller())
        chemin = report.write_markdown(result, self.workdir, source="roue.stl")
        with open(chemin, encoding="utf-8") as flux:
            texte = flux.read()
        self.assertIn("D'ou vient la hauteur", texte)
        self.assertIn("theoreme", texte)
        self.assertIn("pas de la forme des aubes", texte)


if __name__ == "__main__":
    unittest.main()
