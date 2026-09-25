"""Une roue donnee comme on la dessine : plusieurs pieces de corps, toutes les pales, et rien d'autre.

Trois choses que la roue d'essai hel2 a demandees :

- le corps en plusieurs fichiers -- flasques et disque anti-retour --, qui
  peuvent se chevaucher a la soudure : chaque piece garde son propre test de
  volume, un chevauchement faussant la parite d'un rayon qui les traverse ;
- toutes les pales, une par fichier : elles doivent etre les copies tournees
  d'une meme pale. Sur hel2, p2 avait ete exportee 4,19 mm trop haut ;
- ni entree ni sortie fluide : une roue fermee les montre -- l'oeillard d'un
  flasque, la fente ouverte au bord.
"""

import io
import math
import os
import unittest
from contextlib import redirect_stderr, redirect_stdout

from test_composants_declares import ORIGINE_CAO, ComposantsTestCase, decale
import test_pale_isolee as pale

from impeller_analyzer import components, synthetic
from impeller_analyzer.analysis import Options, run_components
from impeller_analyzer.cli import main
from impeller_analyzer.mesh import TriMesh


def tourne(mesh: TriMesh, degres: float, dz: float = 0.0) -> TriMesh:
    """Rotation autour de l'axe de la roue, place en ORIGINE_CAO, et decalage en hauteur."""
    a = math.radians(degres)
    ca, sa = math.cos(a), math.sin(a)
    ox, oy = ORIGINE_CAO[0], ORIGINE_CAO[1]
    return TriMesh(
        [(ox + (v[0] - ox) * ca - (v[1] - oy) * sa, oy + (v[0] - ox) * sa + (v[1] - oy) * ca,
          v[2] + dz) for v in mesh.vertices],
        mesh.faces,
    )


class PiecesTestCase(ComposantsTestCase):
    def corps_en_deux(self):
        """Disque arriere et chambre haute : deux fichiers, comme le corps et l'anti-retour."""
        bas = decale(synthetic.cylinder(0.097, 0.003, z_center=-0.0015))
        haut = decale(synthetic.revolve([
            (0.060, 0.0105), (0.097, 0.0105), (0.097, 0.026), (0.040, 0.026),
            (0.040, 0.023), (0.092, 0.023), (0.092, 0.0115), (0.060, 0.0115),
        ], 180))
        return [self.ecrire(bas, "bas.stl"), self.ecrire(haut, "haut.stl")]

    def pales(self, decalages=None, azimuts=None):
        base = pale.aube(20.0, 30.0, 0.0, 0.010)
        decalages = decalages or [0.0] * 5
        azimuts = azimuts or [72.0 * k for k in range(5)]
        return [self.ecrire(tourne(base, azimuts[k], decalages[k]), f"p{k + 1}.stl") for k in range(5)]


class TestCorpsEnPlusieursPieces(PiecesTestCase):
    def test_pieces_gardees_separees(self):
        piece = components.load_components(components.SLOT_HUB, self.corps_en_deux(), unit="cm")
        self.assertEqual(len(piece.parts), 2)
        self.assertEqual(len(components._testers(piece)), 2)

    def test_un_chevauchement_ne_fait_pas_de_trou(self):
        """Deux disques qui se recouvrent : la zone commune reste une paroi."""
        a = decale(synthetic.cylinder(0.050, 0.004, z_center=0.0))
        b = decale(synthetic.tube(0.040, 0.090, 0.004, z_center=0.0))
        piece = components.load_components(
            components.SLOT_HUB, [self.ecrire(a, "a.stl"), self.ecrire(b, "b.stl")], unit="cm")
        point = (ORIGINE_CAO[0] + 0.045, ORIGINE_CAO[1], 0.0)
        self.assertTrue(any(t.contains(point) for t in components._testers(piece)))
        # Fusionnes en un seul maillage, la parite du rayon y voyait un vide.
        from impeller_analyzer.geometry.inclusion import SolidTester
        self.assertFalse(SolidTester(synthetic.combine([a, b])).contains(point))


class TestPlansDeduits(PiecesTestCase):
    def assembler(self, **paths):
        valeurs = {components.SLOT_HUB: self.corps_en_deux(),
                   components.SLOT_BLADE: self.pales()[0]}
        valeurs.update(paths)
        return components.assemble(valeurs, self.declarations(
            mode="pompe_carenee", blade_topology=components.BLADE_CONVENTIONAL))

    def test_oeillard_et_fente(self):
        assembly = self.assembler()
        self.assertEqual(sorted(assembly.derived), sorted([components.SLOT_INLET, components.SLOT_OUTLET]))
        entree, sortie = assembly.inlet, assembly.outlet
        self.assertFalse(entree.radial)
        self.assertGreater(entree.centroid[2], 0.026, "l'entree est au-dessus du flasque perce")
        self.assertLess(abs(math.sqrt(entree.area / math.pi) - 0.040), 0.002)
        self.assertTrue(sortie.radial)
        self.assertLess(abs(sortie.height - 0.0105), 0.001, "la fente va du disque arriere a la chambre haute")
        controle = next(c for c in assembly.checks if c.name == "entree et sortie")
        self.assertIn("oeillard du flasque superieur", controle.detail)
        self.assertEqual(assembly.confidence.get_level("sections"), "medium")

    def test_passage_d_arbre_n_est_pas_une_entree(self):
        """Flasque arriere perce pour l'arbre : l'entree reste l'oeillard du flasque avant."""
        bas = decale(synthetic.tube(0.006, 0.097, 0.003, z_center=-0.0015))
        haut = decale(synthetic.tube(0.040, 0.097, 0.003, z_center=0.0135))
        corps = [self.ecrire(bas, "bas_perce.stl"), self.ecrire(haut, "haut_perce.stl")]
        assembly = self.assembler(**{components.SLOT_HUB: corps})
        self.assertGreater(assembly.inlet.centroid[2], 0.015)
        controle = next(c for c in assembly.checks if c.name == "entree et sortie")
        self.assertIn("passage de l'arbre", controle.detail)

    def test_sans_corps_les_plans_restent_obligatoires(self):
        with self.assertRaises(ValueError) as ctx:
            components.assemble({components.SLOT_BLADE: self.pales()[0]},
                                self.declarations(mode="pompe_carenee"))
        self.assertIn("moyeu", str(ctx.exception))

    def test_un_plan_fourni_n_est_pas_remplace(self):
        entree = self.ecrire(decale(synthetic.cylinder(0.030, 0.001, z_center=0.040)), "e.stl")
        assembly = self.assembler(**{components.SLOT_INLET: entree})
        self.assertEqual(assembly.derived, [components.SLOT_OUTLET])


class TestPalesDistinctes(PiecesTestCase):
    def controle(self, pales):
        assembly = components.assemble(
            {components.SLOT_HUB: self.corps_en_deux(), components.SLOT_BLADE: pales},
            self.declarations(mode="pompe_carenee", n_blades=5))
        return assembly, next(c for c in assembly.checks if c.name == "pales distinctes")

    def test_copies_conformes(self):
        _, controle = self.controle(self.pales())
        self.assertTrue(controle.passed, controle.detail)

    def test_pale_exportee_trop_haut(self):
        """Le cas de hel2 : une pale deplacee a part, les autres ensemble."""
        assembly, controle = self.controle(self.pales(decalages=[0.0, 0.003, 0.0, 0.0, 0.0]))
        self.assertFalse(controle.passed)
        self.assertIn("p2.stl : decalee de +3.00 mm en hauteur", controle.detail)
        self.assertNotEqual(os.path.basename(assembly.component(components.SLOT_BLADE).path), "p2.stl")

    def test_pas_irregulier(self):
        _, controle = self.controle(self.pales(azimuts=[0.0, 72.0, 152.0, 216.0, 288.0]))
        self.assertIn("p3.stl : a 8.0 deg du pas", controle.detail)

    def test_nombre_de_pales_donne_par_les_fichiers(self):
        result = run_components(Options(
            unit="cm", machine="pompe_carenee", rotation=-1, speeds=(1450.0,),
            blade_topology=components.BLADE_CONVENTIONAL,
            component_paths={components.SLOT_HUB: self.corps_en_deux(),
                             components.SLOT_BLADE: self.pales()},
        ))
        self.assertEqual(result.assembly.declarations.n_blades, 5)
        self.assertTrue(result.curves)

    def test_ligne_de_commande(self):
        """--moyeu et --pale prennent plusieurs fichiers ; ni entree ni sortie a fournir."""
        out, err = io.StringIO(), io.StringIO()
        dossier = self.path("sortie")
        with redirect_stdout(out), redirect_stderr(err):
            code = main(["--machine", "pompe_carenee", "--rotation", "horaire", "--rpm", "1450",
                         "--topologie-pale", "conventionnelle", "--moyeu", *self.corps_en_deux(),
                         "--pale", *self.pales(), "--out", dossier, "--sans-vue3d"])
        self.assertEqual(code, 0, err.getvalue())
        # Declaree conventionnelle : l'outil, qui n'etudie que les toroidales, refuse -- mais
        # apres avoir lu les sept fichiers et deduit l'entree et la sortie.
        self.assertIn("ANALYSE INTERROMPUE : la pale est declaree conventionnelle", out.getvalue())


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
