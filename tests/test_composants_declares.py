"""Import par composants declares : l'outil verifie au lieu de deviner.

Tous les echecs constates sur le fichier reel venaient de l'inference, aucun du
calcul : nombre de pales, type de roue, rayon de moyeu, section d'entree, sens
du fluide, sens de rotation.  Ce mode-ci les remplace par des declarations et
par des mesures faites sur des pieces separees -- et surtout par des controles,
puisqu'une declaration est une entree, jamais une dispense de controle.

Les pieces de synthese sont volontairement placees a plusieurs metres de
l'origine, comme l'est le fichier reel : c'est le cas que le contrat d'export
decrit, et celui ou une erreur de repere se paie sans rien laisser paraitre.
"""

import math
import os
import unittest

from helpers import BaseTestCase

from impeller_analyzer import components, config, synthetic
from impeller_analyzer.analysis import MACHINE_PROPELLER, Options, run_components
from impeller_analyzer.io import writer
from impeller_analyzer.mesh import TriMesh

#: La piece reelle est exportee en (23.8, 17.5) m, soit a 29.5 m de l'origine
#: CAO. On reproduit l'ecart.
ORIGINE_CAO = (23.8, 17.5, 0.0)


def decale(mesh: TriMesh, offset=ORIGINE_CAO) -> TriMesh:
    """Deplace un maillage, comme le ferait un export sans recentrage."""
    return TriMesh(
        [(v[0] + offset[0], v[1] + offset[1], v[2] + offset[2]) for v in mesh.vertices],
        mesh.faces,
    )


def helicoide(sense: int = 1, beta_deg: float = 25.0, wrap_deg: float = 50.0,
              r_hub: float = 0.030, r_tip: float = 0.100, thickness: float = 0.004,
              offset=ORIGINE_CAO) -> TriMesh:
    """Une pale conventionnelle isolee : nappe helicoidale epaissie."""
    nu, nv, vertices, faces = 24, 12, [], []
    wrap = math.radians(wrap_deg)
    for i in range(nu):
        theta = wrap * i / (nu - 1)
        for j in range(nv):
            radius = r_hub + (r_tip - r_hub) * j / (nv - 1)
            z = sense * radius * math.tan(math.radians(beta_deg)) * theta
            for side in (-1, 1):
                vertices.append((
                    radius * math.cos(theta) + offset[0],
                    radius * math.sin(theta) + offset[1],
                    z + side * thickness / 2 + offset[2],
                ))
    index = lambda i, j, s: ((i * nv + j) * 2 + s)
    for i in range(nu - 1):
        for j in range(nv - 1):
            for side in (0, 1):
                a, b = index(i, j, side), index(i + 1, j, side)
                c, d = index(i + 1, j + 1, side), index(i, j + 1, side)
                faces += [(a, b, c), (a, c, d)] if side else [(a, c, b), (a, d, c)]
    return TriMesh(vertices, faces)


def tore(major: float = 0.060, minor: float = 0.012, offset=ORIGINE_CAO) -> TriMesh:
    """Un tore ferme : genre 1, la topologie d'une aube en boucle."""
    nu, nv, vertices, faces = 48, 20, [], []
    for i in range(nu):
        u = 2.0 * math.pi * i / nu
        for j in range(nv):
            v = 2.0 * math.pi * j / nv
            radius = major + minor * math.cos(v)
            vertices.append((radius * math.cos(u) + offset[0],
                             radius * math.sin(u) + offset[1],
                             minor * math.sin(v) + offset[2]))
    index = lambda i, j: (i % nu) * nv + (j % nv)
    for i in range(nu):
        for j in range(nv):
            a, b = index(i, j), index(i + 1, j)
            c, d = index(i + 1, j + 1), index(i, j + 1)
            faces += [(a, b, c), (a, c, d)]
    return TriMesh(vertices, faces)


def secteur_plan(r_in: float, r_out: float, span_deg: float, thickness: float,
                 offset=ORIGINE_CAO, n: int = 24) -> TriMesh:
    """Secteur annulaire plan et epais : une pale sans pas, qui ne s'emboite pas."""
    span = math.radians(span_deg)
    vertices = []
    for z in (-thickness / 2, thickness / 2):
        for k in range(n + 1):
            a = span * k / n
            for r in (r_in, r_out):
                vertices.append((r * math.cos(a) + offset[0], r * math.sin(a) + offset[1], z + offset[2]))
    idx = lambda layer, k, side: layer * 2 * (n + 1) + 2 * k + side
    faces = []
    for k in range(n):
        faces += [(idx(0, k, 0), idx(0, k + 1, 1), idx(0, k, 1)), (idx(0, k, 0), idx(0, k + 1, 0), idx(0, k + 1, 1))]
        faces += [(idx(1, k, 0), idx(1, k, 1), idx(1, k + 1, 1)), (idx(1, k, 0), idx(1, k + 1, 1), idx(1, k + 1, 0))]
        faces += [(idx(0, k, 1), idx(0, k + 1, 1), idx(1, k + 1, 1)), (idx(0, k, 1), idx(1, k + 1, 1), idx(1, k, 1))]
        faces += [(idx(0, k, 0), idx(1, k + 1, 0), idx(0, k + 1, 0)), (idx(0, k, 0), idx(1, k, 0), idx(1, k + 1, 0))]
    for k, sign in ((0, 1), (n, -1)):
        quad = [idx(0, k, 0), idx(0, k, 1), idx(1, k, 1), idx(1, k, 0)]
        tri = [(quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])]
        faces += tri if sign > 0 else [(a, c, b) for a, b, c in tri]
    mesh = TriMesh(vertices, faces)
    if mesh.volume() < 0.0:
        mesh = TriMesh(vertices, [(a, c, b) for a, b, c in faces])
    return mesh


class ComposantsTestCase(BaseTestCase):
    """Ecrit les pieces sur disque et les assemble."""

    def ecrire(self, mesh: TriMesh, name: str) -> str:
        """Ecrit un maillage en STL et renvoie son chemin."""
        path = self.path(name)
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return path

    def pieces(self, blade: TriMesh | None = None, **extra) -> dict:
        """Un jeu d'emplacements complet, dans le repere CAO decale."""
        paths = {
            components.SLOT_INLET: self.ecrire(
                decale(synthetic.tube(0.030, 0.100, 0.003, z_center=+0.060)), "entree.stl"),
            components.SLOT_OUTLET: self.ecrire(
                decale(synthetic.tube(0.030, 0.100, 0.003, z_center=-0.060)), "sortie.stl"),
            components.SLOT_BLADE: self.ecrire(blade or helicoide(), "pale.stl"),
        }
        paths.update(extra)
        return paths

    def declarations(self, **remplacements) -> components.Declarations:
        """Declarations par defaut du banc."""
        valeurs = dict(mode=MACHINE_PROPELLER,
                       blade_topology=components.BLADE_CONVENTIONAL, n_blades=5)
        valeurs.update(remplacements)
        return components.Declarations(**valeurs)


class TestPlansFluide(ComposantsTestCase):
    """Les sections sont mesurees sur des solides, non deduites de rayons."""

    def test_aire_mesuree_sur_la_tranche(self):
        """`aire = volume / epaisseur`, a la precision du maillage.

        C'est tout l'apport : A1 cesse d'etre suspendue a la detection de r1s
        et r1h, dont la fragilite est la cause premiere des ecarts de debit.
        """
        assembly = components.assemble(self.pieces(), self.declarations())
        attendue = math.pi * (0.100 ** 2 - 0.030 ** 2)
        self.assertClose(assembly.inlet.area, attendue, rel=0.01)
        self.assertClose(assembly.outlet.area, attendue, rel=0.01)
        self.assertClose(assembly.inlet.thickness, 0.003, rel=0.02)

    def test_normale_et_sens_debitant(self):
        """La normale sort du tenseur d'inertie, le sens debitant des centroides."""
        assembly = components.assemble(self.pieces(), self.declarations())
        self.assertClose(abs(assembly.inlet.normal[2]), 1.0, abs_=1e-6)
        self.assertSeqClose(assembly.flow_direction, (0.0, 0.0, -1.0), abs_=1e-6)

    def test_l_axe_est_situe_par_les_solides_fluide(self):
        """Le contrat fixe la direction de l'axe, pas sa position.

        Une piece exportee sans recentrage est la ou la CAO l'avait mise -- a
        29.5 m de l'origine dans le cas reel. Mesurer les rayons depuis
        l'origine y donnerait des dizaines de metres.
        """
        assembly = components.assemble(self.pieces(), self.declarations())
        self.assertSeqClose(assembly.axis_origin[:2], ORIGINE_CAO[:2], abs_=1e-6)
        inner, outer = assembly.components[components.SLOT_INLET].radial_extent(
            assembly.axis_origin)
        self.assertClose(inner, 0.030, rel=0.01)
        self.assertClose(outer, 0.100, rel=0.01)

    def test_un_tube_long_n_est_pas_une_tranche(self):
        """Un tube ne dit pas ou se trouve la section de reference."""
        paths = self.pieces()
        paths[components.SLOT_INLET] = self.ecrire(
            decale(synthetic.tube(0.030, 0.100, 0.200, z_center=0.060)), "tube.stl")
        assembly = components.assemble(paths, self.declarations())
        self.assertTrue(any("tranche mince" in w for w in assembly.warnings))


class TestControles(ComposantsTestCase):
    """§6 : une declaration est une entree, jamais une dispense de controle."""

    def resultat(self, assembly, nom: str) -> components.Check:
        """Le controle portant ce nom."""
        return next(c for c in assembly.checks if c.name == nom)

    def test_une_piece_recentree_bloque(self):
        """Le seul controle bloquant, et la raison qui le justifie.

        Sa violation ne se voit sur aucune grandeur publiee : chaque piece se
        lit correctement dans son coin, et seule leur position relative -- donc
        tout ce que le mode apporte -- est fausse.
        """
        paths = self.pieces()
        paths[components.SLOT_BLADE] = self.ecrire(
            helicoide(offset=(0.0, 0.0, 0.0)), "pale_recentree.stl")
        assembly = components.assemble(paths, self.declarations())
        bloquant = assembly.blocked()
        self.assertIsNotNone(bloquant)
        self.assertEqual(bloquant.name, "repere commun")
        self.assertIn("recentree", bloquant.detail)

    def test_un_assemblage_coherent_ne_bloque_pas(self):
        """Toutes les pieces au meme endroit : rien a signaler."""
        assembly = components.assemble(self.pieces(), self.declarations())
        self.assertIsNone(assembly.blocked())
        self.assertTrue(self.resultat(assembly, "repere commun").passed)

    def test_trop_de_pales_declarees(self):
        """Reconstruire par N rotations : un secteur plan de 100 degres ne se repete pas 8 fois."""
        paths = self.pieces(blade=secteur_plan(0.030, 0.100, 100.0, 0.010))
        assembly = components.assemble(paths, self.declarations(n_blades=8))
        check = self.resultat(assembly, "nombre de pales")
        self.assertFalse(check.passed)
        self.assertIn("se recoupe", check.detail)

    def test_pales_helicoidales_entrelacees(self):
        """Huit helicoides de 100 degres s'emboitent comme une vis a huit filets.

        L'ancien critere, par etendue angulaire, les declarait en conflit : une
        copie tournee d'un helicoide en est une copie translatee le long de
        l'axe, qui ne le recoupe pas.
        """
        paths = self.pieces(blade=helicoide(wrap_deg=100.0))
        assembly = components.assemble(paths, self.declarations(n_blades=8))
        self.assertTrue(self.resultat(assembly, "nombre de pales").passed)

    def test_solides_fluide_inverses(self):
        """La pale doit se trouver entre les deux plans, le long du debit."""
        paths = self.pieces()
        paths[components.SLOT_BLADE] = self.ecrire(
            decale(synthetic.cylinder(0.050, 0.010), (23.8, 17.5, 0.400)), "pale_loin.stl")
        assembly = components.assemble(paths, self.declarations())
        check = self.resultat(assembly, "position de la pale")
        self.assertFalse(check.passed)
        self.assertIn("inverses", check.detail)

    def test_topologie_declaree_verifiee(self):
        """Une pale toroidale doit etre de genre 1 : un tore, pas une boule."""
        toroidale = self.pieces(blade=tore())
        conforme = components.assemble(
            toroidale, self.declarations(blade_topology=components.BLADE_TOROIDAL))
        self.assertTrue(self.resultat(conforme, "topologie declaree").passed)

        contredit = components.assemble(
            toroidale, self.declarations(blade_topology=components.BLADE_CONVENTIONAL))
        check = self.resultat(contredit, "topologie declaree")
        self.assertFalse(check.passed)
        self.assertIn("anse", check.detail)

    def test_genre_d_un_tore_et_d_une_sphere(self):
        """Le genre se calcule sur la caracteristique d'Euler."""
        self.assertEqual(components.genus(tore())[0], 1)
        # Un maillage en plusieurs nappes ouvertes n'a pas de genre negatif :
        # la formule doit compter les composantes connexes.
        genre, _, _ = components.genus(helicoide())
        self.assertGreaterEqual(genre, 0)

    def test_declarations_hors_domaine(self):
        """Mode, topologie et nombre de pales sont bornes."""
        for remplacement, motif in (
            (dict(mode="turbine"), "mode inconnu"),
            (dict(blade_topology="spirale"), "topologie de pale inconnue"),
            (dict(n_blades=1), "nombre de pales declare"),
            (dict(n_blades=99), "nombre de pales declare"),
        ):
            with self.subTest(**remplacement):
                with self.assertRaises(ValueError) as capture:
                    self.declarations(**remplacement).check()
                self.assertIn(motif, str(capture.exception))

    def test_emplacement_obligatoire_manquant(self):
        """Entree, sortie et pale sont obligatoires ; coque et moyeu non."""
        paths = self.pieces()
        del paths[components.SLOT_BLADE]
        with self.assertRaises(ValueError) as capture:
            components.assemble(paths, self.declarations())
        self.assertIn("pale", str(capture.exception))


class TestSensDeRotation(ComposantsTestCase):
    """§5.1 : le sens debitant mesure leve l'indetermination."""

    def test_les_deux_sens_d_helice_donnent_des_sens_opposes(self):
        """Le critere `signe(omega) = signe(k . flux_z)` discrimine."""
        sens = []
        for helice in (+1, -1):
            assembly = components.assemble(
                self.pieces(blade=helicoide(sense=helice)), self.declarations())
            signe, _, niveau = components.rotation_from_flow(assembly)
            self.assertEqual(niveau, "high")
            self.assertNotEqual(signe, 0, "une pale helicoidale impose un sens")
            sens.append(signe)
        self.assertEqual(sens[0], -sens[1], "les deux helices doivent s'opposer")

    def test_une_boucle_n_impose_aucun_sens(self):
        """Les deux brins d'une boucle ont des pentes opposees, qui s'annulent.

        Ce n'est pas un defaut de mesure mais une propriete de la forme, et
        l'outil doit le dire au lieu d'inventer un sens.
        """
        assembly = components.assemble(
            self.pieces(blade=tore()),
            self.declarations(blade_topology=components.BLADE_TOROIDAL))
        signe, raison, niveau = components.rotation_from_flow(assembly)
        self.assertEqual(signe, 0)
        self.assertIn("brins", raison)
        self.assertEqual(niveau, "low")

    def test_refoulement_radial_lu_sur_le_recul_des_aubes(self):
        """Sortie en bande : le critere de pente axial ne s'applique pas.

        Deux aubes de recul oppose donnent des sens opposes, en confiance
        moyenne -- le critere suppose des aubes courbees vers l'arriere.
        """
        sens = []
        for recul in (+1, -1):
            aube = TriMesh(
                [(r * math.cos(recul * 4.0 * (r - 0.030)) + ORIGINE_CAO[0],
                  r * math.sin(recul * 4.0 * (r - 0.030)) + ORIGINE_CAO[1], z)
                 for r in (0.030, 0.050, 0.070, 0.090) for z in (0.0, 0.010)],
                [(0, 2, 1), (1, 2, 3)],
            )
            paths = {
                components.SLOT_INLET: self.ecrire(decale(synthetic.cylinder(0.030, 0.001, z_center=0.030)), "e.stl"),
                components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.095, 0.096, 0.012, z_center=0.005)), "s.stl"),
                components.SLOT_BLADE: self.ecrire(aube, "a.stl"),
            }
            assembly = components.assemble(paths, self.declarations(mode="pompe_carenee"))
            signe, raison, niveau = components.rotation_from_flow(assembly)
            self.assertEqual(niveau, "medium")
            self.assertIn("arriere", raison)
            sens.append(signe)
        self.assertEqual(sens[0], -sens[1])

    def test_rotation_facultative_en_mode_composants(self):
        """Sans --rotation, le sens sort de la mesure et non d'une declaration."""
        paths = self.pieces()
        result = run_components(Options(
            unit="cm", machine=MACHINE_PROPELLER, blades=5,
            blade_topology=components.BLADE_CONVENTIONAL, component_paths=paths,
        ))
        self.assertNotEqual(result.blades.rotation_sign, 0)
        self.assertEqual(result.provenance.get_source("sens_de_rotation"), "mesure")
        self.assertEqual(result.confidence.get_level("sens_de_rotation"), "high")


class TestChaineComplete(ComposantsTestCase):
    """Du chargement au rapport, sans jamais inferer ce qui est declare."""

    def test_topologie_batie_sur_les_mesures(self):
        """Rayons et sections viennent des pieces, pas d'une detection."""
        paths = self.pieces(**{components.SLOT_HUB: self.ecrire(
            decale(synthetic.cylinder(0.030, 0.100)), "moyeu.stl")})
        result = run_components(Options(
            unit="cm", machine=MACHINE_PROPELLER, blades=5, component_paths=paths,
        ))
        topology = result.topology
        self.assertClose(topology.r_1s, 0.100, rel=0.01)
        self.assertClose(topology.r_1h, 0.030, rel=0.01)
        self.assertClose(topology.area_1, math.pi * (0.100 ** 2 - 0.030 ** 2), rel=0.01)
        self.assertEqual(topology.blades.n_blades, 5)
        self.assertEqual(topology.machine_type, "axiale")

    def test_le_rapport_dit_quel_mode_a_servi(self):
        """SPEC v2 §8 : les deux modes coexistent, le rapport les distingue."""
        from impeller_analyzer.io import report as report_module

        result = run_components(Options(
            unit="cm", machine=MACHINE_PROPELLER, blades=5,
            component_paths=self.pieces(),
        ))
        self.assertEqual(result.to_dict()["mode_d_import"], "composants")
        path = report_module.write_markdown(result, self.workdir, source="assemblage")
        with open(path, encoding="utf-8") as handle:
            markdown = handle.read()
        self.assertIn("Import par composants declares", markdown)
        self.assertIn("Ce que l'outil a verifie", markdown)

    def test_une_piece_recentree_interrompt_l_analyse(self):
        """Le controle bloquant arrete tout, et le dit en premier."""
        paths = self.pieces()
        paths[components.SLOT_BLADE] = self.ecrire(
            helicoide(offset=(0.0, 0.0, 0.0)), "recentree.stl")
        result = run_components(Options(
            unit="cm", machine=MACHINE_PROPELLER, blades=5, component_paths=paths,
        ))
        self.assertIsNone(result.topology)
        self.assertIn("analyse interrompue", result.warnings[0])
        self.assertEqual(result.overall_confidence(), "low")


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()


class TestPiecesDeplaceesUneParUne(ComposantsTestCase):
    """STLOUT d'AutoCAD exige l'octant positif : pousser les pieces une par une detruit l'assemblage.

    Cas reel : chaque fichier arrivait avec le coin de sa boite a l'origine,
    l'entree a 86 mm de l'axe du corps. Aucune piece n'etant recentree *sur*
    l'origine, l'ancien controle repondait « meme origine ».
    """

    def a_l_octant(self, mesh: TriMesh) -> TriMesh:
        """Deplace la piece seule pour amener le coin de sa boite a l'origine."""
        lo, _ = mesh.bounds()
        return decale(mesh, (-lo[0], -lo[1], 0.0))

    def roue(self):
        return {
            "entree": decale(synthetic.tube(0.0, 0.035, 0.001, z_center=+0.050)),
            "sortie": decale(synthetic.tube(0.030, 0.100, 0.003, z_center=-0.060)),
            "moyeu": decale(synthetic.tube(0.010, 0.030, 0.120, z_center=0.0)),
            "pale": helicoide(),
        }

    def assembler(self, pieces: dict):
        paths = {
            components.SLOT_INLET: self.ecrire(pieces["entree"], "entree.stl"),
            components.SLOT_OUTLET: self.ecrire(pieces["sortie"], "sortie.stl"),
            components.SLOT_HUB: self.ecrire(pieces["moyeu"], "moyeu.stl"),
            components.SLOT_BLADE: self.ecrire(pieces["pale"], "pale.stl"),
        }
        return components.assemble(paths, self.declarations())

    def test_pieces_deplacees_une_par_une_bloquent(self):
        pieces = {nom: self.a_l_octant(mesh) for nom, mesh in self.roue().items()}
        blocked = self.assembler(pieces).blocked()
        self.assertIsNotNone(blocked)
        self.assertEqual(blocked.name, "repere commun")
        self.assertIn("meme axe", blocked.detail)
        self.assertIn("STLOUT", blocked.detail)

    def test_pieces_deplacees_ensemble_passent(self):
        pieces = self.roue()
        lows = [mesh.bounds()[0] for mesh in pieces.values()]
        vecteur = (-min(lo[0] for lo in lows), -min(lo[1] for lo in lows), 0.0)
        ensemble = {nom: decale(mesh, vecteur) for nom, mesh in pieces.items()}
        self.assertIsNone(self.assembler(ensemble).blocked())


class TestRoueCentrifugeFermee(ComposantsTestCase):
    """Refoulement radial et « corps » d'un seul tenant, comme sur la roue d'essai hel1.

    La sortie est une bande cylindrique : lue comme une tranche plate, sa
    hauteur passait pour son epaisseur (6 cm2 au lieu de 71). Et le moyeu,
    exporte avec le flasque, donnait pour r1h le rayon du flasque.
    """

    def pieces(self, blade=None, corps=None):
        entree = decale(synthetic.cylinder(0.035, 0.001, z_center=0.050))
        sortie = decale(synthetic.tube(0.095, 0.096, 0.012, z_center=0.006))
        if corps is None:
            corps = decale(synthetic.combine([
                synthetic.cylinder(0.015, 0.070, z_center=0.035),  # moyeu, traverse le plan d'entree
                synthetic.tube(0.040, 0.096, 0.004, z_center=0.030),  # flasque
            ]))
        return {
            components.SLOT_INLET: self.ecrire(entree, "entree.stl"),
            components.SLOT_OUTLET: self.ecrire(sortie, "sortie.stl"),
            components.SLOT_HUB: self.ecrire(corps, "corps.stl"),
            components.SLOT_BLADE: self.ecrire(
                blade or helicoide(r_hub=0.020, r_tip=0.090, offset=(ORIGINE_CAO[0], ORIGINE_CAO[1], 0.020)),
                "pale.stl"),
        }

    def assembler(self, **kwargs):
        return components.assemble(self.pieces(**kwargs),
                                   self.declarations(mode="pompe_carenee"))

    def test_bande_cylindrique_mesuree(self):
        sortie = self.assembler().outlet
        self.assertTrue(sortie.radial)
        attendu = 2.0 * math.pi * 0.0955 * 0.012
        self.assertLess(abs(sortie.area / attendu - 1.0), 0.03)
        self.assertLess(abs(sortie.height - 0.012), 1e-4)

    def test_topologie_centrifuge_lue_sur_la_sortie(self):
        from impeller_analyzer.analysis import _topology_from_components
        from impeller_analyzer.geometry import topology as topo

        topologie = _topology_from_components(self.assembler())
        self.assertEqual(topologie.machine_type, topo.CENTRIFUGAL)
        self.assertLess(abs(topologie.b_2 - 0.012), 1e-4)
        self.assertLess(abs(topologie.r_2 - 0.0955), 1e-3)

    def test_moyeu_lu_au_plan_d_entree(self):
        """Le corps porte un flasque de 96 mm : r1h reste le moyeu, 15 mm."""
        from impeller_analyzer.analysis import _topology_from_components

        topologie = _topology_from_components(self.assembler())
        self.assertLess(abs(topologie.r_1h - 0.015), 0.001)

    def test_centre_libre_au_plan_d_entree(self):
        """Sans matiere sur l'axe au plan d'entree, rien n'obstrue le centre."""
        from impeller_analyzer.analysis import _topology_from_components

        corps = decale(synthetic.tube(0.040, 0.096, 0.004, z_center=0.030))
        topologie = _topology_from_components(self.assembler(corps=corps))
        self.assertLess(topologie.r_1h, 1e-6)


def disque_sans_centre(rayon: float, epaisseur: float, z: float, n: int = 40) -> TriMesh:
    """Disque plein dont les faces sont triangulees depuis le bord : aucun sommet au centre.

    C'est ainsi qu'AutoCAD exporte un disque : lu sur ses sommets, son rayon
    interieur valait son rayon exterieur.
    """
    vertices, faces = [], []
    for level in (z - epaisseur / 2, z + epaisseur / 2):
        for k in range(n):
            a = 2.0 * math.pi * k / n
            vertices.append((rayon * math.cos(a), rayon * math.sin(a), level))
    for k in range(1, n - 1):
        faces.append((0, k + 1, k))              # dessous
        faces.append((n, n + k, n + k + 1))      # dessus
    for k in range(n):
        a, b = k, (k + 1) % n
        faces += [(a, b, n + b), (a, n + b, n + a)]
    return TriMesh(vertices, faces)


class TestPiegesDeFacettisation(ComposantsTestCase):
    """Ce que les sommets d'un maillage CAO ne disent pas de la piece."""

    def test_disque_sans_sommet_central(self):
        entree = components.load_component(
            components.SLOT_INLET, self.ecrire(decale(disque_sans_centre(0.0356, 0.001, 0.05)), "d.stl"))
        plan = components.fluid_plane(entree)
        self.assertFalse(plan.radial)
        attendu = 40 / 2 * math.sin(2 * math.pi / 40) * 0.0356 ** 2  # aire du polygone
        self.assertLess(abs(plan.area / attendu - 1.0), 0.01)
        axe = (entree.centroid[0], entree.centroid[1], 0.0)
        self.assertLess(components.inner_radius(entree, axe), 1e-9)

    def test_bande_facettisee(self):
        """Paroi d'un millimetre sur 48 facettes : l'aire de passage a 1 % pres."""
        bande = components.load_component(
            components.SLOT_OUTLET,
            self.ecrire(decale(synthetic.tube(0.0955, 0.0965, 0.0118, segments=48)), "b.stl"))
        plan = components.fluid_plane(bande)
        self.assertTrue(plan.radial)
        attendu = 2.0 * math.pi * 0.096 * 0.0118
        self.assertLess(abs(plan.area / attendu - 1.0), 0.01)
        self.assertLess(abs(plan.thickness - 0.001), 5e-5)
