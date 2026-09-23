"""Angles de pale lus sur une pale isolee, et chemin de l'eau entre les parois.

Une pale importee seule se coupe par des plans z = constante, qui sont les
surfaces de courant d'une roue a refoulement radial. Les aubes de synthese ont
une loi beta(r) connue : la lecture doit la retrouver dans une veine plane, et
dire de combien elle s'en ecarte quand la veine s'incline.

Les brins d'une aube en boucle ne se lisent pas comme une aube : on dit lequel
refoule, et si l'autre travaille en parallele ou ne debouche pas -- ce que seul
le trace des parois permet de trancher. C'est la roue d'essai hel1 qui l'a
appris : son brin haut tourne dans une chambre fermee a la peripherie.
"""

import math
import unittest

from helpers import BaseTestCase
from test_composants_declares import ORIGINE_CAO, ComposantsTestCase, decale

from impeller_analyzer import components, synthetic
from impeller_analyzer.analysis import MACHINE_PROPELLER, Options, run_components
from impeller_analyzer.geometry import blade_isolated as iso
from impeller_analyzer.geometry.inclusion import SolidTester
from impeller_analyzer.geometry.meridian import MeridianWalls

AXE = ORIGINE_CAO


def aube(beta1: float, beta2: float, z0: float = 0.0, z1: float = 0.012, sense: int = 1,
         r1: float = 0.035, r2: float = 0.090, thickness: float = 0.003,
         twist: float = 0.0, pitch: float = 0.0, offset=ORIGINE_CAO):
    """Aube radiale extrudee entre z0 et z1, beta variant lineairement de beta1 a beta2.

    `sense` = +1 : l'azimut croit avec le rayon, l'aube est courbee vers
    l'arriere pour une rotation horaire (signe -1). `twist` fait varier beta1
    de 0 a `twist` degres du bas au haut de l'aube. `pitch` (rad/m) fait de
    l'aube une vis : son azimut varie de `pitch` par metre de hauteur.
    """
    nr, ns = 40, 7
    grid_a, grid_b = [], []
    for i in range(nr):
        r = r1 + (r2 - r1) * i / (nr - 1)
        half = 0.5 * thickness / r
        row_a, row_b = [], []
        for j in range(ns):
            z = z0 + (z1 - z0) * j / (ns - 1)
            theta = sense * synthetic._centrifugal_theta(r, r1, r2, beta1 + twist * j / (ns - 1), beta2)
            theta += pitch * z
            row_a.append((r * math.cos(theta - half) + offset[0], r * math.sin(theta - half) + offset[1], z))
            row_b.append((r * math.cos(theta + half) + offset[0], r * math.sin(theta + half) + offset[1], z))
        grid_a.append(row_a)
        grid_b.append(row_b)
    return synthetic._closed_box_from_grids(grid_a, grid_b)


def aube_en_veine_inclinee(beta1: float, beta2: float, r1: float = 0.035, r2: float = 0.090,
                           b1: float = 0.020, b2: float = 0.010, eye: float = 0.035):
    """Aube d'une veine qui descend de 30 degres : beta est defini sur la surface de courant."""
    def z_hi(r):
        return eye * (1.0 - (r - r1) / (r2 - r1))

    def largeur(r):
        return b1 + (b2 - b1) * (r - r1) / (r2 - r1)

    def pente(r):
        h = (r2 - r1) * 1e-4
        return ((z_hi(r + h) - 0.5 * largeur(r + h)) - (z_hi(r - h) - 0.5 * largeur(r - h))) / (2 * h)

    grid_a, grid_b = [], []
    for i in range(60):
        r = r1 + (r2 - r1) * i / 59
        theta = synthetic._centrifugal_theta(r, r1, r2, beta1, beta2, pente)
        half = 0.002 / r
        row_a, row_b = [], []
        for j in range(9):
            z = z_hi(r) - largeur(r) + largeur(r) * j / 8
            row_a.append((r * math.cos(theta - half), r * math.sin(theta - half), z))
            row_b.append((r * math.cos(theta + half), r * math.sin(theta + half), z))
        grid_a.append(row_a)
        grid_b.append(row_b)
    return synthetic._closed_box_from_grids(grid_a, grid_b)


def boucle(pitch: float = 0.0):
    """Deux brins etages de courbures opposees : brin bas vers l'arriere (rotation horaire)."""
    return synthetic.combine([aube(25.0, 35.0, 0.0, 0.010, sense=1, pitch=pitch),
                              aube(25.0, 35.0, 0.012, 0.022, sense=-1, pitch=pitch)])


def lire(mesh, rotation=-1, fente=(0.0, 0.012), entree=0.050, toroidale=False, **kw):
    return iso.read_isolated_blade(mesh, AXE, rotation, fente, entree, toroidale, **kw)


class TestLectureDesAngles(BaseTestCase):
    """La loi beta(r) de l'aube, retrouvee sur ses coupes."""

    def test_angles_connus_en_veine_plane(self):
        """Angles de bord extrapoles au bord : exacts a un demi-degre pres.

        Lus au milieu de la zone de lecture, ils glissaient vers l'angle de
        l'autre bord d'un huitieme de beta2 - beta1.
        """
        for beta1, beta2 in ((20.0, 30.0), (35.0, 50.0), (60.0, 30.0)):
            angles = lire(aube(beta1, beta2))
            self.assertEqual(angles.working.name, "aube")
            self.assertLess(abs(angles.beta1_deg - beta1), 0.5, (beta1, beta2))
            self.assertLess(abs(angles.beta2_deg - beta2), 0.5, (beta1, beta2))
            self.assertEqual(angles.confidence, "medium")

    def test_image_miroir(self):
        """L'aube miroir, tournant en sens inverse, a les memes angles."""
        droite = lire(aube(22.0, 40.0), rotation=-1)
        miroir = lire(aube(22.0, 40.0, sense=-1), rotation=+1)
        self.assertLess(abs(droite.beta1_deg - miroir.beta1_deg), 0.1)
        self.assertLess(abs(droite.beta2_deg - miroir.beta2_deg), 0.1)

    def test_rotation_inverse_courbe_vers_l_avant(self):
        """Pour l'autre sens de rotation, la meme aube est courbee vers l'avant."""
        angles = lire(aube(20.0, 30.0), rotation=+1)
        self.assertFalse(angles.working.backward)
        self.assertLess(abs(angles.beta2_deg - 150.0), 0.5)

    def test_veine_inclinee_lecture_basse(self):
        """Le plan horizontal coupe en biais la surface de courant : lecture basse, bornee."""
        mesh = aube_en_veine_inclinee(20.0, 30.0)
        lo, hi = mesh.bounds()
        angles = iso.read_isolated_blade(mesh, (0.0, 0.0, 0.0), -1, (lo[2], lo[2] + 0.010), hi[2], False)
        for lu, dessine in ((angles.beta1_deg, 20.0), (angles.beta2_deg, 30.0)):
            self.assertLess(lu - dessine, -1.0)
            self.assertGreater(lu - dessine, -5.0)

    def test_pas_de_faux_vrillage(self):
        """Une aube sans vrillage n'est pas dite vrillee.

        Le vrillage se jugeait sur tous les niveaux, y compris ceux qui coupent
        la pale loin de son bord d'attaque : 20 degres d'ecart sur une aube droite.
        """
        mesh = aube_en_veine_inclinee(40.0, 65.0)
        lo, hi = mesh.bounds()
        angles = iso.read_isolated_blade(mesh, (0.0, 0.0, 0.0), -1, (lo[2], lo[2] + 0.010), hi[2], False)
        self.assertFalse(any("vrille" in w for w in angles.warnings), angles.warnings)

    def test_vrillage_signale_et_confiance_basse(self):
        angles = lire(aube(20.0, 35.0, twist=35.0))
        self.assertTrue(any("tres vrille" in w for w in angles.warnings))
        self.assertEqual(angles.confidence, "low")
        self.assertGreater(angles.beta1_spread[1] - angles.beta1_spread[0], 25.0)


class TestBrins(BaseTestCase):
    """Une aube en boucle : deux brins, qui ne se fusionnent pas."""

    def test_deux_brins_de_courbures_opposees(self):
        angles = lire(boucle(), fente=(0.0, 0.010), toroidale=True)
        self.assertEqual([b.name for b in angles.branches], ["brin bas", "brin haut"])
        self.assertEqual(angles.working.name, "brin bas")
        self.assertTrue(angles.branches[0].backward)
        self.assertFalse(angles.branches[1].backward)
        self.assertLess(abs(angles.beta1_deg - 25.0), 0.5)
        self.assertLess(abs(angles.beta2_deg - 35.0), 0.5)
        self.assertFalse(any("toroidale" in w for w in angles.warnings))

    def test_sans_parois_la_disposition_reste_indeterminee(self):
        angles = lire(boucle(), fente=(0.0, 0.010), toroidale=True)
        self.assertEqual(angles.layout, "indeterminee")
        self.assertIn("en serie ou en parallele", angles.layout_detail)

    def test_declaration_contredite(self):
        conventionnelle = lire(boucle(), fente=(0.0, 0.010), toroidale=False)
        self.assertTrue(any("aube en boucle" in w for w in conventionnelle.warnings))
        toroidale = lire(aube(25.0, 35.0), toroidale=True)
        self.assertTrue(any("un seul sens de recul" in w for w in toroidale.warnings))


def parois(chambre_haute_fermee: bool):
    """Disque arriere, disque intermediaire et, au besoin, flasque avant ferme au bord."""
    pieces = [synthetic.cylinder(0.097, 0.003, z_center=-0.0015)]  # disque arriere
    if chambre_haute_fermee:
        pieces += [
            synthetic.tube(0.060, 0.092, 0.001, z_center=0.011),     # disque intermediaire
            synthetic.tube(0.092, 0.097, 0.0155, z_center=0.01825),  # bord ferme, z 0.0105-0.026
            synthetic.cylinder(0.092, 0.003, z_center=0.0245),       # flasque avant
        ]
    return decale(synthetic.combine(pieces))


def murs(mesh, rotating: bool = False) -> MeridianWalls:
    lo, hi = mesh.bounds()
    return MeridianWalls([SolidTester(mesh)], AXE, 0.097, (lo[2], hi[2]), rotating=rotating)


class TestCheminDeLEau(BaseTestCase):
    """Le brin qui refoule, et celui qui ne debouche pas."""

    def test_deux_brins_debouchants_en_parallele(self):
        angles = lire(boucle(), fente=(0.0, 0.010), toroidale=True,
                      walls=murs(parois(False)), outlet_radius=0.0955)
        self.assertEqual(angles.layout, "parallele")
        self.assertTrue(all(b.discharging == len(b.levels) for b in angles.branches))

    def test_brin_haut_enferme(self):
        """Chambre haute fermee au bord : le brin haut ne refoule pas (roue hel1)."""
        angles = lire(boucle(), fente=(0.0, 0.010), toroidale=True,
                      walls=murs(parois(True)), outlet_radius=0.0955)
        self.assertEqual(angles.layout, "un_seul_refoule")
        bas, haut = angles.branches
        self.assertEqual(angles.working, bas)
        self.assertEqual(bas.discharging, len(bas.levels))
        self.assertEqual(haut.discharging, 0)
        self.assertIn("en serie, en amont", angles.layout_detail)
        # Les angles publies restent ceux du brin qui refoule.
        self.assertLess(abs(angles.beta2_deg - 35.0), 0.5)

    def test_l_eau_longe_une_paroi_inclinee(self):
        """Un cone devant le bord de fuite ne ferme pas le chemin : l'eau le longe.

        Une demi-droite radiale tiree du bord de fuite le heurtait, et le brin
        bas de hel1 passait pour ne pas deboucher.
        """
        cone = decale(synthetic.revolve([(0.080, -0.010), (0.097, -0.010), (0.097, 0.001)], 90))
        walls = murs(cone)
        self.assertTrue(walls.reaches(0.085, -0.004, 0.0955, (0.002, 0.010)))
        self.assertFalse(walls.reaches(0.085, -0.004, 0.0955, (-0.012, -0.011)))

    def test_couronne_sans_aube_signalee(self):
        angles = lire(aube(20.0, 30.0), outlet_radius=0.110)
        self.assertTrue(any("couronne sans aube" in w for w in angles.warnings))


class TestAntiRetour(BaseTestCase):
    """Le brin qui ne refoule pas pousse-t-il l'eau vers celui qui refoule ?

    Sur hel1, c'est voulu : la chambre haute, fermee, et le brin qu'elle loge
    empechent l'eau de repartir vers l'oeillard. Une aube dont l'azimut varie
    avec la hauteur est une vis ; selon le sens de rotation, elle pousse l'eau
    vers le brin bas, ou la renvoie d'ou elle vient.
    """

    PAS = -60.0  # rad/m : l'azimut decroit en montant

    def lire_boucle(self, rotation: int, rotating: bool = True):
        return lire(boucle(self.PAS), rotation=rotation, fente=(0.0, 0.010), toroidale=True,
                    walls=murs(parois(True), rotating=rotating), outlet_radius=0.0955)

    def test_pas_de_vis_mesure(self):
        angles = lire(aube(25.0, 35.0, pitch=self.PAS))
        for r in (0.045, 0.060, 0.075):
            self.assertLess(abs(iso.axial_twist(angles.working.levels, r) / self.PAS - 1.0), 0.02)

    def test_pousse_vers_le_brin_bas(self):
        """Rotation horaire, azimut decroissant en montant : l'eau descend vers le brin bas."""
        angles = self.lire_boucle(rotation=-1)
        poussee = angles.axial_push
        self.assertEqual(poussee.branch, "brin haut")
        self.assertEqual(poussee.toward_working, 1.0)
        # Le passage est l'ouverture du disque intermediaire, sous r = 60 mm.
        self.assertLess(poussee.passage[1], 0.0605)
        self.assertIn("s'oppose a son retour vers l'oeillard", angles.layout_detail)
        self.assertIn("tourne en bloc", angles.layout_detail)

    def test_rotation_inverse_renvoie_l_eau(self):
        angles = self.lire_boucle(rotation=+1)
        self.assertEqual(angles.axial_push.toward_working, 0.0)
        self.assertIn("repousse l'eau vers l'oeillard", angles.layout_detail)
        self.assertEqual(angles.confidence, "low")

    def test_coque_fixe_pas_de_rotation_en_bloc(self):
        """Une coque fixe parmi les parois : la chambre ne tourne pas forcement avec la roue."""
        angles = self.lire_boucle(rotation=-1, rotating=False)
        self.assertNotIn("tourne en bloc", angles.layout_detail)


class TestPassageLibre(ComposantsTestCase):
    """Ce qu'une paroi occupe d'un solide fluide n'est pas une section de passage."""

    def pieces(self, levre: bool):
        corps = [synthetic.cylinder(0.015, 0.070, z_center=0.035),        # moyeu, traverse l'entree
                 synthetic.tube(0.040, 0.096, 0.004, z_center=0.030)]     # flasque
        if levre:
            corps.append(synthetic.tube(0.090, 0.097, 0.003, z_center=0.0105))  # levre, z 0.009-0.012
        return {
            components.SLOT_INLET: self.ecrire(decale(synthetic.cylinder(0.035, 0.001, z_center=0.050)), "e.stl"),
            components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.095, 0.096, 0.012, z_center=0.006)), "s.stl"),
            components.SLOT_HUB: self.ecrire(decale(synthetic.combine(corps)), "c.stl"),
            components.SLOT_BLADE: self.ecrire(aube(20.0, 30.0), "p.stl"),
        }

    def test_moyeu_dans_le_disque_d_entree(self):
        assembly = components.assemble(self.pieces(False), self.declarations(mode="pompe_carenee"))
        attendu = 1.0 - (0.015 / 0.035) ** 2
        self.assertLess(abs(assembly.inlet.free - attendu), 0.03)
        controle = next(c for c in assembly.checks if c.name == "passage libre")
        self.assertFalse(controle.passed)
        self.assertFalse(controle.blocking)

    def test_levre_sur_la_bande_de_sortie(self):
        """Une levre couvre 3 mm des 12 de la bande : hauteur libre 9 mm, reprise par b2."""
        from impeller_analyzer.analysis import _topology_from_components

        assembly = components.assemble(self.pieces(True), self.declarations(mode="pompe_carenee"))
        self.assertLess(abs(assembly.outlet.free - 0.75), 0.03)
        topologie = _topology_from_components(assembly)
        self.assertLess(abs(topologie.b_2 - 0.009), 4e-4)
        self.assertLess(abs(topologie.area_2 / assembly.outlet.area - 0.75), 0.03)


class TestChaineComposants(ComposantsTestCase):
    """De l'assemblage aux angles publies."""

    def pieces(self, blade=None):
        corps = decale(synthetic.combine([
            synthetic.cylinder(0.015, 0.070, z_center=0.035),
            synthetic.tube(0.040, 0.096, 0.004, z_center=0.030),
        ]))
        return {
            components.SLOT_INLET: self.ecrire(decale(synthetic.cylinder(0.035, 0.001, z_center=0.050)), "e.stl"),
            components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.095, 0.096, 0.012, z_center=0.006)), "s.stl"),
            components.SLOT_HUB: self.ecrire(corps, "c.stl"),
            components.SLOT_BLADE: self.ecrire(blade or aube(20.0, 30.0), "p.stl"),
        }

    def analyser(self, blade=None, **options):
        valeurs = dict(unit="cm", machine="pompe_carenee", blades=5,
                       component_paths=self.pieces(blade))
        valeurs.update(options)
        return run_components(Options(**valeurs))

    def test_angles_mesures_publies(self):
        result = self.analyser(rotation=-1)
        self.assertLess(abs(result.blades.beta1_deg - 20.0), 0.5)
        self.assertLess(abs(result.blades.beta2_deg - 30.0), 0.5)
        self.assertEqual(result.provenance.get_source("angles_de_pale"), "mesure")
        self.assertEqual(result.confidence.get_level("angles_de_pale"), "medium")
        lecture = result.to_dict()["angles_sur_pale_isolee"]
        self.assertTrue(lecture["parois_consultees"])
        self.assertEqual(lecture["brin_refoulant"], "aube")

    def test_angles_imposes_l_emportent(self):
        result = self.analyser(rotation=-1, beta1_deg=24.0, beta2_deg=28.0)
        self.assertEqual((result.blades.beta1_deg, result.blades.beta2_deg), (24.0, 28.0))
        self.assertEqual(result.provenance.get_source("angles_de_pale"), "declare")
        self.assertEqual(result.confidence.get_level("angles_de_pale"), "high")
        self.assertIsNotNone(result.isolated_angles, "la lecture reste publiee pour comparaison")
        self.assertTrue(any("l'emportant sur la lecture" in n for n in result.blades.notes))

    def test_sans_sens_de_rotation_rien_n_est_lu(self):
        """Boucle declaree, sens non declare : aucun sens ne se lit, aucun angle non plus."""
        result = self.analyser(blade=boucle(), blade_topology=components.BLADE_TOROIDAL)
        self.assertEqual(result.blades.rotation_sign, 0)
        self.assertIsNone(result.isolated_angles)
        self.assertEqual(result.blades.beta2_deg, 0.0)
        self.assertTrue(any("Declarez --rotation" in w for w in result.warnings))

    def test_refoulement_axial_non_lu(self):
        paths = {
            components.SLOT_INLET: self.ecrire(decale(synthetic.tube(0.030, 0.100, 0.003, z_center=0.060)), "e.stl"),
            components.SLOT_OUTLET: self.ecrire(decale(synthetic.tube(0.030, 0.100, 0.003, z_center=-0.060)), "s.stl"),
            components.SLOT_BLADE: self.ecrire(aube(20.0, 30.0, -0.010, 0.010), "p.stl"),
        }
        result = run_components(Options(unit="cm", machine=MACHINE_PROPELLER, blades=5,
                                        component_paths=paths, rotation=1))
        self.assertIsNone(result.isolated_angles)
        self.assertTrue(any("refoulement axial" in w for w in result.warnings))


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
