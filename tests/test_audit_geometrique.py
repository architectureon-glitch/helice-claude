"""Banc d'audit : chaque grandeur relue est confrontee a celle qui a ete dessinee.

Les tests par phase verifient chacun un mecanisme sur une forme.  Celui-ci fait
l'inverse : il balaye des roues **entieres**, de proportions et d'angles varies,
et confronte la sortie de la chaine complete au dessin qui les a produites.  Ce
sont les cas ou les mecanismes se contredisent qui l'interessent, et c'est ainsi
que se sont trouvees les erreurs qu'il verrouille aujourd'hui :

* une roue radialement courte (`r2/r1 = 1.8`, les proportions d'une vraie roue
  de pompe) ou les surfaces de courant n'effleuraient l'aube qu'en fragments :
  beta y montait de trente-cinq degres ;
* une roue fermee a flasque plat annoncee **toroidale**, sur un solide que le
  generateur fermait mal ;
* un maillage dechire, ou les trous donnent la meme signature qu'une aube en
  boucle : le verdict topologique doit y etre suspendu, jamais rendu.

Toutes les roues sont lues avec `rotation=1` : le sens est une entree, pas un
resultat, et un sens impose ne change aucune des grandeurs verifiees ici.
"""

import math
import random
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.io import writer
from impeller_analyzer.mesh import TriMesh

GRID = 120
N_THETA = 240


class AuditTestCase(BaseTestCase):
    """Ecrit une roue de synthese, la relit par la chaine complete."""

    def lire(self, mesh, name="roue.stl", **options):
        """Analyse complete d'un maillage de synthese."""
        path = self.path(name)
        writer.write_stl(mesh, path, unit_factor=config.UNIT_FACTOR)
        return run(path, Options(
            unit="cm", speeds=(1450.0,), grid_nr=GRID, grid_nz=GRID,
            n_theta=N_THETA, symmetry_check=False, rotation=1, **options
        ))

    def assertBeta(self, result, beta1_deg, beta2_deg):
        """Angles relus, a la tolerance de la methode qui les a fournis.

        Quand la cambrure a tenu, c'est la tolerance de la SPEC 8.1 ; quand elle
        a ete ecartee au profit des normales, c'est la bande de biais mesuree de
        cette methode-la, qui est plus large et documentee comme telle.
        """
        par_normales = result.blade_normals is not None
        tolerance = (
            config.VALID_BETA_NORMALS_DEG if par_normales else config.VALID_BETA_DEG
        )
        methode = "normales" if par_normales else "cambrure"
        self.assertClose(
            result.blades.beta1_deg, beta1_deg, abs_=tolerance, msg=f"beta1 ({methode})"
        )
        self.assertClose(
            result.blades.beta2_deg, beta2_deg, abs_=tolerance, msg=f"beta2 ({methode})"
        )


class TestBalayageCentrifuge(AuditTestCase):
    """Nombre d'aubes, rayons et angles sur un eventail de roues centrifuges."""

    def test_angles_de_pale_sur_le_domaine_valide(self):
        """Sur des aubes peu vrillees, les angles se relisent a la tolerance SPEC 8.1."""
        for beta1_deg, beta2_deg in ((15.0, 20.0), (20.0, 25.0), (22.0, 25.0)):
            with self.subTest(beta=(beta1_deg, beta2_deg)):
                result = self.lire(synthetic.centrifugal_impeller(
                    beta1_deg=beta1_deg, beta2_deg=beta2_deg
                ), name=f"b{beta1_deg:.0f}_{beta2_deg:.0f}.stl")
                self.assertBeta(result, beta1_deg, beta2_deg)

    def test_le_vrillage_est_lu_aplati(self):
        """Caracterisation du biais : la lecture rabat beta1 et beta2 vers la moyenne.

        Les angles sont pris comme la moyenne sur les dix premiers et dix
        derniers pour cent de corde ; une moyenne de fenetre rend la valeur au
        **milieu** de la fenetre.  Sur une aube dont beta varie le long de la
        corde, beta1 ressort donc trop grand et beta2 trop petit, et l'ecart
        croit avec le vrillage.  Le biais est mesure, systematique, et non
        corrige : le corriger demanderait un calage sur une loi d'aube
        particuliere, celle du generateur, qui n'est pas celle des roues
        reelles.  Ce test verrouille le sens et l'ordre de grandeur, pour qu'une
        derive future se voie.
        """
        for beta1_deg, beta2_deg in ((20.0, 35.0), (30.0, 40.0), (40.0, 65.0)):
            with self.subTest(beta=(beta1_deg, beta2_deg)):
                result = self.lire(synthetic.centrifugal_impeller(
                    beta1_deg=beta1_deg, beta2_deg=beta2_deg
                ), name=f"v{beta1_deg:.0f}_{beta2_deg:.0f}.stl")
                twist = beta2_deg - beta1_deg
                self.assertGreater(result.blades.beta1_deg, beta1_deg, "beta1 lu trop petit")
                self.assertLess(result.blades.beta2_deg, beta2_deg, "beta2 lu trop grand")
                # Le vrillage relu vaut 65 a 85 % du vrillage dessine.
                read_twist = result.blades.beta2_deg - result.blades.beta1_deg
                self.assertGreater(read_twist, 0.60 * twist)
                self.assertLess(read_twist, 0.90 * twist)

    def test_nombre_d_aubes_et_rayons(self):
        """Quatre, six ou huit aubes : le compte et les rayons ne bougent pas."""
        for n_blades in (4, 6, 8):
            with self.subTest(n_blades=n_blades):
                result = self.lire(synthetic.centrifugal_impeller(n_blades=n_blades))
                self.assertEqual(result.topology.blades.n_blades, n_blades)
                self.assertClose(result.topology.r_1s, 0.035, rel=config.VALID_GEOM_TOL)
                self.assertClose(result.topology.r_2, 0.090, rel=config.VALID_GEOM_TOL)
                self.assertEqual(result.topology.machine_type, "centrifuge")

    def test_lecture_independante_de_l_echelle(self):
        """La meme roue en petit et en grand rend les memes angles."""
        angles = []
        for scale in (0.6, 1.8):
            result = self.lire(synthetic.centrifugal_impeller(
                r1=0.035 * scale, r2=0.090 * scale, b1=0.020 * scale, b2=0.010 * scale,
                eye_height=0.035 * scale, thickness=0.004 * scale,
            ), name=f"echelle{scale}.stl")
            self.assertClose(result.topology.r_1s, 0.035 * scale, rel=config.VALID_GEOM_TOL)
            self.assertClose(result.topology.r_2, 0.090 * scale, rel=config.VALID_GEOM_TOL)
            self.assertBeta(result, 22.0, 25.0)
            angles.append((result.blades.beta1_deg, result.blades.beta2_deg))
        # Un facteur trois sur la taille laisse 1.2 degre d'ecart, du seul
        # decoupage de la grille : la lecture est bien invariante d'echelle.
        self.assertSeqClose(angles[0], angles[1], abs_=1.5)


class TestRoueCourteRadialement(AuditTestCase):
    """La regression qui a motive ce banc : `r2/r1` proche de deux.

    Sur une roue large d'ouie, les surfaces de courant coupent l'aube de biais et
    n'en rendent que des echardes : la cambrure lue sur ces fragments donnait
    quarante-cinq degres pour dix.  Le controle d'enroulement les met en defaut,
    et les angles repassent par les normales.
    """

    def test_grande_ouie(self):
        """Proportions d'une vraie roue de pompe : les angles restent justes."""
        result = self.lire(synthetic.centrifugal_impeller(
            r1=0.0927, r2=0.1670, b1=0.030, b2=0.025, eye_height=0.050,
            thickness=0.012, n_blades=5, beta1_deg=10.0, beta2_deg=20.0,
        ))
        self.assertEqual(result.topology.blades.n_blades, 5)
        self.assertClose(result.topology.r_2, 0.1670, rel=config.VALID_GEOM_TOL)
        self.assertBeta(result, 10.0, 20.0)

    def test_le_controle_d_enroulement_signale_la_contradiction(self):
        """L'enroulement mesure contredit celui qu'impliquent les angles lus."""
        result = self.lire(synthetic.centrifugal_impeller(
            r1=0.0927, r2=0.1670, b1=0.030, b2=0.025, eye_height=0.050,
            thickness=0.012, n_blades=5, beta1_deg=10.0, beta2_deg=20.0,
        ))
        coherence = result.blades.wrap_consistency
        self.assertFalse(
            config.WRAP_CONSISTENCY_MIN <= coherence <= config.WRAP_CONSISTENCY_MAX,
            f"l'enroulement ({coherence:.2f}) aurait du sortir de la bande admise",
        )
        self.assertIsNotNone(result.blade_normals)


class TestRoueFermee(AuditTestCase):
    """Flasque avant conique ou plat : une roue fermee reste une roue simple."""

    def test_flasque_conique(self):
        """Le flasque avant ne deplace ni l'ouie ni les angles."""
        result = self.lire(synthetic.centrifugal_impeller(front_shroud=True))
        self.assertClose(result.topology.r_1s, 0.035, rel=config.VALID_GEOM_TOL)
        self.assertBeta(result, 22.0, 25.0)
        self.assertFalse(result.blade_loops.looped)

    def test_flasque_plat_n_est_pas_une_boucle(self):
        """Une roue fermee a flasque plat n'est pas toroidale.

        Le fond du moyeu se pose sous le point le plus bas de la veine ; pris au
        seul rayon exterieur, il passait au-dessus de la veine a l'ouie, le
        profil se croisait et le solide sortait creux -- la carte d'occupation y
        lisait alors deux troncons en hauteur, la signature d'une aube en boucle.
        """
        result = self.lire(synthetic.centrifugal_impeller(
            front_shroud=True, flat_shroud=True
        ))
        self.assertFalse(
            result.blade_loops.looped,
            "une roue fermee a flasque plat a ete annoncee toroidale",
        )
        self.assertIsNone(result.blade_normals, "la cambrure aurait du suffire")
        self.assertBeta(result, 22.0, 25.0)


class TestMaillageDechire(AuditTestCase):
    """Un maillage troue n'a pas de topologie : aucun verdict topologique.

    Retirer des triangles au hasard n'est pas decimer -- la decimation par
    effondrement d'aretes de la SPEC 8.4 preserve l'etancheite.  C'est le cas du
    fichier abime, et ses trous coupent les troncons en hauteur exactement comme
    le ferait une boucle fermee.  Annoncer "toroidal" y serait la pire des
    sorties, puisque c'est ce mot qui met la cambrure de cote et invalide toute
    la ligne moyenne.
    """

    def dechirer(self, mesh, keep=0.8, seed=0):
        """Copie du maillage privee d'une fraction de ses triangles."""
        generator = random.Random(seed)
        return TriMesh(
            mesh.vertices,
            [face for face in mesh.faces if generator.random() < keep],
        )

    def test_verdict_de_boucle_suspendu(self):
        """La signature est relevee et dite, mais le verdict n'est pas rendu."""
        torn = self.dechirer(synthetic.centrifugal_impeller())
        result = self.lire(torn)

        self.assertFalse(result.import_report.watertight, "le maillage devait etre troue")
        loops = result.blade_loops
        self.assertFalse(loops.looped, "une roue simple dechiree a ete annoncee toroidale")
        if loops.undecided:
            self.assertTrue(loops.warnings)
            self.assertIn("suspendu", " ".join(loops.warnings))

    def test_le_maillage_troue_est_annonce(self):
        """L'utilisateur est prevenu que le fichier n'est pas etanche."""
        result = self.lire(self.dechirer(synthetic.centrifugal_impeller()))
        self.assertTrue(
            any("non etanche" in warning for warning in result.warnings),
            "un maillage troue doit etre signale comme tel",
        )


class TestHeliceAxiale(AuditTestCase):
    """Une helice axiale reste axiale, et ses aubes restent simples."""

    def test_type_et_forme(self):
        """Type de machine, compte des pales et absence de boucle."""
        result = self.lire(synthetic.axial_impeller(n_blades=3, beta_deg=25.0))
        self.assertEqual(result.topology.machine_type, "axiale")
        self.assertEqual(result.topology.blades.n_blades, 3)
        self.assertFalse(result.blade_loops.looped)


class TestDomaineDesEntrees(BaseTestCase):
    """Une entree hors du domaine des modeles doit rendre une phrase, pas une trace.

    Le controle est porte par `Options.check`, appele par `run` : la ligne de
    commande, la page web et l'appel direct passent par la meme definition du
    domaine.  Avant lui, une grille nulle sortait en `ZeroDivisionError` et deux
    secteurs azimutaux en "inf n'est pas serialisable en JSON".
    """

    def refus(self, motif, **kwargs):
        """Verifie qu'une option hors domaine est refusee, et pourquoi."""
        with self.assertRaises(ValueError) as capture:
            Options(**kwargs).check()
        self.assertIn(motif, str(capture.exception))

    def test_regime(self):
        """Un regime nul, negatif ou absurde est refuse."""
        for speeds in ((0.0,), (-1450.0,), (1e9,), (1450.0, -1.0)):
            with self.subTest(speeds=speeds):
                self.refus("regime hors domaine", speeds=speeds)

    def test_grille_et_secteurs(self):
        """Une grille ou une discretisation azimutale trop grossiere est refusee."""
        self.refus("grille trop grossiere", grid_nr=0, grid_nz=0)
        self.refus("grille trop grossiere", grid_nr=config.GRID_MIN - 1)
        self.refus("secteurs azimutaux", n_theta=2)

    def test_conditions_du_site(self):
        """Temperature et altitude restent dans le domaine des correlations."""
        self.refus("temperature hors domaine", temperature_c=500.0)
        self.refus("temperature hors domaine", temperature_c=-20.0)
        self.refus("altitude hors domaine", altitude=-1000.0)
        self.refus("altitude hors domaine", altitude=config.ALTITUDE_MAX + 1.0)

    def test_grandeurs_imposees(self):
        """Pales, angles et rayon d'aspiration imposes sont bornes."""
        self.refus("nombre de pales", blades=1)
        self.refus("beta2 hors domaine", beta2_deg=0.0)
        self.refus("beta1 hors domaine", beta1_deg=95.0)
        self.refus("rayon d'aspiration", r_aspiration_cm=-3.0)

    def test_les_valeurs_par_defaut_passent(self):
        """Les options par defaut sont dans le domaine."""
        Options().check()


if __name__ == "__main__":  # pragma: no cover - execution directe
    unittest.main()
