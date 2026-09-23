"""Stress du domaine : ce qui n'est pas une roue, ce qui y ressemble, ce qui la brouille.

Trois familles de pieces que l'outil lisait de travers :

- des pieces qui ne sont pas une roue -- un corps parasite exporte avec elle,
  deux roues dans le meme fichier, une plaque -- et qui recevaient une hauteur et
  un debit, ou des rayons en confiance haute ;
- des roues legitimes mal orientees par l'inertie : une helice a deux pales
  elancees etait basculee de 90 degres, en confiance haute ;
- des roues brouillees : 0.1 mm de bruit sur les sommets, ce que rend un scan,
  suffisait a faire lire des aubes en boucle et a diviser beta2 par trois.
"""

from __future__ import annotations

import math
import random

from helpers import BaseTestCase

from impeller_analyzer import synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.confidence import MEDIUM
from impeller_analyzer.hydraulics import cavitation
from impeller_analyzer.io import report, writer
from impeller_analyzer.mesh import TriMesh, matvec, rotation_matrix


def shifted(mesh: TriMesh, dx: float = 0.0, dy: float = 0.0, dz: float = 0.0) -> TriMesh:
    return TriMesh([(v[0] + dx, v[1] + dy, v[2] + dz) for v in mesh.vertices], mesh.faces)


def rotated(mesh: TriMesh, axis, angle: float) -> TriMesh:
    matrix = rotation_matrix(axis, angle)
    return TriMesh([matvec(matrix, v) for v in mesh.vertices], mesh.faces)


def noisy(mesh: TriMesh, sigma: float, seed: int = 3) -> TriMesh:
    generator = random.Random(seed)
    return TriMesh(
        [tuple(c + generator.gauss(0.0, sigma) for c in v) for v in mesh.vertices], mesh.faces
    )


class AnalyseCase(BaseTestCase):
    def analyse(self, mesh: TriMesh, name: str = "piece.stl", **kwargs):
        path = writer.write_stl(mesh, self.path(name))
        kwargs.setdefault("speeds", (1450.0,))
        return run(path, Options(**kwargs))


class TestRecevabilite(AnalyseCase):
    """Une piece qui n'est pas une roue ne recoit aucune performance."""

    def test_roue_flanquee_d_un_corps_parasite(self):
        piece = synthetic.combine(
            [synthetic.centrifugal_impeller(), synthetic.cube(0.03, (0.4, 0.0, 0.0))]
        )
        result = self.analyse(piece)
        self.assertIsNotNone(result.rejection)
        self.assertEqual(result.curves, [])
        self.assertEqual(result.overall_confidence(), "low")
        self.assertTrue(result.warnings[0].startswith("analyse interrompue"))
        self.assertEqual(result.to_dict()["analyse_interrompue"], result.rejection)

    def test_deux_roues_cote_a_cote(self):
        """Symetriques d'ordre 2 entre elles : seul l'anneau de matiere les distingue d'un rotor."""
        roue = synthetic.centrifugal_impeller()
        piece = synthetic.combine([shifted(roue, -0.15), shifted(roue, 0.15)])
        result = self.analyse(piece)
        self.assertIsNotNone(result.rejection)
        self.assertIn("anneau de matiere", result.rejection)
        self.assertEqual(result.curves, [])

    def test_le_resume_annonce_l_interruption(self):
        from impeller_analyzer.cli import summarise

        piece = synthetic.combine(
            [synthetic.centrifugal_impeller(), synthetic.cube(0.03, (0.4, 0.0, 0.0))]
        )
        result = self.analyse(piece, grid_nr=100, grid_nz=100)
        self.assertIsNotNone(result.rejection)
        self.assertTrue(summarise(result, {}).startswith("ANALYSE INTERROMPUE"))

    def test_roues_valides_recevables(self):
        for name, roue in (
            ("centrifuge", synthetic.centrifugal_impeller()),
            ("flasque", synthetic.centrifugal_impeller(front_shroud=True)),
            ("axiale", synthetic.axial_impeller()),
            ("toroidale_alesee", synthetic.toroidal_propeller(n_blades=5, bore_radius=0.015)),
        ):
            with self.subTest(roue=name):
                result = self.analyse(roue, name + ".stl", grid_nr=120, grid_nz=120)
                self.assertIsNone(result.rejection)
                self.assertTrue(result.curves)

    def test_controle_levable_a_la_demande(self):
        piece = synthetic.combine(
            [synthetic.centrifugal_impeller(), synthetic.cube(0.03, (0.4, 0.0, 0.0))]
        )
        result = self.analyse(piece, symmetry_check=False, grid_nr=100, grid_nz=100)
        self.assertIsNone(result.rejection)


class TestAxeDesRouesADeuxPales(AnalyseCase):
    """L'inertie designe l'envergure d'une pale elancee ; la periodicite tranche."""

    BIPALE = dict(n_blades=2, hub_margin=0.0, r_hub=0.010, r_tip=0.20, blade_wrap_deg=20.0,
                  beta_deg=25.0, beta2_deg=35.0)

    def test_axe_z_conserve(self):
        """L'inertie designait l'envergure ; la piece reste analysee autour de Z."""
        result = self.analyse(synthetic.axial_impeller(**self.BIPALE))
        self.assertLess(result.axis.angle_to_z_deg, 1.0)
        self.assertFalse(result.axis.realigned)
        self.assertEqual(result.topology.blades.n_blades, 2)
        self.assertLess(abs(result.blades.beta2_deg - 35.0), 1.5)

    def test_bipale_couchee_n_est_pas_affirmee(self):
        """A deux pales helicoidales, trois axes sont d'ordre 2 : l'axe n'est pas verifie."""
        piece = rotated(synthetic.axial_impeller(**self.BIPALE), (1.0, 0.0, 0.0), -math.pi / 2.0)
        result = self.analyse(piece)
        self.assertEqual(result.axis.confidence, MEDIUM)
        self.assertTrue(any("plusieurs axes" in w for w in result.axis.warnings))


class TestRoueCoucheeAPlusDeDeuxPales(AnalyseCase):
    def test_axe_d_ordre_le_plus_eleve(self):
        """Exportee Y en haut, une roue de 5 pales est redressee sur son axe d'ordre 5."""
        roue = synthetic.axial_impeller(n_blades=5, beta_deg=25.0, beta2_deg=35.0)
        result = self.analyse(rotated(roue, (1.0, 0.0, 0.0), -math.pi / 2.0), grid_nr=120, grid_nz=120)
        self.assertTrue(result.axis.realigned)
        self.assertEqual(result.topology.blades.n_blades, 5)


class TestNombreDePales(AnalyseCase):
    def test_harmonique_corrigee_par_la_periodicite(self):
        """Une roue de 3 aubes en boucle rendait un pic en 9 : la piece ne tient qu'en 3."""
        result = self.analyse(synthetic.toroidal_impeller())
        self.assertEqual(result.topology.blades.n_blades, 3)
        self.assertEqual(result.topology.blades.confidence, MEDIUM)
        self.assertIsNone(result.rejection)

    def test_nombre_impose_dementi_par_la_geometrie(self):
        """Le controle porte sur le nombre retenu : impose et faux, il n'est plus en confiance haute."""
        result = self.analyse(synthetic.centrifugal_impeller(n_blades=6), blades=4,
                              grid_nr=120, grid_nz=120)
        self.assertEqual(result.topology.blades.n_blades, 4)
        self.assertEqual(result.topology.blades.confidence, MEDIUM)
        self.assertTrue(any("nombre impose" in w for w in result.topology.blades.warnings))
        # Accompagnee, pas annulee : la roue, reguliere, reste analysee.
        self.assertIsNone(result.rejection)
        self.assertEqual(result.topology.blades.periodic_order, 6)
        self.assertTrue(result.curves)


class TestMaillageBruite(AnalyseCase):
    """0.1 mm de bruit -- un scan -- ne fait plus lire d'aubes en boucle."""

    def test_pas_de_boucle_fantome(self):
        propre = self.analyse(synthetic.centrifugal_impeller(), "propre.stl")
        bruite = self.analyse(noisy(synthetic.centrifugal_impeller(), 1e-4), "bruite.stl")
        self.assertFalse(bruite.blade_loops.looped)
        self.assertLess(abs(bruite.blades.beta2_deg - propre.blades.beta2_deg), 3.0)

    def test_les_vraies_boucles_restent_lues(self):
        result = self.analyse(noisy(synthetic.toroidal_propeller(), 1e-4))
        self.assertTrue(result.blade_loops.looped)


class TestPieceSansAube(AnalyseCase):
    """Une plaque a deux lobes n'a pas d'aube : rien ne doit l'affirmer."""

    def test_plaque(self):
        result = self.analyse(synthetic.box(0.6, 0.77, 0.08), grid_nr=100, grid_nz=100)
        self.assertEqual(result.confidence.get("nombre_de_pales"), MEDIUM)
        self.assertIsNone(result.curves[0].nominal_point())
        self.assertFalse(result.speed_limit.computed)
        self.assertEqual(result.provenance.get("angles_de_pale"), "defaut")
        with open(report.write_markdown(result, self.workdir), encoding="utf-8") as handle:
            texte = handle.read()
        self.assertIn("non lus", texte)
        self.assertIn("Non calculable", texte)
        self.assertNotIn("restent fiables", texte)


class TestVitesseLimite(BaseTestCase):
    def test_zero_peut_etre_un_resultat(self):
        """NPSH disponible negatif : 0 tr/min est une reponse, pas une absence de reponse."""
        limite = cavitation.SpeedLimit(active_limit=cavitation.NPSH_LIMIT, rpm_max_npsh=0.0,
                                       rpm_max_w1s=5000.0)
        self.assertTrue(limite.computed)
        self.assertFalse(cavitation.SpeedLimit().computed)


class TestInvariances(AnalyseCase):
    """La meme roue, exportee autrement, est la meme roue ; son miroir tourne a l'envers."""

    def test_export_y_haut_unite_et_miroir(self):
        roue = synthetic.centrifugal_impeller()
        options = dict(grid_nr=120, grid_nz=120)
        reference = self.analyse(roue, "ref.stl", **options)
        y_haut = self.analyse(rotated(roue, (1.0, 0.0, 0.0), -math.pi / 2.0), "y.stl", **options)
        miroir = self.analyse(
            TriMesh([(-v[0], v[1], v[2]) for v in roue.vertices], [(a, c, b) for a, b, c in roue.faces]),
            "miroir.stl", **options,
        )
        path_mm = writer.write_stl(roue, self.path("mm.stl"), unit_factor=0.001)
        en_mm = run(path_mm, Options(unit="mm", speeds=(1450.0,), **options))
        for nom, autre in (("Y en haut", y_haut), ("mm", en_mm), ("miroir", miroir)):
            with self.subTest(export=nom):
                self.assertEqual(autre.topology.blades.n_blades, reference.topology.blades.n_blades)
                self.assertLess(abs(autre.topology.r_2 / reference.topology.r_2 - 1.0), 0.005)
                self.assertLess(abs(autre.blades.beta2_deg - reference.blades.beta2_deg), 0.5)
        self.assertEqual(
            miroir.blades.observed_rotation_sign, -reference.blades.observed_rotation_sign
        )
        self.assertNotEqual(reference.blades.observed_rotation_sign, 0)
