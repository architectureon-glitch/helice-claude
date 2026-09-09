"""Vue 3D interactive : charge utile, classification des facettes et page produite."""

import base64
import json
import os
import re
import struct
import unittest

from helpers import BaseTestCase

from impeller_analyzer import config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer.geometry import axis as axis_mod
from impeller_analyzer.geometry import occupancy as occ_mod
from impeller_analyzer.io import report, viewer, writer

SMALL = Options(grid_nr=80, grid_nz=80, n_theta=360, symmetry_check=False, speeds=(1450.0,))


class ViewerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.source = self.path("roue.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), self.source, unit_factor=config.UNIT_FACTOR)
        self.result = run(self.source, SMALL)


class TestFaceClassification(BaseTestCase):
    def test_moyeu_et_pales_distingues(self):
        """Le sondage sous la facette separe le moyeu des pales."""
        mesh, _ = axis_mod.align_to_z(synthetic.centrifugal_impeller())
        occupancy = occ_mod.build_occupancy(mesh, nr=100, nz=100)
        classes = viewer.classify_faces(mesh, occupancy)
        self.assertEqual(len(classes), len(mesh.faces))
        counts = {value: classes.count(value) for value in set(classes)}
        self.assertIn(viewer.CLASS_SOLID, counts)
        self.assertIn(viewer.CLASS_BLADE, counts)
        # Le moyeu et le flasque dominent en surface, les pales restent bien presentes.
        self.assertGreater(counts[viewer.CLASS_BLADE], 0.10 * len(classes))

    def test_solide_de_revolution_sans_pale(self):
        """Un cylindre plein n'a que des facettes de plein."""
        mesh = synthetic.cylinder(0.05, 0.10, 120)
        occupancy = occ_mod.build_occupancy(mesh, nr=60, nz=60)
        classes = viewer.classify_faces(mesh, occupancy)
        self.assertEqual(set(classes), {viewer.CLASS_SOLID})


class TestPayload(ViewerTestCase):
    def test_geometrie_encodee(self):
        """Les positions sont en Float32 base64, une classe par facette."""
        payload = viewer.build_payload(self.result)
        raw = base64.b64decode(payload["positions"])
        floats = struct.unpack(f"<{len(raw) // 4}f", raw)
        self.assertEqual(len(floats), payload["faces"] * 9)
        self.assertEqual(len(base64.b64decode(payload["classes"])), payload["faces"])
        self.assertLess(payload["bbox"]["min"][2], payload["bbox"]["max"][2])

    def test_reperes_3d(self):
        """Axe, rayons, sens de rotation et de sortie sont fournis a la page."""
        payload = viewer.build_payload(self.result)
        identifiers = {group["id"] for group in payload["overlays"]}
        self.assertEqual(
            identifiers, {"axe", "aspiration", "sortie", "rotation", "refoulement"}
        )
        for group in payload["overlays"]:
            with self.subTest(group=group["id"]):
                self.assertEqual(len(group["points"]) % 6, 0)  # des segments de ligne
                self.assertTrue(re.fullmatch(r"#[0-9A-Fa-f]{6}", group["color"]))
                self.assertTrue(group["label"])

    def test_sens_de_rotation_transmis(self):
        """Le signe de rotation pilote l'animation de la page."""
        payload = viewer.build_payload(self.result)
        self.assertEqual(payload["resume"]["rotation_signe"], self.result.blades.rotation_sign)
        self.assertEqual(
            payload["resume"]["rotation_suggeree"], self.result.blades.observed_rotation_sign
        )
        self.assertIn(payload["resume"]["confiance"], ("haute", "moyenne", "faible"))

    def test_tableaux_et_carte(self):
        """La page embarque les memes tableaux que le rapport, et la carte meridienne."""
        payload = viewer.build_payload(self.result)
        # Neuf lignes : la nature du centre s'ajoute a r1h.
        self.assertEqual(len(payload["tables"]["geometry"]), 9)
        self.assertEqual(payload["tables"]["speeds"], ["1450 tr/min"])
        self.assertTrue(payload["carte"].startswith("data:image/png;base64,"))
        json.dumps(payload, allow_nan=False)  # doit etre serialisable tel quel

    def test_maillage_trop_lourd_decime(self):
        """Au-dela du plafond, le maillage affiche est decime."""
        original = config.VIEWER_MAX_FACES
        try:
            config.VIEWER_MAX_FACES = 2000
            payload = viewer.build_payload(self.result)
            self.assertLessEqual(payload["faces"], 2400)
        finally:
            config.VIEWER_MAX_FACES = original


class TestPage(ViewerTestCase):
    def test_page_autonome(self):
        """La page se suffit a elle-meme : aucun script externe, tout est embarque."""
        page = viewer.build_page(self.result)
        self.assertTrue(page.startswith("<!doctype html>"))
        self.assertIn("<title>Roue centrifuge a 6 pales</title>", page)
        self.assertIn('<canvas id="gl">', page)
        self.assertIn("webgl2", page)
        # Plus aucune ressource externe. La page chargeait la fonte IBM Plex
        # depuis Google Fonts, ce qui contredisait sa propre garantie
        # d'autonomie et degradait en silence des que le reseau manquait --
        # c'est-a-dire sur le poste d'atelier ou elle sert. Les piles de polices
        # systeme l'ont remplacee.
        self.assertEqual(page.count("<script"), 1)
        self.assertNotIn("<script src", page)
        self.assertEqual(len(re.findall(r'(?:href|src)="https?://', page)), 0)
        self.assertNotIn("fonts.googleapis.com", page)

    def test_page_pour_artefact(self):
        """La variante sans enveloppe garde son titre et n'a ni html ni body."""
        page = viewer.build_page(self.result, standalone=False, title="Inspecteur de roue")
        self.assertTrue(page.startswith("<title>Inspecteur de roue</title>"))
        self.assertNotIn("<!doctype", page)
        self.assertNotIn("<body", page)

    def test_ecriture_du_fichier(self):
        """Le fichier est ecrit et contient bien le maillage."""
        path = viewer.write_page(self.result, self.path("sortie"))
        self.assertTrue(os.path.isfile(path))
        self.assertGreater(os.path.getsize(path), 100_000)

    def test_integration_dans_les_sorties(self):
        """write_all produit la vue 3D, et sait s'en passer."""
        produced = report.write_all(self.result, self.path("avec"), source=self.source)
        self.assertIn("viewer", produced)
        self.assertTrue(produced["viewer"].endswith(report.VIEWER_NAME))
        sans = report.write_all(self.result, self.path("sans"), source=self.source, viewer_page=False)
        self.assertNotIn("viewer", sans)

    def test_page_ouverte_hors_serveur(self):
        """Ouverte depuis le disque, la coquille explique pourquoi elle ne peut rien faire."""
        page = viewer.build_app_page()
        self.assertNotIn("__PORT__", page)
        self.assertIn(str(config.SERVER_PORT), page)
        self.assertIn('id="hors-serveur"', page)
        self.assertIn("python3 -m impeller_analyzer.serve", page)
        # Le bloc est masque par defaut et revele par le protocole file:.
        self.assertIn('id="hors-serveur" hidden', page)
        self.assertIn('location.protocol === "file:"', page)
        self.assertIn("lancer.disabled = horsServeur", page)
        # Un echec de fetch est traduit, pas recopie tel quel.
        self.assertIn("le serveur local ne repond pas", page)

    def test_theme_clair_et_sombre(self):
        """Les trois etats de theme sont couverts par des jetons de couleur."""
        page = viewer.build_page(self.result)
        # Le theme clair est desormais le defaut, et le sombre une bascule
        # explicite : la page ne suit plus la preference du systeme.
        self.assertIn('[data-theme="sombre"]', page)
        self.assertIn('id="bascule-theme"', page)
        self.assertIn("background:var(--paper)", page)
        # La palette vient du style partage avec les figures, pas de valeurs
        # recopiees : une page claire et des figures restees aux reglages
        # d'origine donneraient un resultat incoherent.
        for couleur in (config.COULEUR_FOND, config.COULEUR_MESURE,
                        config.COULEUR_DECLARE, config.COULEUR_LIMITE):
            self.assertIn(couleur, page)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
