"""Application locale : formulaire, analyse d'un fichier depose, telechargements."""

import json
import os
import threading
import unittest
import urllib.error
import urllib.request

from helpers import BaseTestCase

from impeller_analyzer import config, serve, synthetic
from impeller_analyzer.io import writer


class TestFormFields(BaseTestCase):
    def test_regimes(self):
        """Les regimes se saisissent separes par des espaces ou des virgules."""
        self.assertEqual(serve.parse_speeds("1000 2000 3000"), (1000.0, 2000.0, 3000.0))
        self.assertEqual(serve.parse_speeds("1450,2900"), (1450.0, 2900.0))
        self.assertEqual(serve.parse_speeds("  "), config.DEFAULT_RPM)
        with self.assertRaises(serve.BadRequest):
            serve.parse_speeds("abc")
        with self.assertRaises(serve.BadRequest):
            serve.parse_speeds("-500")

    def test_options_depuis_le_formulaire(self):
        """Chaque champ se retrouve dans les options d'analyse."""
        options = serve.options_from_query({
            "unite": ["mm"], "rpm": ["1450"], "r_aspiration": ["4,5"], "pales": ["7"],
            "beta1": ["19"], "beta2": ["28"], "altitude": ["1200"], "temperature": ["45"],
            "hauteur": ["-2.5"], "pertes": ["1.2"], "grille": ["150"],
        })
        self.assertEqual(options.unit, "mm")
        self.assertEqual(options.speeds, (1450.0,))
        self.assertClose(options.r_aspiration_cm, 4.5, rel=1e-12)  # virgule decimale acceptee
        self.assertEqual(options.blades, 7)
        self.assertClose(options.suction_height, -2.5, rel=1e-12)
        self.assertEqual(options.grid_nr, 150)
        self.assertEqual(options.grid_nz, 150)

    def test_valeurs_par_defaut_et_bornes(self):
        """Un formulaire vide reprend les valeurs de config ; la grille est bornee."""
        options = serve.options_from_query({})
        self.assertEqual(options.unit, "cm")
        self.assertEqual(options.speeds, config.DEFAULT_RPM)
        self.assertIsNone(options.blades)
        self.assertClose(options.suction_losses, config.PERTES_ASPIRATION, rel=1e-12)
        self.assertEqual(serve.options_from_query({"grille": ["5"]}).grid_nr, 40)
        self.assertEqual(serve.options_from_query({"grille": ["9000"]}).grid_nr, 400)

    def test_nombre_illisible(self):
        """Un champ numerique invalide donne un message nommant le champ."""
        with self.assertRaises(serve.BadRequest) as ctx:
            serve.options_from_query({"altitude": ["haut"]})
        self.assertIn("altitude", str(ctx.exception))

    def test_extensions(self):
        """Les formats connus passent, le .lisp est refuse avec sa consigne."""
        for name in ("roue.stl", "ROUE.STL", "a/b/roue.obj", "roue.dxf", "roue.step"):
            with self.subTest(name=name):
                self.assertIn(serve.extension_of(name), serve.ACCEPTED)
        with self.assertRaises(serve.BadRequest) as ctx:
            serve.extension_of("roue.lisp")
        self.assertIn("AutoLISP", str(ctx.exception))
        with self.assertRaises(serve.BadRequest):
            serve.extension_of("roue.xyz")


class TestStore(BaseTestCase):
    def test_rotation_des_analyses(self):
        """Seules les dernieres analyses restent sur disque."""
        store = serve.Store(self.path("magasin"), keep=2)
        tokens = [store.create()[0] for _ in range(3)]
        self.assertFalse(os.path.isdir(os.path.join(store.root, tokens[0])))
        self.assertTrue(os.path.isdir(os.path.join(store.root, tokens[2])))

    def test_resolution_stricte(self):
        """Seuls les fichiers produits, et pour un jeton connu, sont servis."""
        store = serve.Store(self.path("magasin"))
        token, directory = store.create()
        target = os.path.join(directory, serve.DOWNLOADABLE[0])
        with open(target, "w", encoding="utf-8") as handle:
            handle.write("x")
        self.assertEqual(store.resolve(token, serve.DOWNLOADABLE[0]), target)
        self.assertIsNone(store.resolve("inconnu", serve.DOWNLOADABLE[0]))
        self.assertIsNone(store.resolve(token, "config.py"))
        self.assertIsNone(store.resolve(token, "../config.py"))
        self.assertIsNone(store.resolve(token, serve.DOWNLOADABLE[1]))  # non produit


class TestAnalyseUpload(BaseTestCase):
    def test_analyse_d_un_fichier_depose(self):
        """Le fichier depose donne la meme charge utile que la vue 3D."""
        store = serve.Store(self.path("magasin"))
        payload = open(self._stl(), "rb").read()
        data = serve.analyse_upload(
            payload, {"nom": ["ma_roue.stl"], "rpm": ["1450"], "grille": ["60"]}, store
        )
        self.assertEqual(data["nom"], "ma_roue.stl")
        self.assertEqual(data["resume"]["pales"], 6)
        self.assertEqual(
            [item["nom"] for item in data["telechargements"]], list(serve.DOWNLOADABLE)
        )
        for item in data["telechargements"]:
            self.assertTrue(item["url"].startswith("/telecharger/"))

    def test_fichier_vide_et_illisible(self):
        """Un corps vide ou un contenu invalide est signale, sans trace."""
        store = serve.Store(self.path("magasin"))
        with self.assertRaises(serve.BadRequest):
            serve.analyse_upload(b"", {"nom": ["roue.stl"]}, store)
        with self.assertRaises(serve.BadRequest):
            serve.analyse_upload(b"pas un stl", {"nom": ["roue.stl"]}, store)

    def _stl(self) -> str:
        path = self.path("roue.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), path, unit_factor=config.UNIT_FACTOR)
        return path


class TestHttp(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.server, self.store = serve.build_server("127.0.0.1", 0)
        self.base = "http://127.0.0.1:%d" % self.server.server_address[1]
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.addCleanup(self._stop)
        self.stl = self.path("roue.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), self.stl, unit_factor=config.UNIT_FACTOR)

    def _stop(self):
        self.server.shutdown()
        self.server.server_close()
        self.store.clear()

    def _post(self, query: str, body: bytes):
        request = urllib.request.Request(self.base + "/analyse?" + query, data=body, method="POST")
        with urllib.request.urlopen(request) as response:
            return response.status, json.loads(response.read())

    def test_page_servie(self):
        """La racine sert la coquille de l'application, en mode serveur."""
        with urllib.request.urlopen(self.base + "/") as response:
            page = response.read().decode("utf-8")
        self.assertEqual(response.status, 200)
        self.assertIn("const SERVEUR = true;", page)
        self.assertIn("const BOOT = null;", page)
        self.assertIn('id="zone"', page)

    def test_analyse_puis_telechargement(self):
        """Un aller-retour complet : depot, analyse, recuperation du rapport."""
        status, data = self._post(
            "nom=roue.stl&rpm=1450&grille=60", open(self.stl, "rb").read()
        )
        self.assertEqual(status, 200)
        self.assertEqual(data["resume"]["type"], "centrifuge")
        url = self.base + data["telechargements"][0]["url"]
        with urllib.request.urlopen(url) as response:
            body = response.read().decode("utf-8")
        self.assertIn("Tableau 1", body)
        self.assertIn("attachment", response.headers["Content-Disposition"])

    def test_erreurs_utilisateur(self):
        """Chaque erreur revient en JSON avec un message lisible."""
        for query, fragment in (
            ("nom=roue.lisp", "AutoLISP"),
            ("nom=roue.xyz", "extension non geree"),
            ("nom=roue.stl&rpm=abc", "regime illisible"),
        ):
            with self.subTest(query=query):
                with self.assertRaises(urllib.error.HTTPError) as ctx:
                    self._post(query, b"contenu")
                self.assertEqual(ctx.exception.code, 400)
                self.assertIn(fragment, json.loads(ctx.exception.read())["erreur"])

    def test_corps_absent(self):
        """Une requete sans fichier est refusee proprement."""
        request = urllib.request.Request(self.base + "/analyse", data=b"", method="POST")
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            urllib.request.urlopen(request)
        self.assertEqual(ctx.exception.code, 400)

    def test_routes_inconnues_et_chemins_refuses(self):
        """Aucune route ne sort du magasin d'analyses."""
        for route in ("/inconnu", "/telecharger/xxx/rapport.md", "/telecharger/x/y/z"):
            with self.subTest(route=route), self.assertRaises(urllib.error.HTTPError) as ctx:
                urllib.request.urlopen(self.base + route)
            self.assertEqual(ctx.exception.code, 404)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
