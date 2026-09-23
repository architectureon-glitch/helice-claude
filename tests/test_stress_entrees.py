"""Stress des entrees : arguments hors domaine, sorties impossibles, page et serveur.

Chaque cas produisait une trace de pile, un processus tue par manque de memoire,
ou une valeur acceptee qui n'avait pas de sens physique. La regle : une phrase
qui nomme l'argument et le domaine attendu, jamais une trace.
"""

from __future__ import annotations

import http.client
import io
import json
import math
import threading
from contextlib import redirect_stderr

from helpers import BaseTestCase

from impeller_analyzer import cli, config, serve, synthetic
from impeller_analyzer.analysis import Options
from impeller_analyzer.io import loader, viewer, writer


class TestBornesDesOptions(BaseTestCase):
    def assertRefused(self, fragment: str, **kwargs) -> None:
        with self.assertRaises(ValueError) as caught:
            Options(**kwargs).check()
        self.assertIn(fragment, str(caught.exception))

    def test_valeurs_non_finies(self):
        for key in ("suction_height", "suction_losses", "r_aspiration_cm"):
            with self.subTest(option=key):
                with self.assertRaises(ValueError):
                    Options(**{key: math.nan}).check()

    def test_pertes_negatives(self):
        self.assertRefused("ajouterait de l'energie", suction_losses=-5.0)

    def test_hauteur_d_aspiration_absurde(self):
        self.assertRefused("faute de saisie", suction_height=1e6)

    def test_trop_de_pales(self):
        self.assertRefused("entre 2 et 24", blades=30)

    def test_memoire(self):
        self.assertRefused("secteurs", n_theta=100_000_000)
        self.assertRefused("trop fine", grid_nr=100_000)

    def test_valeurs_legitimes(self):
        Options(suction_height=-8.0, suction_losses=0.0, blades=24, n_theta=7200).check()


class TestLigneDeCommande(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.roue = writer.write_stl(synthetic.axial_impeller(), self.path("roue.stl"))

    def main(self, *args: str) -> tuple[int, str]:
        buffer = io.StringIO()
        with redirect_stderr(buffer):
            code = cli.main([self.roue, "--quiet", "--sans-vue3d", "--grille", "60", "60", *args])
        return code, buffer.getvalue()

    def test_dossier_de_sortie_impossible(self):
        occupe = self.path("un_fichier")
        with open(occupe, "w", encoding="utf-8") as handle:
            handle.write("x")
        code, message = self.main("--out", occupe)
        self.assertEqual(code, 2)
        self.assertIn("erreur d'ecriture", message)
        self.assertNotIn("Traceback", message)

    def test_rayon_d_aspiration_hors_de_la_roue(self):
        code, message = self.main("--out", self.path("o"), "--r-aspiration", "1000")
        self.assertEqual(code, 2)
        self.assertIn("plus grand que la roue", message)

    def test_rayon_d_aspiration_dans_le_moyeu(self):
        code, message = self.main("--out", self.path("o"), "--r-aspiration", "0.1")
        self.assertEqual(code, 2)
        self.assertIn("dans le moyeu", message)


class TestCoordonneesLointaines(BaseTestCase):
    def test_simple_precision_signalee(self):
        roue = synthetic.centrifugal_impeller()
        loin = type(roue)([(v[0] + 2000.0, v[1] + 1500.0, v[2]) for v in roue.vertices], roue.faces)
        _, report = loader.load_mesh(writer.write_stl(loin, self.path("loin.stl")))
        self.assertTrue(any("simple precision" in w for w in report.warnings))

    def test_piece_proche_non_signalee(self):
        _, report = loader.load_mesh(writer.write_stl(synthetic.centrifugal_impeller(), self.path("p.stl")))
        self.assertFalse(any("simple precision" in w for w in report.warnings))


class TestPageAutonome(BaseTestCase):
    def test_json_inerte_dans_le_script(self):
        charge = json.dumps({"nom": "a</script><!--<b>&c"}, ensure_ascii=False)
        sur = viewer._script_safe(charge)
        self.assertNotIn("<", sur)
        self.assertNotIn(">", sur)
        self.assertEqual(json.loads(sur), json.loads(charge))


class TestServeurLocalDurci(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.server, self.store = serve.build_server("127.0.0.1", 0)
        self.port = self.server.server_address[1]
        thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(self.store.clear)
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)

    def request(self, method: str, path: str, headers: dict, body: bytes = b"") -> int:
        connection = http.client.HTTPConnection("127.0.0.1", self.port, timeout=30)
        connection.putrequest(method, path, skip_host=True)
        for key, value in headers.items():
            connection.putheader(key, value)
        connection.putheader("Content-Length", str(len(body)))
        connection.endheaders(body)
        status = connection.getresponse().status
        connection.close()
        return status

    def test_page_servie_a_elle_meme(self):
        self.assertEqual(self.request("GET", "/", {"Host": f"127.0.0.1:{self.port}"}), 200)

    def test_rebinding_dns_refuse(self):
        self.assertEqual(self.request("GET", "/", {"Host": f"attaquant.example:{self.port}"}), 403)

    def test_envoi_depuis_une_autre_page_refuse(self):
        status = self.request(
            "POST", "/analyse?nom=roue.stl",
            {"Host": f"127.0.0.1:{self.port}", "Origin": "http://attaquant.example"}, b"solid x",
        )
        self.assertEqual(status, 403)

    def test_nombre_non_fini_dans_le_formulaire(self):
        status = self.request(
            "POST", "/analyse?nom=roue.stl&pales=inf",
            {"Host": f"127.0.0.1:{self.port}", "Origin": f"http://127.0.0.1:{self.port}"},
            b"solid x",
        )
        self.assertEqual(status, 400)


class TestFiletDeDerniereInstance(BaseTestCase):
    def test_erreur_interne_annoncee_comme_telle(self):
        buffer = io.StringIO()
        with redirect_stderr(buffer):
            try:
                raise RuntimeError("cas imprevu")
            except RuntimeError:
                code = cli._internal_error()
        self.assertEqual(code, 3)
        self.assertIn("erreur interne de l'outil", buffer.getvalue())
        self.assertIn("cas imprevu", buffer.getvalue())


class TestConfig(BaseTestCase):
    def test_bornes_coherentes(self):
        self.assertLess(config.GRID_MIN, config.GRID_MAX)
        self.assertLess(config.N_THETA_MIN, config.N_THETA_MAX)
        self.assertLess(config.SYM_TOL, config.SYM_REJECT)
        self.assertLess(config.FLOAT32_WARN, config.FLOAT32_LOW)


class TestComposantsPieges(BaseTestCase):
    """Mode composants : les erreurs de saisie qui passaient sans bruit."""

    def setUp(self):
        super().setUp()
        from test_composants_declares import decale, helicoide

        from impeller_analyzer import components

        self.components = components
        ecrire = lambda mesh, name: writer.write_stl(mesh, self.path(name), unit_factor=config.UNIT_FACTOR)
        self.entree = ecrire(decale(synthetic.tube(0.030, 0.100, 0.003, z_center=+0.060)), "entree.stl")
        self.sortie = ecrire(decale(synthetic.tube(0.030, 0.100, 0.003, z_center=-0.060)), "sortie.stl")
        self.pale = ecrire(helicoide(), "pale.stl")

    def assemble(self, inlet: str, outlet: str, blade: str):
        c = self.components
        return c.assemble(
            {c.SLOT_INLET: inlet, c.SLOT_OUTLET: outlet, c.SLOT_BLADE: blade},
            c.Declarations(mode="helice_libre", blade_topology=c.BLADE_CONVENTIONAL, n_blades=5),
        )

    def test_tranche_mince_sans_fausse_alerte(self):
        """Une tranche mince plafonne a 2 : un seuil de 2.0 alertait a chaque analyse."""
        assembly = self.assemble(self.entree, self.sortie, self.pale)
        self.assertGreater(assembly.inlet.anisotropy, config.SLICE_ANISOTROPY_MIN)
        self.assertFalse(any("direction privilegiee" in w for w in assembly.warnings))

    def test_meme_fichier_en_entree_et_en_sortie(self):
        assembly = self.assemble(self.entree, self.entree, self.pale)
        blocked = assembly.blocked()
        self.assertIsNotNone(blocked)
        self.assertEqual(blocked.name, "plans distincts")

    def test_emplacement_nomme_dans_l_erreur(self):
        vide = self.path("vide.stl")
        open(vide, "wb").close()
        with self.assertRaises(loader.ImportError_) as caught:
            self.assemble(self.entree, self.sortie, vide)
        self.assertIn("« pale »", str(caught.exception))
