"""Phase 7 - interface en ligne de commande et sorties (SPEC phase 7)."""

import io
import json
import os
import unittest
from contextlib import redirect_stderr, redirect_stdout

from helpers import BaseTestCase

from impeller_analyzer import __version__, config, synthetic
from impeller_analyzer.analysis import Options, run
from impeller_analyzer import provenance as provenance_module
from impeller_analyzer.cli import build_parser, main, options_from_args
from impeller_analyzer.io import report, writer

#: Grille reduite : les tests verifient le chainage et les sorties, pas la precision.
SMALL_GRID = ["--grille", "80", "80", "--secteurs", "360", "--sans-controle-symetrie"]


class CliTestCase(BaseTestCase):
    """Cas de test disposant d'une roue de reference sur disque."""

    def setUp(self):
        super().setUp()
        self.source = self.path("roue.stl")
        writer.write_stl(
            synthetic.centrifugal_impeller(), self.source, unit_factor=config.UNIT_FACTOR
        )
        self.out = self.path("rapport")

    def invoke(self, *extra):
        """Lance la CLI et renvoie `(code, sortie standard, sortie d'erreur)`."""
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main([self.source, "--out", self.out, *SMALL_GRID, *extra])
        return code, out.getvalue(), err.getvalue()


class TestCommandLine(CliTestCase):
    def test_analyse_complete_produit_les_trois_sorties(self):
        """La SPEC demande resultats.json, rapport.md et courbes.png."""
        code, output, _ = self.invoke()
        self.assertEqual(code, 0)
        for name in (report.JSON_NAME, report.MARKDOWN_NAME, report.CURVES_NAME, report.VIEWER_NAME):
            with self.subTest(name=name):
                self.assertTrue(os.path.isfile(os.path.join(self.out, name)))
        # La carte d'occupation est produite en plus, pour le controle visuel.
        self.assertTrue(os.path.isfile(os.path.join(self.out, report.OCCUPANCY_NAME)))
        self.assertIn("tr/min", output)
        self.assertIn("Confiance globale", output)

    def test_json_en_si_avec_confiances(self):
        """resultats.json contient les grandeurs en SI et les niveaux de confiance."""
        self.invoke()
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertEqual(data["version"], __version__)
        self.assertIn("confiance", data)
        self.assertIn(data["confiance_globale"], ("high", "medium", "low"))
        self.assertEqual(len(data["regimes"]), len(config.DEFAULT_RPM))
        curve = data["regimes"][0]
        self.assertEqual(len(curve["courbe"]), config.Q_SWEEP_POINTS)
        self.assertIn("Q_m3_s", curve["courbe"][0])
        self.assertIn("NPSHa_m", data["installation"])
        self.assertGreater(data["topologie"]["r_1s_m"], 0.0)
        self.assertLess(data["topologie"]["r_1s_m"], 1.0)  # en metres, pas en millimetres
        self.assertEqual(data["incertitude_du_modele"]["hauteur"], config.UNCERTAINTY_H)

    def test_rapport_markdown_complet(self):
        """Le rapport contient les deux tableaux et les deux encadres imposes."""
        self.invoke()
        with open(os.path.join(self.out, report.MARKDOWN_NAME), encoding="utf-8") as handle:
            text = handle.read()
        for fragment in (
            "Tableau 1 - Geometrie extraite",
            "Tableau 2 - Performances",
            "Type de roue",
            "Nombre de pales",
            "Rayon d'aspiration",
            "Sens de rotation",
            "Sens de sortie du liquide",
            "Debit nominal (m3/h)",
            "Vitesse specifique n_q",
            "NPSH disponible (m)",
            "Marge NPSHa / NPSHr",
            "Vitesse maximale sans cavitation",
            "Limite active",
            "Hypotheses d'installation",
            "Niveaux de confiance",
            "essai sur banc",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, text)
        self.assertIn(f"{config.UNCERTAINTY_H * 100.0:.0f} %", text)
        self.assertIn(f"{config.UNCERTAINTY_NPSH * 100.0:.0f} %", text)

    def test_png_bien_formes(self):
        """Les deux images sont des PNG valides."""
        self.invoke()
        for name in (report.CURVES_NAME, report.OCCUPANCY_NAME):
            with self.subTest(name=name), open(os.path.join(self.out, name), "rb") as handle:
                payload = handle.read()
            self.assertEqual(payload[:8], b"\x89PNG\r\n\x1a\n")
            self.assertGreater(len(payload), 1000)

    def test_regimes_personnalises(self):
        """--rpm choisit les regimes analyses."""
        self.invoke("--rpm", "1500", "2500")
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertEqual([c["regime_tr_min"] for c in data["regimes"]], [1500.0, 2500.0])

    def test_rayon_d_aspiration_impose(self):
        """--r-aspiration, en centimetres, prime sur la detection."""
        self.invoke("--r-aspiration", "4.5")
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertClose(data["topologie"]["r_aspiration_m"], 0.045, rel=1e-12)
        self.assertEqual(data["topologie"]["r_aspiration_source"], "utilisateur")
        with open(os.path.join(self.out, report.MARKDOWN_NAME), encoding="utf-8") as handle:
            # Le marquage « (impose) » colle au texte a cede la place a une
            # colonne de provenance, qui vaut pour toutes les grandeurs.
            markdown = handle.read()
            self.assertIn("| Rayon d'aspiration r1s (mm) | 45.00 | declare |", markdown)

    def test_repli_manuel_sur_les_angles_et_les_pales(self):
        """--beta1, --beta2 et --blades remplacent l'extraction."""
        self.invoke("--beta1", "19", "--beta2", "28", "--blades", "7")
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertClose(data["pales"]["beta1_deg"], 19.0, rel=1e-12)
        self.assertClose(data["pales"]["beta2_deg"], 28.0, rel=1e-12)
        self.assertEqual(data["topologie"]["pales"]["nombre_de_pales"], 7)
        self.assertTrue(data["pales"]["angles_imposes"])

    def test_hypotheses_d_installation(self):
        """Altitude, temperature, hauteur et pertes changent le NPSH disponible."""
        self.invoke("--altitude", "1200", "--temperature", "45", "--hauteur-aspiration", "2",
                    "--pertes-aspiration", "1.2")
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        site = data["installation"]
        self.assertEqual(site["altitude_m"], 1200.0)
        self.assertEqual(site["temperature_C"], 45.0)
        self.assertLess(site["NPSHa_m"], 9.61)
        self.assertIn("45", site["hypotheses"])

    def test_unite_du_fichier(self):
        """--unit change l'echelle du modele importe."""
        path = self.path("roue_mm.stl")
        writer.write_stl(synthetic.centrifugal_impeller(), path, unit_factor=0.001)
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main([path, "--out", self.out, "--unit", "mm", *SMALL_GRID])
        self.assertEqual(code, 0)
        with open(os.path.join(self.out, report.JSON_NAME), encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertClose(data["topologie"]["r_2_m"], 0.090, rel=0.02)

    def test_erreur_d_import_signalee_proprement(self):
        """Un fichier non geometrique sort en erreur, avec un message, sans trace."""
        path = self.path("roue.lisp")
        with open(path, "w", encoding="ascii") as handle:
            handle.write("(defun c:roue () (princ))\n")
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main([path, "--out", self.out])
        self.assertEqual(code, 2)
        self.assertIn("AutoLISP", err.getvalue())

    def test_analyseur_d_arguments(self):
        """Les options se convertissent bien en options d'analyse."""
        args = build_parser().parse_args(
            ["roue.stl", "--unit", "mm", "--rpm", "900", "--blades", "5", "--sans-reparation"]
        )
        options = options_from_args(args)
        self.assertEqual(options.unit, "mm")
        self.assertEqual(options.speeds, (900.0,))
        self.assertEqual(options.blades, 5)
        self.assertFalse(options.repair)
        self.assertTrue(options.symmetry_check)

    def test_version_et_aide(self):
        """--version et --help sortent proprement."""
        for flag in ("--version", "--help"):
            with self.subTest(flag=flag), self.assertRaises(SystemExit) as ctx:
                with redirect_stdout(io.StringIO()):
                    main([flag])
            self.assertEqual(ctx.exception.code, 0)


class TestReportPieces(CliTestCase):
    def test_tableaux_construits_depuis_le_resultat(self):
        """Les deux tableaux sont derivables directement du resultat d'analyse."""
        result = run(
            self.source,
            Options(grid_nr=80, grid_nz=80, n_theta=360, symmetry_check=False, speeds=(1450.0,)),
        )
        geometry = report.geometry_table(result)
        # Neuf lignes depuis que la nature du centre -- moyeu plein, alesage
        # traversant, ou ni l'un ni l'autre -- est publiee a cote de r1h.
        self.assertEqual(len(geometry), 9)
        for row in geometry:
            self.assertEqual(len(row), 4, "libelle, valeur, provenance, confiance")
        self.assertEqual(geometry[0][0], "Type de roue")
        for _, _, source, level in geometry:
            self.assertIn(level, report.CONFIDENCE_LABELS.values())
            self.assertIn(source, provenance_module.LABELS.values())

        header, rows = report.performance_table(result)
        self.assertEqual(header, ["1450 tr/min"])
        self.assertEqual(len(rows), 9)
        self.assertEqual(rows[0][0], "Debit nominal (m3/h)")

    def test_ecriture_dans_un_dossier_absent(self):
        """Le dossier de sortie est cree au besoin."""
        result = run(
            self.source,
            Options(grid_nr=60, grid_nz=60, n_theta=360, symmetry_check=False, speeds=(1450.0,)),
        )
        target = self.path("sous/dossier/rapport")
        produced = report.write_all(result, target, source=self.source)
        self.assertTrue(os.path.isdir(target))
        self.assertEqual(set(produced), {"json", "markdown", "curves", "occupancy", "viewer"})


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
