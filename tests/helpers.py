"""Outillage commun aux tests : chemin du paquet, comparaisons approchees, fixtures.

Les tests sont ecrits avec `unittest` de la bibliotheque standard : ils sont donc
executables par `python -m unittest discover -s tests` comme par `pytest`, sans
dependance supplementaire.
"""

from __future__ import annotations

import math
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


class BaseTestCase(unittest.TestCase):
    """Cas de test avec repertoire temporaire et comparaisons tolerantes."""

    def setUp(self) -> None:
        self.workdir = tempfile.mkdtemp(prefix="impeller-test-")
        self.addCleanup(shutil.rmtree, self.workdir, True)

    def path(self, name: str) -> str:
        """Chemin d'un fichier de travail dans le repertoire temporaire."""
        return os.path.join(self.workdir, name)

    def assertClose(self, got, expected, rel=None, abs_=None, msg=""):
        """Egalite approchee, par tolerance relative et/ou absolue."""
        if rel is None and abs_ is None:
            rel = 1e-9
        tolerance = 0.0
        if rel is not None:
            tolerance = max(tolerance, abs(expected) * rel)
        if abs_ is not None:
            tolerance = max(tolerance, abs_)
        if not math.isfinite(got):
            self.fail(f"valeur non finie : {got!r} {msg}")
        self.assertLessEqual(
            abs(got - expected),
            tolerance,
            f"{got!r} != {expected!r} (tolerance {tolerance:g}) {msg}",
        )

    def assertSeqClose(self, got, expected, rel=None, abs_=None, msg=""):
        """Egalite approchee terme a terme de deux sequences."""
        self.assertEqual(len(got), len(expected), f"longueurs differentes {msg}")
        for index, (a, b) in enumerate(zip(got, expected)):
            self.assertClose(a, b, rel=rel, abs_=abs_, msg=f"[{index}] {msg}")
