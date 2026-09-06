"""Permet `python -m impeller_analyzer` (analyse) et `-m impeller_analyzer.validation`."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
