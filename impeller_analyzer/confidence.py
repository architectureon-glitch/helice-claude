"""Niveaux de confiance et leur propagation (SPEC 1, dernier point).

Toute grandeur estimee porte un niveau `high`, `medium` ou `low`.  La regle de
propagation est la plus pessimiste : une grandeur derivee ne peut pas etre plus
sure que la moins sure de ses entrees.
"""

from __future__ import annotations

HIGH = "high"
MEDIUM = "medium"
LOW = "low"

#: Ordre croissant de fiabilite, sert de cle de comparaison.
_ORDER = {LOW: 0, MEDIUM: 1, HIGH: 2}

LEVELS = (HIGH, MEDIUM, LOW)


def is_level(value: object) -> bool:
    """Vrai si `value` est un niveau de confiance connu."""
    return value in _ORDER


def rank(level: str) -> int:
    """Rang numerique du niveau, du moins sur (0) au plus sur (2)."""
    try:
        return _ORDER[level]
    except KeyError:  # pragma: no cover - garde-fou de programmation
        raise ValueError(f"niveau de confiance inconnu : {level!r}") from None


def worst(*levels: str) -> str:
    """Niveau resultant de la combinaison de plusieurs entrees (le plus faible).

    Sans argument, renvoie `high` (element neutre de la combinaison).
    """
    best = HIGH
    for level in levels:
        if rank(level) < rank(best):
            best = level
    return best


def cap(level: str, ceiling: str) -> str:
    """Plafonne `level` a `ceiling` (utilise p.ex. pour forcer medium)."""
    return worst(level, ceiling)


def downgrade(level: str, steps: int = 1) -> str:
    """Abaisse `level` de `steps` crans, sans descendre sous `low`."""
    index = max(0, rank(level) - max(0, steps))
    for name, value in _ORDER.items():
        if value == index:
            return name
    return LOW  # pragma: no cover - inatteignable, _ORDER couvre 0..2


def from_flag(condition: bool, if_true: str = HIGH, if_false: str = LOW) -> str:
    """Niveau choisi par un test booleen, sucre syntaxique lisible."""
    return if_true if condition else if_false


class ConfidenceMap(dict):
    """Table `nom de grandeur -> niveau`, avec propagation integree.

    C'est un `dict` ordinaire (donc serialisable tel quel en JSON) augmente de
    quelques operations de combinaison.
    """

    def set(self, name: str, level: str) -> str:
        """Enregistre `level` pour `name` apres validation, et le renvoie."""
        if not is_level(level):
            raise ValueError(f"niveau de confiance inconnu : {level!r}")
        self[name] = level
        return level

    def get_level(self, name: str, default: str = LOW) -> str:
        """Niveau de `name`, ou `default` si la grandeur est absente."""
        return self.get(name, default)

    def combine(self, *names: str) -> str:
        """Niveau resultant de la combinaison des grandeurs citees."""
        return worst(*(self.get_level(name) for name in names))

    def cap_all(self, ceiling: str, names: object = None) -> None:
        """Plafonne a `ceiling` toutes les grandeurs (ou celles citees)."""
        keys = list(self) if names is None else list(names)
        for key in keys:
            if key in self:
                self[key] = cap(self[key], ceiling)

    def overall(self) -> str:
        """Niveau global : le plus faible de la table (`high` si vide)."""
        return worst(*self.values())
