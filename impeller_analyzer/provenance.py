"""D'ou vient une valeur : mesuree, declaree, ou laissee par defaut.

Le niveau de confiance dit **a quel point** une grandeur est sure.  Il ne dit
pas **d'ou** elle vient, et c'est une autre question : une valeur imposee en
ligne de commande peut etre parfaitement sure sans avoir rien a voir avec le
maillage, et un lecteur qui la prend pour une lecture geometrique se trompe sur
ce que l'outil a fait.

Le cas qui a motive ce module est celui du sens de rotation, deja marque
« (impose) » dans le tableau : la meme distinction vaut pour le nombre de pales,
les angles de pale, le type de roue et le modele hydraulique, tous declarables.
Une valeur declaree n'engage que celui qui la declare ; une valeur mesuree
engage l'outil.
"""

from __future__ import annotations

MESURE = "mesure"  # extraite du maillage par l'outil
DECLARE = "declare"  # imposee par l'utilisateur en ligne de commande
DEFAUT = "defaut"  # ni l'un ni l'autre : valeur de config non contredite

SOURCES = (MESURE, DECLARE, DEFAUT)

#: Libelles courts pour le rapport et la page.
LABELS = {MESURE: "mesure", DECLARE: "declare", DEFAUT: "defaut"}


def is_source(value: object) -> bool:
    """Vrai si `value` est une provenance connue."""
    return value in SOURCES


class ProvenanceMap(dict):
    """Table `nom de grandeur -> provenance`.

    Volontairement un `dict` nu, comme `ConfidenceMap` : serialisable tel quel,
    et lu par le rapport a cote du niveau de confiance.
    """

    def set(self, name: str, source: str) -> str:
        """Enregistre `source` pour `name` apres validation, et la renvoie."""
        if not is_source(source):
            raise ValueError(f"provenance inconnue : {source!r}")
        self[name] = source
        return source

    def get_source(self, name: str, default: str = MESURE) -> str:
        """Provenance de `name`, ou `default` si la grandeur est absente."""
        return self.get(name, default)

    def declare(self, name: str, condition: bool) -> str:
        """Marque `name` declaree si `condition`, mesuree sinon."""
        return self.set(name, DECLARE if condition else MESURE)

    def label(self, name: str, default: str = MESURE) -> str:
        """Libelle affichable de la provenance de `name`."""
        return LABELS[self.get_source(name, default)]
