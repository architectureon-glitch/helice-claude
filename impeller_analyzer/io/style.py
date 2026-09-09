"""Une seule source pour la palette, l'echelle typographique et les tracés.

Le rapport, la page et les deux figures partagent ce module.  Une page claire
avec des figures qui ne partagent ni sa palette ni sa typographie donne un
resultat incoherent, et c'est le defaut qu'on oublie le plus souvent : la page
est refondue, les PNG restent aux reglages d'origine.

**Teal et violet codent, ils ne decorent pas** : le teal designe une grandeur
mesuree sur le maillage, le violet une grandeur declaree par l'utilisateur.  La
convention vaut partout, y compris dans les figures -- une courbe tracee a
partir d'un angle impose est violette, la meme a partir d'un angle lu est teal.

**La couleur ne porte jamais seule.**  Provenance et confiance se lisent aussi
par la forme -- trait plein, tirete, pointille -- et par la graisse.  Une
grandeur en confiance faible ne s'affiche jamais comme une grandeur en confiance
haute, meme quand elle est juste.
"""

from __future__ import annotations

from .. import config
from ..confidence import HIGH, LOW, MEDIUM
from ..provenance import DECLARE, DEFAUT, MESURE

Color = tuple[int, int, int]


def rgb(hexadecimal: str) -> Color:
    """Convertit `#RRGGBB` en triplet 0-255, pour le moteur de rendu PNG."""
    value = hexadecimal.lstrip("#")
    return (int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16))


#: La palette, resolue en composantes pour les figures.
FOND = rgb(config.COULEUR_FOND)
ENCRE = rgb(config.COULEUR_ENCRE)
GRILLE = rgb(config.COULEUR_GRILLE)
MESURE_RGB = rgb(config.COULEUR_MESURE)
DECLARE_RGB = rgb(config.COULEUR_DECLARE)
LIMITE = rgb(config.COULEUR_LIMITE)

#: Couleur de trace par provenance : la meme convention que la page.
PAR_PROVENANCE = {MESURE: MESURE_RGB, DECLARE: DECLARE_RGB, DEFAUT: ENCRE}

#: Style de trait par provenance, pour que la couleur ne porte pas seule.
#: `None` = trait plein ; sinon (longueur du trait, longueur du vide) en pixels.
TIRETS_PAR_PROVENANCE = {MESURE: None, DECLARE: (4, 2), DEFAUT: (1, 3)}

#: Idem pour le CSS de la page.
TIRETS_CSS = {MESURE: "none", DECLARE: "4 2", DEFAUT: "1 3"}


def sequential(t: float) -> Color:
    """Echelle sequentielle a teinte unique, du fond au teal profond.

    Remplace la rampe facon viridis.  Une echelle a teinte unique et a luminance
    monotone ne fabrique pas de frontieres : l'oeil y lit une progression, la ou
    une rampe multicolore invente des paliers que les donnees n'ont pas.  C'est
    le reproche de fond fait a jet, et il vaut aussi, plus discretement, pour
    viridis sur une grandeur qui n'a pas de seuil naturel.
    """
    t = max(0.0, min(1.0, t))
    # Interpolation en luminance perceptuelle approchee : la racine adoucit la
    # marche du clair, ou l'oeil discrimine le mieux.
    weight = t ** 0.85
    return tuple(  # type: ignore[return-value]
        round(FOND[k] + (MESURE_RGB[k] - FOND[k]) * weight) for k in range(3)
    )


#: Echelle typographique : (taille en px, interligne, graisse).
TYPOGRAPHIE = {
    "titre": (32, 1.15, 600),
    "section": (20, 1.25, 600),
    "texte": (15, 1.5, 400),
    "controle": (13, 1.4, 500),
    "nombre": (13, 1.4, 400),
    "note": (11, 1.3, 500),
}

#: Piles de polices. Non embarquees : la page reste legere et s'ouvre partout.
#: Le monospace est reserve aux chiffres, jamais aux libelles.
PILE_TEXTE = (
    '"Public Sans", "Instrument Sans", "Inter", -apple-system, '
    '"Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
)
PILE_NOMBRES = (
    '"JetBrains Mono", "SFMono-Regular", Consolas, "Liberation Mono", '
    'Menlo, monospace'
)


def css_tokens() -> str:
    """Les variables CSS de la page, derivees de la meme palette que les figures.

    Le theme sombre reste disponible en bascule ; il cesse d'etre le defaut.
    """
    taille = lambda nom: TYPOGRAPHIE[nom][0]
    return f""":root {{
  --fond: {config.COULEUR_FOND};
  --encre: {config.COULEUR_ENCRE};
  --grille: {config.COULEUR_GRILLE};
  --mesure: {config.COULEUR_MESURE};
  --declare: {config.COULEUR_DECLARE};
  --limite: {config.COULEUR_LIMITE};
  --pile-texte: {PILE_TEXTE};
  --pile-nombres: {PILE_NOMBRES};
  --t-titre: {taille('titre')}px;
  --t-section: {taille('section')}px;
  --t-texte: {taille('texte')}px;
  --t-controle: {taille('controle')}px;
  --t-nombre: {taille('nombre')}px;
  --t-note: {taille('note')}px;
}}
[data-theme="sombre"] {{
  --fond: {config.COULEUR_FOND_SOMBRE};
  --encre: {config.COULEUR_ENCRE_SOMBRE};
  --grille: {config.COULEUR_GRILLE_SOMBRE};
}}"""


#: Comment une confiance se marque, sans dependre de la seule couleur.
MARQUE_CONFIANCE = {
    HIGH: {"filet": "transparent", "parentheses": False},
    MEDIUM: {"filet": config.COULEUR_DECLARE, "parentheses": False},
    LOW: {"filet": config.COULEUR_LIMITE, "parentheses": True},
}
