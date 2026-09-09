"""Rendu PNG sans dependance externe (zlib de la bibliotheque standard).

Deux sorties sont produites par l'outil :

* la carte d'occupation `f(r, z)`, indispensable pour verifier visuellement tout
  ce que les phases suivantes en deduisent (SPEC 10, point 2) ;
* les courbes H-Q, rendement-Q et NPSHr-Q, une couleur par regime (SPEC 7).

Le texte est trace avec une fonte matricielle 5x7 ; les minuscules sont rendues
par leur capitale, ce qui suffit aux etiquettes d'axes et aux legendes.
"""

from __future__ import annotations

import math
import struct
import unicodedata
import zlib
from typing import Sequence

from .. import config
from . import style

Color = tuple[int, int, int]

# La palette vient de `style`, partagee avec la page : une page claire et des
# figures restees aux reglages d'origine donnent un resultat incoherent, et
# c'est le defaut qu'on oublie le plus souvent.
WHITE: Color = style.FOND
BLACK: Color = style.ENCRE
GREY: Color = style.GRILLE
LIGHT: Color = style.GRILLE

#: Regimes : une seule teinte, celle du mesure, eclaircie du plus lent au plus
#: rapide. Six couleurs categorielles pour une grandeur **ordonnee** feraient
#: lire des familles la ou il n'y a qu'une progression de vitesse.
def series_color(index: int, total: int) -> Color:
    """Teinte du `index`-ieme regime, du plus clair au plus fonce."""
    if total <= 1:
        return style.MESURE_RGB
    return style.sequential(0.42 + 0.58 * index / (total - 1))


SERIES_COLORS: tuple[Color, ...] = tuple(series_color(i, 6) for i in range(6))

#: Fonte 5x7, colonnes de gauche a droite, bit 0 = ligne du haut.
_FONT: dict[str, tuple[int, int, int, int, int]] = {
    " ": (0x00, 0x00, 0x00, 0x00, 0x00), "!": (0x00, 0x00, 0x5F, 0x00, 0x00),
    '"': (0x00, 0x07, 0x00, 0x07, 0x00), "#": (0x14, 0x7F, 0x14, 0x7F, 0x14),
    "$": (0x24, 0x2A, 0x7F, 0x2A, 0x12), "%": (0x23, 0x13, 0x08, 0x64, 0x62),
    "&": (0x36, 0x49, 0x55, 0x22, 0x50), "'": (0x00, 0x05, 0x03, 0x00, 0x00),
    "(": (0x00, 0x1C, 0x22, 0x41, 0x00), ")": (0x00, 0x41, 0x22, 0x1C, 0x00),
    "*": (0x14, 0x08, 0x3E, 0x08, 0x14), "+": (0x08, 0x08, 0x3E, 0x08, 0x08),
    ",": (0x00, 0x50, 0x30, 0x00, 0x00), "-": (0x08, 0x08, 0x08, 0x08, 0x08),
    ".": (0x00, 0x60, 0x60, 0x00, 0x00), "/": (0x20, 0x10, 0x08, 0x04, 0x02),
    "0": (0x3E, 0x51, 0x49, 0x45, 0x3E), "1": (0x00, 0x42, 0x7F, 0x40, 0x00),
    "2": (0x42, 0x61, 0x51, 0x49, 0x46), "3": (0x21, 0x41, 0x45, 0x4B, 0x31),
    "4": (0x18, 0x14, 0x12, 0x7F, 0x10), "5": (0x27, 0x45, 0x45, 0x45, 0x39),
    "6": (0x3C, 0x4A, 0x49, 0x49, 0x30), "7": (0x01, 0x71, 0x09, 0x05, 0x03),
    "8": (0x36, 0x49, 0x49, 0x49, 0x36), "9": (0x06, 0x49, 0x49, 0x29, 0x1E),
    ":": (0x00, 0x36, 0x36, 0x00, 0x00), ";": (0x00, 0x56, 0x36, 0x00, 0x00),
    "<": (0x00, 0x08, 0x14, 0x22, 0x41), "=": (0x14, 0x14, 0x14, 0x14, 0x14),
    ">": (0x41, 0x22, 0x14, 0x08, 0x00), "?": (0x02, 0x01, 0x51, 0x09, 0x06),
    "@": (0x32, 0x49, 0x79, 0x41, 0x3E), "A": (0x7E, 0x11, 0x11, 0x11, 0x7E),
    "B": (0x7F, 0x49, 0x49, 0x49, 0x36), "C": (0x3E, 0x41, 0x41, 0x41, 0x22),
    "D": (0x7F, 0x41, 0x41, 0x22, 0x1C), "E": (0x7F, 0x49, 0x49, 0x49, 0x41),
    "F": (0x7F, 0x09, 0x09, 0x01, 0x01), "G": (0x3E, 0x41, 0x41, 0x51, 0x32),
    "H": (0x7F, 0x08, 0x08, 0x08, 0x7F), "I": (0x00, 0x41, 0x7F, 0x41, 0x00),
    "J": (0x20, 0x40, 0x41, 0x3F, 0x01), "K": (0x7F, 0x08, 0x14, 0x22, 0x41),
    "L": (0x7F, 0x40, 0x40, 0x40, 0x40), "M": (0x7F, 0x02, 0x04, 0x02, 0x7F),
    "N": (0x7F, 0x04, 0x08, 0x10, 0x7F), "O": (0x3E, 0x41, 0x41, 0x41, 0x3E),
    "P": (0x7F, 0x09, 0x09, 0x09, 0x06), "Q": (0x3E, 0x41, 0x51, 0x21, 0x5E),
    "R": (0x7F, 0x09, 0x19, 0x29, 0x46), "S": (0x46, 0x49, 0x49, 0x49, 0x31),
    "T": (0x01, 0x01, 0x7F, 0x01, 0x01), "U": (0x3F, 0x40, 0x40, 0x40, 0x3F),
    "V": (0x1F, 0x20, 0x40, 0x20, 0x1F), "W": (0x7F, 0x20, 0x18, 0x20, 0x7F),
    "X": (0x63, 0x14, 0x08, 0x14, 0x63), "Y": (0x03, 0x04, 0x78, 0x04, 0x03),
    "Z": (0x61, 0x51, 0x49, 0x45, 0x43), "[": (0x00, 0x00, 0x7F, 0x41, 0x41),
    "\\": (0x02, 0x04, 0x08, 0x10, 0x20), "]": (0x41, 0x41, 0x7F, 0x00, 0x00),
    "^": (0x04, 0x02, 0x01, 0x02, 0x04), "_": (0x40, 0x40, 0x40, 0x40, 0x40),
    "|": (0x00, 0x00, 0x77, 0x00, 0x00), "~": (0x02, 0x01, 0x02, 0x04, 0x02),
    "°": (0x00, 0x06, 0x09, 0x09, 0x06),
}

CHAR_WIDTH = 6  # 5 colonnes de glyphe + 1 d'espacement
CHAR_HEIGHT = 7

#: Bas de casse, en lignes lisibles plutot qu'en masques binaires.
#:
#: La fonte d'origine n'avait que des capitales, et `_ascii` remontait donc tout
#: le texte des figures en majuscules.  Un titre de figure en capitales est
#: exactement ce qu'un libelle ne doit pas etre : il se lit moins vite, et il
#: crie.  Sept lignes de cinq colonnes, `#` pour un pixel allume.
_MINUSCULES = {
    "a": (".....", ".....", ".###.", "....#", ".####", "#...#", ".####"),
    "b": ("#....", "#....", "####.", "#...#", "#...#", "#...#", "####."),
    "c": (".....", ".....", ".####", "#....", "#....", "#....", ".####"),
    "d": ("....#", "....#", ".####", "#...#", "#...#", "#...#", ".####"),
    "e": (".....", ".....", ".###.", "#...#", "#####", "#....", ".####"),
    "f": ("..##.", ".#..#", ".#...", "###..", ".#...", ".#...", ".#..."),
    "g": (".....", ".....", ".####", "#...#", ".####", "....#", ".###."),
    "h": ("#....", "#....", "####.", "#...#", "#...#", "#...#", "#...#"),
    "i": ("..#..", ".....", ".##..", "..#..", "..#..", "..#..", ".###."),
    "j": ("...#.", ".....", "..##.", "...#.", "...#.", "#..#.", ".##.."),
    "k": ("#....", "#....", "#..#.", "#.#..", "##...", "#.#..", "#..#."),
    "l": (".##..", "..#..", "..#..", "..#..", "..#..", "..#..", ".###."),
    "m": (".....", ".....", "##.#.", "#.#.#", "#.#.#", "#...#", "#...#"),
    "n": (".....", ".....", "####.", "#...#", "#...#", "#...#", "#...#"),
    "o": (".....", ".....", ".###.", "#...#", "#...#", "#...#", ".###."),
    "p": (".....", ".....", "####.", "#...#", "####.", "#....", "#...."),
    "q": (".....", ".....", ".####", "#...#", ".####", "....#", "....#"),
    "r": (".....", ".....", "#.##.", "##..#", "#....", "#....", "#...."),
    "s": (".....", ".....", ".####", "#....", ".###.", "....#", "####."),
    "t": ("..#..", "..#..", "#####", "..#..", "..#..", "..#.#", "...#."),
    "u": (".....", ".....", "#...#", "#...#", "#...#", "#..##", ".##.#"),
    "v": (".....", ".....", "#...#", "#...#", "#...#", ".#.#.", "..#.."),
    "w": (".....", ".....", "#...#", "#...#", "#.#.#", "#.#.#", ".#.#."),
    "x": (".....", ".....", "#...#", ".#.#.", "..#..", ".#.#.", "#...#"),
    "y": (".....", ".....", "#...#", "#...#", ".####", "....#", ".###."),
    "z": (".....", ".....", "#####", "...#.", "..#..", ".#...", "#####"),
}


def _colonnes(lignes: tuple[str, ...]) -> tuple[int, int, int, int, int]:
    """Convertit sept lignes de cinq caracteres en cinq masques de colonne."""
    return tuple(  # type: ignore[return-value]
        sum(1 << rang for rang, ligne in enumerate(lignes) if ligne[colonne] == "#")
        for colonne in range(5)
    )


_FONT.update({char: _colonnes(lignes) for char, lignes in _MINUSCULES.items()})


def _ascii(text: str) -> str:
    """Retire les accents et remplace les caracteres absents de la fonte."""
    folded = unicodedata.normalize("NFD", text)
    out = []
    for char in folded:
        if unicodedata.combining(char):
            continue
        upper = char.upper()
        if char in _FONT:
            out.append(char)
        elif upper in _FONT:
            out.append(upper)
        else:
            out.append("?")
    return "".join(out)


class Canvas:
    """Image RVB en memoire, avec les primitives de trace necessaires."""

    def __init__(self, width: int, height: int, background: Color = WHITE):
        self.width = width
        self.height = height
        self.pixels = bytearray(bytes(background) * (width * height))

    def set_pixel(self, x: int, y: int, color: Color) -> None:
        """Allume un pixel (hors cadre : ignore)."""
        if 0 <= x < self.width and 0 <= y < self.height:
            offset = (y * self.width + x) * 3
            self.pixels[offset : offset + 3] = bytes(color)

    def fill_rect(self, x0: int, y0: int, x1: int, y1: int, color: Color) -> None:
        """Remplit un rectangle (bornes incluses, ecretees au cadre)."""
        x0, x1 = max(0, min(x0, x1)), min(self.width - 1, max(x0, x1))
        y0, y1 = max(0, min(y0, y1)), min(self.height - 1, max(y0, y1))
        if x1 < x0 or y1 < y0:
            return
        row = bytes(color) * (x1 - x0 + 1)
        for y in range(y0, y1 + 1):
            offset = (y * self.width + x0) * 3
            self.pixels[offset : offset + len(row)] = row

    def line(self, x0: float, y0: float, x1: float, y1: float, color: Color, width: int = 1) -> None:
        """Segment par l'algorithme de Bresenham, epaissi par un carre."""
        xa, ya, xb, yb = int(round(x0)), int(round(y0)), int(round(x1)), int(round(y1))
        dx = abs(xb - xa)
        dy = -abs(yb - ya)
        sx = 1 if xa < xb else -1
        sy = 1 if ya < yb else -1
        error = dx + dy
        half = width // 2
        while True:
            if width <= 1:
                self.set_pixel(xa, ya, color)
            else:
                self.fill_rect(xa - half, ya - half, xa - half + width - 1, ya - half + width - 1, color)
            if xa == xb and ya == yb:
                break
            doubled = 2 * error
            if doubled >= dy:
                error += dy
                xa += sx
            if doubled <= dx:
                error += dx
                ya += sy

    def text(self, x: int, y: int, message: str, color: Color = BLACK, scale: int = 1) -> None:
        """Ecrit `message` avec son coin haut-gauche en (x, y)."""
        cursor = x
        for char in _ascii(message):
            glyph = _FONT.get(char, _FONT[" "])
            for column, bits in enumerate(glyph):
                for row in range(CHAR_HEIGHT):
                    if bits & (1 << row):
                        if scale == 1:
                            self.set_pixel(cursor + column, y + row, color)
                        else:
                            self.fill_rect(
                                cursor + column * scale, y + row * scale,
                                cursor + column * scale + scale - 1, y + row * scale + scale - 1,
                                color,
                            )
            cursor += CHAR_WIDTH * scale

    def text_right(self, x: int, y: int, message: str, color: Color = BLACK, scale: int = 1) -> None:
        """Ecrit `message` cale a droite de l'abscisse `x`."""
        self.text(x - len(_ascii(message)) * CHAR_WIDTH * scale, y, message, color, scale)

    def write_png(self, path: str) -> str:
        """Encode et ecrit l'image en PNG (couleur vraie, 8 bits par canal)."""
        with open(path, "wb") as handle:
            handle.write(self.encode_png())
        return path

    def encode_png(self) -> bytes:
        """Encode l'image en PNG et renvoie les octets."""
        raw = bytearray()
        stride = self.width * 3
        for y in range(self.height):
            raw.append(0)  # filtre "None" pour chaque ligne
            raw.extend(self.pixels[y * stride : (y + 1) * stride])

        def chunk(tag: bytes, payload: bytes) -> bytes:
            return (
                struct.pack(">I", len(payload))
                + tag
                + payload
                + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF)
            )

        pixels_per_metre = int(round(config.PLOT_DPI / 0.0254))
        data = (
            b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", self.width, self.height, 8, 2, 0, 0, 0))
            + chunk(b"pHYs", struct.pack(">IIB", pixels_per_metre, pixels_per_metre, 1))
            + chunk(b"IDAT", zlib.compress(bytes(raw), 6))
            + chunk(b"IEND", b"")
        )
        return data


# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------
def viridis_like(t: float) -> Color:
    """Echelle de la carte d'occupation : teinte unique, du fond au teal.

    Le nom est reste pour ne pas casser les appels ; la rampe, elle, a change.
    Une echelle multicolore fabrique des frontieres que les donnees n'ont pas --
    reproche de fond fait a jet, et vrai plus discretement de viridis sur une
    grandeur sans seuil naturel. `f` est une fraction qui va de zero a un : une
    progression, pas des categories.
    """
    return style.sequential(t)


def occupancy_scale(value: float) -> float:
    """Position d'affichage de `f` sur l'echelle de couleur.

    Echelle affine par morceaux calee sur les deux seuils de la SPEC : la bande
    de pales, dont `f` est souvent tres faible pour des aubes minces tres
    vrillees, occupe l'essentiel de la dynamique de couleur au lieu de se
    confondre avec le vide.
    """
    if value <= config.F_VIDE:
        return config.PLOT_BAND_VIDE * value / config.F_VIDE if config.F_VIDE > 0.0 else 0.0
    if value >= config.F_SOLIDE:
        span = 1.0 - config.F_SOLIDE
        ratio = (value - config.F_SOLIDE) / span if span > 0.0 else 1.0
        return config.PLOT_BAND_PLEIN + (1.0 - config.PLOT_BAND_PLEIN) * ratio
    ratio = (value - config.F_VIDE) / (config.F_SOLIDE - config.F_VIDE)
    return config.PLOT_BAND_VIDE + (config.PLOT_BAND_PLEIN - config.PLOT_BAND_VIDE) * ratio


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _nice_ticks(low: float, high: float, count: int) -> list[float]:
    """Graduations regulieres et lisibles couvrant [low, high]."""
    if high <= low:
        return [low]
    raw = (high - low) / max(1, count)
    magnitude = 10.0 ** math.floor(math.log10(raw))
    for multiple in (1.0, 2.0, 2.5, 5.0, 10.0):
        step = multiple * magnitude
        if step >= raw:
            break
    first = math.ceil(low / step) * step
    ticks = []
    value = first
    while value <= high * (1.0 + 1e-9):
        ticks.append(value)
        value += step
    return ticks or [low, high]


def _format_tick(value: float) -> str:
    """Etiquette courte d'une graduation."""
    if value == 0.0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1000.0 or magnitude < 0.01:
        return f"{value:.1e}"
    if magnitude >= 100.0:
        return f"{value:.0f}"
    if magnitude >= 10.0:
        return f"{value:.1f}"
    return f"{value:.2f}"


def occupancy_canvas(occupancy_map, title: str = "Carte d'occupation f(r, z)") -> "Canvas":
    """Trace la carte `f(r, z)` : abscisse le rayon, ordonnee la hauteur.

    L'echelle de couleur est calee sur les deux seuils de la SPEC (voir
    `occupancy_scale`), de sorte que moyeu, zone de pales et veine fluide se
    distinguent d'un coup d'oeil : c'est l'outil de mise au point visuelle de
    toutes les phases suivantes.  Renvoie le canevas, que l'appelant ecrit ou
    encode selon ses besoins.
    """
    width, height = config.PLOT_WIDTH_PX, config.PLOT_HEIGHT_PX
    margin = config.PLOT_MARGIN_PX
    canvas = Canvas(width, height)
    plot_x0, plot_y0 = margin + 20, margin
    plot_x1, plot_y1 = width - margin - 90, height - margin
    nr, nz = occupancy_map.nr, occupancy_map.nz

    for iz in range(nz):
        # L'axe z est trace vers le haut : la premiere tranche est en bas.
        y_top = plot_y1 - int(round((iz + 1) * (plot_y1 - plot_y0) / nz))
        y_bottom = plot_y1 - int(round(iz * (plot_y1 - plot_y0) / nz)) - 1
        row = occupancy_map.f[iz]
        for ir in range(nr):
            value = row[ir]
            x_left = plot_x0 + int(round(ir * (plot_x1 - plot_x0) / nr))
            x_right = plot_x0 + int(round((ir + 1) * (plot_x1 - plot_x0) / nr)) - 1
            canvas.fill_rect(x_left, y_top, max(x_left, x_right), max(y_top, y_bottom), viridis_like(occupancy_scale(value)))

    canvas.line(plot_x0, plot_y0, plot_x0, plot_y1, BLACK)
    canvas.line(plot_x0, plot_y1, plot_x1, plot_y1, BLACK)
    canvas.text(plot_x0, margin - 30, title, BLACK, 2)
    canvas.text(plot_x0, plot_y1 + 22, "rayon r (mm)", BLACK)
    canvas.text(10, plot_y0, "z (mm)", BLACK)

    for value in _nice_ticks(0.0, occupancy_map.r_max * 1000.0, config.PLOT_TICKS):
        x = plot_x0 + (plot_x1 - plot_x0) * value / (occupancy_map.r_max * config.MM_PER_M)
        canvas.line(x, plot_y1, x, plot_y1 + 4, BLACK)
        canvas.text(int(x) - 12, plot_y1 + 8, _format_tick(value), BLACK)
    z_lo, z_hi = occupancy_map.z_min * config.MM_PER_M, occupancy_map.z_max * config.MM_PER_M
    for value in _nice_ticks(z_lo, z_hi, config.PLOT_TICKS):
        y = plot_y1 - (plot_y1 - plot_y0) * (value - z_lo) / (z_hi - z_lo)
        canvas.line(plot_x0 - 4, y, plot_x0, y, BLACK)
        canvas.text_right(plot_x0 - 6, int(y) - 3, _format_tick(value), BLACK)

    # Barre de couleur et reperage des deux seuils.
    bar_x0, bar_x1 = plot_x1 + 30, plot_x1 + 55
    for pixel in range(plot_y0, plot_y1 + 1):
        t = (plot_y1 - pixel) / max(1, plot_y1 - plot_y0)
        canvas.fill_rect(bar_x0, pixel, bar_x1, pixel, viridis_like(t))
    canvas.text(bar_x0 - 4, plot_y0 - 14, "f", BLACK)
    for level, label in ((0.0, "0"), (config.F_VIDE, f"{config.F_VIDE:g} vide"),
                         (config.F_SOLIDE, f"{config.F_SOLIDE:g} plein"), (1.0, "1")):
        y = plot_y1 - (plot_y1 - plot_y0) * occupancy_scale(level)
        canvas.line(bar_x0 - 5, y, bar_x1, y, BLACK)
        canvas.text(bar_x1 + 4, int(y) - 3, label, BLACK)
    return canvas


def occupancy_png(occupancy_map, path: str, title: str = "Carte d'occupation f(r, z)") -> str:
    """Trace la carte f(r, z) et l'ecrit en PNG."""
    return occupancy_canvas(occupancy_map, title).write_png(path)


def curves_png(
    path: str,
    panels: Sequence[dict],
    title: str = "",
) -> str:
    """Trace plusieurs panneaux empiles, chacun avec ses series.

    Chaque panneau est un dictionnaire `{"titre", "x_label", "y_label",
    "series": [{"label", "x", "y", "color"}], "markers": [(x, y, color)]}`.
    """
    count = max(1, len(panels))
    width = config.PLOT_WIDTH_PX
    height = config.PLOT_HEIGHT_PX * count // 2 + 60
    canvas = Canvas(width, height)
    margin = config.PLOT_MARGIN_PX
    top_offset = 34 if title else 6
    if title:
        canvas.text(margin, 10, title, BLACK, 2)
    panel_height = (height - top_offset - 10) // count

    for index, panel in enumerate(panels):
        base = top_offset + index * panel_height
        x0, y0 = margin + 20, base + 22
        x1, y1 = width - margin - 150, base + panel_height - 40
        series = panel.get("series", [])
        xs = [v for s in series for v in s["x"]]
        ys = [v for s in series for v in s["y"]]
        if not xs or not ys:
            continue
        x_lo, x_hi = min(xs), max(xs)
        y_lo, y_hi = min(min(ys), 0.0), max(ys)
        if x_hi <= x_lo:
            x_hi = x_lo + 1.0
        if y_hi <= y_lo:
            y_hi = y_lo + 1.0
        y_hi += 0.05 * (y_hi - y_lo)

        def to_px(x: float, y: float) -> tuple[float, float]:
            return (
                x0 + (x1 - x0) * (x - x_lo) / (x_hi - x_lo),
                y1 - (y1 - y0) * (y - y_lo) / (y_hi - y_lo),
            )

        for value in _nice_ticks(x_lo, x_hi, config.PLOT_TICKS):
            px, _ = to_px(value, y_lo)
            canvas.line(px, y0, px, y1, LIGHT)
            canvas.line(px, y1, px, y1 + 4, BLACK)
            canvas.text(int(px) - 12, y1 + 8, _format_tick(value), BLACK)
        for value in _nice_ticks(y_lo, y_hi, config.PLOT_TICKS):
            _, py = to_px(x_lo, value)
            canvas.line(x0, py, x1, py, LIGHT)
            canvas.line(x0 - 4, py, x0, py, BLACK)
            canvas.text_right(x0 - 6, int(py) - 3, _format_tick(value), BLACK)

        canvas.line(x0, y0, x0, y1, BLACK)
        canvas.line(x0, y1, x1, y1, BLACK)
        canvas.text(x0, base + 4, panel.get("titre", ""), BLACK)
        canvas.text(x1 - 8 * CHAR_WIDTH, y1 + 22, panel.get("x_label", ""), BLACK)
        # Etiquette d'ordonnee calee a gauche de l'axe pour ne pas recouvrir le titre.
        label = panel.get("y_label", "")
        canvas.text_right(x0 - 8, y0 - 14, label, BLACK)

        for order, serie in enumerate(series):
            color = serie.get("color") or SERIES_COLORS[order % len(SERIES_COLORS)]
            points = [to_px(x, y) for x, y in zip(serie["x"], serie["y"])]
            for a, b in zip(points, points[1:]):
                canvas.line(a[0], a[1], b[0], b[1], color, width=2)
            canvas.fill_rect(x1 + 14, y0 + 4 + order * 16, x1 + 26, y0 + 10 + order * 16, color)
            canvas.text(x1 + 32, y0 + 3 + order * 16, serie.get("label", ""), BLACK)
        for marker in panel.get("markers", []):
            mx, my = to_px(marker[0], marker[1])
            color = marker[2] if len(marker) > 2 else BLACK
            canvas.fill_rect(int(mx) - 3, int(my) - 3, int(mx) + 3, int(my) + 3, color)
    return canvas.write_png(path)
