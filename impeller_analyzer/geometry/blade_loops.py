"""Aubes en boucle fermee : detection, et ce qu'elles invalident.

Le modele de ligne moyenne de la SPEC suppose une aube **simple** : une surface
portant un bord d'attaque et un bord de fuite, dont une coupe sur une surface de
courant donne un profil unique par pale.  Une aube **toroidale** est une boucle
fermee : elle part du moyeu, sort, se retourne au bout et revient.  Une coupe la
traverse donc deux fois, et l'appariement des deux faces dont sort la cambrure
apparie alors la face d'un brin avec celle de l'autre.  Les angles de pale et le
sens de rotation qui en derivent n'ont aucun sens.

La signature est topologique et se lit sur la carte d'occupation sans rien
recouper.  Elle ne se cherche pas en azimut mais **en hauteur** : les deux brins
d'une boucle sont a la meme azimut, separes en z.  A rayon et azimut fixes, une
aube simple donne un tronçon unique le long de z, une boucle en donne deux.  La
grandeur mesuree est donc la fraction des azimuts ou la coupe rencontre deux
tronçons, relevee rayon par rayon.  Elle vaut 0.95 a 1.00 sur toute la portee
d'une aube toroidale et retombe au bout, la ou les brins fusionnent ; sur une
roue centrifuge fermee ordinaire elle ne depasse pas 0.25, sur une helice axiale
elle est nulle.

Cette signature est **topologique**, et un maillage non etanche n'a pas de
topologie : ses trous dedoublent les tronçons exactement comme le ferait une
boucle.  Le verdict est donc suspendu sur un maillage troue -- dire "toroidal"
la ou il n'y a que des trous serait la pire des sorties, puisque c'est ce mot
qui met de cote la cambrure et invalide toute la ligne moyenne.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .. import config
from ..confidence import HIGH, LOW
from .occupancy import OccupancyMap


# Grandeurs que le modele de ligne moyenne ne sait pas etablir sur une aube en
# boucle : toutes celles qui passent par la cambrure d'un profil unique.
INVALIDATED_BY_LOOP = (
    "angles_de_pale",
    "sens_de_rotation",
    "hauteur",
    "debit",
    "puissance",
    "couple",
    "rendement",
    "npshr",
)


@dataclass
class LoopResult:
    """Ce que la lecture des tronçons en hauteur a conclu."""

    looped: bool = False
    undecided: bool = False
    peak_fraction: float = 0.0
    r_inner: float = 0.0
    r_outer: float = 0.0
    merge_radius: float = 0.0
    n_blades: int = 0
    confidence: str = HIGH
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "aubes_en_boucle": self.looped,
            "verdict_suspendu": self.undecided,
            "fraction_dedoublee_max": self.peak_fraction,
            "rayon_interieur_de_boucle_m": self.r_inner or None,
            "rayon_exterieur_de_boucle_m": self.r_outer or None,
            "rayon_de_fusion_m": self.merge_radius or None,
            "confiance": self.confidence,
        }


def z_segments(
    occupancy: OccupancyMap, blade_mask: list[list[bool]], ir: int, itheta: int
) -> int:
    """Nombre de tronçons de pale separes le long de z, a rayon et azimut fixes."""
    previous = False
    count = 0
    for iz in range(occupancy.nz):
        here = blade_mask[iz][ir] and any(
            start <= ir < stop for start, stop in occupancy.intervals[iz][itheta]
        )
        if here and not previous:
            count += 1
        previous = here
    return count


def doubling_profile(
    occupancy: OccupancyMap, blade_mask: list[list[bool]]
) -> list[tuple[float, float]]:
    """Profil radial `(rayon, fraction des azimuts coupant l'aube deux fois)`."""
    profile: list[tuple[float, float]] = []
    for ir in range(occupancy.nr):
        if not any(blade_mask[iz][ir] for iz in range(occupancy.nz)):
            continue
        counts = [z_segments(occupancy, blade_mask, ir, t) for t in range(occupancy.n_theta)]
        seen = [c for c in counts if c >= 1]
        if len(seen) < config.LOOP_MIN_SECTORS:
            continue
        profile.append((
            occupancy.r_centres[ir],
            sum(1 for c in seen if c >= 2) / len(seen),
        ))
    return profile


def detect_looped_blades(
    occupancy: OccupancyMap,
    n_blades: int,
    blade_mask: list[list[bool]] | None = None,
    watertight: bool = True,
) -> LoopResult:
    """Dit si les aubes se referment sur elles-memes (type toroidal).

    `watertight` dit si le maillage a survecu a la reparation sans arete de
    bord.  Faux, la lecture est rendue sans verdict : un maillage troue produit
    la meme signature qu'une boucle.
    """
    from .topology import clean_blade_mask

    result = LoopResult(n_blades=n_blades)
    if blade_mask is None:
        blade_mask = clean_blade_mask(occupancy.blade_mask())

    profile = doubling_profile(occupancy, blade_mask)
    if len(profile) < config.LOOP_MIN_STATIONS:
        result.confidence = LOW
        result.notes.append(
            "zone de pales trop courte pour lire la forme des aubes : le caractere "
            "simple ou boucle n'est pas verifie"
        )
        return result

    result.peak_fraction = max(fraction for _, fraction in profile)
    doubled = [r for r, fraction in profile if fraction >= config.LOOP_DOUBLE_FRACTION]
    if len(doubled) < config.LOOP_MIN_STATIONS:
        result.notes.append(
            f"coupe simple sur {1.0 - result.peak_fraction:.0%} des azimuts au moins : les "
            "aubes sont des surfaces a bord d'attaque et bord de fuite uniques, le modele de "
            "ligne moyenne s'applique"
        )
        return result

    result.r_inner, result.r_outer = min(doubled), max(doubled)
    if not watertight:
        # Les trous d'un maillage non etanche coupent les tronçons en deux
        # exactement comme le ferait une boucle : la mesure ne distingue plus
        # les deux, et l'annoncer toroidal ecarterait la cambrure a tort.
        result.undecided = True
        result.confidence = LOW
        result.warnings.append(
            f"une coupe a azimut fixe traverse la zone de pales deux fois sur "
            f"{result.peak_fraction:.0%} des azimuts, ce qui est la signature d'aubes en boucle "
            "(type toroidal) -- mais le maillage n'est pas etanche, et ses trous donnent la meme "
            "signature. Le verdict est suspendu : les angles restent lus sur la cambrure. "
            "Reparez le maillage pour trancher."
        )
        return result

    result.looped = True
    beyond = [r for r, fraction in profile
              if r > result.r_outer and fraction < config.LOOP_DOUBLE_FRACTION]
    result.merge_radius = min(beyond) if beyond else 0.0
    fusion = (
        f", qui fusionnent a r = {result.merge_radius * config.MM_PER_M:.1f} mm"
        if result.merge_radius else ""
    )
    result.warnings.append(
        f"aubes en boucle fermee (type toroidal) : de r = "
        f"{result.r_inner * config.MM_PER_M:.1f} a {result.r_outer * config.MM_PER_M:.1f} mm, "
        f"une coupe a azimut fixe traverse l'aube deux fois sur "
        f"{result.peak_fraction:.0%} des azimuts{fusion}. Chaque aube a donc deux brins et non "
        "un bord d'attaque et un bord de fuite : la cambrure, qui suppose l'inverse, est mise de "
        "cote, et les angles de pale sont lus sur les normales de la surface (voir plus bas). "
        "Axe, nombre d'aubes, rayons, sections et volumes ne sont pas concernes."
    )
    return result
