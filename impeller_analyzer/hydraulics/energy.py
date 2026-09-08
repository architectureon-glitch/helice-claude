"""Bilan d'energie d'Euler : d'ou vient la hauteur, terme par terme.

La hauteur theorique d'une roue vaut `H = (u2 cu2 - u1 cu1) / g`. Ce n'est pas
un modele mais un theoreme : il sort de la conservation du moment cinetique, et
il vaut pour n'importe quelle forme d'aube. Le couple sur l'arbre est
`rho Q (r2 cu2 - r1 cu1)` ; aucune geometrie ne peut ajouter d'energie autrement
qu'en changeant `cu2`.

La meme quantite se decompose en trois termes, qui disent **par quel mecanisme**
l'energie est transmise :

    H = (u2^2 - u1^2)/2g  +  (w1^2 - w2^2)/2g  +  (c2^2 - c1^2)/2g
        \\____ centrifuge ___/  \\___ diffusion ___/  \\___ cinetique ___/

* le terme **centrifuge** ne depend que des rayons et du regime -- pas de la
  forme des aubes. Deux roues de meme diametre tournant a la meme vitesse en
  tirent exactement la meme chose ;
* la **diffusion** est la conversion de vitesse relative en pression statique
  dans le canal. Elle est positive quand l'ecoulement relatif ralentit,
  `w2 < w1`, ce que fait un canal de pompe bien dessine. Negative, elle
  **detruit** de la hauteur ;
* le terme **cinetique** est l'energie qui sort sous forme de vitesse absolue,
  et que la volute doit recuperer.

La somme des trois retombe exactement sur `u2 cu2 / g` : c'est le controle que
fait ce module. La decomposition ne cree rien, elle repartit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .. import config


@dataclass
class EnergyBudget:
    """Les trois termes du bilan, en metres."""

    centrifugal: float = 0.0
    diffusion: float = 0.0
    kinetic: float = 0.0
    total: float = 0.0  # somme des trois
    euler: float = 0.0  # u2 cu2 / g, pour controle
    w1: float = 0.0
    w2: float = 0.0
    residual: float = 0.0  # ecart relatif entre la somme et Euler

    def to_dict(self) -> dict:
        """Vue serialisable, hauteurs en metres."""
        return {
            "terme_centrifuge_m": self.centrifugal,
            "terme_diffusion_m": self.diffusion,
            "terme_cinetique_m": self.kinetic,
            "somme_m": self.total,
            "euler_u2_cu2_sur_g_m": self.euler,
            "w1_m_par_s": self.w1,
            "w2_m_par_s": self.w2,
            "residu_relatif": self.residual,
        }

    def shares(self) -> tuple[float, float, float]:
        """Part de chaque terme dans la somme, en fraction."""
        if self.total == 0.0:
            return (0.0, 0.0, 0.0)
        return (
            self.centrifugal / self.total,
            self.diffusion / self.total,
            self.kinetic / self.total,
        )


def budget(u1: float, u2: float, cm1: float, cm2: float, cu2: float) -> EnergyBudget:
    """Decompose la hauteur d'Euler, sans prerotation a l'entree (`cu1 = 0`)."""
    result = EnergyBudget()
    result.w1 = math.hypot(cm1, u1)
    result.w2 = math.hypot(cm2, u2 - cu2)
    c1, c2 = cm1, math.hypot(cm2, cu2)

    result.centrifugal = (u2 * u2 - u1 * u1) / (2.0 * config.G)
    result.diffusion = (result.w1 ** 2 - result.w2 ** 2) / (2.0 * config.G)
    result.kinetic = (c2 * c2 - c1 * c1) / (2.0 * config.G)
    result.total = result.centrifugal + result.diffusion + result.kinetic
    result.euler = u2 * cu2 / config.G
    if result.euler != 0.0:
        result.residual = abs(result.total - result.euler) / abs(result.euler)
    return result
