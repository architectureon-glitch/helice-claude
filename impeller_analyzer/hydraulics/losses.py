"""Pertes calculees sur la geometrie du canal, pour comparer deux conceptions.

Le modele de la SPEC (5.3) pose la perte de frottement comme une **fraction du
point nominal** : `K_FROTTEMENT_REL` fois la hauteur theorique. L'incidence y est
nulle au nominal par definition, la fuite et le frottement de disque sont eux
aussi des fractions fixes. Toutes les pertes se retrouvent donc calees sur le
point de fonctionnement, et le rendement qui en sort ne depend pratiquement pas
de la forme des aubes : sur trois roues centrifuges radicalement differentes --
six aubes a 25 degres, trois a 65, huit a 15 -- il varie de sept centiemes de
point. Le modele est bon pour dimensionner une roue donnee ; il ne sait pas dire
si une conception vaut mieux qu'une autre.

Ce module calcule les deux pertes qui, elles, dependent vraiment de la
geometrie, et qui sont les deux criteres classiques du dessin de roue :

* le **frottement de canal**, par Darcy-Weisbach applique a la veine inter-aubes
  -- longueur developpee de l'aube rapportee au diametre hydraulique du canal ;
  c'est lui qui penalise un canal long, etroit, ou une aube dont la surface
  mouillee est doublee ;
* la **diffusion**, par le rapport `w2 / w1` des vitesses relatives. Une roue
  ralentit l'ecoulement relatif ; passe une certaine deceleration la couche
  limite decolle et la perte s'envole. Le seuil classique est celui de de
  Haller, `w2 / w1 >= 0.72`, etabli sur des grilles d'aubes de compresseur et
  repris tel quel en pompe.

Ces deux pertes ne remplacent pas celles de la SPEC : le rendement publie dans
le rapport reste le sien. Elles servent a **comparer** deux roues entre elles,
et le rendement de comparaison est cale pour qu'une roue centrifuge ordinaire
retombe sur la valeur de reference. Un ecart entre deux roues a un sens ; la
valeur absolue n'en a pas davantage que celle de la SPEC.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config


@dataclass
class ChannelLosses:
    """Pertes geometriques d'une roue, et les grandeurs dont elles sortent."""

    length: float = 0.0  # longueur developpee du canal, en m
    hydraulic_diameter: float = 0.0  # diametre hydraulique moyen, en m
    slenderness: float = 0.0  # L / D_h, sans dimension
    wetted_area: float = 0.0  # surface mouillee des aubes et des flasques, en m2
    w1: float = 0.0  # vitesse relative d'entree, en m/s
    w2: float = 0.0  # vitesse relative de sortie, en m/s
    de_haller: float = 0.0  # w2 / w1
    head_friction: float = 0.0  # perte de frottement de canal, en m
    head_diffusion: float = 0.0  # perte de diffusion, en m
    head_theoretical: float = 0.0  # hauteur theorique de reference, en m
    efficiency: float = 0.0  # rendement de comparaison, cale sur une roue ordinaire
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "longueur_de_canal_m": self.length,
            "diametre_hydraulique_m": self.hydraulic_diameter,
            "elancement_L_sur_Dh": self.slenderness,
            "surface_mouillee_m2": self.wetted_area,
            "w1_m_par_s": self.w1,
            "w2_m_par_s": self.w2,
            "de_haller_w2_sur_w1": self.de_haller,
            "perte_frottement_m": self.head_friction,
            "perte_diffusion_m": self.head_diffusion,
            "rendement_de_comparaison": self.efficiency,
        }


def channel_length(r_1: float, r_2: float, beta1_deg: float, beta2_deg: float) -> float:
    """Longueur developpee du canal entre l'entree et la sortie, en m.

    Le long de l'aube, `ds = dm / sin(beta)` puisque `tan(beta) = dm / (r dtheta)`.
    L'integration se fait a `beta` moyen, ce qui suffit a l'ordre de grandeur : ce
    qui compte ici est le rapport entre deux roues, non la valeur absolue. Une
    aube tres couchee -- `beta` petit -- donne un canal tres long, et c'est bien
    le comportement physique.
    """
    beta = math.radians(max(config.BETA_MIN_DEG, 0.5 * (beta1_deg + beta2_deg)))
    return abs(r_2 - r_1) / max(math.sin(beta), 1e-6)


def hydraulic_diameter(width: float, pitch: float) -> float:
    """Diametre hydraulique `4A/P` d'un canal rectangulaire `width` x `pitch`."""
    if width <= 0.0 or pitch <= 0.0:
        return 0.0
    return 2.0 * width * pitch / (width + pitch)


def analyse(
    r_1: float,
    r_2: float,
    b_1: float,
    b_2: float,
    beta1_deg: float,
    beta2_deg: float,
    n_blades: int,
    w1: float,
    w2: float,
    head_theoretical: float,
) -> ChannelLosses:
    """Pertes geometriques et rendement de comparaison d'une roue."""
    result = ChannelLosses(w1=w1, w2=w2, head_theoretical=head_theoretical)
    if min(r_1, r_2, b_1, b_2, w1, w2, head_theoretical) <= 0.0 or n_blades <= 0:
        result.warnings.append(
            "geometrie de canal incomplete : les pertes de comparaison ne sont pas calculees"
        )
        return result

    result.length = channel_length(r_1, r_2, beta1_deg, beta2_deg)
    width = 0.5 * (b_1 + b_2)
    pitch = 2.0 * math.pi * (0.5 * (r_1 + r_2)) / n_blades
    result.hydraulic_diameter = hydraulic_diameter(width, pitch)
    if result.hydraulic_diameter <= 0.0:
        result.warnings.append("diametre hydraulique nul : pertes de comparaison non calculees")
        return result
    result.slenderness = result.length / result.hydraulic_diameter
    # Les deux faces de chaque aube, plus les deux flasques du canal.
    result.wetted_area = n_blades * result.length * (2.0 * width + 2.0 * pitch)

    w_mean = 0.5 * (w1 + w2)
    result.head_friction = (
        4.0 * config.C_F_CANAL * result.slenderness * w_mean ** 2 / (2.0 * config.G)
    )

    result.de_haller = w2 / w1 if w1 > 0.0 else 0.0
    deficit = max(0.0, config.DE_HALLER_MIN - result.de_haller)
    result.head_diffusion = config.XI_DIFFUSION * deficit ** 2 * w1 ** 2 / (2.0 * config.G)
    if result.de_haller < config.DE_HALLER_MIN:
        result.warnings.append(
            f"deceleration relative w2/w1 = {result.de_haller:.2f}, sous le seuil de de Haller "
            f"({config.DE_HALLER_MIN:.2f}) : la couche limite decolle vraisemblablement dans le "
            "canal, et la perte de diffusion domine"
        )

    utile = head_theoretical - result.head_friction - result.head_diffusion
    brut = utile / head_theoretical if head_theoretical > 0.0 else 0.0
    result.efficiency = max(0.0, min(1.0, brut * config.ETA_COMPARAISON))
    return result
