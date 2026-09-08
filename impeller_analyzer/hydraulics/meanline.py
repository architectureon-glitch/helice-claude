"""Modele hydraulique 1D de ligne moyenne (SPEC phase 5).

Euler avec glissement de Wiesner, pertes de frottement et d'incidence, balayage
de la courbe H-Q, puissances et rendements.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import ConfidenceMap, worst
from ..geometry.blade_angles import BladeGeometry
from ..geometry.topology import Topology, specific_speed


@dataclass
class MeanlineInput:
    """Les seules grandeurs geometriques dont le modele 1D a besoin (SI)."""

    r_1: float = 0.0  # rayon quadratique moyen d'entree
    r_1s: float = 0.0  # rayon exterieur d'entree (carter)
    r_1h: float = 0.0  # rayon de moyeu a l'entree
    r_2: float = 0.0  # rayon de sortie
    area_1: float = 0.0  # section d'entree, obstruction comprise
    area_2: float = 0.0  # section de sortie, obstruction comprise
    beta1_deg: float = 0.0
    beta2_deg: float = 0.0
    n_blades: int = 0
    rho: float = config.RHO

    @classmethod
    def from_geometry(
        cls,
        topology: Topology,
        geometry: BladeGeometry,
        rho: float = config.RHO,
    ) -> "MeanlineInput":
        """Assemble les entrees du modele a partir des phases 3 et 4."""
        return cls(
            r_1=topology.r_1,
            r_1s=topology.r_1s,
            r_1h=topology.r_1h,
            r_2=topology.r_2,
            area_1=topology.area_1,
            area_2=topology.area_2,
            beta1_deg=geometry.beta1_deg,
            beta2_deg=geometry.beta2_deg,
            # Le glissement se calcule au refoulement. Sur une aube en boucle les
            # deux brins y ont fusionne en un seul bord de fuite : le nombre de
            # passages a la sortie reste le nombre d'aubes.
            n_blades=topology.blades.n_blades,
            rho=rho,
        )

    def valid(self) -> bool:
        """Vrai si toutes les grandeurs necessaires sont exploitables."""
        return (
            self.r_1 > 0.0
            and self.r_2 > 0.0
            and self.area_1 > 0.0
            and self.area_2 > 0.0
            and self.n_blades >= config.BLADES_MIN
            and config.BETA_MIN_DEG <= self.beta1_deg <= config.BETA_MAX_DEG
            and config.BETA_MIN_DEG <= self.beta2_deg <= config.BETA_MAX_DEG
        )


@dataclass
class OperatingPoint:
    """Un point de la courbe caracteristique, tout en SI."""

    flow: float = 0.0  # m3/s
    head: float = 0.0  # m
    head_theoretical: float = 0.0  # m, Euler avec glissement
    cm1: float = 0.0  # m/s, vitesse meridienne d'entree
    cm2: float = 0.0  # m/s, vitesse meridienne de sortie
    cu2: float = 0.0  # m/s, vitesse tangentielle absolue de sortie
    w1: float = 0.0  # m/s, vitesse relative d'entree au rayon moyen
    w1s: float = 0.0  # m/s, vitesse relative d'entree au carter
    loss_friction: float = 0.0  # m
    loss_incidence: float = 0.0  # m
    hydraulic_power: float = 0.0  # W
    shaft_power: float = 0.0  # W
    torque: float = 0.0  # N.m
    efficiency: float = 0.0  # -
    npshr: float = 0.0  # m, renseigne par la phase 6
    npshr_kinematic: float = 0.0  # m
    npshr_suction_speed: float = 0.0  # m

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "Q_m3_s": self.flow,
            "H_m": self.head,
            "H_theorique_m": self.head_theoretical,
            "cm1_m_s": self.cm1,
            "cm2_m_s": self.cm2,
            "cu2_m_s": self.cu2,
            "w1_m_s": self.w1,
            "w1s_m_s": self.w1s,
            "perte_frottement_m": self.loss_friction,
            "perte_incidence_m": self.loss_incidence,
            "P_hydraulique_W": self.hydraulic_power,
            "P_arbre_W": self.shaft_power,
            "couple_Nm": self.torque,
            "rendement": self.efficiency,
            "NPSHr_m": self.npshr,
            "NPSHr_cinematique_m": self.npshr_kinematic,
            "NPSHr_vitesse_specifique_m": self.npshr_suction_speed,
        }


@dataclass
class PerformanceCurve:
    """Courbe caracteristique complete a une vitesse de rotation donnee."""

    rpm: float = 0.0
    omega: float = 0.0  # rad/s
    u1: float = 0.0  # m/s
    u1s: float = 0.0  # m/s
    u2: float = 0.0  # m/s
    slip: float = 0.0  # coefficient de glissement de Wiesner
    efficiency_raw: float = 0.0  # rendement au BEP des seules pertes modelisees
    flow_nominal: float = 0.0  # m3/s, incidence nulle
    friction_coefficient: float = 0.0  # s2/m5
    points: list[OperatingPoint] = field(default_factory=list)
    nominal_index: int = 0
    bep_index: int = 0
    specific_speed: float = 0.0
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    warnings: list[str] = field(default_factory=list)

    def nominal_point(self) -> OperatingPoint | None:
        """Point de fonctionnement a incidence nulle."""
        if not self.points:
            return None
        return self.points[self.nominal_index]

    def best_efficiency_point(self) -> OperatingPoint | None:
        """Point de meilleur rendement."""
        if not self.points:
            return None
        return self.points[self.bep_index]

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        nominal = self.nominal_point()
        bep = self.best_efficiency_point()
        return {
            "regime_tr_min": self.rpm,
            "omega_rad_s": self.omega,
            "u1_m_s": self.u1,
            "u1s_m_s": self.u1s,
            "u2_m_s": self.u2,
            "glissement_wiesner": self.slip,
            "Q_nominal_m3_s": self.flow_nominal,
            "k_frottement_s2_m5": self.friction_coefficient,
            "vitesse_specifique_nq": self.specific_speed,
            "point_nominal": nominal.to_dict() if nominal else None,
            "point_de_meilleur_rendement": bep.to_dict() if bep else None,
            "courbe": [point.to_dict() for point in self.points],
            "confiance": dict(self.confidence),
            "avertissements": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# Briques du modele
# ---------------------------------------------------------------------------
def angular_velocity(rpm: float) -> float:
    """Vitesse de rotation en rad/s."""
    return rpm * config.RPM_TO_RAD_S


def slip_factor(beta2_deg: float, n_blades: int) -> float:
    """Coefficient de glissement de Wiesner (1967) : `1 - sqrt(sin b2) / N^0.7`."""
    if n_blades <= 0:
        return 1.0
    value = 1.0 - math.sqrt(math.sin(math.radians(beta2_deg))) / n_blades ** config.WIESNER_EXPONENT
    return max(0.0, min(1.0, value))


def nominal_flow(data: MeanlineInput, omega: float) -> tuple[float, float, float]:
    """Point d'incidence nulle : renvoie `(Q_n, u1, cm1)` (SPEC 5.1)."""
    u1 = omega * data.r_1
    cm1 = u1 * math.tan(math.radians(data.beta1_deg))
    return cm1 * data.area_1, u1, cm1


def euler_head(data: MeanlineInput, omega: float, flow: float, slip: float) -> tuple[float, float, float]:
    """Hauteur d'Euler avec glissement : renvoie `(H_th, cm2, cu2)` (SPEC 5.2)."""
    u2 = omega * data.r_2
    cm2 = flow / data.area_2
    cu2 = slip * u2 - cm2 / math.tan(math.radians(data.beta2_deg))
    return u2 * cu2 / config.G, cm2, cu2


def build_curve(
    data: MeanlineInput,
    rpm: float,
    confidence: ConfidenceMap | None = None,
) -> PerformanceCurve:
    """Courbe caracteristique a `rpm` tours par minute (SPEC 5.1 a 5.4).

    Le rendement varie le long de la courbe et vaut `ETA_H * ETA_VOL * ETA_MEC`
    au point de meilleur rendement.  C'est ce qui donne un sens au BEP : la SPEC
    le definit comme le maximum de `rho g Q H / P_arbre` (5.3) tout en posant
    `P_arbre = rho g Q H / eta_global` avec un `eta_global` constant (5.4), ce
    qui rendrait le rendement constant et son maximum indetermine.  Les pertes
    reelles, elles, dependent du debit.
    """
    curve = PerformanceCurve(rpm=rpm)
    curve.confidence = ConfidenceMap(confidence or {})
    if not data.valid():
        curve.warnings.append(
            "geometrie insuffisante pour le modele 1D : rayons, sections, angles de pale "
            "ou nombre de pales manquants ou hors domaine"
        )
        return curve

    omega = angular_velocity(rpm)
    curve.omega = omega
    flow_nominal, u1, cm1_nominal = nominal_flow(data, omega)
    curve.flow_nominal = flow_nominal
    curve.u1 = u1
    curve.u1s = omega * data.r_1s
    curve.u2 = omega * data.r_2
    curve.slip = slip_factor(data.beta2_deg, data.n_blades)

    if flow_nominal <= 0.0 or omega <= 0.0:
        curve.warnings.append("debit nominal nul : vitesse de rotation ou geometrie degeneree")
        return curve

    head_theoretical_nominal, _, _ = euler_head(data, omega, flow_nominal, curve.slip)
    if head_theoretical_nominal <= 0.0:
        curve.warnings.append(
            "hauteur d'Euler negative au point nominal : la roue ne peut pas fonctionner en "
            "pompe avec ces angles de pale (verifiez beta2 et le sens de rotation)"
        )
        return curve

    # k_f est cale pour que la perte de frottement vaille K_FROTTEMENT_REL de la
    # hauteur theorique au point nominal (SPEC 5.3).
    curve.friction_coefficient = config.K_FROTTEMENT_REL * head_theoretical_nominal / flow_nominal ** 2
    w1_nominal = math.hypot(cm1_nominal, u1)

    # Balayage de Q = 0 a Q_SWEEP_MAX * Q_n en Q_SWEEP_POINTS points. Avec 29
    # points et une borne a 1.40, le point nominal tombe exactement sur le
    # vingtieme : il n'y a pas a l'ajouter au balayage.
    flows = [
        max(
            config.Q_MIN_RELATIVE * flow_nominal,
            flow_nominal * config.Q_SWEEP_MAX * index / (config.Q_SWEEP_POINTS - 1),
        )
        for index in range(config.Q_SWEEP_POINTS)
    ]
    curve.nominal_index = min(
        range(len(flows)), key=lambda index: abs(flows[index] - flow_nominal)
    )

    for flow in flows:
        curve.points.append(_operating_point(data, curve, flow, w1_nominal))
    # Les trois rendements de la SPEC 5.4 sont des valeurs **au point nominal** :
    # le rendement volumetrique correspond a une fuite de recirculation a peu
    # pres constante, le rendement mecanique a une puissance de frottement de
    # disque a peu pres constante.  Les faire porter par ces deux constantes,
    # plutot que par un rendement global fixe, laisse leur valeur nominale
    # intacte et donne au rendement une vraie dependance au debit -- sans quoi
    # le BEP, defini en 5.3 comme le maximum de rho g Q H / P_arbre, serait
    # indetermine puisque 5.4 pose P_arbre = rho g Q H / eta_global.
    leakage = (1.0 / config.ETA_VOL - 1.0) * flow_nominal
    disc_power = (1.0 / config.ETA_MEC - 1.0) * data.rho * config.G * (
        flow_nominal + leakage
    ) * head_theoretical_nominal
    for point in curve.points:
        point.hydraulic_power = data.rho * config.G * point.flow * max(point.head, 0.0)
        point.shaft_power = (
            data.rho * config.G * (point.flow + leakage) * max(point.head_theoretical, 0.0)
            + disc_power
        )
        point.efficiency = (
            point.hydraulic_power / point.shaft_power if point.shaft_power > 0.0 else 0.0
        )

    usable = [index for index, point in enumerate(curve.points) if point.head > 0.0 and point.flow > 0.0]
    curve.bep_index = (
        max(usable, key=lambda index: curve.points[index].efficiency) if usable else curve.nominal_index
    )

    # Les pertes que le modele detaille (frottement, incidence, fuite,
    # frottement de disque) restent en deca de celles qu'annonce la SPEC :
    # ETA_H * ETA_VOL * ETA_MEC = 0.803 au point de meilleur rendement. Le
    # facteur ci-dessous porte le complement, suppose independant du debit. Il
    # est borne a 1 pour que le modele n'invente jamais de rendement quand ses
    # pertes propres depassent deja celles de la SPEC.
    eta_target = config.ETA_H * config.ETA_VOL * config.ETA_MEC
    eta_raw = curve.points[curve.bep_index].efficiency
    curve.efficiency_raw = eta_raw
    unmodelled = max(1.0, eta_raw / eta_target) if eta_target > 0.0 else 1.0
    if eta_raw < eta_target:
        curve.warnings.append(
            f"les pertes calculees plafonnent le rendement a {eta_raw * 100.0:.1f} %, "
            f"en deca des {eta_target * 100.0:.1f} % annonces par les rendements de reference : "
            "la valeur calculee est conservee"
        )
    for point in curve.points:
        point.shaft_power *= unmodelled
        point.efficiency = max(0.0, min(1.0, point.efficiency / unmodelled))
        point.torque = point.shaft_power / omega if omega > 0.0 else 0.0
    best = curve.points[curve.bep_index]
    curve.specific_speed = specific_speed(rpm, best.flow, best.head)
    return curve


def _operating_point(
    data: MeanlineInput, curve: PerformanceCurve, flow: float, w1_nominal: float
) -> OperatingPoint:
    """Un point de la courbe, avant calcul des puissances."""
    point = OperatingPoint(flow=flow)
    point.head_theoretical, point.cm2, point.cu2 = euler_head(data, curve.omega, flow, curve.slip)
    point.cm1 = flow / data.area_1
    point.w1 = math.hypot(point.cm1, curve.u1)
    point.w1s = math.hypot(point.cm1, curve.u1s)
    point.loss_friction = curve.friction_coefficient * flow ** 2
    point.loss_incidence = config.XI_INCIDENCE * (point.w1 - w1_nominal) ** 2 / (2.0 * config.G)
    point.head = point.head_theoretical - point.loss_friction - point.loss_incidence
    return point


def curves_for_speeds(
    data: MeanlineInput,
    speeds=config.DEFAULT_RPM,
    confidence: ConfidenceMap | None = None,
) -> list[PerformanceCurve]:
    """Courbes caracteristiques pour une liste de regimes."""
    return [build_curve(data, rpm, confidence) for rpm in speeds]


def head_sensitivity(data: "MeanlineInput", rpm: float) -> float:
    """Variation relative de la hauteur nominale pour +/- 1 degre sur beta2.

    La hauteur d'Euler passe par `cu2 = u2 - cm2 / tan(beta2)`. Aux petits
    angles la tangente varie tres vite : sur une aube tres couchee, un degre
    d'incertitude sur beta2 -- l'ordre de grandeur de ce que sait faire
    n'importe quelle lecture geometrique -- peut deplacer la hauteur bien
    au-dela des 18 % annonces par le modele. Autant le mesurer et le dire.

    Renvoie l'ecart relatif entre les hauteurs nominales obtenues a
    `beta2 - 1 deg` et `beta2 + 1 deg`, rapporte a la hauteur centrale, et
    l'infini si un degre suffit a faire disparaitre le point de fonctionnement.
    """
    import dataclasses

    heads = []
    for delta in (-config.BETA_SENSITIVITY_DEG, 0.0, config.BETA_SENSITIVITY_DEG):
        beta2 = data.beta2_deg + delta
        if not 0.0 < beta2 < 90.0:
            return math.inf
        curve = build_curve(dataclasses.replace(data, beta2_deg=beta2), rpm, ConfidenceMap())
        point = curve.nominal_point()
        if point is None or point.head <= 0.0:
            return math.inf
        heads.append(point.head)
    low, middle, high = heads
    return abs(high - low) / middle if middle > 0.0 else math.inf


def confidence_of(topology: Topology, geometry: BladeGeometry) -> ConfidenceMap:
    """Confiance heritee par les grandeurs hydrauliques (propagation pessimiste)."""
    level = worst(
        topology.confidence.get_level("rayons"),
        topology.confidence.get_level("sections"),
        topology.confidence.get_level("nombre_de_pales"),
        geometry.confidence.get_level("angles_de_pale"),
    )
    result = ConfidenceMap()
    result.set("debit", level)
    result.set("hauteur", level)
    result.set("puissance", level)
    result.set("couple", level)
    result.set("rendement", level)
    return result
