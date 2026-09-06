"""NPSH requis, NPSH disponible et vitesse maximale admissible (SPEC phase 6)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import ConfidenceMap, worst
from ..numeric import linear_interp
from .meanline import PerformanceCurve

NPSH_LIMIT = "NPSH"
SPEED_LIMIT = "vitesse relative w1s"


# ---------------------------------------------------------------------------
# 6.2 Etat du fluide et de l'installation
# ---------------------------------------------------------------------------
def atmospheric_pressure(altitude: float = config.ALTITUDE) -> float:
    """Pression atmospherique a l'altitude donnee, atmosphere standard OACI."""
    factor = 1.0 - config.ATM_LAPSE_COEF * altitude
    if factor <= 0.0:
        return 0.0
    return config.P_ATM_SEA_LEVEL * factor ** config.ATM_EXPONENT


def vapour_pressure(temperature_c: float = config.TEMPERATURE) -> float:
    """Pression de vapeur saturante de l'eau, correlation d'Antoine."""
    denominator = config.ANTOINE_C + temperature_c
    if denominator <= 0.0:
        return 0.0
    exponent = config.ANTOINE_A - config.ANTOINE_B / denominator
    return 10.0 ** exponent * config.ANTOINE_MMHG_TO_PA


def water_density(temperature_c: float = config.TEMPERATURE) -> float:
    """Masse volumique de l'eau, interpolee dans la table de `config`."""
    return linear_interp(config.RHO_TABLE_T_C, config.RHO_TABLE_KG_M3, temperature_c)


@dataclass
class Installation:
    """Hypotheses d'installation et NPSH disponible qui en decoule (SPEC 6.2)."""

    altitude: float = config.ALTITUDE
    temperature_c: float = config.TEMPERATURE
    suction_height: float = config.HAUTEUR_ASPIRATION
    suction_losses: float = config.PERTES_ASPIRATION
    p_atm: float = 0.0
    p_vap: float = 0.0
    rho: float = config.RHO
    npsha: float = 0.0

    def describe(self) -> str:
        """Phrase resumant les hypotheses, pour l'encadre du rapport."""
        charge = "en charge" if self.suction_height >= 0.0 else "en aspiration"
        return (
            f"eau a {self.temperature_c:g} degres C, altitude {self.altitude:g} m, "
            f"hauteur d'aspiration {self.suction_height:+g} m ({charge}), "
            f"pertes de charge d'aspiration {self.suction_losses:g} m"
        )

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "altitude_m": self.altitude,
            "temperature_C": self.temperature_c,
            "hauteur_aspiration_m": self.suction_height,
            "pertes_aspiration_m": self.suction_losses,
            "p_atmospherique_Pa": self.p_atm,
            "p_vapeur_Pa": self.p_vap,
            "rho_kg_m3": self.rho,
            "NPSHa_m": self.npsha,
            "hypotheses": self.describe(),
        }


def installation(
    altitude: float = config.ALTITUDE,
    temperature_c: float = config.TEMPERATURE,
    suction_height: float = config.HAUTEUR_ASPIRATION,
    suction_losses: float = config.PERTES_ASPIRATION,
) -> Installation:
    """NPSH disponible : `(p_atm - p_vap) / (rho g) + z_aspiration - h_pertes`."""
    result = Installation(
        altitude=altitude,
        temperature_c=temperature_c,
        suction_height=suction_height,
        suction_losses=suction_losses,
    )
    result.p_atm = atmospheric_pressure(altitude)
    result.p_vap = vapour_pressure(temperature_c)
    result.rho = water_density(temperature_c)
    result.npsha = (
        (result.p_atm - result.p_vap) / (result.rho * config.G)
        + suction_height
        - suction_losses
    )
    return result


# ---------------------------------------------------------------------------
# 6.1 NPSH requis
# ---------------------------------------------------------------------------
def npshr_kinematic(cm1: float, w1s: float) -> float:
    """Methode A : `lc cm1^2 / 2g + lw w1s^2 / 2g`, au carter d'entree."""
    return (config.LAMBDA_C * cm1 * cm1 + config.LAMBDA_W * w1s * w1s) / (2.0 * config.G)


def npshr_suction_speed(rpm: float, flow: float) -> float:
    """Methode B : `(n sqrt(Q) / n_ss)^(4/3)`, n en tr/min et Q en m3/s."""
    if flow <= 0.0 or rpm <= 0.0 or config.N_SS <= 0.0:
        return 0.0
    return (rpm * math.sqrt(flow) / config.N_SS) ** config.NSS_EXPONENT


def apply_to_curve(curve: PerformanceCurve) -> PerformanceCurve:
    """Renseigne le NPSH requis sur chaque point de la courbe.

    Les deux methodes sont conservees ; la valeur retenue est la plus
    penalisante, comme la SPEC le demande.
    """
    for point in curve.points:
        point.npshr_kinematic = npshr_kinematic(point.cm1, point.w1s)
        point.npshr_suction_speed = npshr_suction_speed(curve.rpm, point.flow)
        point.npshr = max(point.npshr_kinematic, point.npshr_suction_speed)
    return curve


# ---------------------------------------------------------------------------
# 6.3 Vitesse maximale admissible
# ---------------------------------------------------------------------------
@dataclass
class SpeedLimit:
    """Vitesse maximale sans cavitation et limite qui la fixe (SPEC 6.3)."""

    rpm_reference: float = 0.0
    npshr_reference: float = 0.0
    w1s_reference: float = 0.0
    npsha: float = 0.0
    rpm_max_npsh: float = 0.0
    rpm_max_w1s: float = 0.0
    rpm_max: int = 0
    active_limit: str = ""
    assumptions: str = ""
    margin: float = 0.0  # marge NPSHa / NPSHr au regime de reference
    confidence: str = ""
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "regime_de_reference_tr_min": self.rpm_reference,
            "NPSHr_de_reference_m": self.npshr_reference,
            "w1s_de_reference_m_s": self.w1s_reference,
            "NPSHa_m": self.npsha,
            "vitesse_max_npsh_tr_min": self.rpm_max_npsh,
            "vitesse_max_w1s_tr_min": self.rpm_max_w1s,
            "vitesse_max_tr_min": self.rpm_max,
            "limite_active": self.active_limit,
            "hypotheses_installation": self.assumptions,
            "marge_npsh_au_reference": self.margin,
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
        }


def maximum_speed(
    curve: PerformanceCurve,
    site: Installation,
    confidence: ConfidenceMap | None = None,
) -> SpeedLimit:
    """Vitesse maximale admissible, la plus basse des deux limites (SPEC 6.3).

    Le point de reference est le point nominal (incidence nulle) du regime
    fourni : c'est celui que reporte le tableau de performances.
    """
    result = SpeedLimit(rpm_reference=curve.rpm, npsha=site.npsha)
    result.assumptions = site.describe()
    reference = curve.nominal_point()
    if reference is None or curve.rpm <= 0.0:
        result.warnings.append("aucun point nominal : vitesse maximale non calculable")
        result.confidence = "low"
        return result

    result.npshr_reference = reference.npshr
    result.w1s_reference = reference.w1s
    result.margin = site.npsha / reference.npshr if reference.npshr > 0.0 else math.inf

    if site.npsha <= 0.0:
        result.warnings.append(
            f"NPSH disponible negatif ou nul ({site.npsha:.2f} m) : l'installation decrite "
            "ne permet aucun fonctionnement sans cavitation"
        )
        result.rpm_max_npsh = 0.0
    elif reference.npshr > 0.0:
        # NPSHr varie comme n^2 a coefficient de debit constant.
        result.rpm_max_npsh = curve.rpm * math.sqrt(
            site.npsha / (config.MARGE_NPSH * reference.npshr)
        )
    else:
        result.rpm_max_npsh = math.inf

    if reference.w1s > 0.0:
        result.rpm_max_w1s = curve.rpm * config.W1S_MAX / reference.w1s
    else:
        result.rpm_max_w1s = math.inf

    if result.rpm_max_npsh <= result.rpm_max_w1s:
        result.active_limit = NPSH_LIMIT
        limit = result.rpm_max_npsh
    else:
        result.active_limit = SPEED_LIMIT
        limit = result.rpm_max_w1s

    if not math.isfinite(limit):
        result.rpm_max = 0
        result.warnings.append("les deux limites sont indeterminees : geometrie degeneree")
    else:
        result.rpm_max = int(math.floor(limit / config.RPM_ROUNDING) * config.RPM_ROUNDING)

    if result.rpm_max < curve.rpm:
        result.warnings.append(
            f"la vitesse maximale admissible ({result.rpm_max} tr/min) est inferieure au regime "
            f"analyse de {curve.rpm:.0f} tr/min : a ce regime la roue cavite dans l'installation decrite"
        )
    level = (confidence or ConfidenceMap()).overall() if confidence else "medium"
    result.confidence = worst(level, "medium")  # NPSHr : incertitude annoncee de 30 %
    return result
