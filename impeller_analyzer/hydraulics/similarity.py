"""Lois de similitude et vitesse specifique (SPEC 5.5 et 3.3).

Les lois de similitude ne servent pas a produire les resultats -- chaque regime
est calcule directement -- mais a les **verifier** : a geometrie donnee,
`Q` varie comme `n`, `H` comme `n^2` et `P` comme `n^3`.  Un ecart superieur a
`SIMILARITY_TOL` entre le calcul direct et l'extrapolation signale une erreur de
programmation, pas un effet physique.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .. import config
from ..geometry.topology import specific_speed, type_from_specific_speed  # noqa: F401  (re-export)


def scale_flow(flow: float, rpm_from: float, rpm_to: float) -> float:
    """Debit extrapole : `Q` proportionnel a `n`."""
    if rpm_from <= 0.0:
        return 0.0
    return flow * rpm_to / rpm_from


def scale_head(head: float, rpm_from: float, rpm_to: float) -> float:
    """Hauteur extrapolee : `H` proportionnelle a `n^2`."""
    if rpm_from <= 0.0:
        return 0.0
    return head * (rpm_to / rpm_from) ** 2


def scale_power(power: float, rpm_from: float, rpm_to: float) -> float:
    """Puissance extrapolee : `P` proportionnelle a `n^3`."""
    if rpm_from <= 0.0:
        return 0.0
    return power * (rpm_to / rpm_from) ** 3


@dataclass
class SimilarityCheck:
    """Comparaison entre calcul direct et extrapolation depuis un regime de reference."""

    rpm_reference: float = 0.0
    rpm_target: float = 0.0
    flow_direct: float = 0.0
    flow_scaled: float = 0.0
    head_direct: float = 0.0
    head_scaled: float = 0.0
    power_direct: float = 0.0
    power_scaled: float = 0.0
    worst_deviation: float = 0.0
    passed: bool = True
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "regime_de_reference_tr_min": self.rpm_reference,
            "regime_verifie_tr_min": self.rpm_target,
            "debit_direct_m3_s": self.flow_direct,
            "debit_extrapole_m3_s": self.flow_scaled,
            "hauteur_directe_m": self.head_direct,
            "hauteur_extrapolee_m": self.head_scaled,
            "puissance_directe_W": self.power_direct,
            "puissance_extrapolee_W": self.power_scaled,
            "ecart_maximal": self.worst_deviation,
            "conforme": self.passed,
            "avertissements": list(self.warnings),
        }


def _deviation(direct: float, scaled: float) -> float:
    """Ecart relatif entre calcul direct et extrapolation."""
    if direct == 0.0:
        return 0.0 if scaled == 0.0 else 1.0
    return abs(scaled - direct) / abs(direct)


def check(reference, target) -> SimilarityCheck:
    """Verifie la similitude entre deux courbes, a leur point nominal.

    `reference` et `target` sont deux `PerformanceCurve` de la meme roue a des
    vitesses differentes.
    """
    result = SimilarityCheck(rpm_reference=reference.rpm, rpm_target=target.rpm)
    origin = reference.nominal_point()
    goal = target.nominal_point()
    if origin is None or goal is None:
        result.passed = False
        result.warnings.append("point nominal indisponible : similitude non verifiable")
        return result

    result.flow_direct = goal.flow
    result.head_direct = goal.head
    result.power_direct = goal.shaft_power
    result.flow_scaled = scale_flow(origin.flow, reference.rpm, target.rpm)
    result.head_scaled = scale_head(origin.head, reference.rpm, target.rpm)
    result.power_scaled = scale_power(origin.shaft_power, reference.rpm, target.rpm)

    result.worst_deviation = max(
        _deviation(result.flow_direct, result.flow_scaled),
        _deviation(result.head_direct, result.head_scaled),
        _deviation(result.power_direct, result.power_scaled),
    )
    result.passed = result.worst_deviation <= config.SIMILARITY_TOL
    if not result.passed:
        result.warnings.append(
            f"ecart de {result.worst_deviation * 100.0:.1f} % entre le calcul direct a "
            f"{target.rpm:.0f} tr/min et l'extrapolation depuis {reference.rpm:.0f} tr/min "
            f"(seuil {config.SIMILARITY_TOL * 100.0:.0f} %) : le modele n'est pas homogene, "
            "c'est le signe d'une erreur de calcul"
        )
    return result
