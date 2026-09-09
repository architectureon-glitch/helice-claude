"""Ce qu'une erreur d'un degre sur la geometrie fait aux resultats publies.

Le modele de ligne moyenne annonce une incertitude de 18 % sur la hauteur.  Ce
chiffre suppose que la geometrie lue est juste.  Or `cu2 = u2 - cm2 / tan(b2)`
passe par une tangente, et aux petits angles celle-ci varie tres vite : sur une
aube couchee a quatre degres, un degre d'ecart -- l'ordre de grandeur de ce que
sait faire n'importe quelle lecture geometrique -- deplace la hauteur bien
au-dela des 18 %.  L'incertitude affichee devient alors trompeuse, non par ce
qu'elle dit, mais par ce qu'elle laisse croire.

Ce module mesure la chose plutot que de la supposer : hauteur, debit et NPSHr
sont recalcules a `beta1 +/- 1 deg`, `beta2 +/- 1 deg` et diametre `+/- 1 %`, et
l'ecart relatif est publie.  Deux lectures differentes selon la colonne :

* pour les angles, une **sensibilite par degre** -- « un degre de plus sur b2
  et la hauteur bouge de tant » ;
* pour le diametre, une **elasticite** sans dimension, `(dX/X) / (dD/D)`.  Les
  lois de similitude en donnent la valeur attendue : 2 pour la hauteur, 3 pour
  le debit, 2 pour le NPSHr.  La colonne sert donc aussi de controle du modele.

Une sensibilite qui depasse `SENSITIVITY_LOW_PER_DEG` declasse la grandeur en
confiance faible, automatiquement : sur une geometrie pareille, la valeur
publiee n'est plus un resultat, c'est un ordre de grandeur.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import LOW, ConfidenceMap
from . import cavitation as cavitation_module
from .meanline import MeanlineInput, build_curve

#: Grandeurs suivies, et la cle de confiance de chacune.
QUANTITIES = (("hauteur", "hauteur"), ("debit", "debit"), ("npshr", "npshr"))

#: Entrees perturbees, dans l'ordre des colonnes du tableau.
INPUTS = ("beta1", "beta2", "diametre")


@dataclass
class SensitivityRow:
    """Sensibilite d'une grandeur aux trois entrees, en relatif."""

    quantity: str = ""
    beta1: float = 0.0  # variation relative par degre
    beta2: float = 0.0  # variation relative par degre
    diameter: float = 0.0  # elasticite, sans dimension

    def worst_per_degree(self) -> float:
        """La plus grande des deux sensibilites angulaires, en valeur absolue."""
        values = [abs(v) for v in (self.beta1, self.beta2) if math.isfinite(v)]
        if not values:
            return math.inf
        return max(values)

    def alarming(self) -> bool:
        """Vrai si un degre suffit a deplacer la grandeur au-dela du seuil."""
        worst = self.worst_per_degree()
        return not math.isfinite(worst) or worst > config.SENSITIVITY_LOW_PER_DEG

    def to_dict(self) -> dict:
        """Vue serialisable ; l'infini n'etant pas du JSON, il devient `None`."""
        def clean(value: float) -> float | None:
            return value if math.isfinite(value) else None

        return {
            "grandeur": self.quantity,
            "par_degre_de_beta1": clean(self.beta1),
            "par_degre_de_beta2": clean(self.beta2),
            "elasticite_au_diametre": clean(self.diameter),
            "declassante": self.alarming(),
        }


@dataclass
class SensitivityReport:
    """Le tableau complet, et ce qu'il declasse."""

    rows: list[SensitivityRow] = field(default_factory=list)
    downgraded: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def row(self, quantity: str) -> SensitivityRow | None:
        """Ligne d'une grandeur, ou `None` si elle n'a pas ete calculee."""
        return next((r for r in self.rows if r.quantity == quantity), None)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "lignes": [row.to_dict() for row in self.rows],
            "grandeurs_declassees": list(self.downgraded),
            "seuil_par_degre": config.SENSITIVITY_LOW_PER_DEG,
        }


def _operating_values(data: MeanlineInput, rpm: float) -> tuple[float, float, float] | None:
    """Hauteur, debit et NPSHr au point nominal, ou `None` s'il n'existe pas."""
    curve = build_curve(data, rpm, ConfidenceMap())
    cavitation_module.apply_to_curve(curve)
    point = curve.nominal_point()
    if point is None or point.head <= 0.0 or point.flow <= 0.0:
        return None
    return point.head, point.flow, point.npshr


def _scaled(data: MeanlineInput, factor: float) -> MeanlineInput:
    """Copie de la geometrie a l'echelle `factor` : longueurs x f, sections x f2."""
    return dataclasses.replace(
        data,
        r_1=data.r_1 * factor,
        r_1s=data.r_1s * factor,
        r_1h=data.r_1h * factor,
        r_2=data.r_2 * factor,
        area_1=data.area_1 * factor * factor,
        area_2=data.area_2 * factor * factor,
    )


def _relative_spread(low, high, middle: float) -> float:
    """Ecart relatif entre deux perturbations, rapporte a la valeur centrale."""
    if low is None or high is None or middle <= 0.0:
        return math.inf
    return abs(high - low) / middle


def analyse(data: MeanlineInput, rpm: float) -> SensitivityReport:
    """Sensibilite de la hauteur, du debit et du NPSHr aux entrees fragiles.

    Differences centrees : chaque entree est perturbee des deux cotes, et
    l'ecart entre les deux resultats est rapporte a la valeur centrale.  Une
    perturbation qui fait **disparaitre** le point de fonctionnement rend une
    sensibilite infinie, ce qui est l'information juste : la grandeur n'est
    alors pas seulement imprecise, elle n'est pas definie.
    """
    report = SensitivityReport()
    centre = _operating_values(data, rpm)
    if centre is None:
        report.warnings.append(
            "pas de point de fonctionnement au regime de reference : la sensibilite des "
            "resultats a la geometrie n'a pas pu etre chiffree."
        )
        return report

    step = config.BETA_SENSITIVITY_DEG
    perturbations: dict[str, tuple] = {}
    for name, attribute in (("beta1", "beta1_deg"), ("beta2", "beta2_deg")):
        centred = getattr(data, attribute)
        low = high = None
        if 0.0 < centred - step and centred + step < 90.0:
            low = _operating_values(dataclasses.replace(
                data, **{attribute: centred - step}), rpm)
            high = _operating_values(dataclasses.replace(
                data, **{attribute: centred + step}), rpm)
        perturbations[name] = (low, high, 2.0 * step)

    delta = config.SENSITIVITY_DIAMETER_REL
    perturbations["diametre"] = (
        _operating_values(_scaled(data, 1.0 - delta), rpm),
        _operating_values(_scaled(data, 1.0 + delta), rpm),
        2.0 * delta,
    )

    for index, (quantity, key) in enumerate(QUANTITIES):
        row = SensitivityRow(quantity=quantity)
        for name in INPUTS:
            low, high, span = perturbations[name]
            spread = _relative_spread(
                low[index] if low else None,
                high[index] if high else None,
                centre[index],
            )
            setattr(row, "diameter" if name == "diametre" else name,
                    spread / span if math.isfinite(spread) else math.inf)
        report.rows.append(row)
        if row.alarming():
            report.downgraded.append(key)
    return report


def apply(report: SensitivityReport, confidence: ConfidenceMap) -> None:
    """Declasse en confiance faible les grandeurs trop sensibles.

    La puissance et le couple suivent la hauteur : ils en derivent
    directement, et les laisser en confiance haute pendant qu'elle tombe
    ferait mentir la table.
    """
    for key in report.downgraded:
        confidence.set(key, LOW)
        if key == "hauteur":
            for derived in ("puissance", "couple", "rendement"):
                confidence.set(derived, LOW)
