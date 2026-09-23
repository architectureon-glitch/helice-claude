"""Si l'helice etait normale : la meme geometrie, avec des aubes a un seul brin.

L'outil n'etudie que les helices toroidales, et chaque resultat vient avec son
comparatif : la meme piece -- memes rayons, memes sections, memes angles, meme
nombre de pales, meme regime --, dont chaque aube serait un seul brin au lieu
d'une boucle. C'est la forme de l'aube, et elle seule, qui change ; l'ecart
entre les deux colonnes est donc ce que la boucle apporte ou coute, **dans la
mesure ou le modele le voit**.

Ce que le modele voit, et ce qu'il ne voit pas, est publie avec le tableau.
Une ligne de moyenne ne connait pas la forme de l'aube : a geometrie egale elle
rend la meme courbe, et l'ecart d'une pompe se lit sur les pertes de canal --
la surface mouillee d'une aube en boucle compte ses deux brins. Pour une helice
libre, le bilan par element de pale voit ce qui fait l'interet de la boucle :
une pale normale perd de la portance pres de son bout libre, ou le fluide la
contourne ; une boucle n'a pas de bout libre.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import LOW, MEDIUM, worst
from . import cavitation as cavitation_module
from . import losses as losses_module
from . import meanline as meanline_module
from . import propulsion as propulsion_module


@dataclass
class ComparisonRow:
    """Une grandeur, pour la toroidale et pour la normale."""

    quantity: str = ""
    unit: str = ""
    toroidal: float | None = None
    normal: float | None = None
    cause: str = ""  # ce qui fait l'ecart, dans le modele
    toroidal_alt: float | None = None  # seconde estimation de la toroidale (helice : boucle fermee)

    @property
    def absolute(self) -> bool:
        """Un angle ou un rendement se compare par difference, pas par rapport."""
        return self.unit in ("deg", "%")

    @property
    def gap(self) -> float | None:
        """Ecart de la toroidale sur la normale : relatif, ou en unites pour un angle ou un rendement."""
        if self.toroidal is None or self.normal is None:
            return None
        if self.absolute:
            return self.toroidal - self.normal
        if self.normal == 0.0:
            return None
        return (self.toroidal - self.normal) / abs(self.normal)

    def to_dict(self) -> dict:
        return {
            "grandeur": self.quantity,
            "unite": self.unit,
            "toroidale": self.toroidal,
            "toroidale_seconde_estimation": self.toroidal_alt,
            "normale": self.normal,
            "ecart": self.gap,
            "ecart_en_unites": self.absolute,
            "cause": self.cause,
        }


@dataclass
class Comparison:
    """Le comparatif publie avec chaque resultat."""

    machine: str = ""  # "pompe" ou "helice_libre"
    case: str = ""  # le cas calcule, en clair
    equivalent: str = ""  # l'helice normale de reference, en clair
    columns: tuple[str, ...] = ("toroidale", "normale")
    columns_variant: tuple[str, str] = ("toroidale", "variante")
    rows: list[ComparisonRow] = field(default_factory=list)
    variant_title: str = ""  # variante supplementaire (pompe : a diametre egal)
    variant_rows: list[ComparisonRow] = field(default_factory=list)
    modelled: list[str] = field(default_factory=list)
    not_modelled: list[str] = field(default_factory=list)
    confidence: str = LOW
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "machine": self.machine,
            "cas_calcule": self.case,
            "helice_normale_de_reference": self.equivalent,
            "colonnes": list(self.columns),
            "lignes": [row.to_dict() for row in self.rows],
            "variante": self.variant_title or None,
            "colonnes_de_la_variante": list(self.columns_variant),
            "lignes_de_la_variante": [row.to_dict() for row in self.variant_rows],
            "ce_que_le_modele_compte": list(self.modelled),
            "ce_qu_il_ne_compte_pas": list(self.not_modelled),
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# Pompe
# ---------------------------------------------------------------------------
@dataclass
class PumpInlet:
    """Entree de la roue normale quand elle differe de celle de la toroidale.

    Sur une roue en serie (hel1), l'eau entre par la vis du brin amont ; la roue
    normale n'a pas ce brin, et l'eau l'aborde au bord d'attaque de son aube.
    """

    r_1: float = 0.0
    r_1h: float = 0.0
    r_1s: float = 0.0
    area_1: float = 0.0
    beta1_deg: float = 0.0
    description: str = ""


@dataclass
class PumpOutletRim:
    """La fente de sortie, quand le bord de fuite s'arrete avant elle."""

    radius: float = 0.0
    area: float = 0.0  # aire libre de la fente, obstruction des aubes non comprise


def _at(curve: meanline_module.PerformanceCurve, flow: float, attribute: str) -> float | None:
    """Valeur d'une grandeur de la courbe au debit `flow` (interpolation lineaire)."""
    points = [p for p in curve.points if p.head > 0.0]
    for low, high in zip(points, points[1:]):
        if low.flow <= flow <= high.flow and high.flow > low.flow:
            weight = (flow - low.flow) / (high.flow - low.flow)
            return getattr(low, attribute) * (1.0 - weight) + getattr(high, attribute) * weight
    if points and abs(points[0].flow - flow) <= 1e-12:
        return getattr(points[0], attribute)
    return None


def _curve(data: meanline_module.MeanlineInput, rpm: float) -> meanline_module.PerformanceCurve:
    curve = meanline_module.build_curve(data, rpm)
    cavitation_module.apply_to_curve(curve)
    return curve


def _channel(
    data: meanline_module.MeanlineInput, b_2: float, curve: meanline_module.PerformanceCurve,
    strands: int,
) -> losses_module.ChannelLosses | None:
    """Pertes de canal au meilleur rendement de la courbe, comme `_curve_analysis`."""
    best = curve.best_efficiency_point()
    if best is None or data.r_1 <= 0.0:
        return None
    width_1 = data.area_1 / (2.0 * math.pi * data.r_1)
    w_2 = math.hypot(best.cm2, max(0.0, curve.u2 - best.cu2))
    return losses_module.analyse(
        r_1=data.r_1, r_2=data.r_2, b_1=width_1, b_2=b_2,
        beta1_deg=data.beta1_deg, beta2_deg=data.beta2_deg, n_blades=data.n_blades,
        w1=best.w1, w2=w_2, head_theoretical=best.head_theoretical, strands=strands,
    )


def pump(
    data: meanline_module.MeanlineInput,
    curve: meanline_module.PerformanceCurve,
    b_2: float,
    strands: int,
    closed: bool | None,
    inlet: PumpInlet | None = None,
    rim: PumpOutletRim | None = None,
    confidence: str = LOW,
) -> Comparison | None:
    """Comparatif d'une pompe toroidale et de la meme roue a aubes normales.

    `data`, `curve` : le modele et la courbe de la toroidale au regime de
    reference. `strands` : brins par aube de la toroidale. `closed` : roue
    fermee (flasque avant), ouverte, ou inconnu. `inlet` : entree de la roue
    normale si elle differe (roue en serie). `rim` : la fente de sortie, si le
    bord de fuite s'arrete avant elle -- d'ou la variante a diametre egal.
    """
    nominal = curve.nominal_point() or curve.best_efficiency_point()
    if nominal is None or not curve.points:
        return None
    rpm = curve.rpm
    q_case = nominal.flow
    m3h = config.SECONDS_PER_HOUR
    result = Comparison(machine="pompe")
    result.case = (
        f"{rpm:.0f} tr/min, Q = {q_case * m3h:.1f} m3/h : le point d'incidence nulle de la "
        "toroidale"
    )

    normal_data = copy.copy(data)
    if inlet is not None:
        normal_data.r_1, normal_data.r_1h, normal_data.r_1s = inlet.r_1, inlet.r_1h, inlet.r_1s
        normal_data.area_1 = inlet.area_1
        normal_data.beta1_deg = inlet.beta1_deg
    normal = _curve(normal_data, rpm)
    result.equivalent = (
        f"la meme roue -- r1 = {normal_data.r_1 * config.MM_PER_M:.1f} mm, r2 = "
        f"{normal_data.r_2 * config.MM_PER_M:.1f} mm, beta1/beta2 = {normal_data.beta1_deg:.1f}/"
        f"{normal_data.beta2_deg:.1f} deg, {normal_data.n_blades} pales --, chaque aube faite "
        "d'un seul brin"
        + (f" ; {inlet.description}" if inlet is not None else "")
    )

    kw = 1e-3
    omega = meanline_module.angular_velocity(rpm)
    slip = meanline_module.slip_factor(data.beta2_deg, data.n_blades)
    sans_courbe = "" if normal.points else "sans point de fonctionnement (voir la reserve)"
    if sans_courbe:
        result.warnings.append(
            "la roue normale n'a pas de point de fonctionnement : "
            + _why_no_curve(normal_data, normal, omega, slip)
        )
    # Les pertes de la SPEC sont calees sur le point nominal de chaque roue : la
    # normale, dont le debit d'incidence nulle differe, en recevrait un autre
    # frottement, et paraitrait meilleure ou pire pour une raison qui n'est pas
    # physique. Les deux hauteurs sont donc calculees avec le frottement de la
    # toroidale ; l'ecart de frottement du aux brins est lu, lui, sur les pertes
    # de canal plus bas.
    friction = curve.friction_coefficient
    q_zero = curve.points[0].flow
    incidence_n = _incidence(normal_data, omega, q_case)
    has_normal = bool(normal.points)
    h_th_t = meanline_module.euler_head(data, omega, q_case, slip)[0]
    h_th_n = meanline_module.euler_head(normal_data, omega, q_case, slip)[0]
    p_t = _at(curve, q_case, "shaft_power")
    result.rows = [
        ComparisonRow("Debit d'incidence nulle", "m3/h", curve.flow_nominal * m3h,
                      normal.flow_nominal * m3h,
                      "angle et section d'entree" if inlet is not None else "identique par construction"),
        ComparisonRow("Incidence au bord d'attaque, au debit du cas", "deg",
                      _incidence(data, omega, q_case), incidence_n,
                      "beta1 contre angle de l'eau arrivant sans giration"),
        ComparisonRow("Hauteur d'Euler au debit du cas, avant pertes", "m", h_th_t, h_th_n,
                      "meme bord de fuite, meme glissement"),
        ComparisonRow("Hauteur au debit du cas", "m",
                      _head(data, omega, slip, q_case, friction),
                      _head(normal_data, omega, slip, q_case, friction) if has_normal else None,
                      sans_courbe or ("meme frottement ; l'ecart vient de l'incidence"
                                      if inlet is not None
                                      else "identique : meme Euler, memes pertes de la SPEC")),
        ComparisonRow("Hauteur a debit nul", "m",
                      _head(data, omega, slip, q_zero, friction),
                      _head(normal_data, omega, slip, q_zero, friction) if has_normal else None,
                      sans_courbe or "u2 et beta2 identiques"),
        ComparisonRow("Puissance a l'arbre au debit du cas", "kW", _scaled(p_t, kw),
                      _scaled(p_t * h_th_n / h_th_t, kw)
                      if has_normal and p_t is not None and h_th_t > 0.0 else None,
                      sans_courbe or "suit la hauteur d'Euler"),
        ComparisonRow("NPSH requis au debit du cas", "m", _at(curve, q_case, "npshr"),
                      _at(normal, q_case, "npshr"),
                      sans_courbe or ("entree differente" if inlet is not None
                                      else "identique : meme entree")),
    ]
    if has_normal and incidence_n is not None and abs(incidence_n) > config.COMPARISON_INCIDENCE_WARN_DEG:
        result.warnings.append(
            f"la roue normale travaille au debit du cas avec {incidence_n:.0f} deg d'incidence a son "
            "bord d'attaque. Le modele n'en compte que l'ecart de vitesse relative ; un tel angle "
            "fait vraisemblablement decoller l'ecoulement au bord d'attaque, ce qui coute davantage : "
            "sa hauteur au debit du cas est une borne haute."
        )

    toroidal_losses = _channel(data, b_2, curve, strands)
    normal_losses = _channel(normal_data, b_2, normal, 1) if normal.points else None
    if toroidal_losses is not None and normal_losses is not None:
        result.rows += [
            ComparisonRow("Surface mouillee des aubes et flasques", "cm2",
                          toroidal_losses.wetted_area * 1e4, normal_losses.wetted_area * 1e4,
                          f"{strands} brins par aube contre 1"),
            ComparisonRow("Perte de frottement de canal, au meilleur rendement", "m",
                          toroidal_losses.head_friction, normal_losses.head_friction,
                          "surface mouillee"),
            ComparisonRow("Rendement de comparaison, au meilleur rendement", "%",
                          toroidal_losses.efficiency * 100.0, normal_losses.efficiency * 100.0,
                          "frottement et diffusion du canal"),
        ]

    if rim is not None and rim.radius > data.r_2 and rim.area > 0.0:
        # A diametre egal : la question n'est plus la forme de l'aube mais son
        # etendue. La toroidale, bord de fuite prolonge jusqu'a la fente, entree
        # inchangee.
        rim_data = copy.copy(data)
        rim_data.r_2 = rim.radius
        rim_data.area_2 = rim.area * config.TAU_2
        rimmed = _curve(rim_data, rpm)
        result.variant_title = (
            f"a diametre egal : la meme roue, bord de fuite prolonge jusqu'a la fente de sortie "
            f"(r2 = {rim.radius * config.MM_PER_M:.1f} mm au lieu de "
            f"{data.r_2 * config.MM_PER_M:.1f}), entree inchangee"
        )
        result.columns_variant = ("toroidale", "bord de fuite a la fente")
        result.variant_rows = [
            ComparisonRow("Hauteur au debit du cas", "m", _at(curve, q_case, "head"),
                          _at(rimmed, q_case, "head"), "u2 = omega r2"),
            ComparisonRow("Hauteur a debit nul", "m", curve.points[0].head,
                          rimmed.points[0].head if rimmed.points else None, "u2 au carre"),
            ComparisonRow("Puissance a l'arbre au debit du cas", "kW",
                          _scaled(_at(curve, q_case, "shaft_power"), kw),
                          _scaled(_at(rimmed, q_case, "shaft_power"), kw), "suit la hauteur"),
        ]

    result.modelled.append(
        f"la surface mouillee : une aube en boucle presente ses {strands} brins au frottement, "
        "la normale un seul, pour le meme canal (pertes de canal, rendement de comparaison)"
        if strands > 1 else
        "la surface mouillee (pertes de canal, rendement de comparaison)"
    )
    if inlet is not None:
        result.modelled.append(
            f"l'entree : la toroidale aspire par la vis de son brin amont (beta1 = "
            f"{data.beta1_deg:.1f} deg), la normale par le bord d'attaque de son aube (beta1 = "
            f"{inlet.beta1_deg:.1f} deg) -- debit d'incidence nulle et NPSH requis"
        )
    else:
        result.modelled.append(
            "rien d'autre : la ligne moyenne ne connait pas la forme de l'aube, et rend la meme "
            "courbe pour les deux roues ; c'est le rendement de comparaison qui porte l'ecart"
        )
    if closed:
        result.modelled.append(
            "le bout de pale : la roue est fermee, son flasque supprime deja le tourbillon et la "
            "fuite en bout de pale -- pour les deux roues. La boucle n'y apporte rien."
        )
    elif closed is False:
        result.not_modelled.append(
            "la fuite en bout de pale : la roue normale, ouverte, fuit par-dessus ses aubes vers le "
            "carter ; la boucle n'a pas de bout libre. L'ecart, favorable a la toroidale, depend "
            "du jeu au carter et n'est pas chiffre."
        )
    else:
        result.not_modelled.append(
            "la fuite en bout de pale : selon que la roue est ouverte ou fermee, la boucle la "
            "supprime ou n'y change rien ; l'outil ne l'a pas etabli ici."
        )
    result.not_modelled.append(
        "l'interaction entre les brins, les pertes de leur jonction"
        + (" et du coude entre les deux brins de la roue en serie" if inlet is not None else "")
        + ", le bruit, la tenue aux corps etrangers"
    )
    result.confidence = worst(confidence, MEDIUM)
    return result


def _scaled(value: float | None, factor: float) -> float | None:
    return None if value is None else value * factor


def _head(
    data: meanline_module.MeanlineInput, omega: float, slip: float, flow: float, friction: float,
) -> float:
    """Hauteur au debit `flow`, pertes de la SPEC avec un frottement impose.

    Meme formule que `meanline._operating_point` : Euler avec glissement, moins
    le frottement `friction . Q2`, moins l'incidence, comptee depuis le debit
    d'incidence nulle de la roue elle-meme.
    """
    head_theoretical = meanline_module.euler_head(data, omega, flow, slip)[0]
    u1 = omega * data.r_1
    w1 = math.hypot(flow / data.area_1, u1)
    w1_zero = math.hypot(u1 * math.tan(math.radians(data.beta1_deg)), u1)
    return (head_theoretical - friction * flow ** 2
            - config.XI_INCIDENCE * (w1 - w1_zero) ** 2 / (2.0 * config.G))


def _incidence(data: meanline_module.MeanlineInput, omega: float, flow: float) -> float | None:
    """beta1 moins l'angle de l'eau qui arrive sans giration au bord d'attaque, en degres."""
    if data.area_1 <= 0.0 or data.r_1 <= 0.0 or omega <= 0.0:
        return None
    cm1 = flow / data.area_1
    return data.beta1_deg - math.degrees(math.atan2(cm1, omega * data.r_1))


def _why_no_curve(
    data: meanline_module.MeanlineInput, curve: meanline_module.PerformanceCurve,
    omega: float, slip: float,
) -> str:
    """Pourquoi le modele ne trouve pas de point de fonctionnement a une roue."""
    if not data.valid():
        return "geometrie hors du domaine du modele."
    u2 = omega * data.r_2
    zero_head = slip * u2 * math.tan(math.radians(data.beta2_deg)) * data.area_2
    m3h = config.SECONDS_PER_HOUR
    if curve.flow_nominal > zero_head > 0.0:
        return (
            f"son bord d'attaque (beta1 = {data.beta1_deg:.1f} deg) est adapte a "
            f"{curve.flow_nominal * m3h:.0f} m3/h, debit auquel son bord de fuite (beta2 = "
            f"{data.beta2_deg:.1f} deg) ne donne plus de hauteur -- il n'en donne que jusqu'a "
            f"{zero_head * m3h:.0f} m3/h. Entree et sortie ne sont pas accordees."
        )
    return "le modele ne trouve pas de hauteur positive au debit d'incidence nulle."


# ---------------------------------------------------------------------------
# Helice libre
# ---------------------------------------------------------------------------
def propeller(
    prop: propulsion_module.PropulsionResult,
    topology,
    geometry,
    confidence: str = LOW,
) -> Comparison | None:
    """Comparatif d'une helice libre toroidale et de la meme helice a pales normales.

    La toroidale porte `S` brins par boucle, chacun avec la section lue sur la
    piece ; la normale porte une pale par boucle, de la meme section. La
    toroidale est donnee en deux estimations : perte de bout de Prandtl
    conservee, comme si chaque brin avait un bout libre, et supprimee, la boucle
    n'en ayant pas. La realite est entre les deux, et seul un essai la situe.
    """
    if prop.point is None or not prop.stations or prop.diameter <= 0.0:
        return None
    rpm, rho, speed = prop.rpm, prop.rho, prop.point.speed
    sound = propulsion_module.FLUIDS.get(prop.fluid, propulsion_module.FLUIDS["eau"])[1]
    n_loops = max(1, topology.blades.n_blades)
    n_strands = max(1, prop.n_blades)
    strands = max(1, round(n_strands / n_loops))

    def point(stations, v, n_blades, closed_tip=False):
        return propulsion_module.operating_point(
            copy.deepcopy(stations), v, rpm, rho, prop.diameter, prop.hub_diameter, n_blades,
            sound, closed_tip=closed_tip,
        )

    normal_geometry = copy.copy(geometry)
    normal_geometry.n_effective_blades = n_loops
    normal_stations = propulsion_module.blade_stations(topology, normal_geometry)
    if not normal_stations:
        return None
    runs = {
        "ouverte": (prop.stations, n_strands, False),
        "fermee": (prop.stations, n_strands, True),
        "normale": (normal_stations, n_loops, False),
    }
    at_speed = {k: point(st, speed, n, c) for k, (st, n, c) in runs.items()}
    static = {k: point(st, 0.0, n, c) for k, (st, n, c) in runs.items()}

    result = Comparison(machine="helice_libre")
    result.columns = ("toroidale, bout conserve", "toroidale, boucle fermee", "normale")
    result.case = (
        f"{rpm:.0f} tr/min, vitesse d'avance {speed:.2f} m/s"
        + (" (a l'arret)" if speed <= 0.0 else f", J = {at_speed['ouverte'].advance_ratio:.2f}")
        + f", {prop.fluid}"
    )
    result.equivalent = (
        f"la meme helice -- diametre {prop.diameter * config.MM_PER_M:.0f} mm, memes cordes, "
        f"memes calages et cambrures lus sur la piece --, avec {n_loops} pales a bout libre au "
        f"lieu de {n_loops} boucles"
        + (f" de {strands} brins" if strands > 1 else "")
    )

    def row(quantity, unit, getter, cause, source=at_speed):
        return ComparisonRow(quantity, unit, getter(source["ouverte"]), getter(source["normale"]),
                             cause, toroidal_alt=getter(source["fermee"]))

    kw = 1e-3
    result.rows = [
        row("Poussee", "N", lambda p: p.thrust, "surfaces portantes et perte de bout"),
        row("Puissance absorbee", "kW", lambda p: p.power * kw, "surfaces portantes"),
    ]
    if speed > 0.0:
        result.rows += [
            row("Rendement propulsif", "%", lambda p: p.efficiency * 100.0,
                "perte de bout, trainee des brins"),
            row("Plafond de Froude a meme poussee", "%", lambda p: p.froude_efficiency * 100.0,
                "charge du disque"),
        ]
    result.rows += [
        row("Poussee a l'arret", "N", lambda p: p.thrust, "surfaces portantes et perte de bout",
            source=static),
        row("Poussee a l'arret par kilowatt", "N/kW",
            lambda p: p.thrust / (p.power * kw) if p.power > 0.0 else None,
            "ce que coute la poussee", source=static),
    ]

    if strands > 1:
        result.modelled.append(
            f"les surfaces portantes : chaque boucle porte {strands} brins, {n_strands} en tout, "
            f"contre {n_loops} pales pour la normale -- solidite multipliee par {strands}, poussee "
            "et couple plus forts au meme regime"
        )
    else:
        result.not_modelled.append(
            "le second brin de chaque boucle : la lecture ne separe pas les deux brins (une seule "
            f"famille de profils par coupe), et le bilan compte {n_loops} surfaces portantes pour "
            "la toroidale comme pour la normale. L'ecart ne porte que sur le bout de pale"
        )
    result.modelled += [
        "le bout de pale : la normale perd de la portance pres de son bout libre, ou le fluide la "
        "contourne (facteur de Prandtl) ; la boucle n'a pas de bout libre. La colonne « boucle "
        "fermee » supprime cette perte, la colonne « bout conserve » la garde : la realite est "
        "entre les deux",
    ]
    result.not_modelled += [
        "la jonction des brins au bout de la boucle : sa trainee, et la facon dont elle charge "
        "l'ecoulement",
        "l'interaction entre les brins d'une meme boucle, au-dela de ce que le bilan de quantite "
        "de mouvement en dit",
        "le bruit et le tourbillon de bout de pale, dont la suppression est l'argument premier "
        "des helices toroidales",
    ]
    result.confidence = worst(confidence, MEDIUM)
    return result
