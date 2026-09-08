"""Helice libre : poussee, couple et rendement propulsif par element de pale.

Ce module repond a une question que le reste de l'outil ne pose pas.  Les phases
5 et 6 traitent la roue comme une **pompe** : elle est carenee, elle refoule
dans une volute, et ce qu'on lui demande est une hauteur et un debit.  Une meme
geometrie axiale peut aussi tourner en **helice libre** -- bateau, drone, banc
d'essai a l'air libre -- et la question devient alors une poussee et un
rendement propulsif.  Ce sont deux regimes differents, pas deux facons de
regarder le meme, et le module refuse de repondre pour une roue qui n'est pas
axiale.

**Le modele.**  Bilan par element de pale et quantite de mouvement (BEM), dans
la formulation en **vitesses induites** plutot qu'en facteurs d'induction : a
chaque station radiale, la poussee que rend la theorie de la quantite de
mouvement est egalee a celle que rend la portance de l'element de pale, et le
couple de meme.  La formulation en vitesses induites tient a l'arret -- ou
`V = 0` fait s'evanouir le facteur d'induction axial -- ce qui donne la poussee
statique sans formule empirique separee.

**Ce que la geometrie fournit.**  L'angle de calage, la corde, l'epaisseur
relative et la **fleche de cambrure** sont lus sur la piece, coupe par coupe.
C'est la difference avec les modeles d'helice courants, qui demandent a
l'utilisateur de taper un coefficient de portance de dessin et une trainee de
profil qu'il ne connait pas : ici l'angle de portance nulle sort de la ligne de
cambrure mesuree (`alpha_0 = -2 h/c`, profils minces) et la trainee de base de
l'epaisseur relative mesuree.

**Le plafond.**  Chaque rendement calcule est publie a cote du rendement d'un
disque actif ideal de meme poussee (Froude) : c'est le maximum que la quantite
de mouvement autorise, pertes de profil, de bout de pale et de giration mises a
zero.  Un rendement calcule qui le depasserait serait une erreur de programme,
pas une bonne helice -- l'ecart entre les deux dit ou passe l'energie.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import HIGH, LOW, MEDIUM
from ..geometry.blade_angles import BladeGeometry, SectionAngles
from ..geometry.topology import AXIAL, Topology


@dataclass
class BladeStation:
    """Une station radiale de la pale, et ce que le BEM y a converge."""

    radius: float = 0.0  # m
    chord: float = 0.0  # m
    solidity: float = 0.0  # - , sigma = B c / (2 pi r)
    theta_deg: float = 0.0  # calage local, depuis le plan de rotation
    camber_ratio: float = 0.0  # fleche / corde
    thickness_ratio: float = 0.0  # epaisseur max / corde
    phi_deg: float = 0.0  # angle d'ecoulement incident
    alpha_deg: float = 0.0  # incidence
    cl: float = 0.0
    cd: float = 0.0
    prandtl: float = 1.0  # facteur de perte de bout et de pied
    axial_induced: float = 0.0  # m/s
    swirl_induced: float = 0.0  # m/s
    thrust_gradient: float = 0.0  # N/m
    torque_gradient: float = 0.0  # N.m/m
    stalled: bool = False
    converged: bool = True

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "rayon_m": self.radius,
            "corde_m": self.chord,
            "solidite": self.solidity,
            "calage_deg": self.theta_deg,
            "cambrure_relative": self.camber_ratio,
            "epaisseur_relative": self.thickness_ratio,
            "phi_deg": self.phi_deg,
            "incidence_deg": self.alpha_deg,
            "cl": self.cl,
            "cd": self.cd,
            "perte_de_prandtl": self.prandtl,
            "vitesse_induite_axiale_m_s": self.axial_induced,
            "vitesse_induite_tangentielle_m_s": self.swirl_induced,
            "dT_dr_N_m": self.thrust_gradient,
            "dQ_dr_Nm_m": self.torque_gradient,
            "decroche": self.stalled,
            "converge": self.converged,
        }


@dataclass
class PropulsionPoint:
    """Un point de fonctionnement de l'helice, tout en SI."""

    speed: float = 0.0  # m/s, vitesse d'avance
    advance_ratio: float = 0.0  # J = V / (n D)
    thrust: float = 0.0  # N
    power: float = 0.0  # W
    torque: float = 0.0  # N.m
    thrust_coefficient: float = 0.0  # CT = T / (rho n2 D4)
    power_coefficient: float = 0.0  # CP = P / (rho n3 D5)
    torque_coefficient: float = 0.0  # CQ = CP / 2 pi
    efficiency: float = 0.0  # eta = J CT / CP
    froude_efficiency: float = 0.0  # plafond du disque actif de meme poussee
    tip_speed: float = 0.0  # m/s
    tip_mach: float = 0.0  # - , nul en liquide
    stalled_fraction: float = 0.0  # part de l'envergure decrochee

    def ceiling_gap(self) -> float:
        """Points de rendement perdus par rapport au disque actif ideal."""
        if self.froude_efficiency <= 0.0 or self.efficiency <= 0.0:
            return 0.0
        return self.froude_efficiency - self.efficiency

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "vitesse_d_avance_m_s": self.speed,
            "parametre_d_avance_J": self.advance_ratio,
            "poussee_N": self.thrust,
            "puissance_W": self.power,
            "couple_Nm": self.torque,
            "CT": self.thrust_coefficient,
            "CP": self.power_coefficient,
            "CQ": self.torque_coefficient,
            "rendement_propulsif": self.efficiency,
            "rendement_ideal_de_froude": self.froude_efficiency,
            "ecart_au_plafond": self.ceiling_gap(),
            "vitesse_en_bout_de_pale_m_s": self.tip_speed,
            "mach_en_bout_de_pale": self.tip_mach or None,
            "part_d_envergure_decrochee": self.stalled_fraction,
        }


@dataclass
class PropulsionResult:
    """Analyse propulsive complete a un regime donne."""

    rpm: float = 0.0
    rho: float = 0.0
    fluid: str = "eau"
    diameter: float = 0.0
    hub_diameter: float = 0.0
    n_blades: int = 0
    disk_area: float = 0.0
    stations: list[BladeStation] = field(default_factory=list)
    point: PropulsionPoint | None = None
    curve: list[PropulsionPoint] = field(default_factory=list)
    static_thrust: float = 0.0  # N, poussee a l'arret
    zero_thrust_advance: float = 0.0  # J ou la poussee s'annule
    confidence: str = HIGH
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "regime_tr_min": self.rpm,
            "fluide": self.fluid,
            "rho_kg_m3": self.rho,
            "diametre_m": self.diameter,
            "diametre_de_moyeu_m": self.hub_diameter,
            "nombre_de_pales": self.n_blades,
            "surface_de_disque_m2": self.disk_area,
            "poussee_statique_N": self.static_thrust,
            "J_de_poussee_nulle": self.zero_thrust_advance or None,
            "point": self.point.to_dict() if self.point else None,
            "courbe": [point.to_dict() for point in self.curve],
            "stations": [station.to_dict() for station in self.stations],
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# Le profil : ce que la geometrie dit de sa portance et de sa trainee
# ---------------------------------------------------------------------------
def zero_lift_angle(camber_ratio: float) -> float:
    """Angle de portance nulle d'un profil cambre, en radians.

    Theorie des profils minces pour une ligne de cambrure en arc de cercle :
    `alpha_0 = -2 h/c`, ou `h/c` est la fleche relative.  Negatif : un profil
    cambre porte deja a incidence nulle.  La fleche est **mesuree** sur la ligne
    de cambrure extraite, la ou un modele d'helice courant demande un Cl de
    dessin a l'utilisateur.
    """
    return -config.CAMBER_ALPHA0_FACTOR * camber_ratio


def lift_coefficient(alpha: float, camber_ratio: float) -> tuple[float, bool]:
    """Portance du profil a l'incidence `alpha` (rad) ; dit aussi s'il decroche."""
    raw = config.CL_ALPHA * (alpha - zero_lift_angle(camber_ratio))
    clipped = max(config.CL_MIN, min(config.CL_STALL, raw))
    return clipped, clipped != raw


def drag_coefficient(cl: float, thickness_ratio: float) -> float:
    """Trainee du profil : polaire quadratique sur une trainee de base.

    La trainee de base croit avec l'epaisseur relative, elle aussi **mesuree**.
    """
    cd0 = config.CD_BASE + config.CD_THICKNESS_K * thickness_ratio
    return cd0 + config.CD_INDUCED_K * (cl - config.CD_MIN_DRAG_CL) ** 2


def prandtl_loss(radius: float, r_hub: float, r_tip: float, phi: float, n_blades: int) -> float:
    """Facteur de perte de Prandtl, en bout de pale et en pied.

    Une pale n'est pas un disque : pres du bout, le fluide contourne l'extremite
    et la portance s'y annule.  Prandtl en donne le facteur correcteur ; le
    meme raisonnement vaut au pied, contre le moyeu.
    """
    sine = abs(math.sin(phi))
    if sine < config.PRANDTL_LOSS_MIN or radius <= 0.0:
        sine = config.PRANDTL_LOSS_MIN
    factor = 1.0
    for distance in (r_tip - radius, radius - r_hub):
        exponent = 0.5 * n_blades * max(0.0, distance) / (radius * sine)
        factor *= (2.0 / math.pi) * math.acos(min(1.0, math.exp(-exponent)))
    return max(config.PRANDTL_LOSS_MIN, factor)


# ---------------------------------------------------------------------------
# Le bilan a une station
# ---------------------------------------------------------------------------
def _residual(
    phi: float, station: BladeStation, speed: float, omega: float,
    r_hub: float, r_tip: float, n_blades: int,
) -> tuple[float, float, float, float, float, bool]:
    """Residu du bilan a l'angle d'ecoulement `phi`, et l'etat qui va avec.

    Les deux egalites -- quantite de mouvement et element de pale -- se ramenent
    a une seule equation en `phi`.  En posant
    `k = sigma Cn / (4 F sin2 phi)` et `k' = sigma Ct / (4 F sin phi cos phi)`,
    la vitesse axiale au disque vaut `u_a = V / (1 - k)` et la vitesse
    tangentielle `u_t = Omega r / (1 + k')`, tandis que la definition de `phi`
    impose `u_a = u_t tan phi`.  Le residu est la difference des deux, multipliee
    par `(1 - k)` :

        R(phi) = tan(phi) . Omega r (1 - k) / (1 + k') - V

    Cette forme-la vaut aussi **a l'arret** : a `V = 0` elle se reduit a `k = 1`,
    qui est exactement le bilan statique `sin2 phi = sigma Cn / (4 F)`.  Un point
    fixe sur les vitesses induites, lui, y perd sa solution -- et oscille des que
    la pale est solide.
    """
    theta = math.radians(station.theta_deg)
    alpha = theta - phi
    cl, stalled = lift_coefficient(alpha, station.camber_ratio)
    cd = drag_coefficient(cl, station.thickness_ratio)
    loss = prandtl_loss(station.radius, r_hub, r_tip, phi, n_blades)
    normal = cl * math.cos(phi) - cd * math.sin(phi)
    tangential = cl * math.sin(phi) + cd * math.cos(phi)

    sine, cosine = math.sin(phi), math.cos(phi)
    k = station.solidity * normal / (4.0 * loss * sine * sine)
    # (1 + k') est le rapport entre vitesse d'entrainement et vitesse
    # tangentielle vue par la pale. Le laisser passer sous zero ferait tourner
    # la giration induite plus vite que la pale, a l'envers : le residu y partait
    # a l'infini et la bissection ne demarrait pas.
    k_swirl = max(
        station.solidity * tangential / (4.0 * loss * sine * cosine),
        config.BEM_SWIRL_FLOOR - 1.0,
    )
    residual = math.tan(phi) * omega * station.radius * (1.0 - k) / (1.0 + k_swirl) - speed
    return residual, cl, cd, loss, k_swirl, stalled


def solve_station(
    station: BladeStation,
    speed: float,
    omega: float,
    r_hub: float,
    r_tip: float,
    n_blades: int,
    rho: float,
) -> BladeStation:
    """Egalise quantite de mouvement et element de pale a une station radiale.

    L'inconnue est l'**angle d'ecoulement** `phi`, et non le couple de vitesses
    induites : le bilan se ramene a une equation scalaire dont le residu change
    de signe une fois sur `(0, pi/2)`, et une bissection la resout sans jamais
    diverger.  Un point fixe sur les vitesses induites, lui, oscille des que la
    solidite depasse quelques dixiemes -- ce qui est le cas de toute roue de
    pompe axiale.
    """
    low, high = config.BEM_PHI_MIN, config.BEM_PHI_MAX
    r_low = _residual(low, station, speed, omega, r_hub, r_tip, n_blades)[0]
    r_high = _residual(high, station, speed, omega, r_hub, r_tip, n_blades)[0]

    station.converged = False
    phi = high
    if math.isfinite(r_low) and math.isfinite(r_high) and r_low * r_high <= 0.0:
        for _ in range(config.BEM_MAX_ITERATIONS):
            phi = 0.5 * (low + high)
            middle = _residual(phi, station, speed, omega, r_hub, r_tip, n_blades)[0]
            if not math.isfinite(middle):
                break
            if r_low * middle <= 0.0:
                high = phi
            else:
                low, r_low = phi, middle
            if high - low < config.BEM_TOLERANCE:
                station.converged = True
                break

    _, cl, cd, loss, k_swirl, stalled = _residual(
        phi, station, speed, omega, r_hub, r_tip, n_blades
    )
    station.phi_deg = math.degrees(phi)
    station.alpha_deg = math.degrees(math.radians(station.theta_deg) - phi)
    station.cl, station.cd, station.prandtl = cl, cd, loss
    station.stalled = stalled

    u_tangential = omega * station.radius / (1.0 + k_swirl) if k_swirl > -1.0 else 0.0
    u_axial = math.tan(phi) * u_tangential
    station.axial_induced = u_axial - speed
    station.swirl_induced = omega * station.radius - u_tangential
    if u_axial <= 0.0 or u_tangential <= 0.0 or not station.converged:
        station.thrust_gradient = station.torque_gradient = 0.0
        return station

    relative_squared = u_axial ** 2 + u_tangential ** 2
    normal = cl * math.cos(phi) - cd * math.sin(phi)
    tangential = cl * math.sin(phi) + cd * math.cos(phi)
    force = 0.5 * rho * n_blades * station.chord * relative_squared
    station.thrust_gradient = force * normal
    station.torque_gradient = force * tangential * station.radius
    return station



# ---------------------------------------------------------------------------
# De la geometrie lue aux stations radiales
# ---------------------------------------------------------------------------
def _interpolate(sections: list[SectionAngles], radius: float, attribute: str) -> float:
    """Valeur d'une grandeur de coupe au rayon `radius`, par interpolation lineaire.

    Au-dela de la premiere et de la derniere coupe, la valeur est prolongee par
    la plus proche : extrapoler un calage de pale hors de la zone mesuree
    inventerait de la geometrie.
    """
    usable = sorted(
        (s for s in sections if s.radius > 0.0 and getattr(s, attribute) > 0.0),
        key=lambda s: s.radius,
    )
    if not usable:
        return 0.0
    if radius <= usable[0].radius:
        return getattr(usable[0], attribute)
    if radius >= usable[-1].radius:
        return getattr(usable[-1], attribute)
    for low, high in zip(usable, usable[1:]):
        if low.radius <= radius <= high.radius:
            span = high.radius - low.radius
            if span <= 0.0:
                return getattr(low, attribute)
            weight = (radius - low.radius) / span
            return getattr(low, attribute) * (1.0 - weight) + getattr(high, attribute) * weight
    return getattr(usable[-1], attribute)


def blade_stations(
    topology: Topology, geometry: BladeGeometry, count: int = config.BEM_STATIONS
) -> list[BladeStation]:
    """Decoupe la pale en stations radiales, chacune portant sa geometrie lue.

    Les bouts sont retires de `BEM_ROOT_CUTOFF` : au pied et en bout, la corde
    mesuree tend vers zero et la portance avec elle, mais la carte d'occupation
    y est aussi la moins sure.  Mieux vaut ne pas y integrer que d'y integrer du
    bruit.
    """
    r_hub = topology.r_1h
    r_tip = topology.r_blade_tip or topology.r_1s
    if r_tip <= 0.0 or r_tip <= r_hub or not geometry.sections:
        return []
    margin = config.BEM_ROOT_CUTOFF * (r_tip - r_hub)
    low, high = r_hub + margin, r_tip - margin
    n_blades = max(1, geometry.n_effective_blades or topology.blades.n_blades)

    stations = []
    for index in range(count):
        span = (index + 0.5) / count
        radius = low + span * (high - low)
        chord = _interpolate(geometry.sections, radius, "chord")
        if chord <= 0.0:
            continue
        thickness = _interpolate(geometry.sections, radius, "max_thickness")
        stations.append(BladeStation(
            radius=radius,
            chord=chord,
            solidity=n_blades * chord / (2.0 * math.pi * radius),
            theta_deg=_interpolate(geometry.sections, radius, "beta_mean_deg"),
            camber_ratio=_interpolate(geometry.sections, radius, "camber_ratio"),
            thickness_ratio=thickness / chord if chord > 0.0 else 0.0,
        ))
    return stations


# ---------------------------------------------------------------------------
# Le point de fonctionnement, et son plafond
# ---------------------------------------------------------------------------
def froude_efficiency(thrust: float, rho: float, speed: float, disk_area: float) -> float:
    """Rendement d'un disque actif ideal rendant la meme poussee (Froude).

    `eta = 2 / (1 + sqrt(1 + T / (0.5 rho V2 A)))`.  C'est le maximum que la
    conservation de la quantite de mouvement autorise : pertes de profil, de
    bout de pale et de giration toutes mises a zero.  Aucune helice ne peut le
    depasser, et l'ecart entre lui et le rendement calcule dit ce que la pale
    perd -- non pas ce qu'elle vaut dans l'absolu, mais ce qu'il reste a gagner.
    """
    if speed < config.FROUDE_MIN_SPEED or thrust <= 0.0 or disk_area <= 0.0 or rho <= 0.0:
        return 0.0
    loading = thrust / (0.5 * rho * speed ** 2 * disk_area)
    return 2.0 / (1.0 + math.sqrt(1.0 + loading))


def operating_point(
    stations: list[BladeStation],
    speed: float,
    rpm: float,
    rho: float,
    diameter: float,
    hub_diameter: float,
    n_blades: int,
    speed_of_sound: float = 0.0,
) -> PropulsionPoint:
    """Poussee, couple et rendement de l'helice a une vitesse d'avance donnee."""
    point = PropulsionPoint(speed=speed)
    if not stations or rpm <= 0.0 or diameter <= 0.0:
        return point

    omega = rpm * config.RPM_TO_RAD_S
    r_hub, r_tip = 0.5 * hub_diameter, 0.5 * diameter
    width = (r_tip - r_hub) * (1.0 - 2.0 * config.BEM_ROOT_CUTOFF) / len(stations)

    thrust = torque = 0.0
    stalled = 0
    for station in stations:
        solved = solve_station(station, speed, omega, r_hub, r_tip, n_blades, rho)
        thrust += solved.thrust_gradient * width
        torque += solved.torque_gradient * width
        stalled += 1 if solved.stalled else 0

    n = rpm / config.SECONDS_PER_MINUTE
    point.thrust = thrust
    point.torque = torque
    point.power = torque * omega
    point.advance_ratio = speed / (n * diameter) if n > 0.0 else 0.0
    point.thrust_coefficient = thrust / (rho * n ** 2 * diameter ** 4)
    point.power_coefficient = point.power / (rho * n ** 3 * diameter ** 5)
    point.torque_coefficient = point.power_coefficient / (2.0 * math.pi)
    if point.power > 0.0 and thrust > 0.0:
        point.efficiency = thrust * speed / point.power
    disk_area = math.pi * (r_tip ** 2 - r_hub ** 2)
    point.froude_efficiency = froude_efficiency(thrust, rho, speed, disk_area)
    tip_rotational = omega * r_tip
    point.tip_speed = math.hypot(speed, tip_rotational)
    point.tip_mach = point.tip_speed / speed_of_sound if speed_of_sound > 0.0 else 0.0
    point.stalled_fraction = stalled / len(stations)
    return point


# ---------------------------------------------------------------------------
# L'analyse complete
# ---------------------------------------------------------------------------
FLUIDS = {
    # nom -> (masse volumique kg/m3, celerite du son m/s ; 0 si sans objet)
    "eau": (config.RHO, 0.0),
    "air": (config.RHO_AIR, config.SPEED_OF_SOUND_AIR),
}


def analyse(
    topology: Topology,
    geometry: BladeGeometry,
    rpm: float,
    speed: float,
    fluid: str = "eau",
) -> PropulsionResult:
    """Analyse propulsive d'une helice libre a `rpm` et a la vitesse `speed`.

    Rend aussi la courbe complete `poussee(J)`, de l'arret au parametre d'avance
    ou la poussee s'annule : c'est elle qui situe le point demande, et c'est sur
    elle que se lit le rendement maximal atteignable par cette pale.
    """
    rho, speed_of_sound = FLUIDS.get(fluid, FLUIDS["eau"])
    result = PropulsionResult(rpm=rpm, rho=rho, fluid=fluid)

    r_tip = topology.r_blade_tip or topology.r_1s
    result.diameter = 2.0 * r_tip
    result.hub_diameter = 2.0 * topology.r_1h
    result.n_blades = max(1, geometry.n_effective_blades or topology.blades.n_blades)
    result.disk_area = math.pi * (r_tip ** 2 - topology.r_1h ** 2)

    if topology.machine_type != AXIAL:
        result.confidence = LOW
        result.warnings.append(
            f"la roue est de type {topology.machine_type} : le modele d'helice libre ne "
            "s'applique pas. Il suppose un ecoulement axial traversant un disque non carene, "
            "alors qu'une roue centrifuge ou mixte refoule radialement dans une volute. "
            "Poussee et rendement propulsif ne sont pas calcules."
        )
        return result

    result.stations = blade_stations(topology, geometry)
    if not result.stations:
        result.confidence = LOW
        result.warnings.append(
            "aucune station de pale exploitable : corde ou rayons manquants. La poussee n'est "
            "pas calculee."
        )
        return result

    def point_at(v: float) -> PropulsionPoint:
        return operating_point(
            result.stations, v, rpm, rho, result.diameter, result.hub_diameter,
            result.n_blades, speed_of_sound,
        )

    result.static_thrust = point_at(0.0).thrust

    # Balayage en J jusqu'a la poussee nulle. La borne haute vient du calage :
    # une pale calee a theta n'avance pas plus vite que son pas geometrique.
    n = rpm / config.SECONDS_PER_MINUTE
    theta_tip = math.radians(max(s.theta_deg for s in result.stations))
    pitch = 2.0 * math.pi * (0.5 * result.diameter) * math.tan(theta_tip)
    j_max = config.ADVANCE_SWEEP_MARGIN * max(0.2, pitch / result.diameter)
    for index in range(config.ADVANCE_SWEEP_POINTS):
        advance = j_max * index / (config.ADVANCE_SWEEP_POINTS - 1)
        result.curve.append(point_at(advance * n * result.diameter))
    for previous, current in zip(result.curve, result.curve[1:]):
        if previous.thrust > 0.0 >= current.thrust:
            slope = previous.thrust - current.thrust
            weight = previous.thrust / slope if slope > 0.0 else 0.0
            result.zero_thrust_advance = previous.advance_ratio + weight * (
                current.advance_ratio - previous.advance_ratio
            )
            break

    result.point = point_at(speed)
    _add_notes(result)
    return result


def _add_notes(result: PropulsionResult) -> None:
    """Reserves et remarques qu'appelle le point de fonctionnement obtenu."""
    point = result.point
    if point is None:
        return

    if point.thrust <= 0.0 and point.speed > 0.0:
        result.confidence = LOW
        limite = (
            f" La poussee s'annule vers J = {result.zero_thrust_advance:.2f}, soit "
            f"{result.zero_thrust_advance * (result.rpm / config.SECONDS_PER_MINUTE) * result.diameter * 3.6:.0f} km/h "
            f"a ce regime." if result.zero_thrust_advance > 0.0 else ""
        )
        result.warnings.append(
            f"poussee negative a {point.speed:.1f} m/s : a ce regime et a cette vitesse "
            f"d'avance, l'helice freine au lieu de propulser (J = {point.advance_ratio:.2f})."
            + limite
            + " Verifiez le regime, ou le calage si l'helice est a pas variable."
        )

    if point.stalled_fraction > 0.0:
        level = LOW if point.stalled_fraction > 0.5 else MEDIUM
        result.confidence = level if level == LOW else result.confidence
        result.warnings.append(
            f"{point.stalled_fraction:.0%} de l'envergure est au-dela du decrochage "
            f"(Cl plafonne a {config.CL_STALL:g}) : la portance y est bornee et la trainee "
            "sous-estimee, le modele de profil ne decrit plus l'ecoulement. Poussee et "
            "rendement sont des bornes hautes."
        )

    if point.tip_mach > config.TIP_MACH_WARN:
        result.warnings.append(
            f"Mach {point.tip_mach:.2f} en bout de pale, au-dela de {config.TIP_MACH_WARN:g} : "
            "la compressibilite creuse la trainee et rabote la portance, effets qu'un modele "
            "incompressible comme celui-ci ne voit pas."
        )

    # La theorie de la quantite de mouvement cesse de valoir quand la pale est
    # tres chargee : au-dela de a = 0.4, le sillage se retourne (etat de sillage
    # turbulent) et le bilan simple surestime la poussee. Le dire vaut mieux que
    # de rendre un chiffre que la theorie ne porte plus.
    if point.speed > 0.0:
        loaded = sum(
            1 for station in result.stations
            if station.axial_induced / point.speed > config.BEM_INDUCTION_MAX
        )
        if loaded:
            result.confidence = MEDIUM if result.confidence == HIGH else result.confidence
            result.warnings.append(
                f"{loaded} station(s) sur {len(result.stations)} depassent un facteur "
                f"d'induction axial de {config.BEM_INDUCTION_MAX:g} : la pale y est assez "
                "chargee pour que le sillage se retourne, etat que la theorie simple de la "
                "quantite de mouvement ne decrit plus. La poussee y est plutot surestimee."
            )

    unconverged = sum(1 for station in result.stations if not station.converged)
    if unconverged:
        result.confidence = LOW
        result.warnings.append(
            f"{unconverged} station(s) radiale(s) sur {len(result.stations)} n'ont pas converge : "
            "le bilan par element de pale n'y a pas de solution stable, typiquement sur une pale "
            "tres chargee ou decrochee."
        )

    if point.efficiency > point.froude_efficiency > 0.0:
        # Impossible physiquement : le disque actif ideal est un plafond.
        result.confidence = LOW
        result.warnings.append(
            f"rendement calcule ({point.efficiency:.1%}) au-dessus du plafond ideal de Froude "
            f"({point.froude_efficiency:.1%}) : c'est une erreur du modele, pas une bonne helice. "
            "Le resultat n'est pas exploitable."
        )
    elif point.froude_efficiency > 0.0:
        result.notes.append(
            f"rendement propulsif {point.efficiency:.1%} contre un plafond ideal de "
            f"{point.froude_efficiency:.1%} pour la meme poussee : "
            f"{point.ceiling_gap() * 100.0:.1f} points partent en trainee de profil, en perte de "
            "bout de pale et en giration residuelle. C'est la marge que peut reprendre un "
            "meilleur dessin de pale, pas davantage."
        )

    best = max((p for p in result.curve if p.efficiency > 0.0),
               key=lambda p: p.efficiency, default=None)
    if best is not None:
        result.notes.append(
            f"rendement maximal {best.efficiency:.1%} a J = {best.advance_ratio:.2f}, soit "
            f"{best.speed * 3.6:.0f} km/h a {result.rpm:.0f} tr/min."
        )
