"""Enchainement complet des phases 1 a 6 sur un fichier de geometrie.

Ce module ne fait aucun calcul propre : il appelle les phases dans l'ordre,
propage les niveaux de confiance et rassemble les avertissements.  La CLI et les
rapports ne connaissent que sa sortie.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

from . import config
from .confidence import HIGH, LOW, ConfidenceMap, worst
from .geometry import axis as axis_module
from .geometry import blade_angles as blade_module
from .geometry import blade_loops as loops_module
from .geometry import blade_normals as normals_module
from .geometry import occupancy as occupancy_module
from .geometry import sections as sections_module
from .geometry import topology as topology_module
from .hydraulics import cavitation as cavitation_module
from .hydraulics import energy as energy_module
from .hydraulics import losses as losses_module
from .hydraulics import meanline as meanline_module
from .hydraulics import propulsion as propulsion_module
from .hydraulics import similarity as similarity_module
from .io import loader
from .mesh import TriMesh


@dataclass
class Options:
    """Toutes les entrees utilisateur de la CLI (SPEC phase 7)."""

    unit: str | float | None = None
    r_aspiration_cm: float | None = None
    speeds: tuple[float, ...] = config.DEFAULT_RPM
    blades: int | None = None
    beta1_deg: float | None = None
    beta2_deg: float | None = None
    rotation: int | None = None  # +1 anti-horaire, -1 horaire, None : a indiquer
    propulsion_speed: float | None = None  # m/s - vitesse d'avance en helice libre ; None : pas d'analyse propulsive
    fluid: str = "eau"  # fluide de l'analyse propulsive : "eau" ou "air"
    altitude: float = config.ALTITUDE
    temperature_c: float = config.TEMPERATURE
    suction_height: float = config.HAUTEUR_ASPIRATION
    suction_losses: float = config.PERTES_ASPIRATION
    suction: str = axis_module.SUCTION_AUTO
    grid_nr: int = config.GRID_NR
    grid_nz: int = config.GRID_NZ
    n_theta: int = config.N_THETA
    repair: bool = True
    prefer_trimesh: bool = True
    symmetry_check: bool = True

    def check(self) -> None:
        """Refuse les entrees hors du domaine des modeles, avec la raison.

        Chaque borne est celle d'un modele nomme, pas un garde-fou arbitraire :
        le domaine de la correlation d'Antoine pour la temperature, la
        troposphere du modele d'atmosphere OACI pour l'altitude, la finesse sous
        laquelle la carte d'occupation ne resout plus rien pour la grille.  Le
        controle est porte par les options et non par la ligne de commande, de
        sorte que la page web et l'appel direct y passent aussi : sans lui, une
        grille nulle sortait en `ZeroDivisionError` et deux secteurs en "inf
        n'est pas serialisable en JSON" -- des traces de pile la ou il fallait
        une phrase.
        """
        for rpm in self.speeds:
            if not 0.0 < rpm <= config.RPM_MAX:
                raise ValueError(
                    f"regime hors domaine : {rpm:g} tr/min. Attendu entre 0 (exclu) et "
                    f"{config.RPM_MAX:g} tr/min."
                )
        if self.blades is not None and self.blades < config.BLADES_MIN:
            raise ValueError(
                f"nombre de pales impose invalide : {self.blades}. Le modele en demande au "
                f"moins {config.BLADES_MIN}."
            )
        for name, value in (("beta1", self.beta1_deg), ("beta2", self.beta2_deg)):
            if value is not None and not config.BETA_MIN_DEG <= value <= config.BETA_MAX_DEG:
                raise ValueError(
                    f"{name} hors domaine : {value:g} degres. Attendu entre "
                    f"{config.BETA_MIN_DEG:g} et {config.BETA_MAX_DEG:g}, angles mesures depuis "
                    "la direction tangentielle."
                )
        if min(self.grid_nr, self.grid_nz) < config.GRID_MIN:
            raise ValueError(
                f"grille trop grossiere : {self.grid_nr}x{self.grid_nz}. Au moins "
                f"{config.GRID_MIN} cellules dans chaque direction."
            )
        if self.n_theta < config.N_THETA_MIN:
            raise ValueError(
                f"trop peu de secteurs azimutaux : {self.n_theta}. Au moins "
                f"{config.N_THETA_MIN}, sans quoi le nombre de pales ne se lit plus."
            )
        if not config.TEMPERATURE_MIN <= self.temperature_c <= config.TEMPERATURE_MAX:
            raise ValueError(
                f"temperature hors domaine : {self.temperature_c:g} degres C. La correlation "
                f"d'Antoine utilisee pour la pression de vapeur vaut de "
                f"{config.TEMPERATURE_MIN:g} a {config.TEMPERATURE_MAX:g} degres C."
            )
        if not 0.0 <= self.altitude <= config.ALTITUDE_MAX:
            raise ValueError(
                f"altitude hors domaine : {self.altitude:g} m. Le modele d'atmosphere OACI "
                f"utilise pour la pression barometrique vaut de 0 a "
                f"{config.ALTITUDE_MAX:g} m."
            )
        if self.propulsion_speed is not None and not 0.0 <= self.propulsion_speed <= config.PROPULSION_SPEED_MAX:
            raise ValueError(
                f"vitesse d'avance hors domaine : {self.propulsion_speed:g} m/s. Attendu entre 0 "
                f"et {config.PROPULSION_SPEED_MAX:g} m/s."
            )
        if self.fluid not in propulsion_module.FLUIDS:
            raise ValueError(
                f"fluide inconnu : '{self.fluid}'. Attendu : "
                f"{', '.join(sorted(propulsion_module.FLUIDS))}."
            )
        if self.r_aspiration_cm is not None and self.r_aspiration_cm <= 0.0:
            raise ValueError(
                f"rayon d'aspiration impose invalide : {self.r_aspiration_cm:g} cm. Il doit "
                "etre positif."
            )

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "unite": self.unit if self.unit is not None else "cm",
            "r_aspiration_cm": self.r_aspiration_cm,
            "regimes_tr_min": list(self.speeds),
            "pales_imposees": self.blades,
            "beta1_impose_deg": self.beta1_deg,
            "beta2_impose_deg": self.beta2_deg,
            "sens_de_rotation_impose": self.rotation,
            "vitesse_d_avance_m_s": self.propulsion_speed,
            "fluide": self.fluid,
            "cote_aspiration": self.suction,
            "altitude_m": self.altitude,
            "temperature_C": self.temperature_c,
            "hauteur_aspiration_m": self.suction_height,
            "pertes_aspiration_m": self.suction_losses,
            "grille": {"nr": self.grid_nr, "nz": self.grid_nz, "n_theta": self.n_theta},
            "reparation": self.repair,
            "controle_de_symetrie": self.symmetry_check,
        }


@dataclass
class AnalysisResult:
    """Resultat complet d'une analyse, tout en SI."""

    options: Options = field(default_factory=Options)
    import_report: loader.ImportReport | None = None
    axis: axis_module.AxisResult | None = None
    suction: axis_module.SuctionResult | None = None
    occupancy: occupancy_module.OccupancyMap | None = None
    topology: topology_module.Topology | None = None
    blades: blade_module.BladeGeometry | None = None
    blade_loops: loops_module.LoopResult | None = None
    blade_normals: normals_module.NormalAngles | None = None
    propulsion: propulsion_module.PropulsionResult | None = None
    meanline_input: meanline_module.MeanlineInput | None = None
    curves: list[meanline_module.PerformanceCurve] = field(default_factory=list)
    head_sensitivity: float = 0.0  # ecart relatif de hauteur pour +/- 1 deg sur beta2
    channel_losses: losses_module.ChannelLosses | None = None
    energy_budget: energy_module.EnergyBudget | None = None
    installation: cavitation_module.Installation | None = None
    speed_limit: cavitation_module.SpeedLimit | None = None
    similarity: similarity_module.SimilarityCheck | None = None
    discharge: dict = field(default_factory=dict)
    mesh: TriMesh | None = None
    warnings: list[str] = field(default_factory=list)
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    elapsed_s: float = 0.0

    def overall_confidence(self) -> str:
        """Niveau de confiance global de l'analyse."""
        return self.confidence.overall()

    def to_dict(self) -> dict:
        """Vue serialisable en JSON, toutes grandeurs en SI (SPEC 7.1)."""
        return {
            "version": __import__("impeller_analyzer").__version__,
            "entrees": self.options.to_dict(),
            "import": self.import_report.to_dict() if self.import_report else None,
            "axe": self.axis.to_dict() if self.axis else None,
            "cote_aspiration": self.suction.to_dict() if self.suction else None,
            "carte_d_occupation": self.occupancy.summary() if self.occupancy else None,
            "topologie": self.topology.to_dict() if self.topology else None,
            # `None` quand un degre suffit a supprimer le point : JSON ne prend pas l'infini.
            "pertes_de_canal": self.channel_losses.to_dict() if self.channel_losses else None,
            "bilan_d_energie": self.energy_budget.to_dict() if self.energy_budget else None,
            "sensibilite_hauteur_a_beta2": (
                self.head_sensitivity if math.isfinite(self.head_sensitivity) else None
            ),
            "forme_des_aubes": self.blade_loops.to_dict() if self.blade_loops else None,
            "angles_par_normales": self.blade_normals.to_dict() if self.blade_normals else None,
            "propulsion": self.propulsion.to_dict() if self.propulsion else None,
            "pales": self.blades.to_dict() if self.blades else None,
            "sens_de_sortie_du_liquide": self.discharge,
            "installation": self.installation.to_dict() if self.installation else None,
            "regimes": [curve.to_dict() for curve in self.curves],
            "similitude": self.similarity.to_dict() if self.similarity else None,
            "vitesse_maximale": self.speed_limit.to_dict() if self.speed_limit else None,
            "confiance": dict(self.confidence),
            "confiance_globale": self.overall_confidence(),
            "avertissements": list(self.warnings),
            "incertitude_du_modele": {
                "hauteur": config.UNCERTAINTY_H,
                "debit": config.UNCERTAINTY_Q,
                "npshr": config.UNCERTAINTY_NPSH,
                "remarque": (
                    "modele 1D ligne moyenne : a verifier par essai sur banc avant toute "
                    "decision d'achat ou de dimensionnement"
                ),
            },
            "duree_s": self.elapsed_s,
        }


def run(path: str, options: Options | None = None) -> AnalysisResult:
    """Analyse complete d'un fichier de geometrie."""
    options = options or Options()
    options.check()
    started = time.time()
    result = AnalysisResult(options=options)

    # Phase 1 : import.
    mesh, import_report = loader.load_mesh(
        path, unit=options.unit, do_repair=options.repair, prefer_trimesh=options.prefer_trimesh
    )
    result.import_report = import_report
    result.warnings.extend(import_report.warnings)
    result.confidence.update(import_report.confidence)

    # Phase 2 : axe, recentrage, carte d'occupation.
    aligned, axis_result = axis_module.align_to_z(mesh)
    result.mesh = aligned
    result.axis = axis_result
    result.warnings.extend(axis_result.warnings)
    result.confidence.set("axe", axis_result.confidence)
    occupancy = occupancy_module.build_occupancy(
        aligned, nr=options.grid_nr, nz=options.grid_nz, n_theta=options.n_theta
    )

    # La convention impose l'aspiration vers +Z, mais un fichier de CAO n'a
    # aucune raison de la respecter : une roue exportee a l'envers serait lue
    # depuis son cote refoulement, passerait pour axiale et recevrait un sens de
    # sortie faux. Sur une roue a composante radiale la geometrie tranche seule.
    suction = axis_module.detect_suction_side(occupancy, options.suction)
    result.suction = suction
    if suction.flipped:
        aligned = axis_module.flip_axis(aligned)
        result.mesh = aligned
        occupancy = occupancy_module.build_occupancy(
            aligned, nr=options.grid_nr, nz=options.grid_nz, n_theta=options.n_theta
        )
    result.warnings.extend(suction.warnings)
    result.confidence.set("cote_aspiration", suction.confidence)
    result.occupancy = occupancy

    # Phase 3 : topologie.
    topology = topology_module.analyse(
        occupancy,
        aligned if options.symmetry_check else None,
        forced_blades=options.blades,
        r_aspiration_cm=options.r_aspiration_cm,
    )
    result.topology = topology
    result.warnings.extend(topology.warnings)
    result.warnings.extend(topology.blades.warnings)
    result.confidence.update(topology.confidence)

    # Forme des aubes : une aube qui se referme sur elle-meme invalide tout ce
    # que la phase 4 en tirerait, la coupe la traversant deux fois.
    loops = loops_module.detect_looped_blades(
        occupancy, topology.blades.n_blades, watertight=import_report.watertight
    )
    result.blade_loops = loops
    result.warnings.extend(loops.warnings)

    # Phase 4 : coupes, angles de pale, sens de rotation.
    sections = sections_module.extract_sections(aligned, occupancy, topology)
    geometry = blade_module.analyse(
        sections,
        topology,
        forced_beta1_deg=options.beta1_deg,
        forced_beta2_deg=options.beta2_deg,
        forced_rotation=options.rotation,
    )
    result.blades = geometry
    result.confidence.update(geometry.confidence)
    camber_warnings = list(geometry.warnings)

    # Repli sur les normales. Deux cas l'appellent : une aube en boucle, ou la
    # cambrure n'a pas de reponse stable ; et une roue trop courte radialement,
    # ou les surfaces de courant effleurent l'aube au lieu de la traverser et ne
    # rendent que des echardes -- ecartees, il ne reste aucune coupe. Dans les
    # deux cas la lecture par les normales, elle, tient.
    def _ecrete(valeur: float) -> bool:
        """Un angle pose sur une borne du domaine est un ecretage, pas une mesure."""
        return (
            abs(valeur - config.BETA_MIN_DEG) < config.BETA_CLAMP_TOL
            or abs(valeur - config.BETA_MAX_DEG) < config.BETA_CLAMP_TOL
        )

    coherence = geometry.wrap_consistency
    incoherente = coherence > 0.0 and not (
        config.WRAP_CONSISTENCY_MIN <= coherence <= config.WRAP_CONSISTENCY_MAX
    )
    camber_failed = (
        not geometry.sections
        or geometry.beta2_deg <= 0.0
        or _ecrete(geometry.beta1_deg)
        or _ecrete(geometry.beta2_deg)
        or incoherente
    )
    if camber_failed and not loops.looped:
        result.warnings.append(
            "la lecture par la cambrure a echoue : "
            + (f"enroulement mesure {coherence:.2f} fois celui qu'impliquent les angles lus, "
               "les deux se contredisent" if incoherente
               else "aucune coupe exploitable, ou un angle pose sur une borne du domaine")
            + ". Sur une roue courte radialement les surfaces de courant effleurent l'aube au "
            "lieu de la traverser et n'en rendent que des fragments. Les angles sont repris sur "
            "les normales de la surface d'aube."
        )
    if loops.looped or camber_failed:
        normals = normals_module.analyse(aligned, occupancy, topology)
        result.blade_normals = normals
        if normals.beta1_deg > 0.0 and normals.beta2_deg > 0.0:
            # Les normales comblent ce que l'utilisateur n'a pas donne, et rien
            # de plus : imposer beta2 seul ne doit pas faire retomber beta1 sur
            # la cambrure, qui ne s'applique pas a cette forme -- elle y donnait
            # 87 degres la ou les normales en lisent 10.
            if options.beta1_deg is None:
                geometry.beta1_deg = normals.beta1_deg
            if options.beta2_deg is None:
                geometry.beta2_deg = normals.beta2_deg
            # La lecture par les normales alimente la **suggestion**, pas le
            # resultat : le sens retenu reste celui qu'indique l'utilisateur.
            geometry.observed_rotation_sign = blade_module.rotation_sense(
                topology.machine_type, normals.slope_sign
            )
            geometry.observed_rotation_label = blade_module.rotation_label(
                geometry.observed_rotation_sign
            )
            geometry.notes.extend(normals.notes)
            geometry.warnings.append(
                f"beta1 = {geometry.beta1_deg:.1f} deg et beta2 = "
                f"{geometry.beta2_deg:.1f} deg. Ceux que vous n'avez pas imposes sont lus sur les "
                "normales de la surface d'aube, la cambrure n'ayant pas de reponse stable sur "
                "cette forme ; methode basse de 2 a 5 degres sur des roues d'angles connus, biais "
                "non corrige. Imposez --beta1 et --beta2 si vous les connaissez."
            )
            # Les normales fournissent les angles, plus le sens : celui-ci vient
            # de l'utilisateur, et sa confiance ne se lit pas sur une mesure.
            result.confidence.set("angles_de_pale", normals.confidence)
        else:
            result.warnings.append(
                "aubes en boucle et lecture par les normales infructueuse : les angles de pale "
                "restent ceux de la cambrure, qui ne s'applique pas a cette forme. Imposez "
                "--beta1 et --beta2."
            )
            for quantity in loops_module.INVALIDATED_BY_LOOP:
                result.confidence.set(quantity, LOW)

    if result.blade_normals is not None and result.blade_normals.beta2_deg > 0.0 and camber_warnings:
        # Les reserves de la cambrure portent sur une lecture qui n'a pas ete
        # retenue : les garder sans le dire ferait croire a un doute sur les
        # angles publies.  Un paragraphe d'introduction n'y suffit pas -- elles
        # citent des angles chiffres, et le lecteur qui parcourt les puces lit
        # "beta2 = 86 deg" a cote d'un tableau qui publie 3.6 : chaque reserve
        # porte donc sa provenance en tete.
        result.warnings.append(
            "les reserves qui suivent portent sur la lecture par la cambrure, mise de cote au "
            "profit des normales ; elles n'entament pas les angles publies, elles expliquent "
            "pourquoi la cambrure a ete ecartee. Les angles qu'elles citent sont ceux de cette "
            "lecture ecartee, pas ceux du tableau 1."
        )
        camber_warnings = [
            f"[lecture par la cambrure, ecartee] {warning}" for warning in camber_warnings
        ]
    result.warnings.extend(camber_warnings)

    # La confrontation entre le sens impose et celui que suggere la geometrie se
    # fait ici, et pas dans la phase 4 : sur une aube en boucle, la suggestion
    # vient d'etre reprise sur les normales, et la comparer plus tot aurait
    # nomme le sens d'une lecture ecartee.
    contradiction = blade_module.forced_rotation_warning(geometry)
    if contradiction:
        result.warnings.append(contradiction)

    # Phases 5 et 6 : hydraulique et cavitation.
    data = meanline_module.MeanlineInput.from_geometry(topology, geometry)
    result.meanline_input = data
    hydraulic_confidence = meanline_module.confidence_of(topology, geometry)
    result.confidence.update(hydraulic_confidence)

    site = cavitation_module.installation(
        altitude=options.altitude,
        temperature_c=options.temperature_c,
        suction_height=options.suction_height,
        suction_losses=options.suction_losses,
    )
    result.installation = site

    curves = []
    for rpm in options.speeds:
        curve = meanline_module.build_curve(data, rpm, hydraulic_confidence)
        cavitation_module.apply_to_curve(curve)
        curves.append(curve)
        for message in curve.warnings:
            if message not in result.warnings:
                result.warnings.append(message)
    result.curves = curves
    result.confidence.set("npshr", worst(hydraulic_confidence.overall(), "medium"))

    if curves:
        reference = curves[0]
        result.speed_limit = cavitation_module.maximum_speed(reference, site, hydraulic_confidence)
        result.warnings.extend(result.speed_limit.warnings)
        if len(curves) > 1:
            result.similarity = similarity_module.check(curves[0], curves[-1])
            result.warnings.extend(result.similarity.warnings)
        # Controle croise du type de roue par la vitesse specifique (SPEC 3.3).
        message = topology_module.cross_check_type(topology, reference.specific_speed)
        if message:
            result.warnings.insert(0, message)
        best = reference.best_efficiency_point()
        if best is not None:
            result.discharge = blade_module.discharge_direction(
                topology, geometry, cm2=best.cm2, cu2=best.cu2
            )
    if not result.discharge:
        result.discharge = blade_module.discharge_direction(topology, geometry)


    # Garde-fou general : un angle pose sur une borne du domaine n'est pas une
    # mesure, c'est un ecretage. Quelle qu'en soit la cause -- et il en reste
    # forcement que je n'ai pas rencontrees -- il ne doit jamais sortir sans le
    # dire.
    if not geometry.forced_beta:
        for nom, valeur in (("beta1", geometry.beta1_deg), ("beta2", geometry.beta2_deg)):
            if _ecrete(valeur):
                result.warnings.append(
                    f"{nom} = {valeur:.1f} deg est pose sur une borne du domaine "
                    f"({config.BETA_MIN_DEG:.0f} a {config.BETA_MAX_DEG:.0f}) : c'est un "
                    "ecretage, pas une mesure. Les grandeurs qui en derivent ne veulent rien "
                    f"dire ; imposez --{nom} si vous connaissez sa valeur."
                )
                result.confidence.set("angles_de_pale", LOW)

    result.confidence.set(
        "sens_de_rotation", HIGH if geometry.forced_rotation else LOW
    )
    if topology.rotation_ambiguity and not geometry.forced_rotation:
        result.warnings.insert(0, topology.rotation_ambiguity)

    if curves:
        # Pertes calculees sur la geometrie du canal : elles ne remplacent pas
        # celles de la SPEC, elles servent a comparer deux roues entre elles.
        best = curves[0].best_efficiency_point()
        if best is not None:
            # D'ou vient la hauteur : centrifuge, diffusion, energie cinetique.
            result.energy_budget = energy_module.budget(
                u1=curves[0].u1, u2=curves[0].u2,
                cm1=best.cm1, cm2=best.cm2, cu2=best.cu2,
            )
        if best is not None and topology.r_1 > 0.0:
            width_1 = topology.area_1 / (2.0 * math.pi * topology.r_1)
            w_2 = math.hypot(best.cm2, max(0.0, curves[0].u2 - best.cu2))
            result.channel_losses = losses_module.analyse(
                r_1=topology.r_1,
                r_2=topology.r_2,
                b_1=width_1,
                b_2=topology.b_2,
                beta1_deg=geometry.beta1_deg,
                beta2_deg=geometry.beta2_deg,
                n_blades=topology.blades.n_blades,
                w1=best.w1,
                w2=w_2,
                head_theoretical=best.head_theoretical,
                strands=2 if loops.looped else 1,
            )
            result.warnings.extend(result.channel_losses.warnings)

        result.head_sensitivity = meanline_module.head_sensitivity(data, curves[0].rpm)
        if result.head_sensitivity > config.BETA_SENSITIVITY_ALERT:
            # Le texte dit « ordres de grandeur » : la table de confiance doit
            # dire la meme chose, sans quoi le lecteur croit l'une des deux.
            # Debit et NPSHr ne suivent pas : ils ne passent pas par cu2.
            for quantity in ("hauteur", "puissance", "couple", "rendement"):
                result.confidence.set(quantity, LOW)

        if result.head_sensitivity == float("inf"):
            result.warnings.append(
                f"hauteur hypersensible a beta2 : un degre d'ecart suffit a faire disparaitre le "
                f"point de fonctionnement. beta2 vaut {data.beta2_deg:.1f} deg, et aux petits "
                "angles cu2 = u2 - cm2/tan(beta2) varie tres vite. Hauteur et puissance ne sont "
                "pas exploitables ; le debit et le NPSHr, qui n'en dependent pas de la meme "
                "facon, restent fiables."
            )
        elif result.head_sensitivity > config.BETA_SENSITIVITY_ALERT:
            result.warnings.append(
                f"hauteur hypersensible a beta2 : un degre d'ecart la deplace de "
                f"{result.head_sensitivity:.0%}, bien au-dela des 18 % annonces par le modele. "
                f"beta2 vaut {data.beta2_deg:.1f} deg, et aux petits angles cu2 = u2 - cm2/tan(beta2) "
                "varie tres vite. Hauteur et puissance sont a lire comme des ordres de grandeur ; "
                "le debit et le NPSHr, qui n'en dependent pas de la meme facon, restent fiables."
            )

    # Phase propulsive, sur demande. Elle ne remplace pas l'analyse de pompe :
    # elle repond a une autre question -- la meme piece tournant en helice libre,
    # non carenee -- et le module refuse de repondre si la roue n'est pas axiale.
    if options.propulsion_speed is not None:
        result.propulsion = propulsion_module.analyse(
            topology,
            geometry,
            rpm=max(options.speeds) if options.speeds else 0.0,
            speed=options.propulsion_speed,
            fluid=options.fluid,
        )
        result.warnings.extend(result.propulsion.warnings)
        result.confidence.set("propulsion", result.propulsion.confidence)

    result.elapsed_s = time.time() - started
    return result
