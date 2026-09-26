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
from .confidence import HIGH, LOW, MEDIUM, ConfidenceMap, worst
from . import components as components_module
from .provenance import DECLARE, DEFAUT, MESURE, ProvenanceMap
from .geometry import axis as axis_module
from .geometry import blade_angles as blade_module
from .geometry import blade_isolated as isolated_module
from .geometry import blade_loops as loops_module
from .geometry import blade_normals as normals_module
from .geometry import occupancy as occupancy_module
from .geometry import sections as sections_module
from .geometry import topology as topology_module
from .hydraulics import cavitation as cavitation_module
from .hydraulics import comparison as comparison_module
from .hydraulics import energy as energy_module
from .hydraulics import losses as losses_module
from .hydraulics import meanline as meanline_module
from .hydraulics import propulsion as propulsion_module
from .hydraulics import sensitivity as sensitivity_module
from .hydraulics import similarity as similarity_module
from .io import loader
from .mesh import TriMesh, dot


#: Modele hydraulique retenu. `auto` laisse la classification geometrique
#: decider ; les deux autres priment sur elle, parce que l'utilisateur a la
#: piece sous les yeux et la geometrie, elle, peut se tromper.
MACHINE_AUTO = "auto"
MACHINE_PUMP = "pompe_carenee"
MACHINE_PROPELLER = "helice_libre"
MACHINE_MODELS = (MACHINE_AUTO, MACHINE_PUMP, MACHINE_PROPELLER)

#: Famille de roue, declarable de la meme facon.
WHEEL_AUTO = "auto"
WHEEL_TYPES = (WHEEL_AUTO, topology_module.AXIAL, topology_module.MIXED,
               topology_module.CENTRIFUGAL)


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
    component_paths: dict = field(default_factory=dict)  # emplacements du mode composants
    blade_topology: str = components_module.BLADE_CONVENTIONAL  # topologie de pale declaree
    topology_declared: bool = False  # la topologie de pale a-t-elle ete declaree par l'utilisateur
    # L'outil n'etudie que les helices toroidales : il refuse une piece qui n'en est
    # pas une. Le moteur, lui, reste general -- le calcul d'une helice normale sert
    # de reference au comparatif --, d'ou un drapeau que la ligne de commande et la
    # page levent, et que l'appel direct laisse baisse.
    toroidal_only: bool = False
    # Roue classique analysee comme reference (--reference) : le refus est leve,
    # et le rapport le dit en tete.
    reference: bool = False
    machine: str = MACHINE_AUTO  # modele hydraulique : auto, pompe carenee, ou helice libre
    wheel_type: str = WHEEL_AUTO  # famille de roue declaree ; prime sur la classification
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
        if self.blades is not None and not config.BLADES_MIN <= self.blades <= config.DECLARED_BLADES_MAX:
            raise ValueError(
                f"nombre de pales impose invalide : {self.blades}. Attendu entre "
                f"{config.BLADES_MIN} et {config.DECLARED_BLADES_MAX}."
            )
        for name, value in (("beta1", self.beta1_deg), ("beta2", self.beta2_deg)):
            if value is not None and not config.BETA_MIN_DEG <= value <= config.BETA_MAX_DEG:
                raise ValueError(
                    f"{name} hors domaine : {value:g} degres. Attendu entre "
                    f"{config.BETA_MIN_DEG:g} et {config.BETA_MAX_DEG:g}, angles mesures depuis "
                    "la direction tangentielle."
                )
        if max(self.grid_nr, self.grid_nz) > config.GRID_MAX:
            raise ValueError(
                f"grille trop fine : {self.grid_nr}x{self.grid_nz}. Au plus {config.GRID_MAX} "
                "cellules dans chaque direction : au-dela, la memoire et la duree explosent "
                "sans rien resoudre de plus."
            )
        if self.n_theta > config.N_THETA_MAX:
            raise ValueError(
                f"trop de secteurs azimutaux : {self.n_theta}. Au plus {config.N_THETA_MAX}."
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
        if self.machine not in MACHINE_MODELS:
            raise ValueError(
                f"modele hydraulique inconnu : '{self.machine}'. Attendu : "
                f"{', '.join(MACHINE_MODELS)}."
            )
        if self.wheel_type not in WHEEL_TYPES:
            raise ValueError(
                f"type de roue inconnu : '{self.wheel_type}'. Attendu : "
                f"{', '.join(WHEEL_TYPES)}."
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
        if self.r_aspiration_cm is not None and not (
            math.isfinite(self.r_aspiration_cm) and self.r_aspiration_cm > 0.0
        ):
            raise ValueError(
                f"rayon d'aspiration impose invalide : {self.r_aspiration_cm:g} cm. Il doit "
                "etre un nombre positif."
            )
        if not config.HAUTEUR_ASPIRATION_MIN <= self.suction_height <= config.HAUTEUR_ASPIRATION_MAX:
            raise ValueError(
                f"hauteur d'aspiration hors domaine : {self.suction_height:g} m. Attendue entre "
                f"{config.HAUTEUR_ASPIRATION_MIN:g} m (aspiration, pompe au-dessus du plan d'eau) "
                f"et {config.HAUTEUR_ASPIRATION_MAX:g} m (en charge) : au-dela, c'est une faute "
                "de saisie ou d'unite."
            )
        if not 0.0 <= self.suction_losses <= config.PERTES_ASPIRATION_MAX:
            raise ValueError(
                f"pertes de charge d'aspiration hors domaine : {self.suction_losses:g} m. "
                f"Attendues entre 0 et {config.PERTES_ASPIRATION_MAX:g} m : une perte negative "
                "ajouterait de l'energie au liquide."
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
            "modele_hydraulique": self.machine,
            "type_de_roue_impose": self.wheel_type,
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
    assembly: components_module.ComponentAssembly | None = None
    isolated_angles: isolated_module.IsolatedBladeAngles | None = None  # mode composants
    propulsion: propulsion_module.PropulsionResult | None = None
    sensitivity: sensitivity_module.SensitivityReport | None = None
    comparison: comparison_module.Comparison | None = None  # si l'helice etait normale
    meanline_input: meanline_module.MeanlineInput | None = None
    curves: list[meanline_module.PerformanceCurve] = field(default_factory=list)
    head_sensitivity: float = 0.0  # ecart relatif de hauteur pour +/- 1 deg sur beta2
    channel_losses: losses_module.ChannelLosses | None = None
    energy_budget: energy_module.EnergyBudget | None = None
    installation: cavitation_module.Installation | None = None
    speed_limit: cavitation_module.SpeedLimit | None = None
    similarity: similarity_module.SimilarityCheck | None = None
    discharge: dict = field(default_factory=dict)
    rejection: str | None = None  # motif de refus : la piece n'est pas une roue que le modele decrit
    mesh: TriMesh | None = None
    warnings: list[str] = field(default_factory=list)
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    provenance: ProvenanceMap = field(default_factory=ProvenanceMap)
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
            "composants": self.assembly.to_dict() if self.assembly else None,
            "angles_sur_pale_isolee": (
                self.isolated_angles.to_dict() if self.isolated_angles else None
            ),
            "mode_d_import": "composants" if self.assembly else "vrac",
            "propulsion": self.propulsion.to_dict() if self.propulsion else None,
            "sensibilite": self.sensitivity.to_dict() if self.sensitivity else None,
            "si_l_helice_etait_normale": self.comparison.to_dict() if self.comparison else None,
            "pales": self.blades.to_dict() if self.blades else None,
            "sens_de_sortie_du_liquide": self.discharge,
            "installation": self.installation.to_dict() if self.installation else None,
            "regimes": [curve.to_dict() for curve in self.curves],
            "similitude": self.similarity.to_dict() if self.similarity else None,
            "vitesse_maximale": self.speed_limit.to_dict() if self.speed_limit else None,
            "confiance": dict(self.confidence),
            "provenance": dict(self.provenance),
            "confiance_globale": self.overall_confidence(),
            "analyse_interrompue": self.rejection,
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


def _periodicity(aligned: TriMesh) -> tuple[float, int]:
    """Periodicite d'une orientation : `(ecart relatif apres 2*pi/N, N)`.

    Carte grossiere et nombre de pales lu sur elle : le but est de departager
    des axes candidats, pas de mesurer la roue. `(inf, 0)` si aucune pale ne se
    lit.
    """
    occupancy = occupancy_module.build_occupancy(
        aligned, nr=config.AXIS_PROBE_GRID, nz=config.AXIS_PROBE_GRID, n_theta=config.AXIS_PROBE_THETA
    )
    count = topology_module.count_blades(
        occupancy, aligned, samples=config.AXIS_PROBE_SAMPLES, cap=2.0 * config.SYM_REJECT,
        verdict_only=True,
    )
    if count.periodicity is None:
        return math.inf, 0
    return count.periodicity, count.periodic_order


def _parallel(a, b) -> bool:
    return abs(dot(a, b)) >= math.cos(math.radians(config.AXIS_WARN_DEG))


def _choose_axis(mesh: TriMesh, options: Options) -> tuple[TriMesh, axis_module.AxisResult]:
    """Axe de rotation : l'inertie propose, la periodicite tranche quand elle doute.

    Le tenseur d'inertie designe l'axe des deux valeurs propres egales. C'est
    le bon pour une roue de trois pales ou plus ; pour deux pales, une pale
    elancee et un moyeu court donnent l'envergure, en confiance haute, et le
    maillage etait bascule de 90 degres sans que rien ne le signale. Des que
    l'inertie n'est pas nette ou contredit la convention Z, chaque axe principal
    est donc mis a l'epreuve de ce qui definit une roue : se superposer a
    soi-meme apres une rotation de 2*pi/N.

    Plusieurs axes peuvent reussir l'epreuve. Une roue a aubes helicoidales est
    aussi symetrique d'ordre 2 autour d'axes perpendiculaires au sien ; l'axe de
    la roue est alors celui d'ordre le plus eleve. A deux pales, les trois axes
    principaux sont d'ordre 2 et rien dans la geometrie ne les departage : la
    convention Z est gardee, mais l'axe est publie comme non verifie.
    """
    aligned, result = axis_module.align_to_z(mesh)
    nette = result.confidence == HIGH and result.angle_to_z_deg <= config.AXIS_WARN_DEG
    if nette or not options.symmetry_check:
        return aligned, result

    inertia_axis = result.axis
    candidates = []
    for axis in axis_module.candidate_axes(result):
        oriented = axis_module.orient(mesh, result.centre, axis)
        score, order = _periodicity(oriented)
        candidates.append((score, order, axis, oriented))
    passing = [c for c in candidates if c[0] <= config.SYM_TOL and c[1] >= config.BLADES_MIN]
    if not passing:
        # Aucun axe ne rend la piece periodique : on garde la lecture de
        # l'inertie, et le controle de recevabilite dira la suite.
        return aligned, result

    top = max(c[1] for c in passing)
    tied = [c for c in passing if c[1] == top]
    chosen = (
        next((c for c in tied if _parallel(c[2], axis_module.REFERENCE_AXIS)), None)
        or next((c for c in tied if _parallel(c[2], inertia_axis)), None)
        or tied[0]
    )
    score, order, axis, oriented = chosen
    angle = math.degrees(math.acos(max(-1.0, min(1.0, abs(axis[2])))))
    result.realigned = angle > config.AXIS_WARN_DEG
    result.axis = axis if result.realigned else axis_module.REFERENCE_AXIS
    result.angle_to_z_deg = angle if result.realigned else 0.0
    # Les messages de la detection par l'inertie ne valent plus : l'axe est
    # maintenant fixe par la periodicite, et une inertie anisotrope est la
    # signature normale d'une roue a deux pales.
    result.warnings = [
        w for w in result.warnings
        if "le maillage est realigne sur Z" not in w
        and not w.startswith("symetrie de revolution imparfaite")
    ]
    where = (
        f"a {angle:.1f} deg de Z : le maillage est realigne" if result.realigned
        else "soit Z : la convention est confirmee"
    )
    if len(tied) > 1:
        result.confidence = MEDIUM
        result.warnings.append(
            f"plusieurs axes superposent la piece a elle-meme avec le meme ordre ({order}) : c'est "
            "le cas d'une roue a deux pales helicoidales, symetrique autour de ses trois axes "
            f"principaux. Rien dans la geometrie ne designe l'axe de rotation ; celui retenu est {where}"
            ". Il n'est pas verifie : exportez la roue avec son axe de rotation sur Z."
        )
    else:
        result.confidence = HIGH
        if not _parallel(inertia_axis, axis):
            result.warnings.append(
                "l'inertie designait un autre axe que l'axe de rotation. L'axe retenu est le seul "
                f"autour duquel la piece se superpose a elle-meme apres 2*pi/{order} (ecart "
                f"{score:.3f} du rayon exterieur), {where}."
            )
    return oriented, result


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
    aligned, axis_result = _choose_axis(mesh, options)
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
    # Le type de roue declare prime sur la classification geometrique. Celle-ci
    # repose sur le rapport r2/r1s, qui suppose un canal meridien conventionnel
    # entre moyeu et carter : sur une aube en boucle il n'y en a pas, et le
    # critere repond avec aplomb une valeur fausse. L'utilisateur, lui, a la
    # piece sous les yeux.
    if options.wheel_type != WHEEL_AUTO:
        if options.wheel_type != topology.machine_type:
            topology.notes.append(
                f"type de roue impose : {options.wheel_type}. La geometrie, elle, lisait "
                f"{topology.machine_type} (rapport r2/r1s = {topology.ratio_r2_r1s:.2f}). "
                "C'est le type impose qui est retenu."
            )
        topology.machine_type = options.wheel_type
        topology.confidence.set("type_de_roue", HIGH)
        result.confidence.set("type_de_roue", HIGH)

    result.topology = topology
    result.warnings.extend(topology.warnings)
    result.warnings.extend(topology.blades.warnings)
    result.confidence.update(topology.confidence)

    # Recevabilite : tout ce qui suit suppose une roue de N pales identiques
    # autour de l'axe retenu. Si la piece ne se superpose pas a elle-meme apres
    # une rotation de 2*pi/N, les angles, la hauteur et le debit seraient
    # calcules sur une forme que le modele ne decrit pas : l'analyse s'arrete la,
    # avec la carte d'occupation et la topologie pour comprendre pourquoi.
    motif = (
        topology_module.admissibility(topology.blades, occupancy)
        if options.symmetry_check else None
    )
    if motif is not None:
        result.rejection = motif
        result.warnings.insert(0, f"analyse interrompue : {motif}")
        result.confidence.cap_all(LOW)
        _fill_provenance(result, options)
        result.elapsed_s = time.time() - started
        return result

    # Forme des aubes : une aube qui se referme sur elle-meme invalide tout ce
    # que la phase 4 en tirerait, la coupe la traversant deux fois.
    loops = loops_module.detect_looped_blades(
        occupancy, topology.blades.n_blades, watertight=import_report.watertight,
        mesh=aligned,
    )
    result.blade_loops = loops
    result.warnings.extend(loops.warnings)

    if options.toroidal_only:
        refus, reserve = _toroidal_gate_bulk(loops, topology.blades.n_blades, options)
        if refus:
            result.rejection = refus
            result.warnings.insert(0, f"analyse interrompue : {refus}")
            result.confidence.cap_all(LOW)
            _fill_provenance(result, options)
            result.elapsed_s = time.time() - started
            return result
        if reserve:
            result.warnings.insert(0, reserve)
            result.confidence.set("topologie_de_pale", LOW)

    # Les criteres de classification -- rapport r2/r1s, solidite, largeur de
    # sortie -- supposent tous un canal meridien conventionnel, borde par le
    # moyeu et le carter, que le fluide traverse une fois. Une aube en boucle
    # n'en a pas : la meme coupe rencontre deux fois la pale, et les rayons qui
    # nourrissent le rapport ne designent plus ce que le critere croit. La
    # classification garde sa valeur -- il faut bien en publier une -- mais
    # cesse de pouvoir etre affirmee.
    if loops.looped and options.wheel_type == WHEEL_AUTO:
        # Sur le resultat autant que sur la topologie : la table de confiance
        # publiee a deja ete recopiee depuis celle-ci, plus haut.
        topology.confidence.set("type_de_roue", LOW)
        result.confidence.set("type_de_roue", LOW)
        result.warnings.append(
            f"aubes en boucle : la classification geometrique repond « {topology.machine_type} » "
            f"sur un rapport r2/r1s de {topology.ratio_r2_r1s:.2f}, mais ce critere suppose un "
            "canal meridien conventionnel entre moyeu et carter, qu'une aube en boucle n'a pas. "
            "La confiance sur le type de roue est donc abaissee. Si vous connaissez la piece, "
            "declarez-la par --type-de-roue ; et --machine choisit directement le modele "
            "hydraulique, sans passer par cette classification."
        )

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
        elif loops.looped:
            result.warnings.append(
                "aubes en boucle et lecture par les normales infructueuse : les angles de pale "
                "restent ceux de la cambrure, qui ne s'applique pas a cette forme. Imposez "
                "--beta1 et --beta2."
            )
        else:
            result.warnings.append(
                "la lecture par les normales echoue elle aussi : aucune des deux methodes ne "
                "trouve d'angle de pale sur cette piece. Imposez --beta1 et --beta2."
            )
            for quantity in loops_module.INVALIDATED_BY_LOOP:
                result.confidence.set(quantity, LOW)

    if geometry.beta2_deg <= 0.0 and options.beta2_deg is None:
        # Aucune lecture n'a trouve d'aube : ce que la carte compte, ce sont des
        # lobes de matiere -- les coins d'une plaque en font aussi. Le nombre et
        # le type qui en decoulent ne peuvent pas rester affirmes.
        doubtful = []
        if options.blades is None:
            doubtful.append("nombre_de_pales")
        if options.wheel_type == WHEEL_AUTO:
            doubtful.append("type_de_roue")
        for quantity in doubtful:
            for table in (topology.confidence, result.confidence):
                table.set(quantity, worst(table.get(quantity, HIGH), MEDIUM))
        result.warnings.append(
            f"aucune surface d'aube n'a pu etre lue : les {topology.blades.n_blades} « pales » "
            "comptees sont des lobes de matiere sur la carte d'occupation, pas des aubes "
            "confirmees. Le nombre de pales et le type de roue sont donc publies en confiance "
            "moyenne au plus. Si la piece est bien une roue, imposez --beta1, --beta2 et --blades."
        )

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
    _curves(result, options, topology, geometry)

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

    toroidal = loops.looped or (
        options.topology_declared and options.blade_topology == components_module.BLADE_TOROIDAL
    )
    _curve_analysis(result, topology, geometry, strands=2 if toroidal else 1)
    if toroidal and result.curves and options.machine != MACHINE_PROPELLER:
        result.comparison = comparison_module.pump(
            result.meanline_input, result.curves[0], topology.b_2, strands=2,
            closed=topology.closed_impeller,
            confidence=result.confidence.get_level("hauteur"),
        )

    # Phase propulsive, sur demande. Elle ne remplace pas l'analyse de pompe :
    # elle repond a une autre question -- la meme piece tournant en helice libre,
    # non carenee -- et le module refuse de repondre si la roue n'est pas axiale.
    # Le mode helice libre se declare, il ne se deduit pas. Il etait conditionne
    # a la classification geometrique, ce qui le fermait des que celle-ci se
    # trompait -- exactement le cas d'une helice a aubes en boucle, lue
    # « centrifuge » parce que son rapport de rayons y ressemble.
    declared = options.machine == MACHINE_PROPELLER
    refused = options.machine == MACHINE_PUMP
    wanted = options.propulsion_speed is not None or declared
    if wanted and not refused:
        result.propulsion = propulsion_module.analyse(
            topology,
            geometry,
            rpm=max(options.speeds) if options.speeds else 0.0,
            # Sans vitesse d'avance, une helice declaree libre est analysee a
            # l'arret : la poussee statique est une grandeur utile en soi.
            speed=options.propulsion_speed if options.propulsion_speed is not None else 0.0,
            fluid=options.fluid,
            forced=declared,
        )
        result.warnings.extend(result.propulsion.warnings)
        result.confidence.set("propulsion", result.propulsion.confidence)
        if toroidal and declared:
            result.comparison = comparison_module.propeller(
                result.propulsion, topology, geometry,
                confidence=result.propulsion.confidence,
            )
    elif refused and options.propulsion_speed is not None:
        result.warnings.append(
            "une vitesse d'avance est demandee, mais la machine est declaree pompe carenee : "
            "l'analyse en helice libre n'est pas faite. Otez --machine pompe_carenee, ou "
            "declarez --machine helice_libre."
        )

    _fill_provenance(result, options)
    _reference_notice(result, options)
    result.elapsed_s = time.time() - started
    return result


def _reference_notice(result: AnalysisResult, options: Options) -> None:
    """En tete du rapport, quand la roue est analysee comme reference."""
    if options.reference:
        result.warnings.insert(0, (
            "roue analysee comme reference (--reference) : l'outil n'etudie que les helices "
            "toroidales, et celle-ci sert de point de comparaison. Aucun comparatif n'est produit."
        ))


def _curves(
    result: AnalysisResult,
    options: Options,
    topology: topology_module.Topology,
    geometry: blade_module.BladeGeometry,
) -> list:
    """Phases 5 et 6 : courbes caracteristiques, cavitation, vitesse limite."""
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

    return curves


def _curve_analysis(
    result: AnalysisResult,
    topology: topology_module.Topology,
    geometry: blade_module.BladeGeometry,
    strands: int,
) -> None:
    """D'ou vient la hauteur, pertes du canal, sensibilite a la geometrie."""
    curves = result.curves
    data = result.meanline_input
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
                strands=strands,
            )
            result.warnings.extend(result.channel_losses.warnings)

        # Sans point nominal au regime de reference, il n'y a rien dont mesurer la
        # sensibilite : la hauteur n'est pas « hypersensible », elle n'existe pas,
        # et annoncer que le debit « reste fiable » serait parler d'un chiffre
        # qui n'a pas ete calcule.
        has_point = curves[0].nominal_point() is not None
        result.head_sensitivity = (
            meanline_module.head_sensitivity(data, curves[0].rpm) if has_point else 0.0
        )

        # Sensibilite complete : hauteur, debit et NPSHr, contre beta1, beta2 et
        # le diametre. Ce que le modele annonce comme incertitude suppose la
        # geometrie juste ; ce tableau-la mesure ce qu'une erreur de lecture
        # ferait vraiment, et declasse ce qui ne tient pas.
        result.sensitivity = sensitivity_module.analyse(data, curves[0].rpm)
        sensitivity_module.apply(result.sensitivity, result.confidence)
        result.warnings.extend(result.sensitivity.warnings)
        for row in result.sensitivity.rows:
            if not row.alarming():
                continue
            # Pas `worst` : ce nom est celui de la fonction de propagation de
            # confiance importee plus haut, et Python le rendrait local a toute
            # la fonction.
            pire = row.worst_per_degree()
            combien = (
                "suffit a faire disparaitre le point de fonctionnement"
                if not math.isfinite(pire) else f"la deplace de {pire:.0%}"
            )
            result.warnings.append(
                f"{row.quantity} tres sensible a la geometrie : un degre d'ecart sur les angles "
                f"de pale {combien}, au-dela du seuil de "
                f"{config.SENSITIVITY_LOW_PER_DEG:.0%} par degre. La confiance sur cette "
                "grandeur est abaissee a faible : lisez-la comme un ordre de grandeur."
            )
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


def run_components(options: Options) -> AnalysisResult:
    """Analyse en mode composants declares (SPEC v2).

    Le mode d'import en vrac reste disponible et inchange : celui-ci s'ajoute a
    cote.  La difference tient en une phrase -- l'outil cesse de deviner l'axe,
    le nombre de pales, le type de roue et le cote d'aspiration, et se met a les
    verifier.  Les grandeurs qui en decoulent changent donc de provenance, pas
    seulement de valeur : la section d'entree est **mesuree** sur un solide au
    lieu d'etre deduite de rayons fragiles, et le sens debitant est un vecteur
    lu sur l'assemblage au lieu d'une convention.
    """
    started = time.time()
    result = AnalysisResult(options=options)
    blade_files = options.component_paths.get(components_module.SLOT_BLADE)
    fournies = len(blade_files) if isinstance(blade_files, (list, tuple)) else 0
    declarations = components_module.Declarations(
        mode=options.machine,
        blade_topology=options.blade_topology,
        # Les pales fournies une a une donnent leur nombre.
        n_blades=options.blades or (fournies if fournies > 1 else 0),
    )
    assembly = components_module.assemble(
        options.component_paths, declarations, unit=options.unit
    )
    result.assembly = assembly
    result.warnings.extend(assembly.warnings)
    result.confidence.update(assembly.confidence)

    blocked = assembly.blocked()
    refus = _toroidal_gate_components(assembly) if options.toroidal_only else None
    if blocked is not None or refus:
        motif = blocked.detail if blocked is not None else refus
        result.rejection = refus
        result.warnings.insert(0, f"analyse interrompue : {motif}")
        result.confidence.cap_all(LOW)
        result.elapsed_s = time.time() - started
        return result

    # La geometrie de la roue, batie sur ce qui a ete mesure et declare plutot
    # que sur ce qui a ete devine.
    result.topology = _topology_from_components(assembly)
    result.confidence.update(result.topology.confidence)
    result.warnings.extend(result.topology.warnings)

    # §5.1 : le sens debitant etant mesure, le sens de rotation se deduit.
    result.blades = blade_module.BladeGeometry()
    sign, raison, niveau = components_module.rotation_from_flow(assembly)
    if options.rotation is not None:
        result.blades.rotation_sign = 1 if options.rotation > 0 else -1
        result.blades.forced_rotation = True
        result.confidence.set("sens_de_rotation", HIGH)
    elif sign:
        result.blades.rotation_sign = sign
        result.confidence.set("sens_de_rotation", niveau)
        result.provenance.set("sens_de_rotation", MESURE)
    else:
        result.confidence.set("sens_de_rotation", LOW)
        result.warnings.append(raison)
    result.blades.rotation_label = blade_module.rotation_label(result.blades.rotation_sign)
    result.blades.observed_rotation_sign = sign
    result.blades.observed_rotation_label = blade_module.rotation_label(sign)
    if sign:
        result.blades.notes.append(raison)

    _read_isolated_angles(result, assembly)
    _components_hydraulics(result, options)

    _fill_provenance(result, options)
    for quantity in ("nombre_de_pales", "type_de_roue", "cote_aspiration", "axe"):
        result.provenance.set(quantity, DECLARE)
    if options.rotation is not None:
        result.provenance.set("sens_de_rotation", DECLARE)
    elif sign:
        result.provenance.set("sens_de_rotation", MESURE)
    for quantity in ("sections", "rayon_de_moyeu"):
        result.provenance.set(quantity, MESURE)
    _reference_notice(result, options)
    result.elapsed_s = time.time() - started
    return result


#: Ce que l'outil etudie, dit une fois pour toutes les raisons de refus.
TOROIDAL_SCOPE = (
    "L'outil n'etudie que les helices toroidales ; le calcul d'une helice normale n'est "
    "plus publie seul, il sert de reference au comparatif d'une helice toroidale. Pour "
    "analyser une roue classique comme point de comparaison : --reference."
)


def _toroidal_gate_bulk(
    loops: loops_module.LoopResult, n_blades: int, options: Options
) -> tuple[str | None, str | None]:
    """Import global : `(refus, reserve)`. La piece est-elle toroidale ?

    La lecture tranche ; une declaration ne la contredit que pour la completer.
    Une piece que la lecture dit conventionnelle, et que rien ne declare
    toroidale, est refusee. Si l'utilisateur la declare toroidale, l'analyse se
    poursuit en le disant : un detecteur peut se tromper, et la piece est sous
    ses yeux.
    """
    declared = options.topology_declared
    if declared and options.blade_topology == components_module.BLADE_CONVENTIONAL:
        return (
            "la piece est declaree conventionnelle (--topologie-pale conventionnelle). "
            + TOROIDAL_SCOPE, None,
        )
    toroidal_declared = declared and options.blade_topology == components_module.BLADE_TOROIDAL
    if loops.looped:
        return None, None
    if loops.undecided:
        if toroidal_declared:
            return None, (
                "forme des aubes non etablie par la lecture -- le maillage n'est pas etanche, et "
                "un maillage troue produit la meme signature qu'une boucle. La piece est traitee "
                "comme toroidale sur votre declaration."
            )
        return (
            "la forme des aubes n'a pas pu etre etablie : le maillage n'est pas etanche, et un "
            "maillage troue produit la meme signature qu'une boucle. " + TOROIDAL_SCOPE
            + " Reparez le maillage, ou declarez --topologie-pale toroidale si la piece en est une.",
            None,
        )
    if toroidal_declared:
        return None, (
            "la lecture ne voit aucune aube en boucle : une coupe a azimut fixe ne traverse "
            "chaque aube qu'une fois. L'analyse se poursuit parce que la piece est declaree "
            "toroidale ; verifiez-le, car tout ce qui suit en depend."
        )
    return (
        f"la piece n'est pas toroidale : aucune de ses {n_blades} aubes ne se referme sur "
        "elle-meme -- une coupe a azimut fixe ne traverse chaque aube qu'une fois. "
        + TOROIDAL_SCOPE
        + " Si la piece est toroidale et que la lecture se trompe, declarez-la par "
        "--topologie-pale toroidale : l'analyse continuera, en le signalant.",
        None,
    )


def _toroidal_gate_components(assembly: components_module.ComponentAssembly) -> str | None:
    """Mode composants : refus si la pale est declaree, ou mesuree sans ambiguite, conventionnelle.

    Le genre d'un maillage etanche d'un seul tenant ne se discute pas : de genre
    zero, la pale n'a pas d'anse, et ne peut pas etre une boucle fermee. Sur un
    maillage troue ou en plusieurs morceaux, le genre ne tranche pas, et la
    declaration est retenue avec la reserve du controle de topologie.
    """
    if assembly.declarations.blade_topology == components_module.BLADE_CONVENTIONAL:
        return "la pale est declaree conventionnelle. " + TOROIDAL_SCOPE
    blade = assembly.component(components_module.SLOT_BLADE)
    if blade is None or blade.mesh is None or not blade.watertight:
        return None
    g, _, _ = components_module.genus(blade.mesh)
    if g < 1 and components_module._connected_parts(blade.mesh) == 1:
        return (
            "la pale declaree toroidale n'est pas une boucle : son maillage, etanche et d'un "
            "seul tenant, est de genre 0 -- sans anse. " + TOROIDAL_SCOPE
        )
    return None


def _topology_from_components(
    assembly: components_module.ComponentAssembly,
) -> topology_module.Topology:
    """Batit la topologie a partir des pieces, sans rien inferer.

    Les sections d'entree et de sortie viennent des solides fluide, mesurees
    comme `volume / epaisseur` : elles cessent d'etre suspendues a la detection
    de `r_1s` et `r_1h`, dont la fragilite est la cause premiere des ecarts de
    debit.  Les rayons, eux, sont lus sur les memes solides -- une couronne
    donne directement son rayon interieur et son rayon exterieur.
    """
    topology = topology_module.Topology()
    axis = assembly.axis_origin
    inlet = assembly.components[components_module.SLOT_INLET]
    outlet = assembly.components[components_module.SLOT_OUTLET]
    blade = assembly.components[components_module.SLOT_BLADE]

    topology.r_1s = inlet.radial_extent(axis)[1]
    topology.r_2s = outlet.radial_extent(axis)[1]
    topology.r_1h = components_module.inner_radius(inlet, axis)
    topology.r_2h = components_module.inner_radius(outlet, axis)
    topology.r_1 = math.sqrt((topology.r_1s ** 2 + topology.r_1h ** 2) / 2.0)
    topology.r_2 = math.sqrt((topology.r_2s ** 2 + topology.r_2h ** 2) / 2.0)
    topology.r_aspiration = topology.r_1s
    topology.r_tip = max(inlet.radial_extent(axis)[1], outlet.radial_extent(axis)[1],
                         blade.radial_extent(axis)[1])
    topology.r_blade_tip = blade.radial_extent(axis)[1]

    # Mesurees, non deduites : c'est tout l'apport du mode composants. La part
    # d'un solide fluide qu'occupe une paroi n'est pas une section de passage
    # (controle « passage libre ») : seule la part libre est retenue.
    topology.area_1 = assembly.inlet.free_area if assembly.inlet else 0.0
    topology.area_2 = assembly.outlet.free_area if assembly.outlet else 0.0

    # Refoulement radial : la sortie est une bande cylindrique, dont le rayon
    # est r2 et la hauteur b2. Lire b2 comme la largeur radiale d'une couronne
    # rendait l'epaisseur de la paroi du solide, un millimetre.
    if assembly.outlet is not None and assembly.outlet.radial:
        topology.r_2 = assembly.outlet.radius
        topology.r_2h = topology.r_2s = assembly.outlet.radius
        topology.b_2 = assembly.outlet.free_height
    else:
        topology.b_2 = topology.r_2s - topology.r_2h

    hub = assembly.component(components_module.SLOT_HUB)
    hub_radius = (
        components_module.hub_radius_at(hub, assembly.inlet.centroid[2], axis)
        if hub is not None and assembly.inlet is not None else None
    )
    if hub is not None and hub_radius is not None:
        topology.hub_kind = topology_module.HUB_SOLID
        topology.r_1h = max(topology.r_1h, hub_radius)
        topology.confidence.set("rayon_de_moyeu", HIGH)
        topology.notes.append(
            f"rayon de moyeu lu sur le STL de moyeu, coupe au plan d'entree : "
            f"{hub_radius * config.MM_PER_M:.1f} mm."
        )
    elif hub is not None:
        topology.hub_kind = (
            topology_module.HUB_BORE if topology.r_1h > 0.0 else topology_module.HUB_NONE
        )
        topology.confidence.set("rayon_de_moyeu", HIGH)
        topology.notes.append(
            "le STL de moyeu ne porte pas de matiere sur l'axe au plan d'entree : rien n'y "
            "obstrue le centre, et r1h est le rayon interieur du solide d'entree "
            f"({topology.r_1h * config.MM_PER_M:.1f} mm)."
        )
    else:
        topology.hub_kind = (
            topology_module.HUB_BORE if topology.r_1h > 0.0 else topology_module.HUB_NONE
        )
        topology.bore_radius = topology.r_1h
        topology.confidence.set("rayon_de_moyeu", MEDIUM)
        topology.notes.append(
            "aucun STL de moyeu fourni : r1h est le rayon interieur du solide d'entree, ce qui "
            "le confond avec le percement de la veine. Importez le moyeu pour lever le doute."
        )

    topology.blades.n_blades = assembly.declarations.n_blades
    topology.blades.forced = True
    # Le type de roue se lit sur la sortie **mesuree** : une bande cylindrique
    # refoule radialement, une couronne plate axialement. Il etait deduit du
    # seul mode declare, qui faisait de toute pompe carenee une centrifuge.
    if assembly.declarations.mode == MACHINE_PROPELLER:
        topology.machine_type = topology_module.AXIAL
    elif assembly.outlet is not None and assembly.outlet.radial:
        topology.machine_type = topology_module.CENTRIFUGAL
    else:
        topology.machine_type = topology_module.AXIAL
    topology.ratio_r2_r1s = (
        topology.r_blade_tip / topology.r_1s if topology.r_1s > 0.0 else 0.0
    )
    for quantity in ("rayons", "sections", "nombre_de_pales", "type_de_roue", "axe"):
        topology.confidence.set(quantity, HIGH)
    # Une section dont une paroi occupe une part notable a ete corrigee, pas lue.
    planes = [p for p in (assembly.inlet, assembly.outlet) if p is not None]
    if any(1.0 - p.free > config.FREE_PASSAGE_TOLERANCE for p in planes) or \
            assembly.confidence.get_level("sections", HIGH) != HIGH:
        topology.confidence.set("sections", MEDIUM)
    return topology


def _read_isolated_angles(
    result: AnalysisResult, assembly: components_module.ComponentAssembly
) -> None:
    """Angles de pale lus sur la pale isolee (C3.3 et C3.4), roue a refoulement radial.

    L'angle d'une pale se mesure depuis la tangente, **dans le sens oppose a la
    rotation** : sans sens de rotation, on ne sait pas si une aube est courbee
    vers l'arriere ou vers l'avant, et l'on ne publie rien plutot que de choisir.

    Un angle impose par --beta1 ou --beta2 l'emporte sur la lecture, comme en
    import en vrac ; la lecture reste publiee a cote, pour comparaison.
    """
    options = result.options
    declared = (options.beta1_deg, options.beta2_deg)
    blades = result.blades
    if declared[0] is not None:
        blades.beta1_deg = blade_module._clamp_beta(float(declared[0]))
    if declared[1] is not None:
        blades.beta2_deg = blade_module._clamp_beta(float(declared[1]))
    both_declared = None not in declared
    if both_declared:
        blades.forced_beta = True
        result.confidence.set("angles_de_pale", HIGH)

    blade = assembly.component(components_module.SLOT_BLADE)
    outlet = assembly.component(components_module.SLOT_OUTLET)
    if blade is None or blade.mesh is None or assembly.outlet is None or assembly.inlet is None:
        return
    if not assembly.outlet.radial:
        if not both_declared:
            result.warnings.append(
                "refoulement axial : les angles de pale se lisent sur des coupes cylindriques, "
                "que le mode composants ne fait pas encore. Imposez --beta1 et --beta2."
            )
        return
    sign = blades.rotation_sign
    if not sign:
        if not both_declared:
            result.warnings.append(
                "sens de rotation inconnu : l'angle d'une pale se mesure a l'oppose de la "
                "rotation, et sans elle on ne sait pas si une aube est courbee vers l'arriere ou "
                "vers l'avant. Declarez --rotation pour lire les angles."
            )
        return
    toroidal = assembly.declarations.blade_topology == components_module.BLADE_TOROIDAL
    angles = isolated_module.read_isolated_blade(
        blade.mesh,
        assembly.axis_origin,
        sign,
        (outlet.low[2], outlet.high[2]),
        assembly.inlet.centroid[2],
        toroidal,
        walls=components_module.meridian(assembly),
        outlet_radius=assembly.outlet.radius,
    )
    result.isolated_angles = angles
    result.warnings.extend(angles.warnings)
    blades.notes.extend(angles.notes)
    if angles.working is None:
        return
    entree = angles.axial_push if angles.axial_push and angles.axial_push.is_inlet else None
    if declared[0] is None:
        blades.beta1_deg = entree.helix_rms_deg if entree else angles.beta1_deg
    if declared[1] is None:
        blades.beta2_deg = angles.beta2_deg
    if not both_declared:
        result.confidence.set("angles_de_pale", angles.confidence)
    lecture = (
        f"angles lus sur la pale isolee, coupes horizontales {angles.working.du} : "
        + ("angle d'attaque" if entree else "beta1")
        + f" = {angles.beta1_deg:.1f} deg a r = {angles.r_le * config.MM_PER_M:.1f} mm, "
        f"beta2 = {angles.beta2_deg:.1f} deg a r = {angles.r_te * config.MM_PER_M:.1f} mm "
        "(moyennes sur la hauteur des bords). "
        "Un plan horizontal n'est la surface de courant que la ou l'ecoulement est radial : sur "
        "des aubes synthetiques d'angles connus, la lecture est exacte a 0,1 degre dans une veine "
        "plane, et basse de 2 a 4 degres dans une veine inclinee de 30 degres, que le plan coupe "
        "en biais."
    )
    if entree is not None:
        lecture += (
            f" L'eau entre dans la roue le long de l'axe, par le {entree.branch} qu'elle "
            f"rencontre d'abord : beta1 est l'angle de sa vis au rayon moyen du passage, "
            f"{entree.helix_rms_deg:.1f} deg a r = {entree.r_rms * config.MM_PER_M:.1f} mm. Le "
            f"bord d'attaque du {angles.working.name}, en aval, recoit une eau deja entrainee."
        )
    if declared != (None, None):
        lecture += (
            f" Valeurs retenues : beta1 = {blades.beta1_deg:.1f} deg, beta2 = "
            f"{blades.beta2_deg:.1f} deg, les angles imposes l'emportant sur la lecture."
        )
    blades.notes.append(lecture)
    if angles.layout_detail:
        result.warnings.insert(0, angles.layout_detail)


def _blade_edges(topology: topology_module.Topology, angles: isolated_module.IsolatedBladeAngles) -> None:
    """Rayons et sections du modele 1D, pris aux bords de la roue et non aux solides fluide.

    Le modele prend r1, A1 et beta1 la ou l'eau aborde les aubes, r2, A2 et beta2
    la ou elle les quitte. Chaque plan de coupe est une surface de courant
    d'epaisseur `level_step` : la section d'un bord est `2 pi r h`, `h` etant la
    hauteur des niveaux qui le portent, obstruction par l'epaisseur des aubes
    comprise (TAU_1, TAU_2), comme en import global. Quand l'eau entre le long
    de l'axe par un brin amont, l'entree est la couronne de passage qu'elle y
    traverse. La fente de sortie, elle, est au-dela du bord de fuite : l'eau y
    arrive par une couronne sans aube, ou son moment cinetique se conserve.
    """
    mm = config.MM_PER_M
    working = angles.working
    fente_r, fente_a = topology.r_2, topology.area_2
    entree = angles.axial_push if angles.axial_push and angles.axial_push.is_inlet else None
    if entree is not None:
        topology.r_1 = entree.r_rms
        topology.r_1h, topology.r_1s = entree.passage
        topology.area_1 = entree.area * config.TAU_1
        origine = (
            f"entree par la vis du {entree.branch}, couronne de r = {entree.passage[0] * mm:.1f} a "
            f"{entree.passage[1] * mm:.1f} mm"
        )
    else:
        topology.r_1 = angles.r_le
        topology.r_1h, topology.r_1s = angles.r_le_range
        topology.area_1 = 2.0 * math.pi * angles.r_le * angles.le_height * config.TAU_1
        origine = (
            f"entree au bord d'attaque {working.du}, sur {angles.le_height * mm:.1f} mm de hauteur"
        )
    topology.r_2 = angles.r_te
    topology.r_2h, topology.r_2s = angles.r_te_range
    topology.b_2 = angles.te_height
    topology.area_2 = 2.0 * math.pi * angles.r_te * angles.te_height * config.TAU_2
    topology.confidence.set("rayons", MEDIUM)
    topology.confidence.set("sections", MEDIUM)
    topology.notes.append(
        f"modele de ligne moyenne pris aux bords de la roue : {origine} (r1 = "
        f"{topology.r_1 * mm:.1f} mm, A1 = {topology.area_1 * 1e4:.1f} cm2) ; sortie au bord de "
        f"fuite {working.du} (r2 = {topology.r_2 * mm:.1f} mm, b2 = {topology.b_2 * mm:.1f} mm, "
        f"A2 = {topology.area_2 * 1e4:.1f} cm2). La fente de sortie mesuree (r = "
        f"{fente_r * mm:.1f} mm, {fente_a * 1e4:.1f} cm2 libres) est au-dela du bord de fuite."
    )


def _components_hydraulics(result: AnalysisResult, options: Options) -> None:
    """Courbes caracteristiques en mode composants, sur la roue telle que la pale la montre."""
    topology, geometry = result.topology, result.blades
    if options.machine == MACHINE_PROPELLER:
        result.warnings.append(
            "helice libre declaree : la poussee n'est pas encore calculee en mode composants ; "
            "seule la geometrie est publiee."
        )
        return
    angles = result.isolated_angles
    if angles is not None and angles.axial_push is not None and angles.axial_push.upstream \
            and not angles.axial_push.is_inlet:
        result.warnings.append(
            f"courbes non calculees : pour le sens de rotation retenu, le {angles.axial_push.branch}, "
            "que l'eau rencontre en premier, ne la pousse pas vers le brin qui refoule."
        )
        return
    if geometry.beta2_deg <= 0.0 or geometry.beta1_deg <= 0.0:
        result.warnings.append(
            "courbes caracteristiques non calculees : il leur faut beta1 et beta2, lus sur la "
            "pale isolee ou imposes par --beta1 et --beta2."
        )
        return
    if angles is not None and angles.working is not None:
        _blade_edges(topology, angles)
        for quantity in ("rayons", "sections"):
            result.confidence.set(quantity, topology.confidence.get_level(quantity))
    # La confiance des angles et du sens vit dans la table du resultat en mode
    # composants ; le modele la lit sur la geometrie des aubes.
    for quantity in ("angles_de_pale", "sens_de_rotation"):
        geometry.confidence.set(quantity, result.confidence.get_level(quantity))
    assembly = result.assembly
    toroidal = assembly.declarations.blade_topology == components_module.BLADE_TOROIDAL
    strands = 2 if toroidal else 1
    _curves(result, options, topology, geometry)
    _curve_analysis(result, topology, geometry, strands=strands)
    if toroidal and result.curves:
        result.comparison = _components_comparison(result, angles, strands)
    if angles is None or len(angles.branches) < 2:
        return
    for quantity in ("debit", "hauteur", "puissance", "couple", "rendement", "npshr"):
        result.confidence.set(quantity, LOW)
    entree = angles.axial_push if angles.axial_push and angles.axial_push.is_inlet else None
    if entree is not None:
        result.warnings.append(
            f"courbes calculees sur la roue en serie, en confiance basse : entree par la vis du "
            f"{entree.branch} (beta1 = {geometry.beta1_deg:.1f} deg), sortie par le bord de fuite "
            f"du {angles.working.name} (beta2 = {geometry.beta2_deg:.1f} deg). Le travail d'Euler "
            "ne depend que de ces deux bords : l'entrainement que le premier brin donne a l'eau "
            "est interne a la roue. Ne sont pas modelisees l'incidence au bord d'attaque du "
            f"{angles.working.name}, ou l'eau arrive deja entrainee, ni les pertes du coude entre "
            "les deux brins."
        )
    else:
        result.warnings.append(
            f"courbes calculees sur le {angles.working.name} seul, en confiance basse : l'effet "
            "des autres brins n'est pas modelise."
        )


def _components_comparison(
    result: AnalysisResult, angles: isolated_module.IsolatedBladeAngles | None, strands: int
) -> comparison_module.Comparison | None:
    """Comparatif du mode composants : la meme roue, a aubes d'un seul brin."""
    assembly, topology = result.assembly, result.topology
    inlet = rim = None
    closed = None
    if angles is not None and angles.working is not None:
        entree = angles.axial_push if angles.axial_push and angles.axial_push.is_inlet else None
        if entree is not None:
            # Sans son brin amont, la roue normale aspire par le bord d'attaque
            # de son aube.
            inlet = comparison_module.PumpInlet(
                r_1=angles.r_le, r_1h=angles.r_le_range[0], r_1s=angles.r_le_range[1],
                area_1=2.0 * math.pi * angles.r_le * angles.le_height * config.TAU_1,
                beta1_deg=angles.beta1_deg,
                description=(
                    f"sans le {entree.branch}, elle aspire par le bord d'attaque du "
                    f"{angles.working.name}"
                ),
            )
        closed = _closed_wheel(assembly, angles)
    outlet = assembly.outlet
    if outlet is not None and outlet.radial and topology.r_2 < (1.0 - config.EDGE_REACH) * outlet.radius:
        rim = comparison_module.PumpOutletRim(radius=outlet.radius, area=outlet.free_area)
    return comparison_module.pump(
        result.meanline_input, result.curves[0], topology.b_2, strands, closed,
        inlet=inlet, rim=rim, confidence=result.confidence.get_level("hauteur"),
    )


def _closed_wheel(
    assembly: components_module.ComponentAssembly, angles: isolated_module.IsolatedBladeAngles
) -> bool | None:
    """La grille qui refoule est-elle fermee, dessus et dessous, par les parois de la roue ?

    Sonde le plan meridien au rayon moyen du bord de fuite, de part et d'autre de
    la grille, sur trois fois sa hauteur. `None` sans piece de paroi.
    """
    walls = components_module.meridian(assembly)
    if walls is None or not walls.rotating:
        return None
    levels = angles.working.levels
    z_low, z_high = min(l.z for l in levels), max(l.z for l in levels)
    reach = 3.0 * max(angles.te_height, walls.cell)
    steps = max(4, int(reach / walls.cell))

    def wall_within(start: float, direction: int) -> bool:
        return any(walls.wall(angles.r_te, start + direction * reach * k / steps)
                   for k in range(1, steps + 1))

    return wall_within(z_high, +1) and wall_within(z_low, -1)


def _fill_provenance(result: AnalysisResult, options: Options) -> None:
    """Dit, grandeur par grandeur, si la valeur publiee est lue ou declaree.

    Le niveau de confiance repond « a quel point est-ce sur » ; la provenance
    repond « d'ou cela vient ». Une valeur imposee en ligne de commande peut
    etre parfaitement sure sans rien devoir au maillage, et la lire comme une
    mesure serait se tromper sur ce que l'outil a fait.
    """
    provenance = result.provenance
    provenance.declare("type_de_roue", options.wheel_type != WHEEL_AUTO)
    provenance.declare("modele_hydraulique", options.machine != MACHINE_AUTO)
    provenance.declare("nombre_de_pales", options.blades is not None)
    provenance.declare(
        "angles_de_pale", options.beta1_deg is not None and options.beta2_deg is not None
    )
    # Le sens de rotation non impose n'est pas « mesure » : l'outil refuse
    # justement de le trancher, et publie « a indiquer ». Le marquer mesure
    # laisserait croire a une lecture qui n'a pas eu lieu.
    provenance.set(
        "sens_de_rotation", DECLARE if options.rotation is not None else DEFAUT
    )
    provenance.declare("rayon_d_aspiration", options.r_aspiration_cm is not None)
    provenance.declare("cote_aspiration", options.suction != axis_module.SUCTION_AUTO)
    for quantity in ("rayons", "rayon_de_moyeu", "rayon_de_sortie", "sections",
                     "volume", "axe", "maillage"):
        provenance.declare(quantity, False)
    # Un seul des deux angles impose : la valeur publiee melange une declaration
    # et une lecture, ce que ni « declare » ni « mesure » ne decrit. On retient
    # la moins engageante des deux.
    if (options.beta1_deg is None) != (options.beta2_deg is None):
        provenance.set("angles_de_pale", DEFAUT)
    # Des angles que rien n'a lus ne sont pas « mesures » : zero est la valeur
    # d'attente, pas une lecture.
    if result.blades is not None and result.blades.beta2_deg <= 0.0 and options.beta2_deg is None:
        provenance.set("angles_de_pale", DEFAUT)
    # Les conditions de site ne sont jamais mesurees : elles viennent de la
    # ligne de commande ou de config.
    provenance.set(
        "installation",
        DECLARE if any((
            options.altitude != config.ALTITUDE,
            options.temperature_c != config.TEMPERATURE,
            options.suction_height != config.HAUTEUR_ASPIRATION,
            options.suction_losses != config.PERTES_ASPIRATION,
        )) else DEFAUT,
    )
