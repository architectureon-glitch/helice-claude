"""Topologie de la roue : nombre de pales, rayons, type, section de sortie (SPEC phase 3).

Tout est lu dans la carte d'occupation de la phase 2 ; le maillage n'est
re-sollicite que pour le controle croise de symetrie par distance de Hausdorff.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import HIGH, LOW, MEDIUM, ConfidenceMap, worst
from ..mesh import TriMesh, rotation_matrix
from .occupancy import OccupancyMap
from .proximity import hausdorff_distance

#: Ce qu'on trouve au centre de la piece, et qui n'est pas la meme chose.
HUB_SOLID = "moyeu plein"  # de la matiere pleine depuis l'axe : un vrai moyeu
HUB_BORE = "alesage traversant"  # un trou de part en part : pas de moyeu, mais pas de passage non plus
HUB_NONE = "ni moyeu ni alesage"  # rien de continu au centre : le fluide y passe

AXIAL = "axiale"
MIXED = "mixte"
CENTRIFUGAL = "centrifuge"


@dataclass
class BladeCount:
    """Resultat du comptage de pales par analyse spectrale (SPEC 3.1)."""

    n_blades: int = 0
    amplitudes: list[float] = field(default_factory=list)  # amplitudes des harmoniques 0..BLADES_MAX
    ratio_to_runner_up: float = 0.0  # rapport a la deuxieme amplitude, toutes harmoniques
    ratio_to_competing: float = 0.0  # rapport a la plus forte harmonique non multiple de N
    hausdorff_relative: float | None = None  # distance apres rotation de 2*pi/N, rapportee a r_tip
    # Periodicite propre de la piece -- ecart apres 2*pi/N pour le N qu'elle porte,
    # detecte ou corrige -- independante du nombre impose : c'est elle qui dit si
    # la piece est une roue reguliere.
    periodic_order: int = 0
    periodicity: float | None = None
    confidence: str = LOW
    forced: bool = False  # vrai si la valeur vient de --blades
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "nombre_de_pales": self.n_blades,
            "rapport_a_la_suivante": self.ratio_to_runner_up,
            "rapport_aux_harmoniques_concurrentes": self.ratio_to_competing,
            "hausdorff_relatif": self.hausdorff_relative,
            "ordre_de_periodicite_de_la_piece": self.periodic_order or None,
            "ecart_de_periodicite_de_la_piece": self.periodicity,
            "impose_par_l_utilisateur": self.forced,
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
        }


@dataclass
class Topology:
    """Rayons caracteristiques, type de roue et section de sortie (SPEC 3.2 a 3.4)."""

    blades: BladeCount = field(default_factory=BladeCount)
    r_tip: float = 0.0  # rayon exterieur de la matiere (flasque compris)
    r_blade_tip: float = 0.0  # rayon exterieur des pales
    z_1: float = 0.0
    r_1s: float = 0.0
    r_1h: float = 0.0
    hub_kind: str = HUB_NONE  # moyeu plein, alesage traversant, ou ni l'un ni l'autre
    bore_radius: float = 0.0  # rayon interieur de la matiere quand il y a un alesage
    r_1: float = 0.0
    r_aspiration: float = 0.0
    r_aspiration_source: str = "detecte"  # "detecte" ou "utilisateur"
    z_2: float = 0.0
    r_2s: float = 0.0
    r_2h: float = 0.0
    r_2: float = 0.0
    b_2: float = 0.0
    area_1: float = 0.0
    area_2: float = 0.0
    ratio_r2_r1s: float = 0.0
    machine_type: str = AXIAL
    closed_impeller: bool = False
    rotation_ambiguity: str = ""  # reserve sur le sens, a taire si l'utilisateur l'a donne
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)  # remarques qui n'appellent pas d'action

    def to_dict(self) -> dict:
        """Vue serialisable en JSON (grandeurs en SI)."""
        return {
            "pales": self.blades.to_dict(),
            "type_de_roue": self.machine_type,
            "roue_fermee": self.closed_impeller,
            "reserve_sur_le_sens": self.rotation_ambiguity or None,
            "r_tip_m": self.r_tip,
            "r_tip_pales_m": self.r_blade_tip,
            "z_1_m": self.z_1,
            "r_1s_m": self.r_1s,
            "r_1h_m": self.r_1h,
            "nature_du_centre": self.hub_kind,
            "rayon_d_alesage_m": self.bore_radius or None,
            "r_1_m": self.r_1,
            "r_aspiration_m": self.r_aspiration,
            "r_aspiration_source": self.r_aspiration_source,
            "z_2_m": self.z_2,
            "r_2s_m": self.r_2s,
            "r_2h_m": self.r_2h,
            "r_2_m": self.r_2,
            "b_2_m": self.b_2,
            "A_1_m2": self.area_1,
            "A_2_m2": self.area_2,
            "rapport_r2_sur_r1s": self.ratio_r2_r1s,
            "confiance": dict(self.confidence),
            "avertissements": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# 3.1 Nombre de pales
# ---------------------------------------------------------------------------
def count_blades(
    occupancy: OccupancyMap,
    mesh: TriMesh | None = None,
    forced: int | None = None,
    samples: int = config.HAUSDORFF_SAMPLES,
    cap: float = config.HAUSDORFF_CAP,
    verdict_only: bool = False,
) -> BladeCount:
    """Compte les pales par transformee de Fourier du signal `g(theta)`.

    Le nombre de pales est l'indice de l'harmonique dominante dans
    `[BLADES_MIN, BLADES_MAX]`.  Le spectre est celui de
    `OccupancyMap.theta_spectrum`, qui transforme cellule par cellule avant de
    sommer les modules : transformer le signal deja integre sur la zone de pales
    fait disparaitre l'harmonique N des que les pales se recouvrent en
    projection, ce qui est le cas usuel d'une roue a fort enroulement.

    Deux rapports d'amplitude sont calcules et rapportes.  Le critere de
    confiance retient le second : une roue a N pales produit **toujours** des
    harmoniques fortes en 2N, 3N... -- ce sont des consequences mecaniques de la
    periodicite d'ordre N, pas des hypotheses concurrentes.  Comparer l'harmonique
    dominante a la plus forte harmonique **qui n'est pas un multiple de N** est
    donc le seul rapport discriminant ; le rapport brut a la deuxieme amplitude
    est neanmoins reporte, conformement a la lettre de la SPEC.
    """
    result = BladeCount()
    amplitudes = occupancy.theta_spectrum(config.BLADES_MAX)
    result.amplitudes = amplitudes

    candidates = range(config.BLADES_MIN, config.BLADES_MAX + 1)
    dominant = max(candidates, key=lambda k: amplitudes[k])
    peak = amplitudes[dominant]

    runner_up = max((amplitudes[k] for k in candidates if k != dominant), default=0.0)
    competing = max(
        (amplitudes[k] for k in candidates if k % dominant != 0),
        default=0.0,
    )
    result.ratio_to_runner_up = peak / runner_up if runner_up > 0.0 else math.inf
    result.ratio_to_competing = peak / competing if competing > 0.0 else math.inf

    if peak <= 0.0:
        result.n_blades = 0
        result.confidence = LOW
        result.warnings.append(
            "aucune modulation angulaire detectee : la zone de pales est vide ou "
            "la roue est un solide de revolution"
        )
    else:
        result.n_blades = dominant
        if result.ratio_to_competing >= config.FFT_RATIO_MIN:
            result.confidence = HIGH
        elif result.ratio_to_competing >= config.FFT_RATIO_MEDIUM:
            result.confidence = MEDIUM
            result.warnings.append(
                f"harmonique dominante peu marquee (rapport {result.ratio_to_competing:.2f} "
                f"< {config.FFT_RATIO_MIN}) : le nombre de pales est incertain"
            )
        else:
            result.confidence = LOW
            result.warnings.append(
                f"spectre angulaire ambigu (rapport {result.ratio_to_competing:.2f}) : "
                "fournissez --blades si le nombre de pales est connu"
            )

    def ecart(order: int, points: int, verdict: bool) -> float:
        """Hausdorff relatif apres 2*pi/order ; borne inferieure si `verdict`."""
        rotated = mesh.transformed(rotation_matrix((0.0, 0.0, 1.0), 2.0 * math.pi / order))
        stop = config.SYM_TOL * occupancy.r_max if verdict else math.inf
        return hausdorff_distance(
            mesh, rotated, points, cap * occupancy.r_max, stop
        ) / occupancy.r_max

    checkable = mesh is not None and occupancy.r_max > 0.0

    # Periodicite propre de la piece, au nombre que le spectre lui lit. Le
    # spectre a pu retenir une harmonique de la forme d'aube plutot que le
    # nombre d'aubes : une roue de 3 aubes en boucle rend un pic en 9. La
    # periodicite, elle, ne se trompe pas -- la piece se superpose a elle-meme
    # ou non. On essaie les diviseurs du nombre lu et les harmoniques les plus
    # fortes, et l'on retient le plus grand ordre qui tient : une roue de 6
    # aubes tient aussi en 3 et en 2.
    # `verdict_only` : l'appelant ne veut savoir que si la piece tient sous
    # SYM_TOL -- c'est le cas du departage des axes candidats.
    if checkable and result.n_blades >= config.BLADES_MIN:
        lu = result.n_blades
        result.periodic_order = lu
        result.periodicity = ecart(lu, samples, verdict_only)
        if result.periodicity > config.SYM_TOL:
            fortes = sorted(candidates, key=lambda k: amplitudes[k], reverse=True)[: config.PERIOD_CANDIDATES]
            essais = sorted(
                {k for k in fortes if k != lu} | {d for d in range(config.BLADES_MIN, lu) if lu % d == 0},
                reverse=True,
            )
            for order in essais:
                if amplitudes[order] <= 0.0:
                    continue
                if ecart(order, config.AXIS_PROBE_SAMPLES, True) <= config.SYM_TOL:
                    result.periodic_order = order
                    result.periodicity = ecart(order, samples, verdict_only)
                    result.n_blades = order
                    result.confidence = MEDIUM
                    result.warnings.append(
                        f"le spectre angulaire designait {lu} pales, mais la piece ne se superpose "
                        f"pas a elle-meme apres 2*pi/{lu} ; elle le fait apres 2*pi/{order}. "
                        f"C'est {order} qui est retenu : le pic en {lu} est une harmonique de la "
                        "forme des aubes, pas leur nombre. Confiance moyenne ; imposez --blades si "
                        "vous connaissez la roue."
                    )
                    break
        result.hausdorff_relative = result.periodicity

    if forced is not None:
        if result.n_blades and forced != result.n_blades:
            result.warnings.append(
                f"nombre de pales impose a {forced} alors que le spectre en detecte {result.n_blades}"
            )
        result.n_blades = int(forced)
        result.forced = True
        result.confidence = HIGH
        # Le controle porte aussi sur le nombre **retenu** : c'est lui qui
        # nourrit le calcul. Il ne portait que sur le nombre detecte, si bien
        # qu'un nombre impose que la geometrie dementait passait en confiance
        # haute. Une declaration n'est pas annulee pour autant : elle est
        # accompagnee.
        if checkable and result.n_blades >= config.BLADES_MIN:
            result.hausdorff_relative = (
                result.periodicity if result.n_blades == result.periodic_order
                else ecart(result.n_blades, samples, verdict_only)
            )

    if result.hausdorff_relative is not None and result.hausdorff_relative > config.SYM_TOL:
        result.confidence = worst(result.confidence, MEDIUM)
        cause = (
            "le nombre impose ne correspond pas a la periodicite de la piece, ou les pales "
            "sont inegales"
            if result.forced
            else "pales inegales, roue tronquee ou nombre de pales errone"
        )
        result.warnings.append(
            f"la rotation de 2*pi/{result.n_blades} ne superpose pas le maillage a lui-meme "
            f"(Hausdorff {result.hausdorff_relative:.3f} du rayon exterieur, seuil {config.SYM_TOL}) : "
            f"{cause}"
        )
    return result


def has_full_ring(occupancy: OccupancyMap) -> bool:
    """Vrai si au moins une cellule de la carte est pleine sur tout le tour.

    Moyeu plein, moyeu alese, flasque, jante : toute roue a un anneau de matiere
    qui tient ses pales ensemble autour de l'axe. Deux pieces posees cote a cote
    peuvent se superposer par un demi-tour -- deux roues identiques le font --
    sans qu'aucun cercle centre sur l'axe ne soit entierement dans la matiere.
    """
    return any(value >= config.F_SOLIDE for row in occupancy.f for value in row)


def admissibility(blades: "BladeCount", occupancy: OccupancyMap | None = None) -> str | None:
    """Motif de refus si la piece n'est pas une roue reguliere autour de l'axe retenu.

    Une roue de N pales se superpose a elle-meme apres une rotation de 2*pi/N :
    c'est la seule propriete que tout le calcul suppose sans jamais la dire.
    Mesuree sur des roues valides, l'ecart reste sous 2 % du rayon exterieur,
    meme maillage dechire. Au-dela de `SYM_REJECT`, ce n'est plus une roue
    imparfaite, c'est une piece que le modele ne decrit pas : un corps parasite
    exporte avec la roue, deux pieces dans le meme fichier, un carter, un axe
    detecte de travers. Publier une hauteur et un debit, meme en confiance
    faible, laisserait croire a une estimation incertaine la ou il n'y a pas
    d'estimation du tout. `None` si la piece est recevable ou non verifiee.
    """
    if occupancy is not None and occupancy.nr and not has_full_ring(occupancy):
        return (
            "aucun anneau de matiere ne fait le tour de l'axe retenu : ni moyeu, ni alesage, ni "
            "flasque, ni jante. Une roue en a toujours un, qui tient ses pales ensemble ; ce "
            "fichier contient plutot plusieurs pieces disjointes, une pale seule, ou une piece "
            "dont l'axe a ete mal detecte. Aucune performance n'est publiee : exportez la roue "
            "entiere, seule dans le fichier. Pour analyser une pale isolee, utilisez le mode "
            "composants (--pale, --entree-fluide, --sortie-fluide)."
        )
    # La periodicite **propre** de la piece, pas l'ecart au nombre impose : un
    # nombre declare faux est une declaration a accompagner d'un avertissement,
    # pas une raison de refuser une roue reguliere.
    ecart = blades.periodicity
    if ecart is None or ecart <= config.SYM_REJECT:
        return None
    borne = "au moins " if ecart >= config.HAUSDORFF_CAP - 1e-9 else ""
    return (
        f"la piece ne se superpose pas a elle-meme apres une rotation de 2*pi/{blades.periodic_order} "
        f"autour de l'axe retenu : ecart de {borne}{ecart:.0%} du rayon exterieur, pour {config.SYM_REJECT:.0%} "
        "admis. Ce n'est pas une roue reguliere autour de cet axe -- corps parasite exporte avec "
        "la roue, plusieurs pieces dans le meme fichier, carter ou volute inclus, ou axe mal "
        "detecte. Aucune performance n'est publiee : isolez la roue seule dans le fichier. Si le "
        "nombre de pales est faux, imposez-le par --blades ; si les pales sont volontairement "
        "inegales, --sans-controle-symetrie leve ce controle, a vos risques."
    )


# ---------------------------------------------------------------------------
# 3.2 Rayons caracteristiques
# ---------------------------------------------------------------------------
def _hub_index(row: list[float]) -> int:
    """Dernier index radial du moyeu : plage pleine contiguë depuis l'axe."""
    last = -1
    for ir, value in enumerate(row):
        if value >= config.F_SOLIDE:
            last = ir
        else:
            break
    return last


def _inner_material_index(row: list[float]) -> int:
    """Plus petit index radial ou il y a de la matiere ; -1 si la rangee est vide."""
    for ir, value in enumerate(row):
        if value > config.F_MATIERE:
            return ir
    return -1


def hub_nature(occupancy: OccupancyMap, blade_rows: list[int], iz_1: int) -> tuple[str, int]:
    """Ce que la piece porte en son centre, et le rayon qui va avec.

    Trois cas, qu'il ne faut pas confondre parce qu'ils ne se comportent pas
    pareil dans la section d'entree :

    * **moyeu plein** -- de la matiere pleine depuis l'axe.  `r_1h` est le rayon
      du moyeu, et `A1 = pi (r_1s^2 - r_1h^2)` retire bien ce que le moyeu
      occupe.
    * **alesage traversant** -- un trou de part en part, pour l'arbre ou pour le
      montage.  Il n'y a pas de moyeu plein, mais le trou n'est pas non plus une
      section de passage : le fluide n'y circule pas.  Rendre zero serait
      exact au sens du critere du moyeu et **faux** au sens hydraulique, car
      `A1` compterait alors le trou comme du passage et gonflerait le debit.
      `r_1h` vaut donc le rayon interieur de la matiere.
    * **ni l'un ni l'autre** -- rien de continu au centre : le fluide y passe
      vraiment, et `r_1h` vaut zero.

    Renvoie `(nature, index radial)`.  L'alesage n'est reconnu que s'il traverse :
    un evidement borgne, lui, laisse de la matiere en face et se voit sur les
    autres rangees.

    La question porte sur ce qui **obstrue l'entree**, et se tranche donc au plan
    d'entree `iz_1`.  Un plateau arriere de roue centrifuge est bien du plein
    depuis l'axe, mais il est a l'autre bout de la piece et n'obstrue rien a
    l'aspiration : le compter donnerait un rayon de moyeu superieur au rayon
    d'oeillard, et une section `A1 = pi (r_1s^2 - r_1h^2)` negative.  Seule la
    reconnaissance de l'**alesage** regarde toute la hauteur, puisqu'un alesage
    n'en est un que s'il traverse.
    """
    rows = [occupancy.f[iz] for iz in blade_rows] or occupancy.f
    rows = [row for row in rows if _outer_index(row) >= 0]
    if not rows:
        return HUB_NONE, -1

    hub = _hub_index(occupancy.f[iz_1])
    if hub >= 0:
        return HUB_SOLID, hub

    interieurs = [_inner_material_index(row) for row in rows]
    if all(index > 0 for index in interieurs):
        # Aucune rangee ne porte de matiere sur l'axe : le vide central traverse.
        return HUB_BORE, min(interieurs)
    return HUB_NONE, -1


def _outer_index(row: list[float]) -> int:
    """Plus grand index radial ou il y a de la matiere."""
    for ir in range(len(row) - 1, -1, -1):
        if row[ir] > config.F_MATIERE:
            return ir
    return -1


def _blade_outer_index(blade_row: list[bool], row: list[float]) -> int:
    """Plus grand index radial occupe par une pale ; a defaut, par de la matiere."""
    for ir in range(len(blade_row) - 1, -1, -1):
        if blade_row[ir]:
            return ir
    return _outer_index(row)


def _tip_height(
    occupancy: OccupancyMap, blade: list[list[bool]], r_tip: float
) -> float:
    """Hauteur de pale au rayon exterieur, c'est-a-dire la section de refoulement radial."""
    column = occupancy.radial_index(config.B2_RADIUS_FRACTION * r_tip)
    rows = [iz for iz in range(occupancy.nz) if blade[iz][column]]
    if not rows:
        return 0.0
    return occupancy.z_centres[rows[-1]] - occupancy.z_centres[rows[0]] + occupancy.dz


def _shroud_eye(occupancy: OccupancyMap, iz_1: int) -> tuple[int, int]:
    """Index radial de l'oeillard, sur une roue fermee (SPEC 3.2, cas ferme).

    Au-dela du bord d'attaque, cote aspiration, le flasque avant ne laisse
    qu'une couronne pleine percee en son centre.  C'est ce percement qui fixe
    la section d'entree, pas l'extremite des pales : sur une roue fermee les
    aubes courent jusque sous le flasque, et lire `r_1s` sur elles donne le
    rayon exterieur de la roue au lieu du rayon d'aspiration.  Le rayon retenu
    est le plus petit rayon perce des rangees de flasque, c'est-a-dire le col.

    Renvoie `(index de l'oeillard, index du moyeu a ce plan)`, `(-1, -1)` si
    aucune rangee de flasque n'est reconnue.
    """
    best, best_hub = -1, -1
    for iz in range(occupancy.nz - 1, iz_1, -1):
        row = occupancy.f[iz]
        outer = _outer_index(row)
        if outer < 0:  # rangee vide : au-dessus de la piece
            continue
        # Bord interieur de la couronne exterieure : on remonte vers l'axe tant
        # qu'il y a de la matiere.  Le test porte sur la continuite de la
        # couronne, pas sur la vacuite de l'oeillard : une roue fermee peut
        # tres bien porter un bossage d'arbre en son centre, qui donnera r_1h.
        inner = outer
        while inner > 0 and row[inner - 1] > config.F_MATIERE:
            inner -= 1
        if inner == 0:  # la matiere touche l'axe : ce n'est pas un flasque perce
            continue
        if not any(row[ir] >= config.F_SOLIDE for ir in range(inner, outer + 1)):
            continue
        if best < 0 or inner < best:
            best, best_hub = inner, _hub_index(row)
    return best, best_hub


def _shroud_present(row: list[float]) -> bool:
    """Vrai si du plein reapparait au-dela d'une cellule de pale (flasque avant)."""
    seen_blade = False
    for ir in range(_hub_index(row) + 1, len(row)):
        value = row[ir]
        if config.F_MATIERE < value < config.F_SOLIDE:
            seen_blade = True
        elif value >= config.F_SOLIDE and seen_blade:
            return True
    return False


def clean_blade_mask(blade: list[list[bool]]) -> list[list[bool]]:
    """Retire les ilots parasites de la zone de pales.

    Une surface de revolution facettee produit, la ou elle coupe une cellule,
    quelques cellules de `f` intermediaire qui ne sont pas des pales.  Ces
    cellules forment des ilots minuscules, alors que la vraie zone de pales est
    d'un seul tenant : on ne garde donc que les composantes connexes
    suffisamment grosses.  Filtrer par la taille des ilots plutot que par le
    nombre de cellules par ligne laisse le bord d'attaque resolu a la cellule
    pres, ce dont depend la precision de r_1s.
    """
    nz = len(blade)
    nr = len(blade[0]) if nz else 0
    label = [[-1] * nr for _ in range(nz)]
    sizes: list[int] = []
    for iz in range(nz):
        for ir in range(nr):
            if not blade[iz][ir] or label[iz][ir] >= 0:
                continue
            current = len(sizes)
            stack = [(iz, ir)]
            label[iz][ir] = current
            count = 0
            while stack:
                z, r = stack.pop()
                count += 1
                for dz, dr in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    nz_, nr_ = z + dz, r + dr
                    if 0 <= nz_ < nz and 0 <= nr_ < nr and blade[nz_][nr_] and label[nz_][nr_] < 0:
                        label[nz_][nr_] = current
                        stack.append((nz_, nr_))
            sizes.append(count)
    if not sizes:
        return blade
    threshold = max(config.BLADE_BLOB_MIN_CELLS, config.BLADE_BLOB_MIN_FRACTION * max(sizes))
    return [
        [blade[iz][ir] and sizes[label[iz][ir]] >= threshold for ir in range(nr)]
        for iz in range(nz)
    ]


def characteristic_radii(occupancy: OccupancyMap) -> Topology:
    """Rayons, plans d'entree et de sortie, et type de roue (SPEC 3.2 a 3.4)."""
    topology = Topology()
    blade = clean_blade_mask(occupancy.blade_mask())
    nz, nr = occupancy.nz, occupancy.nr

    blade_rows = [
        iz for iz in range(nz)
        if sum(1 for value in blade[iz] if value) >= config.BLADE_ROW_MIN_CELLS
    ]
    if not blade_rows:
        topology.warnings.append(
            "aucune zone de pales identifiee : la geometrie est-elle un solide de "
            "revolution, ou les seuils F_MATIERE / F_SOLIDE sont-ils a revoir ?"
        )
        topology.confidence.set("rayons", LOW)
        topology.confidence.set("type_de_roue", LOW)
        return topology

    # Rayon exterieur de la matiere : c'est le diametre hors tout de la piece.
    r_tip_index = 0
    for ir in range(nr - 1, -1, -1):
        if any(occupancy.f[iz][ir] > config.F_MATIERE for iz in range(nz)):
            r_tip_index = ir
            break
    topology.r_tip = occupancy.r_centres[r_tip_index]

    # Rayon exterieur des **pales**, qui n'est pas le meme : sur une roue semi
    # ouverte le disque arriere deborde souvent les aubes, et sur une roue
    # fermee le flasque avant aussi. C'est le rayon des pales qui fixe u2, donc
    # la hauteur ; le prendre sur la matiere surestimerait les performances.
    blade_tip_index = 0
    for ir in range(nr - 1, -1, -1):
        if any(blade[iz][ir] for iz in range(nz)):
            blade_tip_index = ir
            break
    topology.r_blade_tip = occupancy.r_centres[blade_tip_index] or topology.r_tip

    iz_1 = blade_rows[-1]  # bord d'attaque : z maximal de la zone de pales
    iz_2 = blade_rows[0]  # bord de fuite : z minimal
    topology.z_1 = occupancy.z_centres[iz_1]
    topology.z_2 = occupancy.z_centres[iz_2]

    # r_1s et r_2s sont les rayons exterieurs **des pales**, pas de la matiere :
    # sur une roue fermee, le flasque avant s'etend bien au-dela des pales au
    # plan d'aspiration, et le compter donnerait un rayon d'oeillard trop grand.
    topology.closed_impeller = (
        sum(1 for iz in blade_rows if _shroud_present(occupancy.f[iz])) >= config.MIN_BLADE_SECTIONS
    )

    row_1 = occupancy.f[iz_1]
    outer_1 = _blade_outer_index(blade[iz_1], row_1)
    topology.r_1s = occupancy.r_centres[outer_1] if outer_1 >= 0 else topology.r_tip

    # Ce qu'il y a au centre : un moyeu plein, un alesage traversant, ou rien.
    # Les deux premiers obstruent la section d'entree, le troisieme non -- et un
    # alesage rendu comme un moyeu absent ferait compter le trou comme du
    # passage, ce qui gonfle le debit d'autant.
    topology.hub_kind, hub_1 = hub_nature(occupancy, blade_rows, iz_1)
    topology.r_1h = occupancy.r_centres[hub_1] if hub_1 >= 0 else 0.0
    if topology.hub_kind == HUB_BORE:
        topology.bore_radius = topology.r_1h
        topology.confidence.set("rayon_de_moyeu", MEDIUM)
        topology.notes.append(
            f"aucun moyeu plein : la piece est percee de part en part, la matiere ne commence "
            f"qu'a r = {topology.r_1h * config.MM_PER_M:.1f} mm. C'est ce rayon-la qui est publie "
            "en r1h, et la section d'entree exclut le trou : le fluide n'y circule pas, et l'y "
            "compter gonflerait le debit. Ce n'est pas un rayon de moyeu pour autant, d'ou la "
            "confiance moyenne."
        )
    elif topology.hub_kind == HUB_NONE:
        topology.confidence.set("rayon_de_moyeu", MEDIUM)
        topology.notes.append(
            "aucun moyeu plein ni alesage traversant : le centre est ouvert, r1h vaut zero et "
            "toute la section interieure compte comme passage. Verifiez que la piece n'a pas "
            "de moyeu que le maillage aurait perdu."
        )
    else:
        topology.confidence.set("rayon_de_moyeu", HIGH)

    # Roue fermee : l'entree est le percement du flasque, pas le bout des pales.
    if topology.closed_impeller:
        eye, eye_hub = _shroud_eye(occupancy, iz_1)
        if eye > 0:
            topology.r_1s = occupancy.r_centres[eye]
            topology.r_1h = occupancy.r_centres[eye_hub] if eye_hub >= 0 else 0.0
        else:
            topology.warnings.append(
                "roue fermee mais oeillard introuvable : le rayon d'aspiration est lu sur "
                "les pales, ce qui le surestime. Imposez-le par --r-aspiration."
            )
            topology.confidence.set("rayons", MEDIUM)

    topology.r_1 = math.sqrt((topology.r_1s ** 2 + topology.r_1h ** 2) / 2.0)
    topology.r_aspiration = topology.r_1s

    row_2 = occupancy.f[iz_2]
    outer_2 = _blade_outer_index(blade[iz_2], row_2)
    hub_2 = _hub_index(row_2)
    topology.r_2s = occupancy.r_centres[outer_2] if outer_2 >= 0 else topology.r_tip
    topology.r_2h = occupancy.r_centres[hub_2] if hub_2 >= 0 else 0.0

    # 3.3 Classification. Le rapport est pris sur les rayons exterieurs, seule
    # definition commune aux trois familles : pour une roue axiale le bord de
    # fuite est au meme rayon que le bord d'attaque, pour une centrifuge il est
    # au rayon exterieur de la roue.
    topology.ratio_r2_r1s = topology.r_blade_tip / topology.r_1s if topology.r_1s > 0.0 else 0.0
    if topology.ratio_r2_r1s < config.R_RATIO_AXIAL_MAX:
        topology.machine_type = AXIAL
    elif topology.ratio_r2_r1s < config.R_RATIO_MIXED_MAX:
        topology.machine_type = MIXED
    else:
        topology.machine_type = CENTRIFUGAL

    # 3.4 Sections.
    topology.area_1 = math.pi * (topology.r_1s ** 2 - topology.r_1h ** 2) * config.TAU_1
    if topology.machine_type == CENTRIFUGAL:
        topology.r_2 = topology.r_blade_tip
        column = occupancy.radial_index(config.B2_RADIUS_FRACTION * topology.r_blade_tip)
        heights = [occupancy.z_centres[iz] for iz in range(nz) if blade[iz][column]]
        topology.b_2 = (max(heights) - min(heights) + occupancy.dz) if heights else occupancy.dz
        topology.area_2 = 2.0 * math.pi * topology.r_2 * topology.b_2 * config.TAU_2
    elif topology.r_2h > 0.0 and (topology.r_2s - topology.r_2h) >= _tip_height(
        occupancy, blade, topology.r_blade_tip or topology.r_2s
    ):
        topology.r_2 = math.sqrt((topology.r_2s ** 2 + topology.r_2h ** 2) / 2.0)
        topology.b_2 = topology.r_2s - topology.r_2h
        topology.area_2 = math.pi * (topology.r_2s ** 2 - topology.r_2h ** 2) * config.TAU_2
    else:
        # Le refoulement se fait par la ou la section est ouverte : ici en
        # hauteur, a la peripherie. Deux cas amenent ici. Sans moyeu au plan de
        # sortie, le rayon quadratique moyen suppose une veine annulaire bordee
        # par un moyeu. Sans moyeu au plan de sortie il degenere en r_2s / racine(2),
        # qui n'est pas un rayon de refoulement mais un artefact de formule : sur
        # la roue toroidale de reference il donnait 118 mm pour des aubes qui
        # vont a 167, et faisait basculer u2 d'un tiers selon que la roue etait
        # classee mixte ou centrifuge -- deux familles que son rapport r2/r1s
        # separe justement a un millieme pres.
        topology.r_2 = topology.r_blade_tip or topology.r_2s
        column = occupancy.radial_index(config.B2_RADIUS_FRACTION * topology.r_2)
        heights = [occupancy.z_centres[iz] for iz in range(nz) if blade[iz][column]]
        topology.b_2 = (max(heights) - min(heights) + occupancy.dz) if heights else occupancy.dz
        topology.area_2 = 2.0 * math.pi * topology.r_2 * topology.b_2 * config.TAU_2
        topology.warnings.append(
            "refoulement lu comme radial : la section ouverte au bout des aubes est plus "
            "haute que large. Le rayon de refoulement est pris au bout des aubes, la moyenne "
            "quadratique moyeu-carter ne decrivant pas cette sortie -- soit qu'il n'y ait pas "
            "de moyeu au plan de sortie, soit que le plateau arriere y soit pris pour un moyeu."
        )

    level = HIGH
    if topology.area_1 <= 0.0 or topology.area_2 <= 0.0:
        level = LOW
        topology.warnings.append(
            "section d'entree ou de sortie nulle : les rayons extraits sont incoherents"
        )
    elif len(blade_rows) < config.MIN_BLADE_SECTIONS:
        level = MEDIUM
        topology.warnings.append(
            "zone de pales tres mince dans la direction axiale : rayons peu resolus"
        )
    if abs(topology.ratio_r2_r1s - config.R_RATIO_AXIAL_MAX) < config.VALID_GEOM_TOL or abs(
        topology.ratio_r2_r1s - config.R_RATIO_MIXED_MAX
    ) < config.VALID_GEOM_TOL:
        level = worst(level, MEDIUM)
        if abs(topology.ratio_r2_r1s - config.R_RATIO_MIXED_MAX) < config.VALID_GEOM_TOL:
            # Les deux familles n'appliquent pas la meme regle de sens de rotation :
            # a cheval sur la frontiere, la **suggestion** geometrique est un
            # tirage au sort. La reserve n'est mise en avant que si l'utilisateur
            # n'a pas donne le sens lui-meme -- sinon elle lui demanderait de
            # verifier ce qu'il vient d'affirmer.
            topology.rotation_ambiguity = (
                f"rapport r2/r1s = {topology.ratio_r2_r1s:.3f} a un millieme de la frontiere "
                "mixte / centrifuge, or les deux familles donnent des sens de rotation "
                "**opposes** : la suggestion geometrique n'est pas fiable ici. Donnez le sens "
                "par --rotation, apres l'avoir lu sur la piece -- une aube de pompe fuit le "
                "sens de rotation quand le rayon croit."
            )
        topology.warnings.append(
            f"rapport r2/r1s = {topology.ratio_r2_r1s:.3f} a la frontiere de deux familles : "
            "le type de roue est incertain"
        )
    topology.confidence.set("rayons", level)
    topology.confidence.set("type_de_roue", level)
    topology.confidence.set("sections", level)
    return topology


def apply_user_suction_radius(topology: Topology, r_aspiration_cm: float | None) -> Topology:
    """Applique `--r-aspiration` (en cm) : la valeur utilisateur prime toujours (SPEC 3.2)."""
    if r_aspiration_cm is None:
        if topology.confidence.get_level("rayons") == LOW:
            topology.warnings.append(
                "rayon d'aspiration detecte avec une confiance faible : fournissez "
                "--r-aspiration (en cm) pour fiabiliser tout le calcul hydraulique"
            )
        return topology
    value = float(r_aspiration_cm) * config.UNIT_FACTOR
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError("--r-aspiration doit etre strictement positif")
    # Un rayon impose prime sur la detection, mais pas sur la piece elle-meme :
    # hors de la roue, ou dans son moyeu, ce n'est plus une lecture differente,
    # c'est une faute d'unite (l'option est en centimetres).
    if topology.r_tip > 0.0 and value > topology.r_tip:
        raise ValueError(
            f"rayon d'aspiration impose a {value * config.MM_PER_M:.1f} mm, plus grand que la "
            f"roue elle-meme ({topology.r_tip * config.MM_PER_M:.1f} mm de rayon exterieur). "
            "--r-aspiration s'exprime en centimetres."
        )
    if value <= topology.r_1h:
        raise ValueError(
            f"rayon d'aspiration impose a {value * config.MM_PER_M:.1f} mm, dans le moyeu "
            f"({topology.r_1h * config.MM_PER_M:.1f} mm) : la section d'entree serait negative. "
            "--r-aspiration s'exprime en centimetres."
        )
    if topology.r_1s > 0.0 and abs(value - topology.r_1s) / topology.r_1s > config.VALID_GEOM_TOL:
        topology.warnings.append(
            f"rayon d'aspiration impose a {value * config.MM_PER_M:.1f} mm alors que la detection "
            f"donne {topology.r_1s * config.MM_PER_M:.1f} mm : la valeur utilisateur est retenue"
        )
    topology.r_aspiration = value
    topology.r_aspiration_source = "utilisateur"
    topology.r_1s = value
    topology.r_1 = math.sqrt((topology.r_1s ** 2 + topology.r_1h ** 2) / 2.0)
    topology.area_1 = math.pi * (topology.r_1s ** 2 - topology.r_1h ** 2) * config.TAU_1
    topology.ratio_r2_r1s = topology.r_blade_tip / topology.r_1s
    topology.confidence.set("rayons", HIGH)
    return topology


def analyse(
    occupancy: OccupancyMap,
    mesh: TriMesh | None = None,
    forced_blades: int | None = None,
    r_aspiration_cm: float | None = None,
) -> Topology:
    """Chaine complete de la phase 3."""
    topology = characteristic_radii(occupancy)
    topology.blades = count_blades(occupancy, mesh, forced_blades)
    topology.confidence.set("nombre_de_pales", topology.blades.confidence)
    return apply_user_suction_radius(topology, r_aspiration_cm)


def specific_speed(rpm: float, flow: float, head: float) -> float:
    """Vitesse specifique `n_q = n * sqrt(Q) / H^0.75` (n en tr/min, Q en m3/s, H en m)."""
    if flow <= 0.0 or head <= 0.0:
        return 0.0
    return rpm * math.sqrt(flow) / head ** config.NQ_HEAD_EXPONENT


def type_from_specific_speed(n_q: float) -> str:
    """Famille de roue deduite de la vitesse specifique (SPEC 3.3, controle croise)."""
    if n_q < config.NQ_CENTRIFUGAL_MAX:
        return CENTRIFUGAL
    if n_q <= config.NQ_MIXED_MAX:
        return MIXED
    return AXIAL


def cross_check_type(topology: Topology, n_q: float) -> str | None:
    """Compare la classification geometrique et celle deduite de `n_q`.

    Renvoie l'avertissement a placer en tete de rapport, ou `None` si les deux
    classifications concordent (SPEC 3.3).
    """
    if n_q <= 0.0:
        return None
    by_speed = type_from_specific_speed(n_q)
    if by_speed == topology.machine_type:
        return None
    return (
        f"incoherence de classification : la geometrie donne une roue {topology.machine_type} "
        f"(r2/r1s = {topology.ratio_r2_r1s:.2f}) alors que la vitesse specifique n_q = {n_q:.1f} "
        f"correspond a une roue {by_speed}. Le modele hydraulique applique reste celui de la "
        f"classification geometrique ; verifiez les rayons extraits sur la carte d'occupation."
    )
