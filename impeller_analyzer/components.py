"""Import par composants declares (SPEC v2) : l'outil verifie au lieu de deviner.

Le mode d'import en vrac recoit un seul fichier et doit tout inferer : ou est
l'axe, combien il y a de pales, ou commence le moyeu, de quel cote entre le
fluide.  Chacune de ces inferences peut se tromper, et sur une helice a aubes en
boucle elles se trompent ensemble -- la piece d'essai etait lue « centrifuge » a
cinq pales quand elle est une helice axiale toroidale.

Ce mode-ci renverse la charge.  L'utilisateur declare ce qu'il sait -- le modele
hydraulique, la topologie de pale, le nombre de pales -- et fournit les pieces
separement.  L'outil ne devine plus rien de tout cela : il **verifie**, et chaque
controle produit un avertissement plutot qu'un blocage, sauf celui du repere
commun, qui est le seul dont la violation detruit silencieusement l'assemblage.

Le contrat d'export, impose a l'utilisateur et affiche dans l'aide, est ce qui
rend l'ensemble possible : toutes les pieces sorties du meme repere CAO, sans
recentrage, Z pour axe, en centimetres.  La detection d'axe et le recentrage
disparaissent donc en mode composants -- non par simplification, mais parce
qu'ils sont remplaces par une garantie plus forte.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import config
from .confidence import HIGH, LOW, MEDIUM, ConfidenceMap
from .geometry.inclusion import SolidTester
from .io import loader
from .mesh import TriMesh, rotation_matrix
from .numeric import jacobi_eigen

Vec3 = tuple[float, float, float]

#: Emplacements d'import, et leur statut (SPEC v2 3).
SLOT_SHELL = "coque"
SLOT_HUB = "moyeu"
SLOT_INLET = "entree_fluide"
SLOT_OUTLET = "sortie_fluide"
SLOT_BLADE = "pale"
SLOTS = (SLOT_SHELL, SLOT_HUB, SLOT_INLET, SLOT_OUTLET, SLOT_BLADE)
REQUIRED_SLOTS = (SLOT_INLET, SLOT_OUTLET, SLOT_BLADE)

#: Topologie de pale declaree (SPEC v2 2).
BLADE_CONVENTIONAL = "conventionnelle"
BLADE_TOROIDAL = "toroidale"
BLADE_TOPOLOGIES = (BLADE_CONVENTIONAL, BLADE_TOROIDAL)

#: Le contrat d'export, a afficher au-dessus des emplacements.
EXPORT_CONTRACT = (
    "Toutes les pieces doivent etre exportees depuis le meme repere CAO, sans recentrage, "
    "avec Z pour axe de rotation, en centimetres. C'est ce contrat qui remplace la detection "
    "d'axe et le recentrage : si un exportateur recentre une piece sur l'origine, l'assemblage "
    "est detruit sans que rien n'y paraisse. L'outil le detecte et refuse de continuer."
)


@dataclass
class Declarations:
    """Ce que l'utilisateur affirme, et que l'outil ne cherchera pas a inferer."""

    mode: str = ""  # pompe_carenee ou helice_libre
    blade_topology: str = BLADE_CONVENTIONAL
    n_blades: int = 0

    def check(self) -> None:
        """Refuse une declaration hors domaine, avec la liste des valeurs admises."""
        from .analysis import MACHINE_PROPELLER, MACHINE_PUMP

        if self.mode not in (MACHINE_PUMP, MACHINE_PROPELLER):
            raise ValueError(
                f"mode inconnu : '{self.mode}'. Attendu : {MACHINE_PUMP} ou {MACHINE_PROPELLER}."
            )
        if self.blade_topology not in BLADE_TOPOLOGIES:
            raise ValueError(
                f"topologie de pale inconnue : '{self.blade_topology}'. Attendu : "
                f"{', '.join(BLADE_TOPOLOGIES)}."
            )
        if not config.DECLARED_BLADES_MIN <= self.n_blades <= config.DECLARED_BLADES_MAX:
            raise ValueError(
                f"nombre de pales declare hors domaine : {self.n_blades}. Attendu entre "
                f"{config.DECLARED_BLADES_MIN} et {config.DECLARED_BLADES_MAX}."
            )

    def to_dict(self) -> dict:
        """Vue serialisable."""
        return {
            "mode": self.mode,
            "topologie_pale": self.blade_topology,
            "nombre_de_pales": self.n_blades,
        }


@dataclass
class Component:
    """Une piece importee, gardee dans le repere CAO ou elle a ete exportee."""

    slot: str = ""
    path: str = ""
    mesh: TriMesh | None = None
    watertight: bool = True
    boundary_edges: int = 0
    volume: float = 0.0  # m3
    area: float = 0.0  # m2
    centroid: Vec3 = (0.0, 0.0, 0.0)
    low: Vec3 = (0.0, 0.0, 0.0)
    high: Vec3 = (0.0, 0.0, 0.0)

    def radial_extent(self, axis: Vec3 = (0.0, 0.0, 0.0)) -> tuple[float, float]:
        """Rayons minimal et maximal autour de l'axe, en m.

        Le contrat d'export fixe la **direction** de l'axe -- Z -- mais pas sa
        position : une piece exportee sans recentrage se trouve la ou la CAO
        l'avait mise, souvent a plusieurs metres de l'origine. Mesurer les
        rayons depuis l'origine y donnerait des dizaines de metres.
        """
        if self.mesh is None or not self.mesh.vertices:
            return 0.0, 0.0
        radii = [math.hypot(v[0] - axis[0], v[1] - axis[1]) for v in self.mesh.vertices]
        return min(radii), max(radii)

    def to_dict(self, axis: Vec3 = (0.0, 0.0, 0.0)) -> dict:
        """Vue serialisable, en SI."""
        inner, outer = self.radial_extent(axis)
        return {
            "emplacement": self.slot,
            "fichier": self.path,
            "triangles": len(self.mesh.faces) if self.mesh else 0,
            "volume_m3": self.volume,
            "surface_m2": self.area,
            "centroide_m": list(self.centroid),
            "boite_min_m": list(self.low),
            "boite_max_m": list(self.high),
            "rayon_interieur_m": inner,
            "rayon_exterieur_m": outer,
            "etanche": self.watertight,
            "aretes_de_bord": self.boundary_edges,
        }


@dataclass
class FluidPlane:
    """Section de reference, lue sur une tranche mince (SPEC v2 3.2)."""

    slot: str = ""
    thickness: float = 0.0  # m
    area: float = 0.0  # m2, volume / epaisseur -- mesuree, non deduite de rayons
    normal: Vec3 = (0.0, 0.0, 1.0)
    centroid: Vec3 = (0.0, 0.0, 0.0)
    anisotropy: float = 0.0  # rapport des valeurs propres, dit si la normale est isolee
    radial: bool = False  # bande cylindrique : section traversee radialement (refoulement centrifuge)
    radius: float = 0.0  # m, rayon moyen de la bande quand `radial`
    height: float = 0.0  # m, hauteur axiale de la bande quand `radial`
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "emplacement": self.slot,
            "epaisseur_m": self.thickness,
            "aire_m2": self.area,
            "normale": "radiale" if self.radial else list(self.normal),
            "centroide_m": list(self.centroid),
            "anisotropie": None if self.radial else self.anisotropy,
            "bande_cylindrique": self.radial,
            "rayon_de_bande_m": self.radius or None,
            "hauteur_de_bande_m": self.height or None,
        }


@dataclass
class Check:
    """Un controle du §6 : ce qu'il verifie, et ce qu'il a trouve."""

    name: str = ""
    passed: bool = True
    blocking: bool = False
    detail: str = ""

    def to_dict(self) -> dict:
        """Vue serialisable."""
        return {
            "controle": self.name,
            "conforme": self.passed,
            "bloquant": self.blocking,
            "detail": self.detail,
        }


@dataclass
class ComponentAssembly:
    """Le resultat de l'import par composants : pieces, plans, et controles."""

    declarations: Declarations = field(default_factory=Declarations)
    components: dict = field(default_factory=dict)
    inlet: FluidPlane | None = None
    outlet: FluidPlane | None = None
    flow_direction: Vec3 = (0.0, 0.0, -1.0)
    axis_origin: Vec3 = (0.0, 0.0, 0.0)  # ou passe l'axe de rotation, dans le repere CAO
    checks: list[Check] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    confidence: ConfidenceMap = field(default_factory=ConfidenceMap)

    def blocked(self) -> Check | None:
        """Le premier controle bloquant en echec, s'il y en a un."""
        return next((c for c in self.checks if c.blocking and not c.passed), None)

    def component(self, slot: str) -> Component | None:
        """La piece d'un emplacement, ou `None` s'il est vide."""
        return self.components.get(slot)

    def to_dict(self) -> dict:
        """Vue serialisable, en SI."""
        return {
            "declarations": self.declarations.to_dict(),
            "contrat_d_export": EXPORT_CONTRACT,
            "pieces": {
                slot: c.to_dict(self.axis_origin) for slot, c in self.components.items()
            },
            "plan_d_entree": self.inlet.to_dict() if self.inlet else None,
            "plan_de_sortie": self.outlet.to_dict() if self.outlet else None,
            "sens_debitant": list(self.flow_direction),
            "origine_de_l_axe_m": list(self.axis_origin),
            "controles": [c.to_dict() for c in self.checks],
            "confiance": dict(self.confidence),
            "avertissements": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# Chargement : aucun alignement, aucun recentrage -- c'est le contrat
# ---------------------------------------------------------------------------
def load_component(slot: str, path: str, unit: str | float | None = "cm") -> Component:
    """Charge une piece **telle qu'elle a ete exportee**.

    Ni detection d'axe, ni recentrage : le contrat d'export garantit que Z est
    l'axe et que l'origine est celle de la CAO.  Les recaler ici detruirait
    l'assemblage, puisque c'est leur position relative qui porte l'information.
    """
    try:
        mesh, report = loader.load_mesh(path, unit=unit)
    except loader.ImportError_ as error:
        # Cinq emplacements, cinq fichiers : l'erreur doit dire lequel.
        raise loader.ImportError_(f"emplacement « {slot} » : {error}") from None
    low, high = mesh.bounds()
    return Component(
        slot=slot,
        path=path,
        mesh=mesh,
        watertight=report.watertight,
        boundary_edges=len(mesh.boundary_edges()),
        volume=abs(mesh.volume()),
        area=mesh.area(),
        centroid=mesh.centroid(),
        low=low,
        high=high,
    )


def fluid_plane(component: Component) -> FluidPlane:
    """Aire, normale et position d'une tranche mince (SPEC v2 3.2).

    L'aire est **mesuree** -- `volume / epaisseur` -- et non deduite de rayons,
    ce qui est tout l'interet : elle cesse d'etre suspendue a la detection de
    `r_1s` et `r_1h`, dont la fragilite est la cause premiere des ecarts de
    debit.

    La normale est le vecteur propre **isole** du tenseur d'inertie.  Une
    tranche mince a deux grandes valeurs propres proches, dans son plan, et une
    petite, suivant son epaisseur : c'est cette derniere qui designe la normale.
    Si les trois valeurs se ressemblent, la piece n'est pas une tranche et la
    normale n'a pas de sens -- l'outil le dit plutot que d'en rendre une.
    """
    plane = FluidPlane(slot=component.slot, centroid=component.centroid)
    mesh = component.mesh
    if mesh is None or component.volume <= 0.0:
        plane.warnings.append(
            f"le solide « {component.slot} » n'a pas de volume exploitable : ni aire ni "
            "normale n'en sont tirees."
        )
        return plane

    # Refoulement radial : sur une roue centrifuge, l'eau sort a travers une
    # bande cylindrique, pas a travers un disque. Une tranche plate se lit par
    # son epaisseur la plus faible ; une bande, par son epaisseur **radiale**,
    # plus faible que sa hauteur. La lire comme une tranche prenait sa hauteur
    # pour epaisseur : 6 cm2 au lieu de 71 sur la roue d'essai, et une normale
    # axiale pour un ecoulement radial.
    # Le rayon interieur se lit sur une coupe a mi-hauteur, pas sur les
    # sommets : un disque plein dont les faces sont triangulees depuis le bord
    # n'a aucun sommet au centre, et ses sommets le faisaient passer pour une
    # bande d'epaisseur nulle -- une aire de dix millions de cm2.
    centre = (component.centroid[0], component.centroid[1], 0.0)
    height = component.high[2] - component.low[2]
    cut = section(mesh, 0.5 * (component.low[2] + component.high[2]))
    radii = [math.hypot(v[0] - centre[0], v[1] - centre[1]) for v in mesh.vertices]
    r_out = max(radii) if radii else 0.0
    ring = bool(cut) and not axis_inside(cut, centre)
    # Rayons interieur et exterieur lus tous deux aux sommets de la coupe, pour
    # qu'une facette ne compte pas d'un cote par son plat et de l'autre par son
    # coin : sur une paroi d'un millimetre, l'ecart faisait 20 % d'aire.
    r_in = min(math.hypot(x - centre[0], y - centre[1]) for seg in cut for x, y in seg) if ring else 0.0
    wall = r_out - r_in
    if ring and r_out > 0.0 and r_in > config.BAND_INNER_MIN * r_out and wall < height:
        # Aire de passage : moyenne des surfaces laterales interieure et
        # exterieure, lues sur les facettes dont la normale est radiale. Diviser
        # le volume par une epaisseur de paroi lue sur les sommets heritait de la
        # facettisation -- plat d'un cote, coin de l'autre --, soit 20 % d'erreur
        # sur une paroi d'un millimetre.
        lateral = sum(
            area for area, index in zip(mesh.face_areas(), range(len(mesh.faces)))
            if abs(mesh.face_normal(index)[2]) < config.BAND_LATERAL_NZ
        )
        plane.radial = True
        plane.height = height
        plane.area = 0.5 * lateral
        plane.radius = plane.area / (2.0 * math.pi * height) if height > 0.0 else 0.0
        plane.thickness = component.volume / plane.area if plane.area > 0.0 else 0.0
        wall = plane.thickness
        plane.normal = (0.0, 0.0, 0.0)  # radiale, differente en chaque point de la bande
        if wall > config.SLICE_THICKNESS_MAX:
            plane.warnings.append(
                f"le solide « {component.slot} » est une bande de "
                f"{wall * config.MM_PER_M:.0f} mm d'epaisseur radiale, au-dela des "
                f"{config.SLICE_THICKNESS_MAX * config.MM_PER_M:.0f} mm attendus : l'aire publiee "
                "est une moyenne sur cette epaisseur."
            )
        return plane

    values, vectors = jacobi_eigen(mesh.inertia_tensor(about=component.centroid))
    order = sorted(range(3), key=lambda i: values[i])
    smallest, middle, largest = (values[i] for i in order)
    # Pour une tranche, le moment le plus **grand** est celui autour de la
    # normale : la matiere y est la plus eloignee de l'axe. Les deux autres,
    # dans le plan, sont proches entre eux et plus petits.
    plane.normal = _normalise(vectors[order[2]])
    plane.anisotropy = largest / middle if middle > 0.0 else math.inf

    low, high = component.low, component.high
    extents = [high[k] - low[k] for k in range(3)]
    plane.thickness = min(extents)
    plane.area = component.volume / plane.thickness if plane.thickness > 0.0 else 0.0

    if plane.thickness > config.SLICE_THICKNESS_MAX:
        plane.warnings.append(
            f"le solide « {component.slot} » est epais de "
            f"{plane.thickness * config.MM_PER_M:.0f} mm, au-dela des "
            f"{config.SLICE_THICKNESS_MAX * config.MM_PER_M:.0f} mm attendus d'une tranche "
            "mince. Un tube long ne dit pas ou se trouve la section de reference : l'aire "
            "publiee est une moyenne sur sa longueur."
        )
    if plane.anisotropy < config.SLICE_ANISOTROPY_MIN:
        plane.warnings.append(
            f"le solide « {component.slot} » n'a pas de direction privilegiee nette "
            f"(anisotropie {plane.anisotropy:.2f}) : sa normale est incertaine. Modelisez-le "
            "comme une tranche plate placee dans le plan de reference."
        )
    return plane


#: Direction de la demi-droite du test de parite. Quelconque a dessein : les
#: solides de revolution d'AutoCAD ont leur couture dans le plan y = 0, et une
#: demi-droite le long de +X la suivait exactement -- chaque traversee y etait
#: comptee deux fois ou pas du tout, selon l'arrondi.
_RAY = (math.cos(0.7303), math.sin(0.7303))


def section(mesh, z: float) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Segments de la coupe du maillage par le plan horizontal `z`."""
    segments = []
    for a, b, c in mesh.triangles():
        cut = []
        for p, q in ((a, b), (b, c), (c, a)):
            dp, dq = p[2] - z, q[2] - z
            if (dp > 0.0) != (dq > 0.0):
                t = dp / (dp - dq)
                cut.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
        if len(cut) == 2:
            segments.append((cut[0], cut[1]))
    return segments


def axis_inside(segments, axis: Vec3) -> bool:
    """Vrai si l'axe traverse la matiere dans cette coupe (parite des traversees)."""
    dx, dy = _RAY
    crossings = 0
    for (x1, y1), (x2, y2) in segments:
        ex, ey = x2 - x1, y2 - y1
        det = ex * dy - ey * dx
        if det == 0.0:
            continue
        wx, wy = axis[0] - x1, axis[1] - y1
        s = (wx * dy - wy * dx) / det  # position sur le segment
        t = (wx * ey - wy * ex) / det  # position sur la demi-droite
        if 0.0 <= s < 1.0 and t > 0.0:
            crossings += 1
    return crossings % 2 == 1


def distance_to_axis(segments, axis: Vec3) -> float:
    """Plus courte distance de l'axe au contour de la coupe (pas a ses sommets)."""
    best = math.inf
    for (x1, y1), (x2, y2) in segments:
        ex, ey = x2 - x1, y2 - y1
        length = ex * ex + ey * ey
        t = 0.0 if length == 0.0 else max(0.0, min(1.0, ((axis[0] - x1) * ex + (axis[1] - y1) * ey) / length))
        best = min(best, math.hypot(x1 + t * ex - axis[0], y1 + t * ey - axis[1]))
    return best


def inner_radius(component: Component, axis: Vec3) -> float:
    """Rayon interieur d'une piece autour de l'axe, lu sur une coupe a mi-hauteur.

    Zero si l'axe traverse la matiere. Lu sur les sommets, il valait le rayon
    exterieur d'un disque plein dont les faces sont triangulees depuis le bord,
    et le disque d'entree passait pour une couronne percee jusqu'a son bord.
    """
    if component.mesh is None:
        return 0.0
    cut = section(component.mesh, 0.5 * (component.low[2] + component.high[2]))
    if not cut or axis_inside(cut, axis):
        return 0.0
    return min(math.hypot(x - axis[0], y - axis[1]) for seg in cut for x, y in seg)


def hub_radius_at(component: Component, z: float, axis: Vec3) -> float | None:
    """Rayon de la matiere pleine qui entoure l'axe au plan `z`, ou `None`.

    Le moyeu n'obstrue l'entree que la ou il est : au plan d'entree. Son rayon
    **global** n'y dit rien -- un « corps » exporte d'un seul tenant, moyeu et
    flasque ensemble, a pour rayon exterieur celui du flasque, et le prendre pour
    r1h rendait un moyeu de 96 mm autour d'un oeillard de 36. On coupe donc la
    piece par le plan : si l'axe y est dans la matiere, le rayon du moyeu est la
    distance de l'axe au contour le plus proche ; sinon, rien n'obstrue le
    centre a ce plan.
    """
    if component.mesh is None:
        return None
    segments = section(component.mesh, z)
    if not segments or not axis_inside(segments, axis):
        return None
    return distance_to_axis(segments, axis)


def _normalise(vector) -> Vec3:
    """Vecteur unitaire ; renvoie +Z si le vecteur est nul."""
    norm = math.sqrt(sum(c * c for c in vector))
    if norm <= 0.0:
        return (0.0, 0.0, 1.0)
    return tuple(c / norm for c in vector)  # type: ignore[return-value]


def axis_origin(inlet: FluidPlane, outlet: FluidPlane) -> Vec3:
    """Ou passe l'axe de rotation, dans le repere CAO.

    Le contrat d'export fixe la **direction** de l'axe -- Z -- et rien de plus :
    une piece exportee sans recentrage se trouve la ou la CAO l'avait mise, a
    plusieurs metres de l'origine dans le cas d'essai.  Les deux solides fluide
    sont des couronnes centrees sur l'axe : la moyenne de leurs centroides le
    situe, et le situe **par la mesure** plutot que par une convention.
    """
    return (
        0.5 * (inlet.centroid[0] + outlet.centroid[0]),
        0.5 * (inlet.centroid[1] + outlet.centroid[1]),
        0.0,
    )


def flow_direction(inlet: FluidPlane, outlet: FluidPlane) -> Vec3:
    """Sens debitant : du centroide d'entree vers celui de sortie.

    C'est ce vecteur qui remplace la convention « refoulement vers -Z ».  Il est
    **mesure** sur l'assemblage, ce qui leve l'indetermination du sens de
    rotation : le critere `signe(omega) = signe(dz/dtheta)` se projette dessus.
    """
    delta = tuple(outlet.centroid[k] - inlet.centroid[k] for k in range(3))
    return _normalise(delta)


# ---------------------------------------------------------------------------
# §6 : l'outil verifie, il n'accepte pas
# ---------------------------------------------------------------------------
def check_distinct_planes(assembly: ComponentAssembly) -> Check:
    """Entree et sortie distinctes : sans elles, le sens debitant n'existe pas.

    Le meme fichier donne pour les deux emplacements -- ou deux tranches
    superposees -- rend un vecteur nul, que la normalisation changeait sans
    rien dire en +Z. Le sens de rotation en etait ensuite deduit. Le controle
    est **bloquant** : tout ce que le mode apporte en depend.
    """
    check = Check(name="plans distincts", blocking=True)
    inlet, outlet = assembly.inlet, assembly.outlet
    gap = math.dist(inlet.centroid, outlet.centroid)
    if gap < config.PLANES_DISTINCT_MIN:
        check.passed = False
        check.detail = (
            f"les solides d'entree et de sortie sont confondus (centroides a "
            f"{gap * config.MM_PER_M:.2f} mm l'un de l'autre) : le sens debitant n'est pas "
            "defini, ni rien de ce qui en decoule. Verifiez que les deux emplacements ne "
            "designent pas le meme fichier."
        )
    else:
        check.detail = f"entree et sortie distantes de {gap * config.MM_PER_M:.1f} mm"
    return check


def check_common_frame(assembly: ComponentAssembly) -> Check:
    """1. Repere commun : une piece recentree a l'export detruit l'assemblage.

    C'est l'un des deux controles **bloquants**, avec les plans distincts, et
    pour une raison precise : sa violation ne se voit sur aucune grandeur publiee.  Les pieces se lisent
    toutes correctement, chacune dans son coin, et seule leur position relative
    -- donc tout ce que le mode composants apporte -- est fausse.
    """
    pieces = list(assembly.components.values())
    if len(pieces) < 2:
        return Check("repere commun", True, True, "moins de deux pieces : rien a confronter")

    axis = assembly.axis_origin
    # L'echelle est le rayon de la piece, mesure autour de **l'axe** ; la
    # distance testee est celle a **l'origine CAO**. Les confondre ferait
    # signaler comme recentree toute piece axisymetrique, qui est sur l'axe par
    # construction -- moyeu, coque et solides fluide en tete.
    rayon = max(max(c.radial_extent(axis)) for c in pieces) or 1.0
    distances = {
        c.slot: math.sqrt(sum(v * v for v in c.centroid)) for c in pieces
    }
    seuil = config.RECENTRE_TOLERANCE * rayon
    colles = [slot for slot, d in distances.items() if d < seuil]
    loin = [slot for slot, d in distances.items() if d >= seuil]
    if colles and loin:
        return Check(
            "repere commun", False, True,
            f"le centroide de « {', '.join(colles)} » est a moins de "
            f"{config.RECENTRE_TOLERANCE:.0%} du rayon exterieur de l'origine, alors que "
            f"« {', '.join(loin)} » en est loin. Cette piece a probablement ete recentree a "
            "l'export : l'assemblage serait faux sans que rien n'y paraisse. Reexportez-la "
            "depuis le meme repere CAO, sans recentrage."
        )

    # Pieces deplacees chacune de son cote. `STLOUT` d'AutoCAD n'exporte que
    # des objets situes dans l'octant positif ; on les y pousse volontiers une
    # par une, et chaque piece arrive avec le coin de sa boite a l'origine.
    # Aucune n'est recentree sur l'origine -- le controle precedent passe --
    # mais l'assemblage est detruit. Deux signatures le trahissent.
    corner = config.CORNER_TOLERANCE
    paires = []
    for index, first in enumerate(pieces):
        for second in pieces[index + 1:]:
            same_corner = (
                abs(first.low[0] - second.low[0]) < corner
                and abs(first.low[1] - second.low[1]) < corner
            )
            sizes = [
                abs((first.high[k] - first.low[k]) - (second.high[k] - second.low[k]))
                for k in (0, 1)
            ]
            if same_corner and max(sizes) > 10.0 * corner:
                paires.append(f"{first.slot} / {second.slot}")
    # Les pieces de revolution -- solides fluide et moyeu -- sont centrees sur
    # l'axe par construction : leurs centres doivent coincider dans le plan XY.
    # La coque en est exclue, une volute n'etant pas de revolution.
    revolution = [
        c for c in pieces if c.slot in (SLOT_INLET, SLOT_OUTLET, SLOT_HUB)
    ]
    ecart, pire = 0.0, ""
    for index, first in enumerate(revolution):
        for second in revolution[index + 1:]:
            d = math.hypot(
                first.centroid[0] - second.centroid[0], first.centroid[1] - second.centroid[1]
            )
            if d > ecart:
                ecart, pire = d, f"{first.slot} / {second.slot}"
    decentre = ecart > config.COAXIAL_TOLERANCE * rayon
    if paires or decentre:
        causes = []
        if decentre:
            causes.append(
                f"les pieces de revolution ne sont pas sur le meme axe ({pire} : centres a "
                f"{ecart * config.MM_PER_M:.1f} mm l'un de l'autre dans le plan XY)"
            )
        if paires:
            causes.append(
                f"des pieces de tailles differentes ont exactement le meme coin de boite "
                f"englobante ({'; '.join(paires)}), signature de pieces deplacees chacune de "
                "son cote"
            )
        return Check(
            "repere commun", False, True,
            "; ".join(causes) + ". L'assemblage est detruit : la position de la pale par "
            "rapport a l'axe, dont dependent angles, sens de rotation et sections, n'est plus "
            "connue. Sous AutoCAD, STLOUT exige l'octant positif : selectionnez toutes les "
            "pieces ensemble, deplacez-les en une seule fois du meme vecteur, puis exportez-les "
            "une a une sans plus rien deplacer."
        )
    return Check("repere commun", True, True,
                 "toutes les pieces partagent le meme repere et le meme axe")


def _inside_fraction(outer: Component, inner: Component) -> float:
    """Part de la surface de `inner` situee dans le volume de `outer`."""
    if outer.mesh is None or inner.mesh is None:
        return 0.0
    tester = SolidTester(outer.mesh)
    return tester.fraction_inside(inner.mesh.sample_surface(config.INTERSECTION_SAMPLES))


def check_no_interpenetration(assembly: ComponentAssembly) -> Check:
    """2. Non-interpenetration de moyeu, pale et coque, **en volume**.

    Le test portait sur les boites englobantes. Sur une roue fermee, celle du
    corps -- moyeu et flasque -- contient toutes les pales par construction :
    le controle echouait sur toute roue fermee, et un avertissement qui sonne
    toujours n'avertit plus de rien. On mesure maintenant la part de la surface
    de chaque piece qui se trouve dans le volume de l'autre.
    """
    solides = [assembly.component(s) for s in (SLOT_HUB, SLOT_BLADE, SLOT_SHELL)]
    solides = [c for c in solides if c is not None and c.mesh is not None]
    pires, mesures = [], []
    for index, first in enumerate(solides):
        for second in solides[index + 1:]:
            part = max(_inside_fraction(first, second), _inside_fraction(second, first))
            mesures.append(f"{first.slot} / {second.slot} : {part:.1%}")
            if part > config.OVERLAP_TOLERANCE:
                pires.append(f"{first.slot} / {second.slot} : {part:.0%}")
    if pires:
        return Check(
            "non-interpenetration", False, False,
            f"des pieces entrent l'une dans l'autre au-dela de {config.OVERLAP_TOLERANCE:.0%} de "
            f"leur surface ({'; '.join(pires)}). Une pale soudee a son moyeu y plonge un peu ; "
            "au-dela, verifiez que les pieces sont bien celles du meme assemblage, et placees."
        )
    return Check(
        "non-interpenetration", True, False,
        "aucune piece n'entre dans une autre au-dela de la tolerance"
        + (f" ({'; '.join(mesures)})" if mesures else "")
    )


def check_blade_between_planes(assembly: ComponentAssembly) -> Check:
    """3. La pale se trouve entre le plan d'entree et celui de sortie.

    Sinon, les deux solides fluide sont probablement inversés -- et le sens
    debitant avec eux, donc le sens de rotation qui s'y ancre.
    """
    blade = assembly.component(SLOT_BLADE)
    if blade is None or assembly.inlet is None or assembly.outlet is None:
        return Check("position de la pale", True, False, "pieces manquantes : non verifie")

    direction = assembly.flow_direction
    def projection(point) -> float:
        return sum(point[k] * direction[k] for k in range(3))

    entree, sortie = projection(assembly.inlet.centroid), projection(assembly.outlet.centroid)
    pale = projection(blade.centroid)
    if entree <= pale <= sortie:
        return Check("position de la pale", True, False,
                     "le centroide de la pale est bien entre les deux plans fluide")
    return Check(
        "position de la pale", False, False,
        f"le centroide de la pale se projette a {(pale - entree) * config.MM_PER_M:+.0f} mm de "
        f"l'entree le long du sens debitant, hors de l'intervalle "
        f"[0, {(sortie - entree) * config.MM_PER_M:.0f}] mm. Les solides d'entree et de sortie "
        "sont probablement inverses : le sens debitant, et donc le sens de rotation qui s'y "
        "ancre, seraient pris a l'envers."
    )


def check_blade_count(assembly: ComponentAssembly) -> Check:
    """4. Reconstruire l'assemblage par N rotations, et verifier qu'il tient.

    Deux copies qui se recoupent signifient que le nombre declare est trop
    grand pour cette pale. Le test portait sur l'etendue **angulaire** de la
    pale -- plus de 2 pi / N, et les copies etaient declarees en conflit. Une
    aube en boucle ou tres enroulee fait presque le tour de l'axe tout en
    laissant la place a ses voisines, a un autre rayon ou a une autre hauteur :
    sur la roue d'essai, 356 degres pour cinq pales qui existent bel et bien.
    Le test porte maintenant sur les volumes : la surface de chaque copie
    tournee ne doit pas entrer dans la pale d'origine.
    """
    blade = assembly.component(SLOT_BLADE)
    n_blades = assembly.declarations.n_blades
    if blade is None or blade.mesh is None or n_blades <= 0:
        return Check("nombre de pales", True, False, "pale absente : non verifie")

    axis = assembly.axis_origin
    tester = SolidTester(blade.mesh)
    samples = blade.mesh.sample_surface(config.INTERSECTION_SAMPLES)
    pire, pire_k = 0.0, 0
    for k in range(1, n_blades // 2 + 1):
        matrix = rotation_matrix((0.0, 0.0, 1.0), 2.0 * math.pi * k / n_blades)
        tournes = []
        for p in samples:
            local = (p[0] - axis[0], p[1] - axis[1], p[2])
            x = matrix[0][0] * local[0] + matrix[0][1] * local[1]
            y = matrix[1][0] * local[0] + matrix[1][1] * local[1]
            tournes.append((x + axis[0], y + axis[1], p[2]))
        part = tester.fraction_inside(tournes)
        if part > pire:
            pire, pire_k = part, k
    angle = 360.0 * pire_k / n_blades
    if pire > config.COPY_OVERLAP_TOLERANCE:
        return Check(
            "nombre de pales", False, False,
            f"reconstruite par {n_blades} rotations, la pale se recoupe : {pire:.0%} de la "
            f"surface de sa copie tournee de {angle:.0f} degres est dans son volume. Le nombre "
            "declare est trop grand pour cette pale, ou le fichier importe contient deja "
            "plusieurs pales."
        )
    return Check(
        "nombre de pales", True, False,
        f"les {n_blades} copies de la pale tiennent : au plus {pire:.1%} de la surface d'une "
        f"copie entre dans sa voisine (tolerance {config.COPY_OVERLAP_TOLERANCE:.0%}, pour les "
        "pales soudees entre elles). Ce controle ecarte une erreur grossiere ; il ne departage "
        "pas N et N+1."
    )


def genus(mesh: TriMesh) -> tuple[int, int, int]:
    """Genre topologique du maillage, par la caracteristique d'Euler.

    Renvoie `(genre, composantes de bord, caracteristique)`.  La relation
    generale d'une surface orientable est `chi = 2 c - 2 g - b`, ou `c` est le
    nombre de composantes connexes et `b` celui des composantes de bord.  La
    version simplifiee `2 - 2g - b`, qui suppose une seule composante, rend un
    genre **negatif** des que le maillage en compte plusieurs -- une pale
    exportee en deux nappes non cousues, par exemple -- ce qui n'existe pas et
    trahirait le calcul plutot que la piece.

    Une boule est de genre 0, un tore de genre 1 ; une aube en boucle est
    topologiquement un tore, et c'est ce qui la distingue mecaniquement.
    """
    vertices = len(mesh.vertices)
    faces = len(mesh.faces)
    edges = len(mesh.edge_map())
    boundary = _boundary_loops(mesh)
    euler = vertices - edges + faces
    parts = _connected_parts(mesh)
    return max(0, (2 * parts - euler - boundary) // 2), boundary, euler


def _connected_parts(mesh: TriMesh) -> int:
    """Nombre de composantes connexes du maillage, par les aretes partagees."""
    voisins: dict[int, set] = {}
    for a, b, c in mesh.faces:
        for x, y in ((a, b), (b, c), (c, a)):
            voisins.setdefault(x, set()).add(y)
            voisins.setdefault(y, set()).add(x)
    vus: set = set()
    parts = 0
    for depart in voisins:
        if depart in vus:
            continue
        parts += 1
        pile = [depart]
        vus.add(depart)
        while pile:
            for suivant in voisins[pile.pop()]:
                if suivant not in vus:
                    vus.add(suivant)
                    pile.append(suivant)
    return max(1, parts)


def _boundary_loops(mesh: TriMesh) -> int:
    """Nombre de composantes connexes du bord (chaines fermees d'aretes libres)."""
    edges = mesh.boundary_edges()
    if not edges:
        return 0
    voisins: dict[int, set] = {}
    for a, b in edges:
        voisins.setdefault(a, set()).add(b)
        voisins.setdefault(b, set()).add(a)
    vus: set = set()
    boucles = 0
    for depart in voisins:
        if depart in vus:
            continue
        boucles += 1
        pile = [depart]
        vus.add(depart)
        while pile:
            courant = pile.pop()
            for suivant in voisins[courant]:
                if suivant not in vus:
                    vus.add(suivant)
                    pile.append(suivant)
    return boucles


def check_declared_topology(assembly: ComponentAssembly) -> Check:
    """5. Une pale declaree toroidale doit etre de genre 1.

    Une boucle fermee a une anse : c'est un tore, pas une boule.  Si la pale
    importee est simplement connexe, la declaration contredit la geometrie -- et
    c'est la declaration qui pilote la logique de coupe, d'ou l'avertissement.
    """
    blade = assembly.component(SLOT_BLADE)
    declaree = assembly.declarations.blade_topology
    if blade is None or blade.mesh is None:
        return Check("topologie declaree", True, False, "pale absente : non verifie")

    g, boundary, euler = genus(blade.mesh)
    detail = f"genre {g}, {boundary} composante(s) de bord, caracteristique d'Euler {euler}"
    if declaree == BLADE_TOROIDAL and g < 1:
        return Check(
            "topologie declaree", False, False,
            f"la pale est declaree toroidale, mais le maillage est de {detail} : une boucle "
            "fermee serait de genre 1 au moins. La declaration contredit la geometrie, et "
            "c'est elle qui pilote la logique de coupe."
        )
    if declaree == BLADE_CONVENTIONAL and g >= 1:
        return Check(
            "topologie declaree", False, False,
            f"la pale est declaree conventionnelle, mais le maillage est de {detail} : il "
            "porte une anse, signature d'une boucle fermee. Une coupe la traversera deux fois, "
            "et la cambrure lue sur un seul profil sera fausse."
        )
    return Check("topologie declaree", True, False, f"conforme a la declaration ({detail})")


def check_watertight(assembly: ComponentAssembly) -> Check:
    """6. Etancheite : inchange, mais piece par piece."""
    troues = [
        f"{c.slot} ({c.boundary_edges} aretes de bord)"
        for c in assembly.components.values() if not c.watertight
    ]
    if troues:
        return Check(
            "etancheite", False, False,
            f"maillage non etanche : {', '.join(troues)}. Les grandeurs volumiques de ces "
            "pieces sont plafonnees a la confiance moyenne."
        )
    return Check("etancheite", True, False, "toutes les pieces sont fermees")


# ---------------------------------------------------------------------------
# L'assemblage
# ---------------------------------------------------------------------------
def assemble(
    paths: dict, declarations: Declarations, unit: str | float | None = "cm"
) -> ComponentAssembly:
    """Charge les pieces declarees, en tire les plans fluide, et les controle.

    Ne calcule rien d'hydraulique : c'est l'etape 1 de l'ordre demande, et elle
    s'arrete la ou commence le modele.  Un emplacement vide n'est pas une erreur
    -- l'outil signale les grandeurs qu'il ne pourra pas produire.
    """
    declarations.check()
    assembly = ComponentAssembly(declarations=declarations)

    manquants = [slot for slot in REQUIRED_SLOTS if not paths.get(slot)]
    if manquants:
        raise ValueError(
            f"emplacement(s) obligatoire(s) vide(s) : {', '.join(manquants)}. "
            f"L'entree fluide, la sortie fluide et la pale sont necessaires ; la coque et le "
            "moyeu sont facultatifs."
        )

    for slot in SLOTS:
        path = paths.get(slot)
        if path:
            assembly.components[slot] = load_component(slot, path, unit=unit)

    for slot in (SLOT_SHELL, SLOT_HUB):
        if slot not in assembly.components:
            assembly.warnings.append(
                f"emplacement « {slot} » vide : "
                + ("rayon de carter et jeu en bout de pale ne seront pas produits."
                   if slot == SLOT_SHELL else
                   "r1h et la longueur de moyeu seront lus sur la pale et les plans fluide, "
                   "avec la confiance moindre que cela suppose.")
            )

    assembly.inlet = fluid_plane(assembly.components[SLOT_INLET])
    assembly.outlet = fluid_plane(assembly.components[SLOT_OUTLET])
    assembly.warnings.extend(assembly.inlet.warnings)
    assembly.warnings.extend(assembly.outlet.warnings)
    assembly.flow_direction = flow_direction(assembly.inlet, assembly.outlet)
    assembly.axis_origin = axis_origin(assembly.inlet, assembly.outlet)

    assembly.checks = [
        check_distinct_planes(assembly),
        check_common_frame(assembly),
        check_no_interpenetration(assembly),
        check_blade_between_planes(assembly),
        check_blade_count(assembly),
        check_declared_topology(assembly),
        check_watertight(assembly),
    ]
    for check in assembly.checks:
        if not check.passed:
            assembly.warnings.append(f"[{check.name}] {check.detail}")

    # La confiance suit les controles : un echec non bloquant n'annule pas la
    # declaration, il abaisse ce qui en depend.
    rate = sum(1 for c in assembly.checks if not c.passed)
    assembly.confidence.set("assemblage", HIGH if rate == 0 else (MEDIUM if rate == 1 else LOW))
    assembly.confidence.set(
        "sections",
        MEDIUM if any(assembly.inlet.warnings + assembly.outlet.warnings) else HIGH,
    )
    volumique = all(c.watertight for c in assembly.components.values())
    assembly.confidence.set("volume", HIGH if volumique else MEDIUM)
    return assembly


# ---------------------------------------------------------------------------
# §5.1 : le sens de rotation, ancre sur le sens debitant
# ---------------------------------------------------------------------------
def helical_slope(mesh: TriMesh, axis: Vec3) -> float:
    """Pente helicoidale moyenne `dz/dtheta` de la surface de pale, en m/rad.

    Une pale est une portion d'helicoide : en la parcourant en azimut, sa
    hauteur derive.  La pente de cette derive est ce que le critere de sens de
    rotation utilise.  Elle est obtenue par une regression de `z` sur `theta`,
    l'azimut etant deroule autour de sa moyenne pour qu'une pale a cheval sur
    l'origine des angles ne casse pas la regression.
    """
    if mesh is None or len(mesh.vertices) < 3:
        return 0.0
    points = [
        (math.atan2(v[1] - axis[1], v[0] - axis[0]), v[2])
        for v in mesh.vertices
    ]
    reference = points[0][0]
    deroules = []
    for theta, z in points:
        ecart = (theta - reference + math.pi) % (2.0 * math.pi) - math.pi
        deroules.append((reference + ecart, z))

    n = float(len(deroules))
    mean_t = sum(t for t, _ in deroules) / n
    mean_z = sum(z for _, z in deroules) / n
    variance = sum((t - mean_t) ** 2 for t, _ in deroules)
    if variance <= 0.0:
        return 0.0
    return sum((t - mean_t) * (z - mean_z) for t, z in deroules) / variance


def sweep_slope(mesh: TriMesh, axis: Vec3) -> float:
    """Pente `dtheta/dr` de la pale, en rad/m : de combien elle recule en s'eloignant de l'axe.

    Regression de l'azimut deroule sur le rayon. Sur une roue centrifuge, les
    aubes sont presque toujours courbees vers l'arriere : leur bout trane
    derriere leur pied, dans le sens oppose a la rotation.
    """
    if mesh is None or len(mesh.vertices) < 3:
        return 0.0
    points = [
        (math.hypot(v[0] - axis[0], v[1] - axis[1]), math.atan2(v[1] - axis[1], v[0] - axis[0]))
        for v in mesh.vertices
    ]
    reference = points[0][1]
    deroules = [(r, reference + (t - reference + math.pi) % (2.0 * math.pi) - math.pi) for r, t in points]
    n = float(len(deroules))
    mean_r = sum(r for r, _ in deroules) / n
    mean_t = sum(t for _, t in deroules) / n
    variance = sum((r - mean_r) ** 2 for r, _ in deroules)
    if variance <= 0.0:
        return 0.0
    return sum((r - mean_r) * (t - mean_t) for r, t in deroules) / variance


def rotation_from_flow(assembly: ComponentAssembly) -> tuple[int, str, str]:
    """Sens de rotation deduit de l'assemblage : `(signe, justification, confiance)`.

    Deux criteres, selon ce que la sortie **mesuree** dit de la machine :

    * refoulement axial -- `signe(omega) = signe(dz/dtheta . flux_z)` : une pale
      de pente `k` tournant a `omega` pousse le fluide a `k omega`, comme un
      filet de vis. Le sens debitant etant mesure et non suppose, le critere
      tranche, en confiance haute ;
    * refoulement radial -- le critere de pente est celui d'une machine axiale
      et ne s'applique pas : il etait pourtant publie en confiance haute. On lit
      le recul des aubes, `signe(omega) = -signe(dtheta/dr)`, qui suppose des
      aubes courbees vers l'arriere -- le cas de presque toutes les pompes --,
      d'ou la confiance moyenne.

    Une aube **en boucle** n'a ni pente ni recul nets -- ses deux brins vont en
    sens opposes, c'est ce qui la distingue -- et une regression faite sur la
    boucle entiere ne mesure rien : le sens reste a declarer.
    """
    blade = assembly.component(SLOT_BLADE)
    if blade is None or blade.mesh is None:
        return 0, "pale absente : le sens de rotation reste a declarer", LOW
    if assembly.declarations.blade_topology == BLADE_TOROIDAL:
        return 0, (
            "aube en boucle : ses deux brins ont des pentes et des reculs opposes. Une "
            "regression menee sur la boucle entiere ne mesure rien -- elle rendait pourtant un "
            "sens en confiance haute. Le sens de rotation reste a declarer par --rotation, "
            "jusqu'a la lecture brin par brin."
        ), LOW

    if assembly.outlet is not None and assembly.outlet.radial:
        recul = sweep_slope(blade.mesh, assembly.axis_origin)
        if abs(recul) < 1e-6:
            return 0, (
                "refoulement radial et aubes sans recul mesurable (aubes droites) : aucun sens "
                "ne s'impose. Declarez-le par --rotation."
            ), LOW
        sign = -1 if recul > 0.0 else 1
        return sign, (
            f"refoulement radial : les aubes reculent de {recul:+.2f} rad/m en s'eloignant de "
            "l'axe. En supposant des aubes courbees vers l'arriere -- le cas de presque toutes "
            "les pompes centrifuges --, la roue tourne dans ce sens. Confiance moyenne : "
            "declarez --rotation si les aubes sont courbees vers l'avant."
        ), MEDIUM

    slope = helical_slope(blade.mesh, assembly.axis_origin)
    flow_z = assembly.flow_direction[2]
    produit = slope * flow_z
    if abs(slope) < 1e-9 or abs(flow_z) < 1e-9:
        return 0, (
            "la pale n'a pas de pente helicoidale mesurable, ou le sens debitant est "
            "perpendiculaire a l'axe : aucun sens de rotation ne s'impose. Declarez-le."
        ), LOW
    sign = 1 if produit > 0.0 else -1
    return sign, (
        f"pente helicoidale de la pale dz/dtheta = {slope * config.MM_PER_M:+.1f} mm/rad, "
        f"sens debitant projete sur l'axe {flow_z:+.2f} : pour pousser le fluide de l'entree "
        "vers la sortie, la roue doit tourner dans ce sens. Le sens debitant etant mesure sur "
        "l'assemblage et non suppose, l'indetermination disparait."
    ), HIGH
