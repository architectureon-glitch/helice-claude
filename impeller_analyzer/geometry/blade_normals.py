"""Angles de pale lus sur les normales de la surface, sans passer par la cambrure.

La phase 4 de la SPEC etablit beta a partir de la ligne de cambrure d'un profil
de pale. C'est la bonne methode -- tant que la coupe donne un profil. Sur une
aube en boucle elle n'en donne pas : selon l'envergure ou l'on coupe, on obtient
le contour de la boucle entiere, aller et retour en un seul tour ferme, ou un
seul brin, ou le bourrelet ou les deux brins fusionnent. Sur la roue toroidale
de reference, beta2 vaut ainsi 8 degres pres du plateau et 89 au milieu de la
veine : le modele n'a pas de reponse stable, et ce n'est pas un defaut de mise
en oeuvre mais la limite de son domaine.

Ce module mesure autrement, localement, sans jamais apparier deux faces. Une
surface de pale ne contient pas sa propre normale : si l'aube fait l'angle beta
avec la direction tangentielle, sa direction dans la surface de courant vaut
`cos(beta) e_u + sin(beta) e_m`, et la normale, qui lui est perpendiculaire,
vaut `-sin(beta) e_u + cos(beta) e_m`. D'ou

    tan(beta) = |N_u| / |N_m|

face par face, pondere par les aires. Les faces d'intrados et d'extrados
portent des normales opposees, mais leurs deux composantes changent de signe
ensemble : le **rapport** garde son signe, qui donne le sens d'enroulement de
l'aube et donc le sens de rotation.

Restent a ecarter les **chants** de l'aube, ces bandes etroites ou elle vient
mourir contre le moyeu et le flasque : elles ne portent aucun angle et tirent la
moyenne vers le bas. Le tri est geometrique et non directionnel -- on ecarte les
faces trop proches des parois de la veine. Filtrer sur la direction de la
normale serait plus simple mais ne marche pas : sur une aube helicoidale la
surface elle-meme porte une grande composante d'envergure, et le meme filtre qui
nettoie une roue centrifuge ordinaire supprime alors toutes les faces.

Portee et limites, mesurees sur des roues synthetiques d'angles imposes
(beta1/beta2 de 15/20 a 40/65 degres) : la lecture est basse de 2 a 5 degres,
d'autant plus que l'angle est grand, et le sens d'enroulement est toujours
juste. Ce biais n'est pas corrige -- il n'a pas ete explique, et le corriger
d'apres le seul generateur interne reviendrait a caler l'instrument sur
lui-meme. La methode ne remplace donc pas la cambrure la ou celle-ci
s'applique : elle prend le relais la ou celle-ci n'a pas de reponse.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import LOW, MEDIUM
from ..mesh import TriMesh
from .occupancy import OccupancyMap
from .topology import AXIAL, Topology


@dataclass
class NormalAngles:
    """Angles de pale et sens d'enroulement lus sur les normales."""

    beta1_deg: float = 0.0
    beta2_deg: float = 0.0
    slope_sign: int = 0
    area_1: float = 0.0  # aire de pale exploitee a l'entree, en m2
    area_2: float = 0.0  # aire de pale exploitee a la sortie, en m2
    stations: list[tuple[float, float]] = field(default_factory=list)  # (rayon, beta)
    confidence: str = LOW
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable, angles en degres et aires en SI."""
        return {
            "beta1_deg": self.beta1_deg,
            "beta2_deg": self.beta2_deg,
            "signe_d_enroulement": self.slope_sign,
            "aire_entree_m2": self.area_1,
            "aire_sortie_m2": self.area_2,
            "stations": [{"rayon_m": r, "beta_deg": b} for r, b in self.stations],
            "confiance": self.confidence,
        }


def passage_runs(
    occupancy: OccupancyMap, blade_mask: list[list[bool]]
) -> list[list[tuple[int, int]]]:
    """Plages axiales **contigues** de pale, colonne radiale par colonne radiale.

    La marge de paroi se prend sur la plage qui contient la face, non sur
    l'etendue totale de la colonne : sur une aube en boucle les deux brins sont
    separes en hauteur, et une etendue prise du plus bas au plus haut enjamberait
    le vide entre eux -- la marge supprimerait alors les deux brins au lieu de
    leurs seuls chants.
    """
    runs: list[list[tuple[int, int]]] = []
    for ir in range(occupancy.nr):
        column: list[tuple[int, int]] = []
        start = -1
        for iz in range(occupancy.nz):
            if blade_mask[iz][ir]:
                if start < 0:
                    start = iz
            elif start >= 0:
                column.append((start, iz - 1))
                start = -1
        if start >= 0:
            column.append((start, occupancy.nz - 1))
        runs.append(column)
    return runs


def _run_at(column: list[tuple[int, int]], iz: int) -> tuple[int, int] | None:
    """Plage contigue contenant la rangee `iz`, s'il y en a une."""
    for low, high in column:
        if low <= iz <= high:
            return (low, high)
    return None


def band_beta(
    mesh: TriMesh,
    occupancy: OccupancyMap,
    blade_mask: list[list[bool]],
    runs: list[list[tuple[int, int]]],
    r_low: float,
    r_high: float,
    axial_span: bool = True,
) -> tuple[float, int, float]:
    """`(beta moyen en degres, signe d'enroulement, aire retenue)` sur une couronne.

    Seules comptent les faces de la zone de pales assez loin des deux parois de
    la veine : les autres sont les chants de l'aube, qui ne portent pas d'angle.
    """
    weighted = 0.0
    area_total = 0.0
    signed = 0.0
    for i, j, k in mesh.faces:
        a, b, c = mesh.vertices[i], mesh.vertices[j], mesh.vertices[k]
        cx = (a[0] + b[0] + c[0]) / 3.0
        cy = (a[1] + b[1] + c[1]) / 3.0
        cz = (a[2] + b[2] + c[2]) / 3.0
        radius = math.hypot(cx, cy)
        if radius <= 0.0 or not (r_low <= radius <= r_high):
            continue
        ir = occupancy.radial_index(radius)
        iz = occupancy.axial_index(cz)
        if not blade_mask[iz][ir]:
            continue
        span = _run_at(runs[ir], iz)
        if span is None:
            continue
        low, high = span
        margin = config.NORMAL_WALL_MARGIN * (high - low)
        if iz < low + margin or iz > high - margin:
            continue

        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        norm = math.sqrt(nx * nx + ny * ny + nz * nz)
        if norm <= 0.0:
            continue
        area = 0.5 * norm
        nx, ny, nz = nx / norm, ny / norm, nz / norm

        n_u = nx * (-cy / radius) + ny * (cx / radius)
        n_r = nx * (cx / radius) + ny * (cy / radius)
        n_m = math.sqrt(max(0.0, 1.0 - n_u * n_u))

        weighted += area * math.degrees(math.atan2(abs(n_u), n_m))
        area_total += area
        # Le signe se lit sur le rapport a la composante meridienne : intrados et
        # extrados portent des normales opposees, mais leurs deux composantes
        # changent de signe ensemble, donc le rapport garde le sien.
        reference = n_r if axial_span else nz
        if abs(reference) > config.NORMAL_SIGN_MIN:
            signed += area * (1.0 if n_u / reference > 0.0 else -1.0)

    if area_total <= 0.0:
        return 0.0, 0, 0.0
    return weighted / area_total, (1 if signed >= 0.0 else -1), area_total


def analyse(
    mesh: TriMesh,
    occupancy: OccupancyMap,
    topology: Topology,
    blade_mask: list[list[bool]] | None = None,
) -> NormalAngles:
    """beta1, beta2 et sens d'enroulement d'une roue dont la cambrure est hors jeu."""
    from .topology import clean_blade_mask

    result = NormalAngles()
    if blade_mask is None:
        blade_mask = clean_blade_mask(occupancy.blade_mask())

    if topology.machine_type == AXIAL:
        # Sur une helice axiale l'entree et la sortie sont separees en z, non en
        # rayon : le decoupage en couronnes n'a pas de sens. La cambrure y donne
        # de toute facon un profil propre, coupe par coupe.
        result.notes.append(
            "roue axiale : les angles se lisent sur la cambrure, pas sur les couronnes de rayon"
        )
        return result
    axial = False
    r_low, r_high = topology.r_1s, topology.r_2
    if r_high <= r_low:
        result.notes.append("rayons d'entree et de sortie confondus : angles non mesurables")
        return result

    runs = passage_runs(occupancy, blade_mask)
    width = config.NORMAL_BAND_FRACTION * (r_high - r_low)
    steps = max(2, int(round(1.0 / config.NORMAL_BAND_FRACTION)))
    for step in range(steps):
        low = r_low + (r_high - r_low) * step / steps
        beta, _, area = band_beta(
            mesh, occupancy, blade_mask, runs, low, low + width, not axial
        )
        if area > 0.0:
            result.stations.append((low + 0.5 * width, beta))

    beta1, sign1, result.area_1 = band_beta(
        mesh, occupancy, blade_mask, runs, r_low, r_low + width, not axial
    )
    beta2, sign2, result.area_2 = band_beta(
        mesh, occupancy, blade_mask, runs, r_high - width, r_high, not axial
    )
    result.beta1_deg, result.beta2_deg = beta1, beta2
    result.slope_sign = sign2 if sign2 else sign1

    if result.area_1 <= 0.0 or result.area_2 <= 0.0:
        result.notes.append(
            "aucune face de pale dans la couronne d'entree ou de sortie : angles non mesurables"
        )
        return result
    result.confidence = MEDIUM if sign1 == sign2 else LOW
    result.notes.append(
        f"beta lus sur les normales de la surface d'aube, sur {len(result.stations)} couronnes "
        f"entre {r_low * config.MM_PER_M:.1f} et {r_high * config.MM_PER_M:.1f} mm ; methode "
        "basse de 2 a 5 degres sur des roues d'angles connus, biais non corrige"
    )
    if sign1 != sign2:
        result.notes.append(
            "le sens d'enroulement differe entre l'entree et la sortie : aube fortement "
            "vrillee, ou lecture perturbee"
        )
    return result
