"""Detection de l'axe de revolution, realignement et recentrage (SPEC 2.1).

La convention du projet impose l'axe Z et l'aspiration vers +Z ; l'outil le
verifie tout de meme a partir du tenseur d'inertie et realigne le maillage si
necessaire, en signalant l'ecart.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .. import config
from ..confidence import HIGH, LOW, MEDIUM
from ..mesh import TriMesh, dot, normalize, rotation_between, rotation_matrix
from ..numeric import jacobi_eigen

Vec3 = tuple[float, float, float]

#: Axe de reference du projet : la rotation se fait autour de Z, l'aspiration
#: vient de +Z (SPEC 0).
REFERENCE_AXIS: Vec3 = (0.0, 0.0, 1.0)


@dataclass
class AxisResult:
    """Resultat de la detection d'axe."""

    axis: Vec3 = REFERENCE_AXIS  # axe detecte, unitaire, dans le repere d'entree
    centre: Vec3 = (0.0, 0.0, 0.0)  # barycentre du volume, dans le repere d'entree
    eigenvalues: tuple[float, float, float] = (0.0, 0.0, 0.0)
    eigenvectors: tuple[Vec3, Vec3, Vec3] = (REFERENCE_AXIS,) * 3
    pair_ratio: float = 0.0  # |lambda_a - lambda_b| / max(lambda) du couple retenu
    second_pair_ratio: float = 0.0  # meme grandeur pour le couple concurrent
    angle_to_z_deg: float = 0.0
    realigned: bool = False
    confidence: str = LOW
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "axe_detecte": list(self.axis),
            "centre_m": list(self.centre),
            "valeurs_propres": list(self.eigenvalues),
            "rapport_couple_egal": self.pair_ratio,
            "rapport_couple_concurrent": self.second_pair_ratio,
            "ecart_a_z_deg": self.angle_to_z_deg,
            "realigne": self.realigned,
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
        }


def detect_axis(mesh: TriMesh) -> AxisResult:
    """Detecte l'axe de revolution par diagonalisation du tenseur d'inertie.

    Une roue a N pales (N >= 3) a deux valeurs propres quasi egales ; l'axe est
    le vecteur propre de la valeur propre isolee.  Le critere d'egalite est
    `|lambda_a - lambda_b| / max(lambda) < AXIS_TOL`.
    """
    result = AxisResult()
    centre = mesh.centroid()
    result.centre = centre
    tensor = mesh.inertia_tensor(about=centre)
    values, vectors = jacobi_eigen(tensor)
    result.eigenvalues = tuple(values)  # type: ignore[assignment]
    result.eigenvectors = tuple(tuple(v) for v in vectors)  # type: ignore[assignment]

    largest = max(abs(v) for v in values) or 1.0
    # Les trois couples possibles, avec l'indice de la valeur propre isolee.
    pairs = [
        (abs(values[0] - values[1]) / largest, 2),
        (abs(values[0] - values[2]) / largest, 1),
        (abs(values[1] - values[2]) / largest, 0),
    ]
    pairs.sort()
    result.pair_ratio = pairs[0][0]
    result.second_pair_ratio = pairs[1][0]
    axis = normalize(vectors[pairs[0][1]])

    if result.pair_ratio < config.AXIS_TOL and result.second_pair_ratio >= config.AXIS_TOL:
        result.confidence = HIGH
    elif result.second_pair_ratio < config.AXIS_TOL:
        # Les trois valeurs propres sont quasi egales (corps quasi spherique ou
        # cubique) : aucun axe n'emerge, on garde la convention.
        result.confidence = LOW
        axis = REFERENCE_AXIS
        result.warnings.append(
            "inertie quasi isotrope : aucun axe de revolution ne se degage, "
            "l'axe Z de la convention est conserve"
        )
    else:
        result.confidence = MEDIUM
        result.warnings.append(
            f"symetrie de revolution imparfaite (ecart relatif des valeurs propres "
            f"{result.pair_ratio:.3f} > {config.AXIS_TOL})"
        )

    # L'axe propre est defini au signe pres : on choisit celui qui pointe vers +Z.
    if dot(axis, REFERENCE_AXIS) < 0.0:
        axis = (-axis[0], -axis[1], -axis[2])
    result.axis = axis
    result.angle_to_z_deg = math.degrees(math.acos(max(-1.0, min(1.0, dot(axis, REFERENCE_AXIS)))))
    if result.angle_to_z_deg > config.AXIS_WARN_DEG and result.confidence != LOW:
        result.warnings.append(
            f"l'axe de revolution detecte s'ecarte de Z de {result.angle_to_z_deg:.1f} deg "
            f"(> {config.AXIS_WARN_DEG} deg) : le maillage est realigne sur Z"
        )
    return result


def align_to_z(mesh: TriMesh, result: AxisResult | None = None) -> tuple[TriMesh, AxisResult]:
    """Renvoie une copie du maillage d'axe Z, recentree sur le barycentre.

    Le maillage d'entree n'est pas modifie.  Le realignement n'est applique que
    si l'ecart depasse `AXIS_WARN_DEG` ; en deca, la convention Z est appliquee
    telle quelle pour ne pas introduire de rotation parasite.
    """
    if result is None:
        result = detect_axis(mesh)
    rotate = result.angle_to_z_deg > config.AXIS_WARN_DEG and result.confidence != LOW
    aligned = orient(mesh, result.centre, result.axis if rotate else REFERENCE_AXIS)
    result.realigned = rotate
    return aligned, result


def orient(mesh: TriMesh, centre: Vec3, axis: Vec3) -> TriMesh:
    """Copie du maillage recentree sur `centre` et tournee pour amener `axis` sur Z.

    L'origine finale est sur l'axe, a la hauteur du barycentre du volume.
    """
    aligned = mesh.copy()
    aligned.apply_translation((-centre[0], -centre[1], -centre[2]))
    if dot(normalize(axis), REFERENCE_AXIS) < 1.0 - 1e-12:
        aligned.apply_rotation(rotation_between(normalize(axis), REFERENCE_AXIS))
    centre = aligned.centroid()
    aligned.apply_translation((-centre[0], -centre[1], -centre[2]))
    return aligned


def candidate_axes(result: AxisResult) -> list[Vec3]:
    """Axes a departager : Z de la convention, puis les trois axes principaux.

    L'axe de rotation d'une roue a N pales est toujours un axe principal
    d'inertie -- la periodicite l'impose -- mais pas forcement celui des deux
    valeurs propres egales : pour N = 2, une pale elancee est quasi symetrique
    autour de sa propre envergure, et c'est l'envergure que ce critere designe.
    """
    axes: list[Vec3] = [REFERENCE_AXIS]
    for vector in result.eigenvectors:
        axis = normalize(vector)
        if dot(axis, REFERENCE_AXIS) < 0.0:
            axis = (-axis[0], -axis[1], -axis[2])
        if all(abs(dot(axis, known)) < math.cos(math.radians(config.AXIS_WARN_DEG)) for known in axes):
            axes.append(axis)
    return axes


#: Reponses possibles a la question « de quel cote aspire cette roue ? ».
SUCTION_AUTO = "auto"
SUCTION_PLUS_Z = "+z"
SUCTION_MINUS_Z = "-z"


@dataclass
class SuctionResult:
    """Cote par lequel la roue aspire (SPEC 2.1, verification de la convention)."""

    sign: int = 0  # +1 conforme a la convention, -1 maillage a retourner, 0 indetermine
    asymmetry: float = 0.0  # (rayon moyen en haut - en bas) / rayon exterieur
    flipped: bool = False
    forced: bool = False
    confidence: str = LOW
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Vue serialisable en JSON."""
        return {
            "signe": self.sign,
            "asymetrie_radiale": self.asymmetry,
            "maillage_retourne": self.flipped,
            "impose_par_l_utilisateur": self.forced,
            "confiance": self.confidence,
            "avertissements": list(self.warnings),
            "remarques": list(self.notes),
        }


def meridional_asymmetry(occupancy, blade_mask=None) -> float:
    """Asymetrie radiale de la veine entre le haut et le bas de la zone de pales.

    Le rayon moyen des cellules de pales est calcule separement au-dessus et
    au-dessous du milieu de la zone, pondere par l'element de volume `r dr dz`,
    puis leur ecart est rapporte au rayon exterieur.  Une roue a composante
    radiale voit sa veine partir du petit rayon (l'oeillard) vers le grand : la
    grandeur est nettement negative quand la convention est respectee.  Sur une
    roue purement axiale elle reste proche de zero, la veine gardant le meme
    rayon d'un bout a l'autre.
    """
    from .topology import clean_blade_mask

    if blade_mask is None:
        blade_mask = clean_blade_mask(occupancy.blade_mask())
    rows = [iz for iz in range(occupancy.nz) if any(blade_mask[iz])]
    if not rows or occupancy.r_max <= 0.0:
        return 0.0
    middle = 0.5 * (occupancy.z_centres[rows[0]] + occupancy.z_centres[rows[-1]])
    sums = [0.0, 0.0]
    weights = [0.0, 0.0]
    for iz in rows:
        side = 0 if occupancy.z_centres[iz] > middle else 1
        for ir in range(occupancy.nr):
            if not blade_mask[iz][ir]:
                continue
            radius = occupancy.r_centres[ir]
            sums[side] += radius * radius
            weights[side] += radius
    if weights[0] <= 0.0 or weights[1] <= 0.0:
        return 0.0
    return (sums[0] / weights[0] - sums[1] / weights[1]) / occupancy.r_max


def detect_suction_side(occupancy, forced: str = SUCTION_AUTO) -> SuctionResult:
    """Determine si la roue aspire bien vers +Z, comme la convention l'impose.

    Un fichier issu d'un logiciel de CAO n'a aucune raison de respecter cette
    convention : une roue exportee a l'envers voit son plan d'aspiration lu du
    cote du refoulement, ce qui la fait passer pour axiale et lui donne un sens
    de sortie faux.  Sur une roue centrifuge ou mixte l'ambiguite se leve sans
    rien demander a l'utilisateur ; sur une roue axiale elle ne se leve pas, et
    la convention est alors conservee telle quelle.
    """
    result = SuctionResult()
    result.asymmetry = meridional_asymmetry(occupancy)

    if forced == SUCTION_PLUS_Z:
        result.sign, result.forced, result.confidence = 1, True, HIGH
        result.notes.append("cote aspiration impose vers +Z par l'utilisateur")
        return result
    if forced == SUCTION_MINUS_Z:
        result.sign, result.forced, result.confidence = -1, True, HIGH
        result.flipped = True
        result.notes.append("cote aspiration impose vers -Z : le maillage est retourne")
        return result

    if abs(result.asymmetry) < config.SUCTION_ASYMMETRY_MIN:
        result.sign = 0
        result.confidence = LOW
        result.notes.append(
            f"veine de rayon constant (asymetrie {result.asymmetry:+.3f}) : le cote aspiration "
            "ne se deduit pas de la geometrie, la convention +Z est conservee. C'est le cas "
            "normal d'une helice axiale ; si la roue est montee a l'envers, utilisez "
            "--aspiration -z"
        )
        return result

    result.confidence = HIGH
    if result.asymmetry < 0.0:
        result.sign = 1
        result.notes.append(
            f"la veine s'ecarte de l'axe vers -Z (asymetrie {result.asymmetry:+.3f}) : "
            "l'aspiration est bien du cote +Z, conforme a la convention"
        )
    else:
        result.sign = -1
        result.flipped = True
        result.warnings.append(
            f"la veine s'ecarte de l'axe vers +Z (asymetrie {result.asymmetry:+.3f}) : la roue "
            "est fournie a l'envers, aspiration du cote -Z. Le maillage est retourne pour "
            "respecter la convention ; sens de rotation et sens de sortie sont donnes dans "
            "le repere corrige"
        )
    return result


def flip_axis(mesh: TriMesh) -> TriMesh:
    """Retourne le maillage bout pour bout : demi-tour autour de X."""
    return mesh.transformed(matrix=rotation_matrix((1.0, 0.0, 0.0), math.pi))
