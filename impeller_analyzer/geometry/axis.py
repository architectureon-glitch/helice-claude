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
from ..mesh import TriMesh, dot, normalize, rotation_between
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
    aligned = mesh.copy()
    aligned.apply_translation((-result.centre[0], -result.centre[1], -result.centre[2]))
    if result.angle_to_z_deg > config.AXIS_WARN_DEG and result.confidence != LOW:
        aligned.apply_rotation(rotation_between(result.axis, REFERENCE_AXIS))
        result.realigned = True
    # Recentrage final : l'origine est sur l'axe, a la hauteur du barycentre.
    centre = aligned.centroid()
    aligned.apply_translation((-centre[0], -centre[1], -centre[2]))
    return aligned, result
