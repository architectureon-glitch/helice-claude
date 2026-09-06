"""impeller-analyzer : analyse hydraulique d'une roue de pompe a partir d'un maillage 3D.

Convention d'orientation du paquet (SPEC 0) : axe de rotation = Z, aspiration
vers +Z, refoulement vers -Z (axial) ou radial (centrifuge).  Toutes les
grandeurs manipulees et stockees sont en SI.
"""

__all__ = ["config", "confidence"]
__version__ = "0.8.0"
