"""Interface en ligne de commande (SPEC phase 7).

    python -m impeller_analyzer roue.stl \\
        --unit cm --r-aspiration 4.5 --rpm 1000 2000 3000 \\
        --temperature 20 --altitude 0 --hauteur-aspiration 0 --out rapport/
"""

from __future__ import annotations

import argparse
import sys

from . import __version__, config
from .analysis import Options, run
from .io import loader, report

EPILOG = """\
Toutes les grandeurs de sortie sont en SI dans resultats.json ; le rapport
Markdown les convertit en unites d'atelier (mm, m3/h, kW).

L'unite d'import par defaut est le centimetre. Le rayon d'aspiration se
fournit en centimetres et prime toujours sur la valeur detectee.

Si l'extraction des angles de pale est peu sure, --beta1, --beta2 et --blades
permettent de les imposer et le reste du calcul se poursuit normalement.

Pour travailler dans le navigateur plutot qu'en ligne de commande :

    python -m impeller_analyzer.serve

deposez la roue dans la page, l'analyse s'execute sur cette machine.
"""


def build_parser() -> argparse.ArgumentParser:
    """Construit l'analyseur d'arguments."""
    parser = argparse.ArgumentParser(
        prog="impeller-analyzer",
        description=(
            "Analyse hydraulique d'une helice ou d'une roue de pompe a partir d'un fichier 3D : "
            "geometrie extraite, sens de rotation, performances et cavitation."
        ),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("fichier", help="fichier de geometrie (.stl, .obj, .ply, .off, .step, .dxf...)")
    parser.add_argument(
        "--unit",
        default="cm",
        help="unite du fichier d'entree ou facteur vers le metre (defaut : cm)",
    )
    parser.add_argument(
        "--r-aspiration",
        type=float,
        default=None,
        metavar="CM",
        help="rayon d'aspiration impose, en centimetres ; prime sur la detection",
    )
    parser.add_argument(
        "--rpm",
        type=float,
        nargs="+",
        default=list(config.DEFAULT_RPM),
        metavar="N",
        help="regimes analyses, en tr/min (defaut : %(default)s)",
    )
    parser.add_argument("--blades", type=int, default=None, help="nombre de pales impose")
    parser.add_argument("--beta1", type=float, default=None, metavar="DEG", help="angle de pale d'entree impose")
    parser.add_argument("--beta2", type=float, default=None, metavar="DEG", help="angle de pale de sortie impose")
    parser.add_argument(
        "--aspiration", choices=["auto", "+z", "-z"], default="auto",
        help=(
            "cote par lequel la roue aspire (defaut : auto, deduit de la geometrie ; "
            "sur une roue axiale la deduction est impossible et la convention +Z est gardee)"
        ),
    )
    parser.add_argument(
        "--temperature", type=float, default=config.TEMPERATURE, metavar="C",
        help="temperature du liquide, en degres Celsius (defaut : %(default)s)",
    )
    parser.add_argument(
        "--altitude", type=float, default=config.ALTITUDE, metavar="M",
        help="altitude du site, en metres (defaut : %(default)s)",
    )
    parser.add_argument(
        "--hauteur-aspiration", type=float, default=config.HAUTEUR_ASPIRATION, metavar="M",
        help="hauteur d'aspiration, positive en charge (defaut : %(default)s)",
    )
    parser.add_argument(
        "--pertes-aspiration", type=float, default=config.PERTES_ASPIRATION, metavar="M",
        help="pertes de charge de la conduite d'aspiration (defaut : %(default)s)",
    )
    parser.add_argument("--out", default="rapport", metavar="DOSSIER", help="dossier de sortie (defaut : %(default)s)")
    parser.add_argument(
        "--grille", type=int, nargs=2, default=[config.GRID_NR, config.GRID_NZ], metavar=("NR", "NZ"),
        help="taille de la carte d'occupation (defaut : %(default)s)",
    )
    parser.add_argument(
        "--secteurs", type=int, default=config.N_THETA, metavar="N",
        help="nombre de secteurs angulaires de la carte (defaut : %(default)s)",
    )
    parser.add_argument("--sans-reparation", action="store_true", help="ne pas reparer le maillage a l'import")
    parser.add_argument(
        "--sans-controle-symetrie", action="store_true",
        help="ne pas verifier la periodicite par distance de Hausdorff (plus rapide)",
    )
    parser.add_argument("--sans-trimesh", action="store_true", help="forcer le lecteur interne pour les maillages")
    parser.add_argument("--sans-vue3d", action="store_true", help="ne pas produire la vue 3D interactive")
    parser.add_argument("--quiet", action="store_true", help="n'ecrire que les fichiers, sans resume console")
    parser.add_argument("--version", action="version", version=f"impeller-analyzer {__version__}")
    return parser


def options_from_args(args: argparse.Namespace) -> Options:
    """Convertit les arguments de la ligne de commande en options d'analyse."""
    return Options(
        unit=args.unit,
        r_aspiration_cm=args.r_aspiration,
        speeds=tuple(args.rpm),
        blades=args.blades,
        beta1_deg=args.beta1,
        beta2_deg=args.beta2,
        suction=args.aspiration,
        altitude=args.altitude,
        temperature_c=args.temperature,
        suction_height=args.hauteur_aspiration,
        suction_losses=args.pertes_aspiration,
        grid_nr=args.grille[0],
        grid_nz=args.grille[1],
        n_theta=args.secteurs,
        repair=not args.sans_reparation,
        prefer_trimesh=not args.sans_trimesh,
        symmetry_check=not args.sans_controle_symetrie,
    )


def summarise(result, produced: dict[str, str]) -> str:
    """Resume console de l'analyse."""
    lines = []
    topology = result.topology
    blades = result.blades
    looped = result.blade_loops is not None and result.blade_loops.looped
    if topology is not None and blades is not None:
        lines.append(
            f"Roue {topology.machine_type}, {topology.blades.n_blades} pales, "
            f"r1s = {topology.r_1s * config.MM_PER_M:.1f} mm, r2 = {topology.r_2 * config.MM_PER_M:.1f} mm, "
            f"beta1/beta2 = {blades.beta1_deg:.1f}/{blades.beta2_deg:.1f} deg"
        )
        lines.append(f"Sens de rotation : {blades.rotation_label}")
    if looped:
        # Le premier chiffre lu est celui qu'on croit : autant dire tout de suite
        # que la suite du resume ne decrit pas cette roue.
        lines.append(
            "Aubes en boucle fermee (type toroidal) : beta1 et beta2 sont lus sur les normales "
            "de la surface, la cambrure ne s'appliquant pas. Voir les reserves du rapport."
        )
    for curve in result.curves:
        point = curve.nominal_point()
        if point is None:
            continue
        lines.append(
            f"{curve.rpm:6.0f} tr/min : Q = {point.flow * config.SECONDS_PER_HOUR:8.1f} m3/h, "
            f"H = {point.head:7.2f} m, P = {point.shaft_power / config.W_PER_KW:7.2f} kW, "
            f"NPSHr = {point.npshr:5.2f} m"
        )
    if result.speed_limit is not None:
        lines.append(
            f"Vitesse maximale sans cavitation : {result.speed_limit.rpm_max} tr/min "
            f"(limite {result.speed_limit.active_limit})"
        )
    lines.append(f"Confiance globale : {result.overall_confidence()}")
    if result.warnings:
        lines.append(f"{len(result.warnings)} avertissement(s), voir le rapport")
    for key in ("viewer", "markdown", "json", "curves", "occupancy"):
        if key in produced:
            lines.append(f"  ecrit : {produced[key]}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """Point d'entree ; renvoie le code de sortie du processus."""
    args = build_parser().parse_args(argv)
    try:
        result = run(args.fichier, options_from_args(args))
    except loader.ImportError_ as error:
        print(f"erreur d'import : {error}", file=sys.stderr)
        return 2
    except ValueError as error:
        print(f"erreur : {error}", file=sys.stderr)
        return 2

    produced = report.write_all(result, args.out, source=args.fichier, viewer_page=not args.sans_vue3d)
    if not args.quiet:
        print(summarise(result, produced))
    return 0
