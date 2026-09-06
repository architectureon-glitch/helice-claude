"""Fichier de calage unique de l'outil.

Toutes les constantes physiques, tous les coefficients empiriques et tous les
seuils numériques de `impeller_analyzer` sont déclarés ici, une constante par
ligne, commentée avec sa source.  Aucune autre valeur numérique significative ne
doit apparaître ailleurs dans le paquet : pour recaler l'outil on ne modifie que
ce fichier.

Convention d'unités : tout est en SI (m, m^3/s, rad/s, Pa, kg/m^3), sauf les
grandeurs explicitement suffixées (`_RPM`, `_DEG`, `_MM`, `_C`).
"""

# ---------------------------------------------------------------------------
# Fluide (eau, 20 degres C)
# ---------------------------------------------------------------------------
RHO = 998.2  # kg/m3 - masse volumique de l'eau a 20 C (SPEC 9 ; NIST/IAPWS-95)
P_VAP = 2339.0  # Pa - pression de vapeur saturante de l'eau a 20 C (SPEC 9)
G = 9.80665  # m/s2 - acceleration normale de la pesanteur (CGPM 1901, SPEC 9)

# Correlation d'Antoine pour l'eau, forme log10(P_mmHg) = A - B / (C + T_C),
# domaine 1-100 C (Antoine 1888 ; constantes NIST Webbook, jeu "Stull 1947").
ANTOINE_A = 8.07131  # - (log10, P en mmHg, T en degres C)
ANTOINE_B = 1730.63  # degres C
ANTOINE_C = 233.426  # degres C
ANTOINE_MMHG_TO_PA = 133.322368421  # Pa/mmHg - definition du millimetre de mercure

# Masse volumique de l'eau : table (T en degres C -> kg/m3), interpolation lineaire.
# Source : NIST/IAPWS-95, eau liquide a pression atmospherique.
RHO_TABLE_T_C = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0)
RHO_TABLE_KG_M3 = (999.84, 999.96, 999.70, 999.10, 998.20, 997.05, 995.65, 992.22, 988.04, 983.20, 977.76, 971.79, 965.31, 958.35)

# Atmosphere standard OACI : p_atm = P_ATM_SEA_LEVEL * (1 - K * h)^EXP.
P_ATM_SEA_LEVEL = 101325.0  # Pa - pression atmospherique normale au niveau de la mer (ISO 2533)
ATM_LAPSE_COEF = 2.25577e-5  # 1/m - coefficient de decroissance barometrique (ISO 2533, SPEC 6.2)
ATM_EXPONENT = 5.2559  # - - exposant barometrique (ISO 2533, SPEC 6.2)

# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------
UNIT_FACTOR = 0.01  # m/unite - facteur par defaut cm -> m (SPEC 1.1)
MERGE_TOL = 1e-6  # m - tolerance de fusion des sommets dupliques (SPEC 1.2)
STEP_TESSELLATION = 2e-4  # m - tolerance de tessellation STEP/IGES, 0.2 mm (SPEC 1)

# Facteurs des unites acceptees par --unit, vers le metre.
UNIT_FACTORS = {
    "m": 1.0,  # metre
    "dm": 0.1,  # decimetre
    "cm": 0.01,  # centimetre - unite d'import par defaut du projet
    "mm": 0.001,  # millimetre
    "um": 1e-6,  # micrometre
    "in": 0.0254,  # pouce international (1959)
    "ft": 0.3048,  # pied international (1959)
}

# Extensions reconnues par le lecteur, groupees par famille de traitement.
EXT_MESH_NATIVE = (".stl", ".obj", ".ply", ".off")  # lus par le parseur interne
EXT_MESH_TRIMESH_ONLY = (".3ds",)  # necessitent trimesh (SPEC phase 1, priorite 1)
EXT_CAD_BREP = (".step", ".stp", ".iges", ".igs")  # tesselles par cadquery/OCC (priorite 2)
EXT_DXF = (".dxf",)  # maillages POLYFACE / 3DSOLID tesselle, via ezdxf (priorite 3)
EXT_DWG = (".dwg",)  # converti en DXF par ODA File Converter (priorite 4)
EXT_REFUSED = (".lisp", ".lsp")  # AutoLISP : pas un format geometrique (SPEC phase 1)

FILL_HOLES_MAX_EDGES = 64  # aretes - taille maximale d'un trou rebouche par triangulation en eventail
REPAIR_MAX_PASSES = 4  # - - nombre maximal de passes d'orientation/rebouchage

# ---------------------------------------------------------------------------
# Grille d'occupation
# ---------------------------------------------------------------------------
GRID_NR = 200  # cellules radiales de la carte f(r, z) (SPEC 2.2)
GRID_NZ = 200  # cellules axiales de la carte f(r, z) (SPEC 2.2)
N_THETA = 720  # echantillons angulaires du signal g(theta) (SPEC 2.2 et 3.1)
F_SOLIDE = 0.98  # - - seuil de fraction angulaire au-dela duquel la cellule est du moyeu/flasque
F_VIDE = 0.02  # - - seuil en deca duquel la cellule est de la veine fluide ou l'exterieur

GRID_R_MARGIN = 1e-9  # m - retrait applique aux bornes de la grille pour eviter les tangences exactes
SLICE_EPS_REL = 1e-9  # - - tolerance relative de coupe plan/triangle (fraction de la taille de bbox)
OCCUPANCY_MIN_CELL = 1e-12  # m - taille minimale d'une cellule, garde-fou contre les bbox degenerees

# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------
BLADES_MIN = 2  # nombre de pales minimal cherche (SPEC 3.1)
BLADES_MAX = 12  # nombre de pales maximal cherche (SPEC 3.1)
FFT_RATIO_MIN = 3.0  # - - confiance high si l'harmonique dominante depasse 3x la suivante (SPEC 3.1)
FFT_RATIO_MEDIUM = 1.5  # - - en deca de ce rapport la detection tombe en confiance low
AXIS_TOL = 0.05  # - - tolerance d'egalite relative des valeurs propres d'inertie (SPEC 2.1)
SYM_TOL = 0.02  # - - distance de Hausdorff relative admise apres rotation de 2*pi/N (SPEC 3.1)
AXIS_WARN_DEG = 5.0  # degres - ecart axe detecte / Z au-dela duquel on avertit et on realigne (SPEC 2.1)
HAUSDORFF_SAMPLES = 4000  # points - echantillonnage du maillage pour la distance de Hausdorff
PROXIMITY_CELLS = 32  # cellules - resolution du hachage spatial des triangles (plus grande dimension)
BLADE_ROW_MIN_CELLS = 1  # cellules - occupation minimale d'une ligne pour compter dans la zone de pales
BLADE_BLOB_MIN_CELLS = 20  # cellules - taille minimale absolue d'un ilot de zone de pales retenu
BLADE_BLOB_MIN_FRACTION = 0.05  # - - taille minimale d'un ilot, en fraction du plus gros ilot
JACOBI_MAX_SWEEPS = 64  # - - nombre maximal de balayages de la diagonalisation de Jacobi
JACOBI_TOL = 1e-14  # - - seuil d'arret relatif de la diagonalisation de Jacobi

# ---------------------------------------------------------------------------
# Sections et cambrure
# ---------------------------------------------------------------------------
N_SECTIONS = 11  # coupes cylindriques reparties du moyeu au tip (SPEC 4.1)
N_STATIONS = 25  # stations de decoupe de la corde pour la ligne de cambrure (SPEC 4.2)
SPLINE_SMOOTH = 0.001  # - - lissage relatif : residu RMS vise = SPLINE_SMOOTH * corde (SPEC 4.2)
SECTION_MARGIN = 0.02  # - - retrait relatif aux deux bouts de la plage [r_moyeu, r_tip] (SPEC 4.1)
BETA_CHORD_FRACTION = 0.10  # - - fraction de corde sur laquelle sont moyennes beta1 et beta2 (SPEC 4.3)
SPLINE_LAMBDA_MIN = 1e-18  # - - borne basse de la recherche du parametre de penalisation
SPLINE_LAMBDA_MAX = 1e12  # - - borne haute de la recherche du parametre de penalisation
SPLINE_BISECTION_STEPS = 80  # - - iterations de bissection sur log(lambda)
MIN_PROFILE_POINTS = 8  # points - taille minimale d'une polyligne pour etre exploitee comme profil
MIN_BLADE_SECTIONS = 3  # coupes - nombre minimal de coupes exploitables pour une confiance > low

# ---------------------------------------------------------------------------
# Classification de la roue
# ---------------------------------------------------------------------------
R_RATIO_AXIAL_MAX = 1.15  # - - r2/r1s en deca duquel la roue est axiale (SPEC 3.3)
R_RATIO_MIXED_MAX = 1.80  # - - r2/r1s en deca duquel la roue est mixte, au-dela centrifuge (SPEC 3.3)
NQ_CENTRIFUGAL_MAX = 35.0  # tr/min - vitesse specifique n_q maximale d'une roue centrifuge (SPEC 3.3)
NQ_MIXED_MAX = 80.0  # tr/min - n_q maximale d'une roue mixte, au-dela axiale (SPEC 3.3)
B2_RADIUS_FRACTION = 0.98  # - - rayon relatif ou est mesuree la largeur de sortie b2 (SPEC 3.4)
HELIX_PITCH_RATIO_MAX = 8.0  # - - pas helicoidal / rayon en deca duquel le sens axial est fiable (SPEC 4.4)
BETA2_RADIAL_DEG = 85.0  # degres - au-dela, aubes radiales : sens de rotation ambigu (SPEC 4.4)

# ---------------------------------------------------------------------------
# Obstruction
# ---------------------------------------------------------------------------
TAU_1 = 0.90  # - - coefficient d'obstruction de la section d'entree (SPEC 9)
TAU_2 = 0.92  # - - coefficient d'obstruction de la section de sortie (SPEC 9)

# ---------------------------------------------------------------------------
# Hydraulique
# ---------------------------------------------------------------------------
ETA_H = 0.88  # - - rendement hydraulique (SPEC 5.4)
ETA_VOL = 0.96  # - - rendement volumetrique (SPEC 5.4)
ETA_MEC = 0.95  # - - rendement mecanique (SPEC 5.4)
K_FROTTEMENT_REL = 0.06  # - - perte de frottement au nominal, en fraction de H_th (SPEC 5.3)
XI_INCIDENCE = 0.75  # - - coefficient de perte par incidence (SPEC 5.3)
Q_SWEEP_MAX = 1.40  # - - borne haute du balayage de debit, en fraction de Q_n (SPEC 5.3)
Q_SWEEP_POINTS = 29  # points du balayage H-Q (SPEC 5.3)

WIESNER_EXPONENT = 0.7  # - - exposant du nombre de pales dans la correlation de Wiesner (1967)
RPM_TO_RAD_S = 2.0 * 3.141592653589793 / 60.0  # rad/s par tr/min - conversion de vitesse de rotation
SIMILARITY_TOL = 0.02  # - - ecart maximal admis entre calcul direct et similitude (SPEC 5.5)
DEFAULT_RPM = (1000.0, 2000.0, 3000.0)  # tr/min - regimes analyses par defaut (SPEC 0.4)
BETA_MIN_DEG = 1.0  # degres - garde-fou bas sur beta1/beta2 (evite tan(beta) -> 0)
BETA_MAX_DEG = 89.0  # degres - garde-fou haut sur beta1/beta2 (evite tan(beta) -> inf)
Q_MIN_RELATIVE = 1e-6  # - - debit plancher du balayage, en fraction de Q_n (evite la division par zero)

# ---------------------------------------------------------------------------
# Cavitation
# ---------------------------------------------------------------------------
LAMBDA_C = 1.10  # - - coefficient de la composante meridienne du NPSHr (SPEC 6.1, methode A)
LAMBDA_W = 0.28  # - - coefficient de la composante relative du NPSHr (SPEC 6.1, methode A)
N_SS = 180.0  # tr/min - vitesse specifique d'aspiration (SPEC 6.1, methode B ; unites SI+tr/min)
NSS_EXPONENT = 4.0 / 3.0  # - - exposant de la methode B du NPSHr (SPEC 6.1)
MARGE_NPSH = 1.30  # - - marge de securite NPSHa/NPSHr imposee (SPEC 6.3)
W1S_MAX = 30.0  # m/s - vitesse relative maximale admise en entree, eau/roue metallique (SPEC 6.3)
RPM_ROUNDING = 10.0  # tr/min - pas d'arrondi (a la dizaine inferieure) de la vitesse maximale (SPEC 6.3)

# ---------------------------------------------------------------------------
# Installation par defaut
# ---------------------------------------------------------------------------
ALTITUDE = 0.0  # m - altitude du site (SPEC 6.2)
TEMPERATURE = 20.0  # degres C - temperature du liquide (SPEC 6.2)
HAUTEUR_ASPIRATION = 0.0  # m - hauteur d'aspiration, positif = en charge (SPEC 6.2)
PERTES_ASPIRATION = 0.5  # m - pertes de charge de la conduite d'aspiration (SPEC 6.2)

# ---------------------------------------------------------------------------
# Incertitude annoncee du modele 1D (SPEC 7, encadre systematique)
# ---------------------------------------------------------------------------
UNCERTAINTY_H = 0.18  # - - incertitude relative annoncee sur la hauteur
UNCERTAINTY_Q = 0.25  # - - incertitude relative annoncee sur le debit
UNCERTAINTY_NPSH = 0.30  # - - incertitude relative annoncee sur le NPSHr

# ---------------------------------------------------------------------------
# Rendu graphique (PNG produits sans dependance externe)
# ---------------------------------------------------------------------------
PLOT_WIDTH_PX = 1200  # px - largeur des figures produites
PLOT_HEIGHT_PX = 800  # px - hauteur des figures produites
PLOT_MARGIN_PX = 70  # px - marge autour de la zone tracee
PLOT_TICKS = 6  # - - nombre de graduations par axe
PLOT_DPI = 96  # px/pouce - resolution declaree dans l'entete PNG
PLOT_BAND_VIDE = 0.15  # - - position d'affichage du seuil F_VIDE sur l'echelle de couleur
PLOT_BAND_PLEIN = 0.85  # - - position d'affichage du seuil F_SOLIDE sur l'echelle de couleur

# ---------------------------------------------------------------------------
# Tolerances de la campagne de validation (SPEC phase 8)
# ---------------------------------------------------------------------------
VALID_GEOM_TOL = 0.02  # - - ecart relatif admis sur N, beta1, beta2, r1s, r2 (SPEC 8.1)
VALID_REFERENCE_TOL = 0.18  # - - ecart relatif admis sur H au BEP du cas de reference (SPEC 8.2)
VALID_INVARIANCE_TOL = 0.005  # - - ecart relatif admis apres rotation/translation (SPEC 8.3)
VALID_DECIMATION_DEG = 3.0  # degres - ecart admis sur beta2 apres decimation a 20 % (SPEC 8.4)
VALID_DECIMATION_RATIO = 0.20  # - - fraction de triangles conservee par la decimation (SPEC 8.4)
