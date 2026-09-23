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
FLOAT32_EPS = 2.0 ** -24  # - - demi-ecart relatif entre deux flottants simple precision voisins (IEEE 754) : resolution d'un STL binaire
FLOAT32_WARN = 2e-4  # - - resolution simple precision rapportee a la taille de la piece au-dela de laquelle on avertit (piece a ~2 km de l'origine pour 18 cm : ecarts de 0.1 % mesures)
FLOAT32_LOW = 5e-3  # - - meme rapport au-dela duquel la geometrie importee n'est plus fiable : confiance du maillage abaissee a faible
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
F_MATIERE = 1e-4  # - - seuil de **presence de matiere** dans une cellule, distinct de F_VIDE : l'occupation angulaire d'une pale vaut N.e / (2.pi.r), soit 0.017 pour quatre pales de 4 mm a r = 90 mm -- sous F_VIDE, alors que la matiere y est bien

SLICE_EPS_REL = 1e-9  # - - tolerance relative de coupe plan/triangle (fraction de la taille de bbox)
OCCUPANCY_MIN_CELL = 1e-12  # m - taille minimale d'une cellule, garde-fou contre les bbox degenerees

# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------
BLADES_MIN = 2  # nombre de pales minimal cherche (SPEC 3.1)
BLADES_MAX = 12  # nombre de pales maximal cherche (SPEC 3.1)
FFT_RATIO_MIN = 3.0  # - - confiance high si l'harmonique dominante depasse 3x la suivante (SPEC 3.1)
SPECTRUM_MAX_ROWS = 60  # tranches - hauteurs echantillonnees pour le spectre angulaire cellule par cellule
FFT_RATIO_MEDIUM = 1.5  # - - en deca de ce rapport la detection tombe en confiance low
AXIS_TOL = 0.05  # - - tolerance d'egalite relative des valeurs propres d'inertie (SPEC 2.1)
SYM_TOL = 0.02  # - - distance de Hausdorff relative admise apres rotation de 2*pi/N (SPEC 3.1)
AXIS_PROBE_GRID = 60  # cellules - carte grossiere servant a departager les axes candidats quand l'inertie doute
AXIS_PROBE_THETA = 360  # secteurs - resolution azimutale de cette carte grossiere
SYM_REJECT = 0.10  # - - Hausdorff relatif au-dela duquel la piece n'est pas une roue reguliere autour de l'axe : aucune performance publiee. Roues valides mesurees <= 0.018 (dont une roue privee de 5 % de ses triangles), roue flanquee d'un corps parasite ou deux roues cote a cote 0.29-0.30
AXIS_WARN_DEG = 5.0  # degres - ecart axe detecte / Z au-dela duquel on avertit et on realigne (SPEC 2.1)
SUCTION_ASYMMETRY_MIN = 0.05  # - - asymetrie radiale minimale de la veine pour trancher le cote aspiration
HAUSDORFF_SAMPLES = 4000  # points - echantillonnage du maillage pour la distance de Hausdorff
HAUSDORFF_CAP = 1.0  # - - ecretage de la distance de Hausdorff, rapporte au rayon exterieur : au-dela, la piece est de toute facon refusee et la recherche exacte ne coute que du temps
AXIS_PROBE_SAMPLES = 1000  # points - echantillonnage reduit pour departager les axes candidats
PERIOD_CANDIDATES = 5  # harmoniques - nombre de pics du spectre essayes quand la rotation de 2*pi/N ne superpose pas la piece
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
SPLINE_SMOOTH = 0.001  # - - lissage relatif : residu RMS vise = SPLINE_SMOOTH * amplitude de cambrure (SPEC 4.2)
SECTION_MARGIN = 0.02  # - - retrait relatif aux deux bouts de la plage [r_moyeu, r_tip] (SPEC 4.1)
BETA_CHORD_FRACTION = 0.10  # - - fraction de corde sur laquelle sont moyennes beta1 et beta2 (SPEC 4.3)
CAMBER_END_TRIM = 1.0  # epaisseurs - longueur ecartee a chaque bout de corde (zone des faces de bout)
CAMBER_TRIM_TAPER = 2.0  # - - vitesse d'ouverture d'epaisseur, en multiples de e_max/corde, marquant une face de bout
CAMBER_MIN_STATIONS = 6  # stations - nombre minimal conserve apres ecretage des bouts
SPLINE_LAMBDA_MIN = 1e-18  # - - borne basse de la recherche du parametre de penalisation
SPLINE_LAMBDA_MAX = 1e12  # - - borne haute de la recherche du parametre de penalisation
SPLINE_BISECTION_STEPS = 80  # - - iterations de bissection sur log(lambda)
MIN_PROFILE_POINTS = 8  # points - taille minimale d'une polyligne pour etre exploitee comme profil
MIN_BLADE_SECTIONS = 3  # coupes - nombre minimal de coupes exploitables pour une confiance > low
LOOP_DOUBLE_FRACTION = 0.50  # - - fraction des azimuts coupant l'aube deux fois au-dela de laquelle elle est une boucle
LOOP_FRINGE_ROWS = 2  # cellules - un troncon de pale d'au plus cette hauteur, colle a une cellule pleine, est la frange d'un moyeu ou d'un flasque bossele (maillage scanne), pas un brin d'aube
LOOP_MIN_STATIONS = 3  # stations - nombre minimal de rayons dedoubles pour conclure a une boucle
LOOP_MIN_SECTORS = 3  # secteurs - azimuts voyant de la pale, en deca desquels le rayon n'est pas exploite
NORMAL_WALL_MARGIN = 0.15  # - - part de la hauteur de veine ecartee contre chaque paroi, ou l'aube n'a que des chants
NORMAL_BAND_FRACTION = 0.10  # - - largeur relative des couronnes ou beta est moyenne, rapportee a r2-r1s
WRAP_SENSE_BANDS = 12  # couronnes - decoupage radial ou est suivie la derive azimutale de l'aube
WRAP_SENSE_MIN_BANDS = 3  # couronnes - minimum portant de la pale pour que la derive soit exploitable
WRAP_SENSE_MIN_DEG = 20.0  # deg - derive azimutale en deca de laquelle le sens d'enroulement n'est pas tranche
BETA_SENSITIVITY_DEG = 1.0  # deg - perturbation de beta2 servant a chiffrer la sensibilite de la hauteur
SENSITIVITY_DIAMETER_REL = 0.01  # - - perturbation relative du diametre servant a chiffrer l'elasticite (differences centrees)
SENSITIVITY_LOW_PER_DEG = 0.30  # - - sensibilite relative par degre au-dela de laquelle la grandeur est declassee en confiance faible
BETA_SENSITIVITY_ALERT = 0.09  # - - variation relative de hauteur PAR DEGRE de beta2 au-dela de laquelle la sensibilite est signalee ; vaut la moitie des 18 % d'incertitude annonces par le modele, qui portent sur deux degres d'ecart
SECTION_REFINE_STEPS = 6  # - - bissections de recalage d'un point de coupe sur la surface exacte
SECTION_MAX_WRAP = 5.235987755982989  # rad - enroulement maximal d'un profil de pale, 300 deg (au-dela : contour de revolution, moyeu ou flasque, qui fait le tour complet)
WRAP_CONSISTENCY_MIN = 0.40  # - - enroulement mesure / enroulement implique par beta, en deca duquel la cambrure se contredit
WRAP_CONSISTENCY_MAX = 2.50  # - - meme rapport, au-dela duquel elle se contredit aussi
TWIST_NOTE_DEG = 8.0  # degres - vrillage beta2 - beta1 au-dela duquel l'aplatissement de la lecture est signale a l'utilisateur
TWIST_RECOVERY_MIN = 0.60  # - - part basse du vrillage reel que la lecture restitue (mesure sur roues de synthese 20/35 a 40/65, voir tests/test_audit_geometrique.py)
TWIST_RECOVERY_MAX = 0.90  # - - part haute du vrillage reel que la lecture restitue (meme source)
CHORD_THICKNESS_MIN = 2.0  # - - rapport corde/epaisseur en deca duquel un profil est juge mal conditionne

# ---------------------------------------------------------------------------
# Classification de la roue
# ---------------------------------------------------------------------------
R_RATIO_AXIAL_MAX = 1.15  # - - r2/r1s en deca duquel la roue est axiale (SPEC 3.3)
R_RATIO_MIXED_MAX = 1.80  # - - r2/r1s en deca duquel la roue est mixte, au-dela centrifuge (SPEC 3.3)
NQ_CENTRIFUGAL_MAX = 35.0  # tr/min - vitesse specifique n_q maximale d'une roue centrifuge (SPEC 3.3)
NQ_MIXED_MAX = 80.0  # tr/min - n_q maximale d'une roue mixte, au-dela axiale (SPEC 3.3)
NQ_HEAD_EXPONENT = 0.75  # - - exposant de la hauteur dans n_q = n sqrt(Q) / H^0.75 (SPEC 3.3)
MIXED_DISCHARGE_DEG = 45.0  # degres - inclinaison meridienne de repli d'une roue mixte degeneree
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
C_F_CANAL = 0.005  # - - coefficient de frottement de peau du canal inter-aubes, valeur usuelle en pompe (Darcy-Weisbach)
DE_HALLER_MIN = 0.72  # - - rapport w2/w1 minimal avant decollement (critere de de Haller, grilles d'aubes)
XI_DIFFUSION = 2.0  # - - coefficient de la perte de diffusion sous le seuil de de Haller
ETA_COMPARAISON = 0.8101  # - - part des pertes que le modele 1D ne voit pas (volute, ecoulements secondaires, rugosite) ; calee pour qu'une roue centrifuge ordinaire retombe sur ETA_H x ETA_VOL x ETA_MEC
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
BETA_CLAMP_TOL = 0.05  # degres - distance a une borne en deca de laquelle un angle est declare ecrete, non mesure
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
# Propulsion : helice libre, element de pale et quantite de mouvement (BEM)
# ---------------------------------------------------------------------------
BEM_STATIONS = 24  # stations radiales du bilan par element de pale
BEM_MAX_ITERATIONS = 200  # bissections maximales de la resolution en phi
BEM_TOLERANCE = 1e-9  # rad - largeur d'intervalle sous laquelle la bissection sur phi s'arrete
BEM_PHI_MIN = 0.0020  # rad - borne basse de la bissection sur l'angle d'ecoulement
BEM_PHI_MAX = 1.5688  # rad - borne haute de la bissection, pi/2 moins deux milliradians
BEM_SWIRL_FLOOR = 0.05  # - - plancher de (1 + k') dans le bilan de couple : sous zero, la giration induite depasserait a l'envers la vitesse d'entrainement, ce qui n'a pas de sens
BEM_INDUCTION_MAX = 0.40  # - - facteur d'induction axial au-dela duquel la theorie de la quantite de mouvement cesse de valoir (etat de sillage turbulent, Glauert 1926)
BEM_ROOT_CUTOFF = 0.02  # - - fraction de rayon retiree en pied et en bout, ou la portance s'annule
CL_ALPHA = 6.10  # rad^-1 - pente de portance d'un profil mince reel (2 pi corrige de la viscosite, Abbott & von Doenhoff)
CL_STALL = 1.40  # - - portance maximale avant decrochage d'un profil d'helice usuel
CL_MIN = -0.80  # - - portance minimale (decrochage negatif, fonctionnement en moulinet)
CD_INDUCED_K = 0.020  # - - coefficient du terme quadratique de la polaire, cd = cd0 + k (cl - cl_min_drag)^2
CD_MIN_DRAG_CL = 0.20  # - - portance au minimum de trainee du profil
CD_BASE = 0.0090  # - - trainee de profil a portance nulle, aube usinee lisse (Re ~ 1e6)
CD_THICKNESS_K = 0.060  # - - part de trainee ajoutee par l'epaisseur relative, cd0 = CD_BASE + k (e/c)
CAMBER_ALPHA0_FACTOR = 2.0  # - - angle de portance nulle d'un arc de cercle : alpha_0 = -2 h/c (theorie des profils minces, Glauert 1926)
PRANDTL_LOSS_MIN = 0.05  # - - plancher du facteur de perte de Prandtl, evite la division par zero en bout de pale
FROUDE_MIN_SPEED = 0.10  # m/s - vitesse d'avance sous laquelle le rendement de Froude n'est pas defini
SPEED_OF_SOUND_AIR = 340.3  # m/s - celerite du son dans l'air standard a 15 C (atmosphere OACI)
RHO_AIR = 1.225  # kg/m3 - masse volumique de l'air standard au niveau de la mer (atmosphere OACI)
TIP_MACH_WARN = 0.78  # - - nombre de Mach en bout de pale au-dela duquel la compressibilite degrade la portance
PROPULSION_SPEED_MAX = 350.0  # m/s - vitesse d'avance au-dela de laquelle l'entree releve de la faute de frappe (transsonique dans l'air, impossible en eau)
ADVANCE_SWEEP_POINTS = 41  # points du balayage du parametre d'avance J
ADVANCE_SWEEP_MARGIN = 1.20  # - - borne haute du balayage, en fraction du J de poussee nulle

GENUS_DAMAGE_MAX = 0.01  # - - part d'aretes de bord au-dela de laquelle le genre topologique n'est plus lisible : chaque dechirure fabrique une anse, et un maillage troue a 5 % rend deja un genre de 185 la ou une piece a peine ouverte en rend 6
GENUS_PER_BLADE_MAX = 3.0  # - - genre maximal admis par pale pour qu'une lecture toroidale soit plausible : une boucle vaut une anse, trois laissent la marge des conges et des raccords

# ---------------------------------------------------------------------------
# Import par composants declares (SPEC v2)
# ---------------------------------------------------------------------------
SLICE_THICKNESS_MAX = 0.05  # m - epaisseur au-dela de laquelle un solide fluide n'est plus une tranche mince mais un tube, qui ne dit plus ou est la section de reference
SLICE_ANISOTROPY_MIN = 1.5  # - - rapport minimal entre le plus grand moment d'inertie (autour de la normale) et le suivant, pour qu'une tranche ait une normale isolee. Une tranche mince tend vers 2 sans l'atteindre (axes perpendiculaires : I_n = I_1 + I_2) ; 1.5 est un disque aussi epais que son rayon. Le seuil valait 2.0, qu'aucune tranche ne passait
PLANES_DISTINCT_MIN = 0.001  # m - distance minimale entre les centroides des solides d'entree et de sortie : en deca, le sens debitant n'est pas defini
BAND_INNER_MIN = 0.5  # - - rayon interieur minimal d'une bande cylindrique, rapporte a son rayon exterieur : en deca, le solide est un disque ou une couronne plate, pas une bande de refoulement radial
BAND_LATERAL_NZ = 0.5  # - - composante axiale maximale de la normale d'une facette laterale de bande (faces de la paroi, par opposition aux faces de dessus et de dessous)
RECENTRE_TOLERANCE = 0.01  # - - fraction du rayon exterieur en deca de laquelle un centroide colle a l'origine trahit une piece recentree a l'export (SPEC v2 6.1)
COAXIAL_TOLERANCE = 0.02  # - - ecart admis entre les centres, dans le plan XY, des pieces de revolution (entree, sortie, moyeu), rapporte au rayon exterieur : au-dela, elles ne sont pas sur le meme axe
CORNER_TOLERANCE = 0.0005  # m - deux pieces de tailles differentes dont les boites ont le meme coin a cette tolerance pres ont ete deplacees chacune de son cote (STLOUT d'AutoCAD exige l'octant positif)
OVERLAP_TOLERANCE = 0.05  # - - recouvrement admis entre boites englobantes de composants, en fraction du plus petit volume (SPEC v2 6.2)
BLADE_COPY_CLEARANCE = 0.0  # m - jeu minimal exige entre deux copies de pale voisines ; zero, seul le recoupement compte
DECLARED_BLADES_MIN = 2  # nombre de pales declarable minimal (SPEC v2 2)
DECLARED_BLADES_MAX = 24  # nombre de pales declarable maximal (SPEC v2 2)

# ---------------------------------------------------------------------------
# Domaine des entrees : bornes au-dela desquelles les modeles ne valent plus
# ---------------------------------------------------------------------------
RPM_MAX = 100000.0  # tr/min - regime au-dela duquel l'entree releve de la faute de frappe, pas de la pompe
GRID_MIN = 20  # cellules - grille minimale en r et en z sous laquelle plus rien ne se resout
N_THETA_MIN = 24  # secteurs - discretisation azimutale minimale : quatre points par pale a six pales
GRID_MAX = 2000  # cellules - grille maximale en r et en z : au-dela, memoire et duree explosent sans rien resoudre de plus (100 microns sur une roue de 20 cm)
N_THETA_MAX = 7200  # secteurs - discretisation azimutale maximale (0.05 degre) ; 1e8 secteurs epuisaient la memoire
TEMPERATURE_MIN = 1.0  # degres C - borne basse du domaine de la correlation d'Antoine
TEMPERATURE_MAX = 100.0  # degres C - borne haute du domaine d'Antoine et de la table de masse volumique
ALTITUDE_MAX = 11000.0  # m - plafond de la troposphere, borne haute du modele d'atmosphere OACI

# ---------------------------------------------------------------------------
# Installation par defaut
# ---------------------------------------------------------------------------
ALTITUDE = 0.0  # m - altitude du site (SPEC 6.2)
TEMPERATURE = 20.0  # degres C - temperature du liquide (SPEC 6.2)
HAUTEUR_ASPIRATION = 0.0  # m - hauteur d'aspiration, positif = en charge (SPEC 6.2)
PERTES_ASPIRATION = 0.5  # m - pertes de charge de la conduite d'aspiration (SPEC 6.2)
HAUTEUR_ASPIRATION_MIN = -30.0  # m - aspiration la plus haute acceptee : trois fois la colonne d'eau atmospherique, au-dela c'est une faute de saisie
HAUTEUR_ASPIRATION_MAX = 300.0  # m - charge a l'aspiration la plus haute acceptee (reservoir en hauteur) ; au-dela, faute de saisie
PERTES_ASPIRATION_MAX = 50.0  # m - pertes de charge d'aspiration les plus fortes acceptees ; negatives, elles creeraient de l'energie

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
VIEWER_MAX_FACES = 120000  # triangles - au-dela, le maillage est decime pour la vue 3D
VIEWER_CIRCLE_SEGMENTS = 128  # segments - finesse des cercles de reperage de la vue 3D
VIEWER_CLASS_DEPTH = 2.0  # cellules - profondeur de sondage sous une facette pour la classer

# ---------------------------------------------------------------------------
# Application locale (python -m impeller_analyzer.serve)
# ---------------------------------------------------------------------------
SERVER_HOST = "127.0.0.1"  # adresse d'ecoute par defaut : la machine locale seule
SERVER_PORT = 8765  # port d'ecoute par defaut
SERVER_MAX_UPLOAD = 250 * 1024 * 1024  # octets - taille maximale d'un fichier depose
SERVER_HISTORY = 8  # analyses conservees sur disque pour le telechargement des rapports
# ---------------------------------------------------------------------------
# Palette et typographie des sorties visuelles (une seule source, C6)
# ---------------------------------------------------------------------------
# Teal et violet ne sont pas decoratifs : ils codent mesure contre declare, et
# la convention tient partout, page comme figures.
RPM_SLIDER_MARGIN = 1.15  # - - marge au-dela de la vitesse limite couverte par le curseur de regime, pour que le franchissement soit atteignable
COULEUR_FOND = "#FBFBFA"  # fond de page et de figure
COULEUR_ENCRE = "#16232B"  # texte, axes, trace neutre
COULEUR_GRILLE = "#DDE3E0"  # filets de tableau, grille des graphiques
COULEUR_MESURE = "#1D6F6A"  # teal - grandeur mesuree, trace principal
COULEUR_DECLARE = "#5B4B8A"  # violet - grandeur declaree ou imposee
COULEUR_LIMITE = "#9B1D20"  # rouge - cavitation, confiance faible, franchissement
COULEUR_FOND_SOMBRE = "#12181C"  # fond du theme sombre, en bascule
COULEUR_ENCRE_SOMBRE = "#E8EDEA"  # encre du theme sombre
COULEUR_GRILLE_SOMBRE = "#2A353B"  # filets du theme sombre

PLOT_BAND_VIDE = 0.15  # - - position d'affichage du seuil F_VIDE sur l'echelle de couleur
PLOT_BAND_PLEIN = 0.85  # - - position d'affichage du seuil F_SOLIDE sur l'echelle de couleur

# ---------------------------------------------------------------------------
# Conversions d'unites, pour la presentation seule (les calculs restent en SI)
# ---------------------------------------------------------------------------
SECONDS_PER_HOUR = 3600.0  # s/h - conversion m3/s -> m3/h des tableaux et des courbes
SECONDS_PER_MINUTE = 60.0  # s/min - conversion tr/min -> tr/s
MM_PER_M = 1000.0  # mm/m - conversion des longueurs pour l'affichage
CM3_PER_M3 = 1.0e6  # cm3/m3 - conversion du volume pour le rapport d'import
W_PER_KW = 1000.0  # W/kW - conversion des puissances pour l'affichage

# ---------------------------------------------------------------------------
# Tolerances de la campagne de validation (SPEC phase 8)
# ---------------------------------------------------------------------------
VALID_GEOM_TOL = 0.02  # - - ecart relatif admis sur N, r1s, r2 (SPEC 8.1)
VALID_BETA_DEG = 2.0  # degres - ecart absolu admis sur beta1 et beta2 (SPEC 8.1, voir README)
VALID_REFERENCE_TOL = 0.18  # - - ecart relatif admis sur H au BEP du cas de reference (SPEC 8.2)
VALID_INVARIANCE_TOL = 0.005  # - - ecart relatif admis apres rotation/translation (SPEC 8.3)
VALID_DECIMATION_DEG = 3.0  # degres - ecart admis sur beta2 apres decimation a 20 % (SPEC 8.4)
VALID_BETA_NORMALS_DEG = 5.0  # degres - ecart absolu admis quand les angles sont lus sur les normales et non sur la cambrure : la methode est basse de 2 a 5 degres sur des roues d'angles connus, biais mesure et non corrige (voir blade_normals)
VALID_DECIMATION_RATIO = 0.20  # - - fraction de triangles conservee par la decimation (SPEC 8.4)
