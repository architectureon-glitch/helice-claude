# SPEC — `impeller-analyzer`

## 0. Objectif

Outil en ligne de commande (Python) qui :

1. importe un fichier 3D d'hélice / roue de pompe à eau (unité d'import : **cm**),
2. extrait automatiquement la géométrie hydraulique (axe, nombre de pales, moyeu, rayons, angles de pale),
3. détermine le **sens de rotation requis** et le **sens de sortie du liquide**,
4. calcule les performances à **1000, 2000 et 3000 tr/min**,
5. calcule le **NPSH requis** et la **vitesse de rotation maximale** avant cavitation.

Convention d'orientation imposée : **axe = Z, aspiration vers +Z**. Le fluide entre par le haut
et est refoulé vers −Z (axial) ou radialement (centrifuge).

## 1. Contraintes de développement

- Python 3.11, aucune dépendance propriétaire.
- **Toutes les constantes physiques et coefficients empiriques dans un seul fichier `config.py`**,
  une constante par ligne, commentée, avec sa source. Aucune valeur codée en dur ailleurs.
- Développement **test-driven** : chaque phase livre ses tests avant d'être déclarée finie.
- Pas de refonte globale entre les phases : on ajoute des modules, on ne réécrit pas les précédents.
- Chaque grandeur calculée est stockée en **SI** (m, m³/s, rad/s, Pa). La conversion cm→m se fait
  une seule fois, à l'import.
- Toute grandeur estimée porte un **niveau de confiance** (`high` / `medium` / `low`) propagé
  jusqu'au rapport final.

## 2. Architecture

```
impeller_analyzer/
├── config.py            # toutes les constantes (§9)
├── io/
│   ├── loader.py        # lecture multi-format → maillage triangulaire unifié
│   └── report.py        # sortie JSON + Markdown + PNG
├── geometry/
│   ├── axis.py          # détection de l'axe et recentrage
│   ├── occupancy.py     # carte méridienne d'occupation angulaire  ← cœur du système
│   ├── topology.py      # moyeu, carter, rayons, nombre de pales, type de roue
│   ├── sections.py      # coupes cylindriques → profils déroulés
│   └── blade_angles.py  # ligne de cambrure, β1, β2, corde, épaisseur
├── hydraulics/
│   ├── meanline.py      # Euler + glissement + pertes → courbe H-Q
│   ├── cavitation.py    # NPSHr, NPSHa, vitesse max
│   └── similarity.py    # lois de similitude, vitesse spécifique
├── cli.py
└── tests/
```

---

## Phase 1 — Import géométrique

**Formats à supporter, par ordre de priorité :**

| Format | Bibliothèque | Priorité |
|---|---|---|
| `.stl`, `.obj`, `.ply`, `.3ds`, `.off` | `trimesh` | 1 — implémenter d'abord |
| `.step`, `.stp`, `.iges` | `cadquery` (OCC) → tessellation, tolérance 0.2 mm | 2 |
| `.dxf` (maillages POLYFACE / 3DSOLID tesselé) | `ezdxf` | 3 |
| `.dwg` | conversion externe via ODA File Converter → `.dxf`, puis chemin ci-dessus | 4 |

**`.lisp` n'est pas un format géométrique** : si l'utilisateur en fournit un, refuser avec un message
explicite lui demandant d'exporter en STL depuis son logiciel CAO.

Traitements à l'import :

1. Conversion d'unité : facteur par défaut **0.01** (cm → m), surchargeable par `--unit`.
2. Fusion des sommets dupliqués, tolérance 1e-6 m.
3. Réparation : `trimesh.repair.fix_normals`, `fill_holes`. Si le maillage reste non étanche,
   continuer mais marquer `watertight=False` et forcer la confiance à `medium` sur tous les volumes.
4. Rapport d'import : nombre de triangles, bbox en mm, volume, étanchéité.

**Test :** un cube de 10 cm importé doit donner un volume de 1.0e-3 m³.

---

## Phase 2 — Axe, recentrage, carte d'occupation

### 2.1 Détection de l'axe

Même si la convention impose Z, l'outil **vérifie** :

1. Calculer le tenseur d'inertie du maillage, diagonaliser.
2. Une roue à N pales (N ≥ 3) possède deux valeurs propres quasi égales ; l'axe de révolution est
   le vecteur propre dont la valeur propre est isolée. Critère : `|λa − λb| / max(λ) < 0.05`.
3. Comparer l'axe trouvé à Z. Écart > 5° → avertissement, et l'outil réaligne le maillage sur Z.
4. Recentrer : origine sur l'axe, à la hauteur du barycentre.

Sens de l'aspiration : imposé à +Z par la convention. Vérification de cohérence en 3.3.

### 2.2 Carte d'occupation angulaire — structure centrale

Construire une grille méridienne `(r, z)` de **200 × 200 cellules** sur la bbox cylindrique.
Pour chaque cellule, calculer la **fraction angulaire occupée par la matière**
`f(r, z) ∈ [0, 1]`, par lancer de 720 rayons en θ et test d'inclusion dans le maillage.

Cette carte unique donne tout le reste :

- `f ≥ 0.98` → **moyeu** (ou flasque, si la zone est extérieure aux pales)
- `0.02 < f < 0.98` → **zone de pales**
- `f ≤ 0.02` → **veine fluide ou extérieur**

C'est la structure sur laquelle s'appuient toutes les phases suivantes. À optimiser
(vectorisation numpy + `trimesh.proximity`), c'est le point chaud en temps de calcul.

**Test :** un cylindre plein doit donner f = 1 partout dans son emprise, 0 ailleurs.

---

## Phase 3 — Topologie de la roue

### 3.1 Nombre de pales N

Extraire le signal d'occupation angulaire `g(θ)` intégré sur toute la zone de pales
(720 échantillons). Appliquer une FFT. **N = indice de l'harmonique dominante.**
Chercher N dans [2, 12]. Confiance `high` si l'amplitude de l'harmonique dominante
dépasse **3×** celle de la deuxième.

Contrôle croisé : rotation du maillage de 2π/N et distance de Hausdorff au maillage
d'origine < 2 % du rayon extérieur.

### 3.2 Rayons caractéristiques

- `r_moyeu(z)` = plus grand r tel que `f(r,z) ≥ 0.98` (côté intérieur)
- `r_tip` = plus grand r tel que `f > 0.02`
- Roue **fermée** (avec flasque avant) si `f ≥ 0.98` réapparaît au-delà de la zone de pales ;
  sinon roue **ouverte / semi-ouverte**.
- Plan d'aspiration `z_1` = z max de la zone de pales (bord d'attaque).
  - `r_1s` = rayon extérieur des pales en `z_1`
  - `r_1h` = `r_moyeu(z_1)`
  - `r_1` = √((r_1s² + r_1h²)/2)  ← rayon quadratique moyen, à utiliser pour u1
- **Rayon d'aspiration** : `r_asp = r_1s`. Si la détection échoue (confiance `low`) ou si
  `--r-aspiration` est fourni en cm, la valeur utilisateur **prime toujours**.

### 3.3 Classification du type de roue

Calculer `R = r_2 / r_1s` où `r_2` = rayon en sortie (voir 3.4) :

| Condition | Type | Modèle hydraulique |
|---|---|---|
| R < 1.15 | axiale | déviation en cascade |
| 1.15 ≤ R < 1.80 | mixte (hélico-centrifuge) | Euler, ligne de courant méridienne |
| R ≥ 1.80 | centrifuge | Euler radial |

Contrôle croisé obligatoire après la phase 5 : la vitesse spécifique
`n_q = n·√Q / H^0.75` (n en tr/min, Q en m³/s, H en m) doit tomber dans
`n_q < 35` (centrifuge), `35–80` (mixte), `> 80` (axiale). **Incohérence entre les deux
classifications → avertissement en tête de rapport.**

### 3.4 Section de sortie

- Axiale/mixte : plan `z_2` = z min de la zone de pales ; `A_2 = π(r_2s² − r_2h²)·τ2`
- Centrifuge : cylindre `r_2 = r_tip` ; largeur `b_2` = hauteur en z de la zone de pales
  à `r = 0.98·r_tip` ; `A_2 = 2π·r_2·b_2·τ2`

---

## Phase 4 — Angles de pale

### 4.1 Coupes

**11 rayons** régulièrement répartis de `r_moyeu + 0.02·(r_tip − r_moyeu)` à
`r_tip − 0.02·(r_tip − r_moyeu)`. Pour chaque rayon, intersecter le maillage avec le cylindre
correspondant → polylignes fermées, une par pale.

Dérouler : `x = r·θ`, `y = z`. On obtient N profils de cascade dans un plan 2D.

### 4.2 Ligne de cambrure

Pour un profil déroulé :

1. Bord d'attaque / bord de fuite = les deux points les plus éloignés du profil (corde maximale).
   Le bord d'attaque est celui de **z le plus grand** (côté aspiration).
2. Découper la corde en **25 stations**. À chaque station, prendre le milieu du segment
   intrados/extrados perpendiculaire à la corde.
3. Lisser la ligne de cambrure par spline cubique, paramètre de lissage 0.001.

### 4.3 Angles

Angles mesurés **depuis la direction tangentielle** (convention pompe) :

- Machine axiale/mixte : `tan β = dz / (r·dθ)` le long de la cambrure
- Machine centrifuge : `tan β = dr / (r·dθ)`

Sortir :

- `β1(r)` : angle sur les 10 premiers % de la corde
- `β2(r)` : angle sur les 10 derniers % de la corde
- corde, épaisseur max, pas hélicoïdal `p = 2πr·tan β`
- Pour le modèle 1D, retenir β1 et β2 au **rayon quadratique moyen**.

**Repli manuel** : si la confiance d'extraction est `low`, l'outil demande β1, β2 et N en entrée
directe (`--beta1`, `--beta2`, `--blades`) et le reste du calcul continue normalement.

### 4.4 Sens de rotation et sens de sortie ← exigence explicite

**Cas axial / mixte.** La surface moyenne de pale est localement hélicoïdale : `z ≈ z₀ + k·θ`,
avec `k = ∂z/∂θ` mesuré sur la ligne de cambrure au rayon moyen. Une surface hélicoïdale
tournant à ω se translate axialement à `v_z = −k·ω`. On veut refouler vers −Z, donc `v_z < 0`,
donc :

> **signe(ω) = signe(k)**

Si `k > 0` : rotation **anti-horaire vue de dessus** (depuis +Z, côté aspiration).
Si `k < 0` : rotation **horaire vue de dessus**. Confiance `high` si `|k|` est significatif
(pas hélicoïdal < 8× le rayon).

**Cas centrifuge.** Les aubes de pompe sont quasi systématiquement **incurvées vers l'arrière** :
l'aube fuit le sens de rotation quand r augmente. Donc :

> **signe(ω) = −signe(dθ/dr)** le long de la cambrure

Confiance `low` si `β2 > 85°` (aubes radiales, sens ambigu) ou si `dθ/dr` change de signe.
Signaler explicitement dans le rapport que la règle s'inverse pour des aubes incurvées vers
l'avant (rare en pompe, courant en ventilateur à cage d'écureuil).

**Sens de sortie du liquide** à reporter :

- composante méridienne : axiale (−Z) pour une roue axiale, radiale (+r) pour une centrifuge,
  angle `atan(vitesse_axiale / vitesse_radiale)` pour une mixte ;
- composante tangentielle : **toujours dans le sens de rotation**, valeur `cu2` calculée en phase 5 ;
- angle absolu de sortie `α2 = atan(cm2 / cu2)` en degrés.

**Test :** générer une hélice synthétique à 4 pales planes calées à 30°, moyeu 30 mm, tip 100 mm.
L'extraction doit retrouver N = 4, β = 30° ± 1° à tous les rayons, et le sens de rotation attendu.

---

## Phase 5 — Modèle hydraulique (ligne moyenne)

### 5.1 Point de fonctionnement nominal (incidence nulle)

```
ω    = 2π·n/60
u1   = ω·r1
A1   = π(r1s² − r1h²)·τ1
cm1  = u1·tan(β1)
Q_n  = cm1·A1
```

### 5.2 Hauteur — équation d'Euler avec glissement

```
u2   = ω·r2
cm2  = Q / A2
σ    = 1 − √(sin β2) / N^0.7                    (Wiesner)
cu2  = σ·u2 − cm2 / tan(β2)
H_th = u2·cu2 / g
```

### 5.3 Pertes et courbe H–Q

```
h_frottement = k_f · Q²                  avec k_f calé pour donner 0.06·H_th au point nominal
h_incidence  = ξ · (w1 − w1_nominal)² / (2g)     ξ = 0.75
H(Q)         = H_th(Q) − h_frottement − h_incidence
```

Balayage : **29 points** de Q = 0 à Q = 1.40·Q_n.
Le BEP est le maximum de `η(Q) = ρ·g·Q·H / P_arbre`.

### 5.4 Puissances et rendements

```
P_hydraulique = ρ·g·Q·H
η_global      = η_h · η_vol · η_méc = 0.88 · 0.96 · 0.95 = 0.803
P_arbre       = P_hydraulique / η_global
Couple        = P_arbre / ω
```

### 5.5 Similitude

Fournir aussi les lois de similitude comme **contrôle de cohérence** :
Q ∝ n, H ∝ n², P ∝ n³. L'écart entre le calcul direct à 3000 tr/min et
l'extrapolation depuis 1000 tr/min doit être < 2 %. Sinon, il y a un bug.

---

## Phase 6 — Cavitation

### 6.1 NPSH requis — deux méthodes, on retient la plus pénalisante

**Méthode A (cinématique, au carter d'entrée — point le plus critique) :**

```
u1s   = ω·r1s
w1s   = √(cm1² + u1s²)
NPSHr = λc·cm1²/(2g) + λw·w1s²/(2g)          λc = 1.10 ; λw = 0.28
```

**Méthode B (vitesse spécifique d'aspiration) :**

```
NPSHr = (n·√Q / n_ss)^(4/3)                   n_ss = 180  (n en tr/min, Q en m³/s)
```

`NPSHr_retenu = max(A, B)`. Reporter les deux valeurs dans le rapport.

### 6.2 NPSH disponible — dépend de l'installation, pas de l'hélice

```
NPSHa = (p_atm − p_vap)/(ρ·g) + z_aspiration − h_pertes
```

Entrées CLI avec valeurs par défaut :

- `--altitude` 0 m → `p_atm = 101325·(1 − 2.25577e-5·altitude)^5.2559`
- `--temperature` 20 °C → `p_vap` par la corrélation d'Antoine, `ρ` par table
- `--hauteur-aspiration` 0 m (en charge ; négatif si la pompe aspire au-dessus du plan d'eau)
- `--pertes-aspiration` 0.5 m

Valeur par défaut résultante : **NPSHa = 9.61 m**.

### 6.3 Vitesse maximale admissible

Deux limites, on retient la plus basse :

**Limite 1 — NPSH.** `NPSHr ∝ n²` à coefficient de débit constant, donc :

```
n_max_npsh = n_ref · √( NPSHa / (marge · NPSHr(n_ref)) )      marge = 1.30
```

**Limite 2 — vitesse relative en entrée.** `w1s ≤ 30 m/s` (eau, roue métallique) :

```
n_max_w = n_ref · 30 / w1s(n_ref)
```

Sortie : **une seule valeur en tr/min**, entière, arrondie à la dizaine inférieure,
avec mention de la limite active et de l'hypothèse d'installation utilisée.

---

## Phase 7 — Interface et sortie

```bash
python -m impeller_analyzer roue.stl \
    --unit cm \
    --r-aspiration 4.5 \
    --rpm 1000 2000 3000 \
    --temperature 20 --altitude 0 --hauteur-aspiration 0 \
    --out rapport/
```

Produire trois choses :

1. `resultats.json` — toutes les grandeurs, en SI, avec les niveaux de confiance.
2. `rapport.md` — lisible, avec les tableaux ci-dessous.
3. `courbes.png` — H–Q, η–Q, NPSHr–Q, une couleur par régime.

**Tableau 1 — Géométrie extraite**

| Grandeur | Valeur | Confiance |
|---|---|---|
| Type de roue | axiale / mixte / centrifuge | |
| Nombre de pales | | |
| Rayon d'aspiration r1s | mm | |
| Rayon moyeu r1h | mm | |
| Rayon de sortie r2 | mm | |
| β1 / β2 au rayon moyen | ° | |
| Sens de rotation | horaire / anti-horaire vu du côté aspiration | |
| Sens de sortie du liquide | axial −Z / radial / mixte, α2 = ...° | |

**Tableau 2 — Performances**

| | 1000 tr/min | 2000 tr/min | 3000 tr/min |
|---|---|---|---|
| Débit nominal (m³/h) | | | |
| Hauteur (m) | | | |
| Puissance arbre (kW) | | | |
| Couple (N·m) | | | |
| Rendement estimé (%) | | | |
| Vitesse spécifique n_q | | | |
| NPSH requis (m) | | | |
| NPSH disponible (m) | | | |
| Marge NPSHa/NPSHr | | | |

**Encadré final :** vitesse maximale sans cavitation = **____ tr/min**, limite active = ____,
hypothèses d'installation = ____.

**Encadré d'incertitude, à imprimer systématiquement :**
modèle 1D ligne moyenne — hauteur ±18 %, débit ±25 %, NPSHr ±30 %. À vérifier par essai
sur banc avant toute décision d'achat ou de dimensionnement.

---

## Phase 8 — Validation

1. **Géométries synthétiques** générées par script, dont la réponse est connue analytiquement :
   hélice à pales planes, roue centrifuge à aubes en arc de cercle. Vérifier N, β1, β2, r1s, r2
   à moins de 2 % de la valeur théorique.
2. **Cas de référence** : une pompe dont la courbe constructeur est publiée. Écart sur H au BEP
   attendu sous 18 %.
3. **Invariance** : le même maillage tourné de 37° autour de Z et translaté doit donner des
   résultats identiques à 0.5 % près.
4. **Robustesse** : maillage décimé à 20 % des triangles → écart sur β2 sous 3°.

---

## 9. `config.py` — valeurs initiales

Toutes modifiables une par une, c'est le fichier de calage de l'outil.

```python
# Fluide (eau, 20 °C)
RHO = 998.2            # kg/m³
P_VAP = 2339.0         # Pa
G = 9.80665            # m/s²

# Import
UNIT_FACTOR = 0.01     # cm → m
MERGE_TOL = 1e-6       # m
STEP_TESSELLATION = 2e-4   # m

# Grille d'occupation
GRID_NR = 200
GRID_NZ = 200
N_THETA = 720
F_SOLIDE = 0.98        # seuil moyeu/flasque
F_VIDE = 0.02          # seuil veine fluide

# Détection
BLADES_MIN = 2
BLADES_MAX = 12
FFT_RATIO_MIN = 3.0    # confiance high si harmonique dominante > 3× la suivante
AXIS_TOL = 0.05        # tolérance d'égalité des valeurs propres
SYM_TOL = 0.02         # Hausdorff relatif après rotation de 2π/N

# Sections et cambrure
N_SECTIONS = 11
N_STATIONS = 25
SPLINE_SMOOTH = 0.001

# Obstruction
TAU_1 = 0.90           # entrée
TAU_2 = 0.92           # sortie

# Hydraulique
ETA_H = 0.88
ETA_VOL = 0.96
ETA_MEC = 0.95
K_FROTTEMENT_REL = 0.06    # perte de frottement au nominal, en fraction de H_th
XI_INCIDENCE = 0.75
Q_SWEEP_MAX = 1.40
Q_SWEEP_POINTS = 29

# Cavitation
LAMBDA_C = 1.10
LAMBDA_W = 0.28
N_SS = 180.0               # vitesse spécifique d'aspiration
MARGE_NPSH = 1.30
W1S_MAX = 30.0             # m/s

# Installation par défaut
ALTITUDE = 0.0             # m
TEMPERATURE = 20.0         # °C
HAUTEUR_ASPIRATION = 0.0   # m (positif = en charge)
PERTES_ASPIRATION = 0.5    # m
```

---

## 10. Ordre d'exécution demandé à Claude Code

Ne pas tout construire d'un coup. Livrer et faire valider dans cet ordre :

1. Phase 1 + tests d'import (STL/OBJ seulement)
2. Phase 2 (axe + carte d'occupation) + visualisation PNG de la carte `f(r,z)` — **indispensable
   pour déboguer visuellement tout le reste**
3. Phase 3 (topologie)
4. Phase 4 (angles + sens de rotation)
5. Phase 5 (hydraulique)
6. Phase 6 (cavitation)
7. Phase 7 (CLI + rapport)
8. Phase 8 (validation) puis élargissement des formats d'import (STEP, DXF)
