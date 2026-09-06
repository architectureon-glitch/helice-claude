# impeller-analyzer

Analyse hydraulique d'une hélice ou d'une roue de pompe à eau à partir d'un
fichier 3D : géométrie extraite, sens de rotation requis, sens de sortie du
liquide, performances à plusieurs régimes, NPSH requis et vitesse maximale
avant cavitation.

Le cahier des charges complet est dans [`SPEC.md`](SPEC.md).

## Installation

Aucune installation n'est nécessaire : l'outil ne tient qu'à la bibliothèque
standard de **Python 3.11**.

```bash
git clone <ce dépôt> && cd helice
python3 -m impeller_analyzer examples/roue_centrifuge.stl --out rapport/
```

Les paquets de [`requirements.txt`](requirements.txt) sont **optionnels** : ils
élargissent les formats d'import (`trimesh` pour `.3ds`, `cadquery` pour
`.step`/`.iges`, `ezdxf` pour les `.dxf` à entités `MESH`), sans rien changer au
calcul. Les formats `.stl`, `.obj`, `.ply`, `.off` et les `.dxf` à `3DFACE` ou
`POLYFACE` sont lus par un parseur interne.

## Utilisation

### Dans le navigateur

```bash
python3 -m impeller_analyzer.serve
```

Le navigateur s'ouvre sur l'application : **déposez votre fichier 3D dans la
page**, réglez l'unité et les régimes, lancez l'analyse. La roue apparaît en 3D
et le rapport complet se télécharge depuis le panneau de droite.

Le calcul est en Python : il ne peut pas s'exécuter dans le navigateur. Plutôt
que d'entretenir une seconde implémentation de la physique en JavaScript — deux
versions à garder d'accord, dont une sans tests — l'application est servie par
un petit serveur HTTP de la bibliothèque standard qui appelle **exactement le
même code que la ligne de commande**. Il n'écoute que sur `127.0.0.1` : le
fichier déposé ne quitte pas votre poste. `--port` et `--host` si besoin.

### En ligne de commande

```bash
python -m impeller_analyzer roue.stl \
    --unit cm \
    --r-aspiration 4.5 \
    --rpm 1000 2000 3000 \
    --temperature 20 --altitude 0 --hauteur-aspiration 0 \
    --out rapport/
```

`python -m impeller_analyzer --help` détaille toutes les options. Les plus
utiles au quotidien :

| Option | Rôle |
|---|---|
| `--unit` | unité du fichier (défaut `cm`), ou un facteur vers le mètre |
| `--r-aspiration` | rayon d'aspiration imposé, en cm ; **prime toujours** sur la détection |
| `--blades`, `--beta1`, `--beta2` | repli manuel quand l'extraction est peu sûre |
| `--altitude`, `--temperature`, `--hauteur-aspiration`, `--pertes-aspiration` | hypothèses d'installation, qui fixent le NPSH disponible |
| `--grille NR NZ`, `--secteurs N` | finesse de la carte d'occupation |
| `--sans-controle-symetrie` | saute la vérification de périodicité par distance de Hausdorff (plus rapide) |

### Sorties

Dans le dossier `--out` :

- `resultats.json` — toutes les grandeurs, **en SI**, avec leur niveau de confiance ;
- `rapport.md` — les deux tableaux de la SPEC, l'encadré de vitesse maximale et
  l'encadré d'incertitude ;
- `courbes.png` — H–Q, rendement–Q et NPSHr–Q, une couleur par régime ;
- `carte_occupation.png` — la carte `f(r, z)`, l'outil de contrôle visuel de
  toute la chaîne : on y lit d'un coup d'œil le moyeu, la zone de pales et la
  veine fluide ;
- `vue3d.html` — la roue en 3D, manipulable.

### La vue 3D

Ouvrez `vue3d.html` dans n'importe quel navigateur : un seul fichier, aucune
dépendance, le maillage y est embarqué. Il fonctionne hors ligne et s'envoie par
courriel tel quel.

Le maillage y est **colorié selon ce que l'outil a compris** — acier pour le
moyeu et le flasque, sarcelle pour les pales — ce qui vérifie la segmentation
bien plus directement que la carte en coupe. S'y ajoutent les cercles `r1s`,
`r1h` et `r2` tracés à leurs vrais plans `z1` et `z2`, la flèche de sens de
sortie du liquide, et l'arc de sens de rotation.

Le bouton **Faire tourner** met la roue en rotation dans le sens calculé, et
**Vue de dessus (+Z)** vous place exactement au point de vue depuis lequel ce
sens est énoncé (« vu de +Z, côté aspiration ») : c'est là qu'on vérifie le
résultat d'un coup d'œil. **Demi-coupe** tranche la roue par un plan passant par
l'axe et découvre la veine méridienne.

La même vue sert d'interface à l'application locale
(`python3 -m impeller_analyzer.serve`) : elle s'y ouvre vide, en attente d'un
fichier, et se remplit sans rechargement à chaque analyse.

`--sans-vue3d` s'en passe si le fichier de sortie vous paraît trop lourd (il
pèse à peu près la taille du maillage : environ 1 Mo pour 17 000 triangles ;
au-delà de 120 000 triangles le maillage est décimé pour l'affichage seul, le
calcul restant fait sur le maillage complet).

## Comment ça marche

L'outil ne travaille pas directement sur le maillage mais sur une **carte
méridienne d'occupation angulaire** `f(r, z)` : pour chaque cellule d'une grille
(r, z), la fraction du tour occupée par de la matière. Le maillage est coupé
plan par plan, puis le nombre d'enroulement est accumulé depuis l'extérieur le
long de chaque rayon. Ce comptage signé — plutôt qu'une simple parité — traite
correctement une roue livrée comme **union de solides qui s'interpénètrent**
(moyeu + pales), cas le plus fréquent des exports CAO.

Tout le reste s'en déduit : `f ≥ 0.98` est du moyeu ou du flasque, `f ≤ 0.02` de
la veine fluide, l'entre-deux est la zone de pales.

## Écarts assumés par rapport à la SPEC

Cinq points où l'application littérale de la SPEC ne donnait pas le résultat
attendu. Chacun est commenté à l'endroit du code concerné.

1. **Critère de confiance du nombre de pales** (3.1). Une roue à N pales produit
   toujours des harmoniques fortes en 2N, 3N : ce sont des conséquences de la
   périodicité d'ordre N, pas des hypothèses concurrentes. Comparée à « la
   deuxième amplitude » du spectre, l'harmonique dominante ne dépasse jamais un
   rapport de 2 et la confiance ne serait jamais `high`. Le critère retenu la
   compare à la plus forte harmonique **qui n'est pas un multiple de N** ; le
   rapport brut à la deuxième amplitude reste calculé et reporté.

2. **Spectre calculé cellule par cellule** (3.1). Transformer le signal `g(θ)`
   déjà intégré sur toute la zone de pales fait disparaître l'harmonique N dès
   que les pales se recouvrent en projection : une hélice à 6 pales calée à 15°
   était comptée 12. L'analyse se fait donc avant sommation, les modules étant
   additionnés.

3. **Surface de coupe des pales** (4.1). La coupe cylindrique est la bonne pour
   une hélice axiale, mais dégénère sur une roue centrifuge : à rayon constant
   une aube radiale se réduit à un rectangle sans corde, et
   `tan β = dr / (r dθ)` — la formule que la SPEC donne elle-même pour ce cas
   (4.3) — n'est pas calculable puisque `r` y est constant. Les coupes se font
   donc sur la **surface de courant méridienne**, dont le cylindre est le cas
   particulier : `tan β = dm / (r dθ)` redonne `dz / (r dθ)` pour une courbe
   verticale et `dr / (r dθ)` pour une courbe horizontale.

4. **Appariement des deux faces d'un profil** (4.2). Le milieu de la coupe
   perpendiculaire à la corde n'est sur la cambrure que si les deux faces y sont
   parallèles. Vrai pour une pale d'épaisseur constante en z, faux pour une aube
   d'épaisseur **tangentielle** constante — la règle de dessin la plus courante
   en centrifuge, dont l'épaisseur perpendiculaire croît avec le rayon. Les
   faces sont donc appariées par abscisse curviligne normalisée, exact pour les
   deux conventions ; l'épaisseur reste mesurée par la coupe perpendiculaire de
   la SPEC.

5. **Définition du BEP** (5.3 et 5.4). La SPEC définit le BEP comme le maximum
   de `ρgQH / P_arbre` tout en posant `P_arbre = ρgQH / η_global` avec un
   `η_global` constant : le rendement serait constant et son maximum
   indéterminé. Les trois rendements sont donc traités comme des valeurs **au
   point de meilleur rendement**, portées par les deux mécanismes physiques
   correspondants — une fuite de recirculation à peu près constante pour
   `η_vol`, une puissance de frottement de disque à peu près constante pour
   `η_méc`. Le rendement dépend alors du débit, vaut exactement 0.8026 au BEP,
   et le maximum est bien défini.

Deux précisions de moindre portée :

- le paramètre de lissage `SPLINE_SMOOTH` est un résidu **relatif à l'amplitude
  de la cambrure**, échelle propre de la grandeur lissée : calé sur la corde, il
  rabattait la spline vers la droite et faussait de plusieurs degrés la pente
  aux deux stations où β1 et β2 sont justement mesurés ;
- la tolérance de validation sur β est exprimée en **degrés** (`VALID_BETA_DEG`)
  et non en pour-cent : 2 % n'a pas de sens sur un angle qui varie fortement le
  long de la pale.

## Validation

```bash
python -m impeller_analyzer.validation
python -m impeller_analyzer.validation --cas-de-reference examples/cas_de_reference.json
```

État des quatre contrôles de la phase 8 :

| Contrôle | Résultat |
|---|---|
| 8.1 géométries synthétiques | **conforme** — N exact, r1s et r2 à mieux que 0.4 %, β à 0.01° sur les hélices (15 à 45°) et à moins de 0.6° sur la roue centrifuge |
| 8.2 cas de référence | **à compléter** — voir ci-dessous |
| 8.3 invariance | **conforme** — rotation de 37° et translation : écart nul à la précision machine |
| 8.4 robustesse | **conforme** — maillage décimé à 20 % des triangles : β2 varie de 1.3° (limite 3°) |

### Cas de référence (8.2) : ce qui manque

Le mécanisme est en place (`examples/cas_de_reference.json` décrit la géométrie
et le point publié, `check_reference_case` mesure l'écart), mais **aucune courbe
constructeur n'a pu être vérifiée** : l'environnement de développement n'a pas
d'accès réseau. Le fichier livré est explicitement marqué « À REMPLACER ».

Sur ce jeu de valeurs provisoire (pompe normalisée EN 733 / ISO 2858 calibre
32-200 à 2900 tr/min, géométrie de roue typique), le modèle donne :

- **débit au BEP : 12.4 m³/h** contre 12.5 m³/h attendu — écart de 1 % ;
- **hauteur au BEP : 68.6 m** contre 50 m attendu — écart de **+37 %**, au-delà
  des 18 % visés.

Le débit tombe juste, la hauteur non : le modèle 1D avec la perte de frottement
de la SPEC (`K_FROTTEMENT_REL = 0.06`, soit 94 % de la hauteur d'Euler
conservée) est trop optimiste pour une petite pompe réelle, dont le rendement
hydraulique est plutôt de 0.70 à 0.80. C'est exactement ce à quoi sert
`config.py` : porter `K_FROTTEMENT_REL` autour de **0.30** ramènerait ce cas
dans les 18 %. Le calage n'a **pas** été appliqué, faute de données vérifiées :
il ne faut pas caler un modèle sur un chiffre non sourcé.

**À faire pour compléter la phase 8.2** : relever la courbe d'une pompe dont la
fiche technique est disponible, ainsi que les cotes de sa roue, les mettre dans
`examples/cas_de_reference.json`, relancer la validation, puis recaler
`K_FROTTEMENT_REL` (et au besoin `ETA_H`) sur cet écart.

## Calage

Toutes les constantes physiques, tous les coefficients empiriques et tous les
seuils numériques sont dans [`impeller_analyzer/config.py`](impeller_analyzer/config.py),
une par ligne, commentée avec sa source. Aucune valeur numérique significative
n'apparaît ailleurs. Pour recaler l'outil, on ne modifie que ce fichier.

## Tests

```bash
python -m unittest discover -s tests -t tests
```

149 tests, une phase par module. Ils passent aussi sous `pytest` si vous
l'avez : ce sont des `unittest.TestCase`.

## Architecture

```
impeller_analyzer/
├── config.py            # toutes les constantes
├── confidence.py        # niveaux high/medium/low et leur propagation
├── mesh.py              # maillage triangulaire unifié, en mètres
├── numeric.py           # Jacobi, DFT, spline de lissage, solveur bande
├── synthetic.py         # géométries à réponse analytique connue
├── analysis.py          # enchaînement des phases 1 à 6
├── validation.py        # campagne de la phase 8
├── cli.py
├── serve.py             # application locale : dépôt de fichier dans le navigateur
├── io/
│   ├── loader.py        # lecture multi-format → maillage unifié
│   ├── repair.py        # orientation des normales, rebouchage
│   ├── writer.py        # écriture STL/OBJ/PLY/OFF/DXF
│   ├── plot.py          # rendu PNG par zlib seul
│   ├── viewer.py        # vue 3D WebGL autonome, en un fichier
│   └── report.py        # resultats.json, rapport.md, courbes.png
├── geometry/
│   ├── axis.py          # détection de l'axe, recentrage
│   ├── occupancy.py     # carte f(r, z) — cœur du système
│   ├── topology.py      # pales, rayons, type de roue, sections
│   ├── proximity.py     # distance point-maillage, Hausdorff
│   ├── sections.py      # coupes de pale, profils déroulés
│   └── blade_angles.py  # cambrure, β1, β2, sens de rotation
└── hydraulics/
    ├── meanline.py      # Euler + glissement + pertes → H-Q
    ├── cavitation.py    # NPSHr, NPSHa, vitesse maximale
    └── similarity.py    # lois de similitude, vitesse spécifique
```

## Limites

Modèle 1D ligne moyenne : **hauteur ±18 %, débit ±25 %, NPSHr ±30 %** — et le
cas de référence ci-dessus suggère que l'écart sur la hauteur peut être plus
grand encore sur une pompe réelle tant que le calage n'a pas été fait.

**À vérifier par essai sur banc avant toute décision d'achat ou de
dimensionnement.**
