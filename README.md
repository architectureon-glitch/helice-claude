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
| `--aspiration` | `auto` (défaut), `+z` ou `-z` : quel bout du maillage est le côté aspiration ; force la détection décrite plus bas |
| `--rotation` | `horaire` ou `antihoraire`, vu de +Z. **Sans lui le sens reste non renseigné** : la géométrie le suggère, elle ne le tranche pas |
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

### Quel bout du maillage est l'aspiration

La SPEC pose la convention « axe = Z, aspiration vers +Z », mais un fichier
sorti d'un logiciel de CAO n'a aucune raison de la respecter. Une roue exportée
à l'envers est alors lue depuis son refoulement : sa veine paraît converger vers
l'axe, l'outil la classe *axiale* et lui donne un sens de sortie faux — sur une
roue de pompe, qui aspire par le centre et refoule **latéralement**, l'erreur
est franche.

L'orientation se tranche donc sur la géométrie, sans rien demander. Le rayon
moyen des cellules de pales est calculé séparément dans la moitié haute et la
moitié basse de la zone, pondéré par `r dr dz`, et leur écart rapporté au rayon
extérieur donne une **asymétrie méridienne** : négative quand la veine part du
petit rayon (l'ouïe) vers le grand, c'est-à-dire quand la convention est
respectée. Au-delà de `SUCTION_ASYMMETRY_MIN` la conclusion est `high` ; si
l'asymétrie est positive, le maillage est retourné d'un demi-tour autour de X et
toute l'analyse reprend sur le maillage corrigé — sens de rotation et sens de
sortie sont donnés dans le repère corrigé, et un avertissement le signale.

Sur une hélice axiale la veine garde le même rayon d'un bout à l'autre :
l'asymétrie reste proche de zéro, la question ne se tranche pas, la convention
+Z est conservée et la confiance descend à `low`. `--aspiration +z|-z` impose
alors la réponse.

### Les aubes qui se referment sur elles-mêmes

Le modèle de ligne moyenne de la SPEC suppose une aube **simple** : une surface à
bord d'attaque et bord de fuite uniques, dont une coupe sur une surface de
courant donne un profil par pale. Une aube **toroïdale** est une boucle : elle
part du moyeu, sort, se retourne au bout et revient. Une coupe la traverse deux
fois, et l'appariement des deux faces dont sort la cambrure apparie alors la
face d'un brin avec celle de l'autre. β1, β2, le sens de rotation et toute
l'hydraulique qui en découle ne veulent rien dire.

La signature est topologique et se lit sur la carte d'occupation, sans rien
recouper. Elle ne se cherche pas en azimut mais **en hauteur** : les deux brins
d'une boucle sont au même azimut, séparés en z. À rayon et azimut fixés, une aube
simple donne un tronçon unique le long de z, une boucle en donne deux. La
grandeur mesurée est donc la fraction des azimuts où la coupe rencontre deux
tronçons, relevée rayon par rayon. Sur une hélice toroïdale réelle de Ø 335 elle
vaut **1,00 de r = 74 à 161 mm** puis retombe à 0 au-delà de 162, là où les brins
fusionnent ; sur une roue centrifuge fermée ordinaire elle plafonne à 0,23, sur
une hélice axiale elle est nulle. Le seuil est à 0,50, au milieu d'un fossé.

Reconnue, la boucle ne dégrade pas la confiance : elle **retire** les grandeurs
concernées. Le résumé console et le tableau des performances portent une mention
« non applicable », et β1, β2, sens de rotation, hauteur, débit, puissance,
couple, rendement et NPSHr passent en confiance basse. Axe, nombre d'aubes,
rayons, sections et volumes restent valables : ils ne passent pas par la
cambrure.

### Le sens de rotation est une entrée, pas un résultat

L'outil le mesurait et l'annonçait. Deux cas montrent qu'il n'en a pas le droit.
Quand le rapport `r2/r1s` tombe à cheval sur la frontière mixte / centrifuge, les
deux familles appliquent des règles **opposées** (`ω = +signe(k)` d'un côté,
`−signe(dθ/dr)` de l'autre) : le sens annoncé bascule pour un millième d'écart.
Et sur une aube quasi radiale, la lecture n'a plus aucune marge. Dans les deux
cas l'outil affichait un sens avec l'aplomb d'un résultat mesuré.

Il le demande désormais. `--rotation horaire|antihoraire` (ou le menu du même nom
dans l'application locale) le renseigne ; sans lui, la case du tableau porte
`à indiquer (--rotation)` et la confiance reste au plus bas. La lecture
géométrique n'est pas perdue pour autant : elle est conservée et **présentée comme
une suggestion**, dans le tableau, dans le résumé console et sur la vue 3D — où la
flèche s'affiche en gris tant qu'elle n'est pas confirmée. Un sens indiqué qui
contredit la suggestion est retenu quand même, avec un avertissement.

Corollaire du même principe : la réserve « la suggestion n'est pas fiable ici »,
émise quand le rapport `r2/r1s` tombe sur la frontière, ne sort plus dès lors que
le sens est fourni — elle demanderait de vérifier ce que l'utilisateur vient
d'affirmer.

Le calcul n'en dépend pas : le modèle de ligne moyenne ne connaît que `|ω|` et
les angles de pale, et les courbes sont bit à bit identiques dans les deux cas.
Ce que le sens change, c'est ce que l'outil **affirme** — et c'est vous qui avez
la pièce sous les yeux.

### Calculer une roue dont les aubes sont des boucles

Reconnue toroïdale, la roue n'est pas abandonnée : la lecture des angles bascule
sur une autre méthode. La cambrure n'a pas de réponse stable sur cette forme —
selon l'envergure où l'on coupe, on obtient le contour de la boucle entière
(aller et retour en un seul tour fermé), un seul brin, ou le bourrelet où les
deux brins fusionnent ; sur la roue de référence β2 vaut 8° près du plateau et
89° au milieu de la veine. Ce n'est pas un défaut de mise en œuvre, c'est la
limite du domaine du modèle.

`geometry/blade_normals.py` mesure autrement, localement, sans jamais apparier
deux faces. Une surface de pale ne contient pas sa propre normale : si l'aube
fait l'angle β avec la direction tangentielle, sa normale vaut
`−sin β · e_u + cos β · e_m`, d'où **tan β = |N_u| / |N_m|**, face par face,
pondéré par les aires. Intrados et extrados portent des normales opposées, mais
leurs deux composantes changent de signe ensemble : le *rapport* garde le sien,
ce qui donne le sens d'enroulement et donc le sens de rotation.

Restent à écarter les **chants**, ces bandes étroites où l'aube meurt contre le
moyeu et le flasque : elles ne portent aucun angle et tirent la moyenne vers le
bas. Le tri est géométrique, pas directionnel — on écarte les faces trop proches
des parois, la marge étant prise sur la plage *contiguë* qui contient la face et
non sur l'étendue totale de la colonne, sans quoi elle enjamberait le vide entre
les deux brins et les supprimerait tous les deux. Filtrer sur la direction de la
normale serait plus simple mais ne marche pas : sur une aube hélicoïdale la
surface porte elle-même une grande composante d'envergure, et le filtre qui
nettoie une roue centrifuge ordinaire supprime alors *toutes* les faces.

Portée et limites, mesurées sur des roues synthétiques d'angles imposés (β1/β2
de 15/20 à 40/65) : la lecture est **basse de 2 à 5 degrés**, le sens
d'enroulement est toujours juste, et le résultat ne bouge pas d'un degré entre
une grille de 100 et une de 150. Le biais n'est **pas corrigé** : il n'a pas été
expliqué, et le corriger d'après le seul générateur interne reviendrait à caler
l'instrument sur lui-même. β2 reste le point faible sur une boucle — il est
mesuré là où les deux brins fusionnent en un bout massif, qui bloque plus qu'il
ne guide. `--beta1` et `--beta2` court-circuitent toute la lecture.

Deux corrections d'appoint accompagnent ce chemin. Les profils d'une coupe sont
désormais **regroupés en familles** de N copies périodiques, triées par
enroulement : les mélanger donnait un angle qui ne décrivait aucune famille et
un sens de cambrure qui basculait d'une coupe à la suivante. Et le plafond
d'enroulement d'un profil passe de 180° à 300° — il est là pour rejeter les
contours de révolution, qui font le tour complet, et il rejetait au passage des
aubes réelles enroulées de 190°.

### Sur une roue fermée, l'entrée est le percement du flasque

Le rayon d'aspiration se lit normalement sur l'extrémité des pales au plan
d'entrée. Sur une roue **fermée** cette règle se retourne : les aubes courent
jusque sous le flasque avant, et les lire donne le rayon extérieur de la roue au
lieu de l'ouïe. C'est le **percement du flasque** qui fait l'entrée. L'ouïe est
donc cherchée au-dessus du bord d'attaque, côté aspiration : dans chaque rangée
on remonte vers l'axe depuis le bord extérieur tant qu'il y a de la matière, et
le bord intérieur de cette couronne est l'ouïe — le plus petit trouvé, c'est-à-dire
le col. Le test porte sur la continuité de la couronne, pas sur la vacuité de
l'ouïe : une roue fermée peut porter un bossage d'arbre en son centre, qui donne
alors `r1h`. Sur un flasque conique la couronne se referme dès le bord d'attaque
et les deux règles coïncident.

L'écart n'est pas marginal. Sur une roue de Ø 335 à 3 pales très enveloppées,
lire `r1s` sur les aubes donnait 160,9 mm au lieu de 92,7 : le rapport `r2/r1s`
tombait à 1,04, la roue était classée *axiale*, β2 valait 9,7°, la hauteur
d'Euler devenait négative et l'outil ne trouvait aucun point de fonctionnement.
Avec l'ouïe, le rapport vaut 1,80 et la sortie redevient radiale.

### Sans moyeu, la moyenne quadratique moyeu-carter n'a pas de sens

Pour une roue mixte ou axiale, la SPEC prend `r2 = √((r2s² + r2h²) / 2)` — le
rayon quadratique moyen d'une veine annulaire bordée par un moyeu. Sans moyeu au
plan de sortie, `r2h = 0` et la formule dégénère en `r2s / √2`, qui n'est pas un
rayon de refoulement mais un artefact. Sur une roue d'essai à aubes courant
jusqu'à 89,6 mm, elle donnait 39,5. `r2` est donc pris au bout des aubes dans ce
cas, avec un avertissement.

L'enjeu dépassait la définition. Sur l'hélice toroïdale de référence, le rapport
`r2/r1s` vaut **1,800** pour un seuil mixte/centrifuge à **1,80** : selon le pas
de grille, la roue basculait d'une famille à l'autre, et avec elle `r2` de 167 à
118 mm — un tiers sur `u2`, deux tiers sur la hauteur. Le rapport annonçait
sereinement deux résultats contradictoires pour le même fichier. La correction
stabilise `r2` à 166,7–167,1 mm sur des grilles de 100 à 180.

Reste que les deux familles appliquent des règles de **sens de rotation
opposées** (`ω = +signe(k)` pour l'axial et le mixte, `−signe(dθ/dr)` pour le
centrifuge). À cheval sur la frontière, le sens annoncé est un tirage au sort :
il est désormais signalé comme tel, et sa confiance forcée au plus bas.

### Comparer deux conceptions : ce que le rendement de la SPEC ne sait pas faire

Le rendement du tableau 2 ne juge pas la forme des aubes, et il ne l'a jamais
prétendu — mais il est affiché à côté de grandeurs calculées, ce qui prête à
confusion. La raison est structurelle : dans la SPEC, **chaque perte est une
fraction du point nominal**. Le frottement vaut `K_FROTTEMENT_REL` fois la
hauteur théorique, l'incidence est nulle au nominal par définition, fuite et
frottement de disque sont des fractions fixes. Passées à la moulinette, trois
roues centrifuges radicalement différentes — six aubes à 25°, trois à 65°, huit
à 15° — rendent 86,25 %, 86,18 % et 86,18 % de pertes brutes : **sept centièmes
de point d'écart**. Le modèle dimensionne une roue donnée ; il ne compare pas.

`hydraulics/losses.py` calcule les deux pertes qui, elles, dépendent vraiment de
la géométrie, et qui sont les deux critères classiques du dessin de roue :

- le **frottement de canal**, par Darcy-Weisbach sur la veine inter-aubes —
  longueur développée `L = (r2 − r1) / sin β` rapportée au diamètre hydraulique
  `4A/P` du canal. C'est lui qui pénalise un canal long, étroit, ou une aube dont
  la surface mouillée est doublée ;
- la **diffusion**, par le rapport `w2/w1`. Passé une certaine décélération
  relative la couche limite décolle ; le seuil retenu est celui de **de Haller**,
  `w2/w1 ≥ 0,72`.

Le rendement qui en sort est **calé** par `ETA_COMPARAISON` pour qu'une roue
centrifuge ordinaire retombe sur `ETA_H × ETA_VOL × ETA_MEC`. Ce facteur porte
tout ce que le modèle 1D ne voit pas — volute, écoulements secondaires,
rugosité. **C'est l'écart entre deux roues qui a un sens, pas la valeur absolue**,
et le rapport le dit à l'endroit où il l'affiche.

Deux propriétés le rendent utilisable comme critère de conception : il est
indépendant du **diamètre** et du **régime**. Vérifié sur les deux échelles d'une
même hélice réelle (Ø 335 et Ø 184, soit un rapport 0,55) et à deux régimes —
`L/Dh = 9,92`, `w2/w1 = 1,50` et 75,96 % dans les trois cas.

### La hauteur est-elle seulement calculable ?

`cu2 = u2 − cm2 / tan β2`. Aux petits angles la tangente varie très vite, et un
degré d'incertitude sur β2 — l'ordre de grandeur de ce que sait faire n'importe
quelle lecture géométrique — peut déplacer la hauteur bien au-delà des ±18 %
annoncés. Plutôt que de laisser croire à une précision qu'il n'a pas, l'outil
**mesure** cette sensibilité : il recalcule le point nominal à β2 ± 1° et compare.
Au-delà de 18 % d'écart, il le dit avec le chiffre ; si un degré suffit à faire
disparaître le point de fonctionnement, il le dit aussi. Sur l'hélice toroïdale,
β2 vaut 4° et l'écart atteint **89 %** — hauteur et puissance n'y sont que des
ordres de grandeur, alors que le débit et le NPSHr, qui n'en dépendent pas de la
même façon, restent stables à 1 % près.

Deux conséquences de cette mesure. La table de confiance **suit le texte** :
au-delà du seuil, `hauteur`, `puissance`, `couple` et `rendement` passent en
`faible`, tandis que `debit` et `npshr`, qui ne passent pas par `cu2`, gardent
leur niveau — annoncer « ordres de grandeur » et coter « moyenne » se
contredisait. Et `--beta1` / `--beta2` se posent **indépendamment** : imposer
seulement β2 laisse β1 à la lecture par les normales, au lieu de le faire
retomber sur la cambrure qui, sur une aube en boucle, donnait 87° là où les
normales en lisent 10 — et faisait disparaître le point de fonctionnement.

### Le rayon de sortie est celui des pales, pas de la matière

Une roue de pompe semi-ouverte est ouverte à l'avant et **fermée au dos** : les
pales sont portées par un disque arrière qui les déborde souvent de plusieurs
millimètres. Prendre `r2` sur le rayon de matière compte ce débord comme de la
pale — sur la roue d'essai, 104,65 mm au lieu de 89,95 mm, soit 16 % sur `u2` et
près de 35 % sur la hauteur. `r_tip` (matière) et `r_blade_tip` (pales) sont donc
deux grandeurs distinctes ; `r1s`, `r2s`, `b2` et le rapport `r2/r1s` se lisent
tous sur le masque de pales, et le rapport de matière n'est plus qu'une
indication reportée telle quelle.

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

196 tests, une phase par module. Ils passent aussi sous `pytest` si vous
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
│   ├── axis.py          # détection de l'axe, recentrage, côté aspiration
│   ├── blade_loops.py   # aubes en boucle fermée (type toroïdal)
│   ├── blade_normals.py # angles lus sur les normales, quand la cambrure ne s'applique pas
│   ├── occupancy.py     # carte f(r, z) — cœur du système
│   ├── topology.py      # pales, rayons, type de roue, sections
│   ├── proximity.py     # distance point-maillage, Hausdorff
│   ├── sections.py      # coupes de pale, profils déroulés
│   └── blade_angles.py  # cambrure, β1, β2, sens de rotation
└── hydraulics/
    ├── meanline.py      # Euler + glissement + pertes → H-Q
    ├── losses.py        # pertes de canal, pour comparer deux conceptions
    ├── cavitation.py    # NPSHr, NPSHa, vitesse maximale
    └── similarity.py    # lois de similitude, vitesse spécifique
```

## Limites

Modèle 1D ligne moyenne : **hauteur ±18 %, débit ±25 %, NPSHr ±30 %** — et le
cas de référence ci-dessus suggère que l'écart sur la hauteur peut être plus
grand encore sur une pompe réelle tant que le calage n'a pas été fait.

**À vérifier par essai sur banc avant toute décision d'achat ou de
dimensionnement.**
