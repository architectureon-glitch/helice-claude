# impeller-analyzer

Étude des hélices et roues de pompe **toroïdales** — aubes en boucle fermée — à
partir d'un fichier 3D : géométrie des boucles, sens de rotation, performances
(pompe : courbes hauteur-débit, NPSH requis, vitesse maximale avant cavitation ;
hélice libre : poussée et rendement propulsif), et pour chaque résultat le
**comparatif de la même géométrie à aubes normales**, avec ce que le modèle
compte de l'écart et ce qu'il n'en compte pas.

Une pièce qui n'est pas toroïdale est refusée, avec la raison. Le calcul d'une
hélice normale reste dans le moteur : il sert de référence au comparatif.

Le cahier des charges complet est dans [`SPEC.md`](SPEC.md) ; son extension,
l'import par composants déclarés, dans [`SPEC_V2.md`](SPEC_V2.md).

## Installation

Aucune installation n'est nécessaire : l'outil ne tient qu'à la bibliothèque
standard de **Python 3.11**.

```bash
git clone https://github.com/architectureon-glitch/helice-claude && cd helice-claude
python3 -m impeller_analyzer examples/helice_toroidale.stl --machine helice_libre \
    --vitesse-avance 2 --rpm 1000 --out rapport/
```

Sous Windows, la commande s'appelle `py` (ou `python`) plutôt que `python3`, et
l'on peut aussi télécharger le dépôt en ZIP depuis GitHub (bouton **Code →
Download ZIP**) au lieu de le cloner.

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
et le rapport complet se télécharge depuis le panneau de droite. Le champ
**Forme des aubes** laisse la lecture trancher ; « déclarée toroïdale » poursuit
l'analyse quand la lecture ne voit pas de boucle, comme `--topologie-pale`.

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

## Si l'hélice était normale : le comparatif

Chaque résultat vient avec son jumeau : la **même géométrie** — mêmes rayons,
mêmes sections, mêmes angles, même nombre de pales, même régime — dont chaque
aube serait un seul brin au lieu d'une boucle. Seule la forme de l'aube change ;
l'écart entre les deux colonnes est ce que la boucle apporte ou coûte, **dans la
mesure où le modèle le voit**, et le rapport dit à chaque fois ce qu'il compte et
ce qu'il ne compte pas.

**Pompe.** Une ligne moyenne ne connaît pas la forme de l'aube : à géométrie
égale, elle rend la même courbe, et le rapport le dit plutôt que d'afficher un
écart fabriqué. L'écart se lit sur les **pertes de canal** : une aube en boucle
présente ses deux brins au frottement, pour le même canal — surface mouillée,
perte de frottement, rendement de comparaison. Sur une roue **fermée**, le
flasque supprime déjà le tourbillon et la fuite en bout de pale, pour les deux
roues : la boucle n'y apporte rien, et le rapport le dit. Sur une roue
**ouverte**, la normale fuit par-dessus ses aubes, la boucle non ; l'écart,
favorable à la toroïdale, dépend du jeu au carter et n'est pas chiffré.

Sur une roue **en série** comme hel1, la jumelle normale n'a pas le brin amont :
elle aspire par le bord d'attaque de son aube. Le comparatif chiffre alors
l'incidence de chacune au débit du cas, leurs hauteurs d'Euler, et dit quand la
normale n'a pas de point de fonctionnement. Sur hel1 : son bord d'attaque, à 47°,
est adapté à 101 m³/h, débit auquel son bord de fuite, à 23°, ne donne plus de
hauteur. C'est le brin haut qui rend cette roue cohérente à 32 m³/h.

Quand le bord de fuite s'arrête avant la fente de sortie, une **variante à
diamètre égal** prolonge l'aube jusqu'à la fente, entrée inchangée. Sur hel1, à
1450 tr/min et 32 m³/h : 10,6 m au lieu de 8,6, et 17,3 m au lieu de 13,5 à débit
nul.

**Hélice libre.** Le bilan par élément de pale voit ce qui fait l'intérêt de la
boucle : une pale normale perd de la portance près de son bout libre, où le
fluide la contourne (facteur de Prandtl) ; une boucle n'a pas de bout libre. La
toroïdale est donnée en **deux estimations** — perte de bout conservée, comme si
chaque brin avait un bout libre, et supprimée — et la réalité est entre les
deux : la jonction des brins au bout de la boucle perturbe elle aussi
l'écoulement, d'une façon que le modèle ne décrit pas. Seul un essai la situe.
Quand la lecture sépare les deux brins de chaque boucle, la toroïdale compte deux
fois plus de surfaces portantes que la normale ; quand elle ne les sépare pas, le
rapport le dit, et l'écart ne porte que sur le bout de pale. Le bruit et le
tourbillon de bout, argument premier des hélices toroïdales, sont hors du modèle.

**Le refus.** Une pièce que la lecture dit conventionnelle est refusée. La
lecture peut se tromper : `--topologie-pale toroidale` poursuit alors l'analyse,
en le signalant en tête du rapport. Un maillage troué laisse le verdict
suspendu, et la pièce est refusée sauf déclaration. En mode composants, la pale
est toroïdale par défaut ; une pale étanche d'un seul tenant, de genre 0, n'a pas
d'anse et ne peut pas être une boucle : elle est refusée. Le drapeau est levé par
la ligne de commande et la page ; l'appel direct au moteur le laisse baissé.

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

### Optimiser une roue à aubes en boucle plutôt que la remplacer

`synthetic.toroidal_impeller` engendre une roue de pompe à aubes toroïdales —
distincte de l'hélice nue de `toroidal_propeller`. Chaque aube y est faite de
deux **rubans** suivant la même loi de cambrure, l'un contre le dessus de la
veine, l'autre contre le dessous, qui se rejoignent avant le rayon extérieur
pour ne faire qu'une aube pleine hauteur au refoulement : la topologie mesurée
sur les roues réelles, où une coupe à azimut fixe traverse l'aube deux fois au
milieu et une seule au bout. Les rubans sont **verticaux**, comme sur une roue
coulée ; les prendre normaux à la veine les ferait déborder en rayon là où elle
descend le plus fort.

La roue est **ouverte** : plateau arrière, pas de flasque. C'est la seule
configuration où la boucle sert à quelque chose — sur une roue fermée le flasque
supprime déjà le tourbillon de bout de pale, et la boucle n'apporte que sa
surface mouillée.

Ce que le balayage montre sur cette forme : le rendement monte de 79,0 à 80,1 %
quand β2 passe de 8 à 38°, et le moteur en est le **frottement**, pas la
diffusion. `L = (r2−r1)/sin β` fait passer l'élancement de 8,1 à 3,4 et la
surface mouillée de 14 062 à 4 991 cm², tandis que le terme de diffusion reste
proche de zéro d'un bout à l'autre. Une aube en boucle présente deux fois ses
faces sur toute la longueur du canal : c'est là que se joue son rendement.

### D'où vient la hauteur

`H = u2·cu2 / g` n'est pas un modèle mais un **théorème** : il sort de la
conservation du moment cinétique et vaut pour n'importe quelle forme d'aube. Le
couple sur l'arbre est `ρQ(r2·cu2 − r1·cu1)` ; aucune géométrie ne peut ajouter
d'énergie autrement qu'en changeant `cu2`. `hydraulics/energy.py` en donne la
décomposition classique en trois termes :

```
H = (u2² − u1²)/2g  +  (w1² − w2²)/2g  +  (c2² − c1²)/2g
    \___ centrifuge __/   \___ diffusion __/   \___ cinétique __/
```

Le terme **centrifuge** ne dépend que des rayons et du régime — deux roues de
même diamètre à la même vitesse en tirent exactement la même chose, quelle que
soit la forme des aubes. La **diffusion** est la conversion de vitesse relative
en pression statique dans le canal : positive quand l'écoulement relatif
ralentit (`w2 < w1`), elle **détruit** de la hauteur quand il accélère. Le terme
**cinétique** sort en vitesse absolue et c'est à la volute de le récupérer.

La somme retombe sur `u2·cu2/g` à la précision machine, et le module le vérifie :
la décomposition ne crée rien, elle répartit. C'est ce qui la rend utile pour
diagnostiquer — sur l'hélice toroïdale de référence le terme de diffusion est
**négatif** et retire 5,7 m des 27,7 que l'effet centrifuge apporte ; le rapport
le signale explicitement.

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

Une aube **en boucle** présente deux fois ses faces pour un même canal : la
surface mouillée compte ses deux brins, et le frottement suit. Sans cette
correction le modèle aurait traité une hélice toroïdale comme une roue à aubes
simples et sous-estimé son coût — sur l'hélice réelle, de six dixièmes de point.

Deux propriétés le rendent utilisable comme critère de conception : il est
indépendant du **diamètre** et du **régime**. Vérifié sur les deux échelles d'une
même hélice réelle (Ø 335 et Ø 184, soit un rapport 0,55) et à deux régimes —
`L/Dh = 9,92`, `w2/w1 = 1,50` et le même rendement dans les trois cas.

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

## Les sorties visuelles : un instrument, pas un document

Ce n'est pas un rapport imprimé. On tourne un bouton, une aiguille bouge, une
limite s'allume.

### La carte d'occupation est l'élément principal

Elle est le premier élément de la page, en grand, et **interrogeable** : au
survol on lit `r`, `z` et la fraction angulaire occupée. L'image seule donne une
impression ; la grille se lit. Elle est embarquée quantifiée sur un octet par
cellule — une grille 200×200 tient en 53 Ko, dix fois moins qu'en JSON.

Surimpressions commutables : axe, `r1h` `r1s` `r2`, plans d'entrée et de sortie
avec le sens débitant, et — quand une topologie en boucle est détectée — **les
deux brins et le rayon de fusion**, puisque c'est ce qui invalide la lecture par
la cambrure. Sur la pièce d'essai réelle, la carte *montre* la boucle : les deux
branches se séparent puis se rejoignent.

### Les deux manipulations

Aucune ne réimplémente le modèle. Un second modèle en JavaScript aurait été libre
de diverger du Python — exactement l'incohérence que cet outil passe son temps à
retirer.

- **Régime** : les lois de similitude, que le contrôle du rapport vérifie déjà à
  0,00 %. Mesuré dans un vrai navigateur : ×2 sur le régime donne **Q ×2, H ×4,
  P ×8, couple ×4**. C'est le modèle lui-même. La borne du curseur va au-delà de
  la vitesse limite pour que le franchissement soit **atteignable** : un curseur
  qui s'arrête avant la réponse ne pose pas la question.
- **β2 ± 1°** : trois courbes réellement calculées en Python, la page interpole
  entre elles. L'avertissement du rapport devient manipulable.

### Une seule source de palette

```
#FBFBFA fond   #16232B encre   #DDE3E0 grille
#1D6F6A teal   → MESURÉ        #5B4B8A violet → DÉCLARÉ    #9B1D20 rouge → LIMITE
```

`io/style.py` sert la page **et** les deux figures. Une page claire à côté de
figures restées aux réglages d'origine donne un résultat incohérent, et c'est le
défaut qu'on oublie le plus souvent.

La couleur ne porte jamais seule : mesuré en graisse 600, déclaré avec un filet
violet, défaut en italique atténué, confiance faible avec un filet rouge **et**
la valeur entre parenthèses.

L'échelle de la carte est séquentielle à teinte unique, du fond au teal, à
luminance monotone (vérifié par test). Ni viridis ni jet : une échelle
multicolore fabrique des frontières que les données n'ont pas. Les régimes
suivent la même logique — une seule teinte éclaircie, parce que la vitesse est
une grandeur *ordonnée* et non des catégories.

### Deux défauts trouvés en chemin

**La page chargeait IBM Plex depuis Google Fonts.** Elle affirmait son autonomie
et ne l'avait pas : hors ligne — c'est-à-dire sur le poste d'atelier où elle sert
— elle dégradait en silence. Un test *affirmait* même cette dépendance. Piles
système désormais, et **zéro ressource externe**, ce que le test vérifie.

**Le moteur de rendu des figures n'avait pas de bas de casse.** Tout le texte
des PNG remontait en capitales, ce qu'un libellé ne doit pas être. Vingt-six
glyphes ajoutés à la fonte 5×7 — avec la hauteur d'x commençant à la même ligne
pour toutes, jambages compris : démarrer `p` une ligne plus haut lui donne la
taille d'une capitale, et « disponible » se lit « disPonible ».

Le thème sombre reste disponible en bascule ; il cesse d'être le défaut.

## Import par composants déclarés

Le mode d'import global reste disponible, inchangé. Celui-ci s'ajoute à côté, et
renverse la charge : l'utilisateur déclare ce qu'il sait, fournit les pièces
séparément, et **l'outil vérifie au lieu de deviner**.

```bash
python -m impeller_analyzer --machine helice_libre --blades 5 \
    --topologie-pale toroidale \
    --entree-fluide entree.stl --sortie-fluide sortie.stl \
    --moyeu moyeu.stl --pale pale.stl --coque coque.stl
```

Trois déclarations — modèle hydraulique, topologie de pale, nombre de pales —
qui ne sont jamais inférées en mode composants. Cinq emplacements ; la pale est
obligatoire, l'entrée et la sortie fluide aussi sauf si le corps de la roue est
fourni, qui les donne. Une seule pale suffit : les N−1 autres sont reconstruites
par rotation.

Une roue se donne aussi comme on la dessine — le corps en plusieurs pièces,
toutes les pales, et rien d'autre :

```bash
python -m impeller_analyzer --machine pompe_carenee --rotation horaire \
    --moyeu corps.stl anti_retour.stl --pale p1.stl p2.stl p3.stl p4.stl p5.stl
```

### Le contrat d'export

> Toutes les pièces exportées depuis le même repère CAO, **sans recentrage**,
> Z pour axe, en centimètres.

C'est ce contrat qui remplace la détection d'axe et le recentrage. Le fichier
d'essai réel, [`examples/Drawing1.stl`](examples/Drawing1.stl), se trouve à
**29,5 m de l'origine CAO** (x = 23,8 m, y = 17,5 m) — et c'est normal, c'est
là que la CAO l'avait mis. Le contrat fixe la *direction* de l'axe, pas sa position :
celle-ci est mesurée sur les solides fluide, qui sont des couronnes centrées
dessus.

### Ce que le mode apporte

| | Import global | Composants déclarés |
|---|---|---|
| Section d'entrée A1 | π(r1s²−r1h²), suspendue à deux détections | **volume/épaisseur** du solide d'entrée, à 0,02 % |
| Sens du fluide | convention « refoulement vers −Z » | vecteur mesuré entrée→sortie |
| Sens de rotation | indéterminé, à déclarer | déduit : `signe(ω) = signe(k·flux_z)` |
| Nombre de pales | analyse de Fourier | déclaré |
| Type de roue | classification géométrique | déclaré |

### Les huit contrôles

Une déclaration est une entrée, jamais une dispense de contrôle. Deux sont
**bloquants**. Le repère commun, parce que sa violation ne se voit sur aucune
grandeur publiée : chaque pièce se lit correctement dans son coin, et seule leur
position relative, donc tout ce que le mode apporte, est fausse. Et les plans
distincts : le même fichier donné en entrée et en sortie rendait un sens débitant
nul, que la normalisation changeait sans bruit en +Z, et dont le sens de rotation
était ensuite déduit.

Les six autres avertissent sans annuler : interpénétration des pièces (mesurée en
volume), position de la pale entre les deux plans, recouvrement des copies de la
pale tournées de 2π/N, topologie déclarée confrontée au **genre topologique**,
étanchéité, et **passage libre** : la part d'un solide fluide qu'occupe une paroi
n'est pas une section de passage. L'aire retenue pour la suite est toujours la
part libre ; le contrôle signale l'écart au-delà de 5 %. Deux s'ajoutent selon ce
qui est fourni : **pales distinctes**, quand chaque pale a son fichier, et
**entrée et sortie**, quand elles sont déduites du corps.

### Une roue donnée comme on la dessine

La roue d'essai hel2 est arrivée en sept fichiers : le corps, le disque
anti-retour, et chacune des cinq pales. Aucun plan fluide. Trois choses l'ont
rendue calculable telle quelle :

- **Le corps en plusieurs pièces.** `--moyeu` prend plusieurs fichiers. Ils sont
  réunis pour les rayons et les coupes, mais chaque pièce garde son propre test
  de volume : le disque anti-retour de hel2 entre de 0,4 mm dans l'anneau du
  corps, et dans un maillage fusionné la zone commune, traversée deux fois par
  le rayon du test de parité, passait pour vide — juste là où la chambre haute se
  ferme.
- **Toutes les pales.** `--pale` prend plusieurs fichiers, et leur nombre donne N.
  Un contrôle vérifie qu'elles sont les copies tournées d'une même pale : même
  hauteur, même volume, un pas de 360/N degrés. Sur hel2, p2 était 4,18 mm plus
  haute que les autres — exportée avant d'avoir été déplacée avec elles. Le
  calcul prend pour modèle la pale qui s'accorde avec le plus d'autres.
- **Ni entrée ni sortie.** Une roue fermée les montre. L'entrée est l'**œillard**,
  le percement d'un flasque d'extrémité autour de l'axe ; la sortie, la **fente**
  ouverte au bord de la roue entre deux parois. Les deux se lisent sur le plan
  méridien du corps, et le contrôle « entrée et sortie » dit où elles ont été
  posées. Une roue percée aux deux bouts, ou sans fente au bord, ne se laisse pas
  deviner : l'outil demande alors les plans.

### Ce que la première roue réelle a appris au mode composants

La roue d'essai hel1 — pompe fermée, cinq aubes en boucle, refoulement radial —
est arrivée en cinq STL exportés d'AutoCAD. Quatre défauts sont apparus, tous
corrigés et couverts par des tests :

- **Pièces déplacées une à une.** `STLOUT` n'exporte que dans l'octant positif ;
  chaque pièce avait été poussée de son côté, coin de boîte à l'origine, l'entrée
  à 86 mm de l'axe du corps. Aucune n'étant recentrée *sur* l'origine, le contrôle
  répondait « même origine ». Il vérifie maintenant que les pièces de révolution
  (entrée, sortie, moyeu) sont sur le même axe, et qu'aucune paire de pièces de
  tailles différentes ne partage le même coin de boîte. Le contrôle est bloquant, et
  son message dit comment exporter : tout sélectionner, déplacer une seule fois,
  exporter sans plus rien bouger.
- **Refoulement radial.** La sortie d'une roue centrifuge est une bande
  cylindrique, pas une tranche plate. Lue comme une tranche, sa hauteur passait
  pour son épaisseur : 6 cm² au lieu de 71. La bande est reconnue (paroi radiale
  plus mince que sa hauteur), son aire de passage est la moyenne de ses surfaces
  latérales, et elle fixe r2, b2 et le type de roue — mesuré, et non plus déduit
  du mode déclaré. Sur hel1, la paroi relue vaut 0,997 mm.
- **Moyeu et flasque d'un seul tenant.** Le « corps » exporté portait le flasque :
  son rayon extérieur, 96 mm, devenait r1h autour d'un œillard de 36. Le moyeu est
  maintenant lu par une coupe au plan d'entrée.
- **Disque sans sommet central.** AutoCAD triangule un disque depuis son bord :
  lu sur ses sommets, son rayon intérieur valait son rayon extérieur, et l'aire
  d'entrée dix millions de cm². Les rayons intérieurs se lisent désormais sur une
  coupe, avec un test de parité dont la demi-droite évite la couture des solides
  de révolution d'AutoCAD, placée à y = 0.

Replacée dans son vrai repère, la même roue a fait tomber trois contrôles qui
affirmaient trop :

- **Copies de pale.** Le contrôle comparait l'étendue angulaire de la pale au
  secteur de 2π/N : une aube en boucle en couvre 356 degrés, et cinq pales qui
  existent bel et bien étaient déclarées en conflit. Huit hélicoïdes de 100 degrés
  l'étaient aussi, alors qu'ils s'emboîtent comme une vis à huit filets. Le test
  porte maintenant sur les volumes : la part de la surface d'une copie tournée qui
  entre dans la pale d'origine. Sur hel1, 4,5 % — les pales réelles sont soudées
  entre elles. Le seuil (10 %) écarte l'erreur grossière ; le contrôle ne départage
  pas N et N+1, et le dit.
- **Interpénétration.** Le recouvrement de boîtes échouait sur toute roue fermée,
  dont le corps contient les pales par construction. La mesure porte sur la part
  de surface d'une pièce dans le volume de l'autre : 9 % de la pale plonge dans
  moyeu et flasque, pour la soudure.
- **Sens de rotation.** La pente hélicoïdale est un critère de machine axiale.
  Sur une aube en boucle, une régression sur la boucle entière ne mesure rien ; sur
  un refoulement radial, le critère ne s'applique pas. Il rendait pourtant un sens
  en confiance haute. Aube en boucle : le sens reste à déclarer. Refoulement radial :
  il se lit sur le recul des aubes, en confiance moyenne, l'hypothèse d'aubes
  courbées vers l'arrière étant dite.

### Les angles d'une pale importée seule

En import global, la cambrure se lit mal : la coupe traverse plusieurs pales et
une aube en boucle la traverse deux fois. Importée seule, la pale lève les deux
difficultés. Sur une roue à refoulement radial, le plan z = constante est la
surface de courant ; chaque plan coupe la pale selon un profil fermé, coupé en
deux faces entre son point le plus proche de l'axe et le plus éloigné, et la
ligne moyenne est prise à mi-chemin des deux faces. L'angle de pale se mesure
depuis la tangente, **à l'opposé de la rotation** : sans sens de rotation, rien
n'est publié, plutôt que de choisir entre une aube courbée vers l'arrière et la
même courbée vers l'avant.

- **Angles de bord extrapolés au bord.** Lus sur la zone de 5 à 20 % de l'étendue
  radiale qui suit chaque bord — l'arrondi du bord lui-même n'en dit rien —, ils
  glissaient vers l'angle de l'autre bord d'un huitième de β2 − β1 : 5 degrés sur
  le β2 de hel1. L'angle local y est maintenant ajusté en droite et extrapolé au
  bord. Sur des aubes de synthèse d'angles connus, la lecture est exacte à 0,1 degré
  dans une veine plane ; dans une veine inclinée de 30 degrés, que le plan coupe en
  biais, elle est basse de 2 à 4 degrés, et la note le dit.
- **Vrillage jugé sur le seul bord.** Il se jugeait sur tous les niveaux, y compris
  ceux qui coupent la pale loin de son bord d'attaque : une aube droite passait pour
  vrillée de 20 degrés. Au-delà de 20 degrés d'écart sur un bord, un β moyen ne
  représente plus l'aube, et la confiance des angles tombe à basse.
- **Les angles imposés l'emportent.** `--beta1` et `--beta2` étaient ignorés en
  mode composants ; ils priment maintenant sur la lecture, qui reste publiée à côté.

### Par où l'eau traverse une aube en boucle

Les plans de coupe regroupent les niveaux d'une aube en boucle en deux **brins**
de recul opposé. Ils ne se fusionnent pas : le modèle de ligne moyenne ne décrit
qu'une grille d'aubes, et il faut savoir laquelle refoule. Cela, seules les parois
le disent. Leur tracé dans le plan méridien — une case est paroi si la matière
l'occupe à deux azimuts sur trois — permet de chercher, depuis le bord de fuite de
chaque niveau, un **chemin** vers la fente de sortie : l'eau peut longer une paroi,
monter ou descendre, mais ne revient jamais vers l'axe. Une demi-droite radiale ne
suffisait pas : sur hel1, elle heurtait le cône du fond, que l'eau longe jusqu'à la
fente.

Sur hel1, ce tracé a renversé une conclusion. Le corps porte un disque
intermédiaire, ouvert seulement entre r = 35 et 62 mm, et la chambre du dessus est
fermée à sa périphérie. Le brin bas débouche vers la fente par ses huit niveaux ; le
brin haut par aucun. L'eau venue de l'œillard le traverse pour gagner le brin bas :
les deux brins sont en série, et non en parallèle. Sans pièce de paroi, l'outil
déclare la disposition indéterminée.

C'est voulu : la chambre fermée et le brin qu'elle loge doivent empêcher l'eau de
repartir vers l'œillard. L'outil le vérifie. Une aube dont l'azimut varie avec la
hauteur est une vis : en tournant dans le sens s, elle pousse l'eau le long de
l'axe dans le sens −s · signe(dθ/dz). Le passage où l'eau traverse le brin haut se
lit dans le plan qui sépare les deux brins — les rayons qu'aucune paroi n'y
occupe, l'ouverture du disque intermédiaire. Sur hel1, en rotation horaire vue de
l'entrée, le brin haut y pousse l'eau vers le brin bas sur tous les rayons sondés,
de 39 à 61 mm, comme une vis inclinée de 10 degrés sur la tangente. En rotation
anti-horaire, il la renverrait vers l'œillard, et l'outil le signale en confiance
basse. Au-delà de l'ouverture, la chambre, bordée par le seul corps de la roue,
tourne en bloc avec l'eau qu'elle contient : le brin n'y fait pas travailler
l'eau, et il n'y a pas de brassage à compter. Le modèle de ligne moyenne ne décrit
que le brin bas ; la prérotation que le brin haut donne à l'eau n'y est pas
comptée.

Le même tracé a mesuré deux écarts sur hel1 :

- **La fente de sortie chevauche le corps.** La bande de sortie couvre les deux
  lèvres du corps au bord de la roue : 71 cm² déclarés, 56 de libres, pour 9,3 mm
  de hauteur libre au lieu de 11,8.
- **La pale s'arrête avant la fente.** Le bord de fuite du brin bas est incliné,
  de r = 76 à 94 mm selon la hauteur, et vaut 85 mm en moyenne quand la fente est
  à 95,9 mm. Entre les deux, une couronne sans aube, où l'eau conserve son moment
  cinétique : pour la ligne moyenne, r2 est le rayon du bord de fuite.

Les angles lus sur le brin bas sont β1 ≈ 47° (de 19 à 71° le long d'un bord
d'attaque incliné) et β2 ≈ 23° (de 14 à 39°), en confiance basse. Les premiers
chiffres publiés, 67,5 et 27°, venaient des seuls niveaux du haut et d'une lecture
décalée vers l'autre bord.

### Les courbes en mode composants

Le modèle de ligne moyenne prend r1, A1 et β1 là où l'eau aborde les aubes, r2, A2
et β2 là où elle les quitte. En mode composants, ces bords se lisent sur la pale
elle-même, et non sur les solides fluide :

- **Sortie : le bord de fuite de la grille qui refoule.** r2 est son rayon moyen,
  b2 sa hauteur, A2 = 2π·r2·b2·τ2. Chaque plan de coupe est une surface de courant
  d'épaisseur connue ; la hauteur d'un bord est celle des niveaux qui le portent.
  La fente de sortie, mesurée et publiée, est au-delà : l'eau y arrive par une
  couronne sans aube.
- **Entrée : le bord d'attaque de cette même grille**, A1 = 2π·r1·h1·τ1, sauf quand
  l'eau rencontre d'abord un autre brin en amont, qu'elle traverse le long de l'axe
  et qui la pousse vers l'aval. La roue est alors une seule roue à deux grilles en
  série : l'entrée est la couronne de passage, β1 l'angle de la vis de ce brin au
  rayon moyen. Le travail d'Euler ne dépend que du premier bord et du dernier ;
  l'entraînement que le premier brin donne à l'eau est interne à la roue.
  L'incidence au bord d'attaque de la seconde grille et les pertes du coude entre
  les deux ne sont pas modélisées, et les courbes sortent en confiance basse.
- **Si ce brin amont renvoie l'eau vers l'œillard**, pour le sens de rotation
  déclaré, aucune courbe n'est publiée.

Sur hel1, en rotation horaire, l'entrée est l'ouverture du disque intermédiaire
(r = 39 à 61 mm, A1 = 62,5 cm² obstruction comprise, β1 = 10,3°) et la sortie le
bord de fuite du brin bas (r2 = 85 mm, b2 = 12,7 mm, β2 = 23°). À 1450 tr/min :
32 m³/h au point d'incidence nulle pour 8,6 m, 13,5 m à débit nul, NPSHr de
1,4 m ; à 2900 tr/min, 64 m³/h pour 34 m. Tout est en confiance basse : les deux
bords du brin bas sont très vrillés, et la seconde grille est traitée comme la
suite de la première. Le rendement affiché au meilleur point n'est pas une
prédiction : c'est la valeur de référence sur laquelle le modèle cale ses pertes.

Lire β1 comme la moyenne des angles d'attaque du brin bas, 47°, plaçait le débit
d'incidence nulle au-delà du débit de hauteur nulle, et aucune courbe ne sortait :
dans l'ouverture, l'eau descend le long de l'axe, et c'est l'angle de la vis
qu'elle voit, pas celui d'un profil coupé à plat.

### Le genre topologique tranche ce que l'occupation ne peut pas

Une aube en boucle est un tore : une anse par pale. Le fichier d'essai réel rend
un **genre de 6 pour 5 pales** — et le genre ne dépend pas de l'étanchéité, là où
la signature par carte d'occupation devait suspendre son verdict sur un maillage
à 60 arêtes de bord. Les deux lectures concordent, et le verdict cesse d'être
indécidable.

Avec un garde-fou, mesuré et nécessaire : **chaque déchirure fabrique une anse**.
La même roue centrifuge conventionnelle privée de 5 % de ses triangles rend un
genre de **185**. Le genre n'est donc lu que sous deux conditions cumulées — moins
de 1 % d'arêtes ouvertes, et un genre compris entre une et trois anses par pale.
Sans elles, on échangerait un verdict faux contre un autre.

## Ce que l'outil déclare, et ce qu'il mesure

Quatre corrections de la même famille : l'outil affirmait, à quelques endroits,
plus qu'il ne savait.

### Le type de machine se déclare

Sur une hélice à aubes en boucle, la classification géométrique répondait
« centrifuge » **en confiance haute**. C'était le seul endroit du programme où
une valeur fausse était affirmée sans réserve — et cette erreur fermait le mode
dont l'utilisateur avait besoin, l'analyse en hélice libre étant conditionnée à
elle.

Les critères de classification — rapport r2/r1s, solidité, largeur de sortie —
supposent tous un canal méridien conventionnel, bordé par le moyeu et le carter,
que le fluide traverse une fois. Une aube en boucle n'en a pas. La classification
garde donc une valeur, il faut bien en publier une, mais **une topologie en
boucle la fait passer en confiance faible d'elle-même**, avec l'invitation à
déclarer.

```bash
--machine {auto, pompe_carenee, helice_libre}    # choisit le modèle hydraulique
--type-de-roue {auto, axiale, mixte, centrifuge} # prime sur la classification
```

Le mode hélice libre ne dépend plus de la classification : `--machine
helice_libre` l'ouvre quoi qu'en dise la géométrie, et `--machine pompe_carenee`
le ferme même sur une roue lue axiale.

### Alésage et moyeu ne sont pas la même chose

Une pièce percée de part en part n'a pas de moyeu plein. Rendre `r1h = 0` était
exact au sens du critère et **faux** au sens hydraulique : la valeur part dans
`A1 = π(r1s² − r1h²)`, qui comptait alors le trou central comme section de
passage. Sur la pièce d'essai, cela gonflait le débit de **38 %**.

Trois cas sont maintenant distingués et nommés dans le tableau :

| Nature du centre | r1h publié | Confiance |
|---|---|---|
| moyeu plein | rayon du moyeu, au plan d'entrée | haute |
| alésage traversant | rayon intérieur de la matière | moyenne |
| ni moyeu ni alésage | 0 | moyenne |

La confiance tombe à *moyenne* dès qu'aucun moyeu plein n'est trouvé : un rayon
intérieur de matière n'est pas un rayon de moyeu, et le dire serait abusif.

La nature se tranche **au plan d'entrée**, pas sur toute la hauteur : le plateau
arrière d'une roue centrifuge est bien du plein depuis l'axe, mais il n'obstrue
rien à l'aspiration — le compter donnait un rayon intérieur supérieur au rayon
d'œillard, et une section d'entrée négative.

### Provenance : mesuré, déclaré, ou par défaut

Le niveau de confiance dit *à quel point* une grandeur est sûre. Il ne dit pas
*d'où* elle vient — et une valeur imposée en ligne de commande peut être
parfaitement sûre sans rien devoir au maillage. Le tableau 1 porte donc une
colonne de plus :

| Grandeur | Valeur | Provenance | Confiance |
|---|---|---|---|
| Type de roue | axiale | déclaré | haute |
| Rayon intérieur de matière r1h (mm) | 26.73 | mesuré | moyenne |
| Nature du centre | alésage traversant | mesuré | moyenne |
| Sens de rotation | horaire | déclaré | haute |

Le `(imposé)` collé au sens de rotation a disparu de la valeur : la colonne le
porte désormais pour toutes les grandeurs.

### Sensibilité publiée, et non plus seulement signalée

L'incertitude annoncée par le modèle — ±18 % sur la hauteur — suppose la
géométrie juste. Un tableau mesure l'autre moitié de la question :

| Grandeur | par degré de β1 | par degré de β2 | élasticité au diamètre |
|---|---|---|---|
| hauteur | 1.2 % | 0.6 % | 2.00 |
| débit | 5.0 % | 0.0 % | 3.00 |
| NPSHr | 2.9 % | 0.0 % | 2.00 |

Au-delà de **30 % par degré**, la grandeur passe automatiquement en confiance
faible — puissance et couple avec elle, puisqu'ils en dérivent.

La colonne de droite est une élasticité sans dimension, `(dX/X)/(dD/D)`. Sa
valeur est **connue d'avance** : les lois de similitude donnent 2 pour la
hauteur, 3 pour le débit, 2 pour le NPSHr. Elle sert donc aussi de contrôle du
modèle, et elle retombe dessus à 0.01 près.

Corrigé au passage : `head_sensitivity` divisait par un intervalle de **deux**
degrés tout en annonçant l'effet d'**un seul**. Elle surestimait donc du facteur
deux ce qu'elle décrivait. Le seuil d'alerte a été ramené de 18 % à 9 % pour que
le déclenchement reste identique.

## En hélice libre : poussée, et le plafond qui juge le rendement

Les phases 5 et 6 traitent la roue en **pompe** : carénée, refoulant dans une
volute, jugée sur une hauteur et un débit. La même pièce axiale peut tourner en
**hélice libre** — bateau, drone, banc d'essai — et la question devient alors une
poussée et un rendement propulsif. Deux régimes distincts, pas deux façons de
regarder le même :

```bash
python -m impeller_analyzer helice.stl --rotation horaire --rpm 1450 \
    --vitesse-avance 3.0 --fluide eau
```

Sans `--vitesse-avance`, rien n'est calculé. Sur une roue centrifuge ou mixte, le
module **refuse de répondre** et dit pourquoi : son modèle suppose un disque non
caréné traversé axialement.

### Le rendement ne se juge pas dans l'absolu

C'est l'idée reprise d'un simulateur d'hélice tiers, et la seule qui méritait de
l'être : un rendement propulsif ne veut rien dire seul, il se juge contre le
maximum que la conservation de la quantité de mouvement autorise **pour cette
poussée-là** — le disque actif idéal de Froude, pertes de profil, de bout de pale
et de giration toutes mises à zéro.

| | Valeur |
|---|---|
| Rendement propulsif calculé | 58.4 % |
| **Plafond idéal (Froude)** | **78.2 %** |
| Écart au plafond | 19.8 points |

L'écart dit ce qu'un meilleur dessin de pale peut reprendre, et rien de plus.
Aucune hélice ne franchit ce plafond : si le calcul le dépasse, c'est une erreur
de programme, et l'analyse le déclare comme telle plutôt que d'imprimer le
chiffre.

### La géométrie remplace les curseurs

Le modèle est un bilan **par élément de pale et quantité de mouvement** (BEM).
Là où un simulateur d'hélice demande à l'utilisateur de taper un coefficient de
portance de dessin et une traînée de profil qu'il ne connaît pas, ici :

- le **calage** vient de β mesuré coupe par coupe ;
- la **corde** et l'**épaisseur relative** sont lues sur les profils extraits ;
- l'**angle de portance nulle** sort de la flèche de la ligne de cambrure
  mesurée (`α₀ = −2 h/c`, théorie des profils minces), et c'est lui qui fait
  qu'un profil cambré porte déjà à incidence nulle.

Il reste deux entrées non géométriques — la pente de portance et la traînée de
base d'un profil lisse — toutes deux dans `config.py` avec leur source.

### Comment il est vérifié

Aucune courbe d'essai n'étant disponible hors ligne, les tests vérifient les
**invariants** plutôt qu'un chiffre absolu :

- **le plafond n'est jamais franchi**, sur toute la courbe ;
- **traînée de profil mise à zéro, l'écart au plafond se referme** de 25 à
  17 points — sans s'annuler, puisque la giration et le bout de pale subsistent.
  Un modèle dont l'écart ne bougerait pas ne ferait pas passer la traînée par où
  il faut ;
- la poussée part d'un maximum à l'arrêt, décroît, change de signe : la signature
  d'une hélice à pas fixe, qui freine au-delà de son avance de poussée nulle ;
- **CT et CP sont invariants par le régime** à J égal, à 10⁻⁶ près.

Le bilan est résolu par **bissection sur l'angle d'écoulement** et non par point
fixe sur les vitesses induites : le résidu change de signe une fois sur (0, π/2)
et la bissection ne diverge jamais, là où le point fixe oscille dès que la
solidité dépasse quelques dixièmes — ce qui est le cas de toute roue de pompe
axiale. Cette forme du résidu vaut aussi **à l'arrêt**, où elle se réduit à
`sin²φ = σ·Cn/(4F)` : la poussée statique sort du même calcul, sans formule
empirique séparée.

Ce que le modèle ne sait pas : les pertes d'interaction entre pales, le nombre de
Reynolds, la compressibilité (signalée au-delà de Mach 0.78 en bout de pale), et
l'état de sillage turbulent au-delà d'un facteur d'induction de 0.4 — signalé
lui aussi, station par station.

## Audit : ce qu'une relecture complète de l'outil a corrigé

### Ce que le banc d'audit a trouvé

Les tests par phase vérifient chacun un mécanisme sur une forme choisie pour
lui. Un banc distinct ([`tests/test_audit_geometrique.py`](tests/test_audit_geometrique.py))
fait l'inverse : il balaie des roues **entières**, de proportions, d'angles et
de nombres d'aubes variés, et confronte chaque grandeur relue à celle qui a été
dessinée. Ce sont les cas où les mécanismes se contredisent qu'il cherche, et
quatre défauts y sont apparus.

**Le seuil de présence de matière dépendait du nombre d'aubes.** L'occupation
angulaire d'une pale vaut `N·e / (2πr)` : quatre pales de 4 mm à r = 90 mm
donnent 0.017, sous le seuil `F_VIDE = 0.02` qui séparait le vide de la
matière. Une couronne entière passait pour du vide, le filtre d'îlots emportait
avec elle toute la bande extérieure devenue isolée, et **r2 se raccourcissait de
4 % — donc la hauteur, qui va comme `u2²`, de 8 %**. Toute roue à pales peu
nombreuses ou minces au grand rayon était concernée. `F_MATIERE = 1e-4` sépare
désormais « il y a de la matière » de « la cellule est de la veine fluide », et
les îlots parasites que ce seuil bas laisse passer sont écartés par le filtre
dont c'est le rôle. Sur six aubes, rien ne change ; sur quatre, l'écart sur r2
passe de −3.7 % à −0.4 %.

**Un verdict topologique était rendu sur un maillage qui n'a pas de topologie.**
Les trous d'un fichier déchiré coupent les tronçons en hauteur exactement comme
le ferait une aube en boucle : une roue centrifuge ordinaire, privée de 20 % de
ses triangles, était annoncée **toroïdale**. C'est la pire des sorties, puisque
c'est ce mot qui met la cambrure de côté et invalide toute la ligne moyenne. Le
verdict est maintenant suspendu quand le maillage n'est pas étanche : la
signature est relevée et dite, la conclusion ne l'est pas, et le message demande
de réparer le fichier. (Ceci ne concerne pas la décimation de la phase 8.4, qui
procède par effondrement d'arêtes et préserve l'étanchéité.)

**Une roue fermée à flasque plat était annoncée toroïdale elle aussi**, mais
pour une autre raison : le générateur posait le fond du moyeu sous la veine au
seul rayon extérieur. Sur un flasque plat le dessous de veine monte avec le
rayon, ce fond passait donc au-dessus de la veine à l'ouïe, le profil méridien
se croisait et le solide sortait creux. Défaut du générateur, pas de l'analyse —
mais il produisait un faux positif parfaitement crédible.

**Le même bloc de 86 lignes était défini deux fois** dans `blade_angles.py`
(`_local_derivatives`, `_taper_rates`, `_trim_ends`, `beta_windows`), la seconde
définition masquant silencieusement la première. Les corps étaient identiques à
un texte de docstring près, donc sans effet sur les résultats ; le doublon est
supprimé.

Deux constantes de `config.py` ne servaient plus (`GRID_R_MARGIN`,
`CAMBER_ITERATIONS`, cette dernière héritée d'une méthode de cambrure
abandonnée) : une constante que personne ne lit est un mensonge sur le code,
elles sont retirées.

### Le vrillage est lu aplati

β1 et β2 sont pris comme la moyenne sur les dix premiers et dix derniers pour
cent de corde. Une moyenne de fenêtre rend la valeur au **milieu** de la
fenêtre, pas à son bord : la lecture rabat donc les deux extrémités vers la
moyenne. Mesuré sur des roues d'angles connus :

| dessiné | relu | écart |
|---|---|---|
| 15 / 20 | 14.5 / 19.5 | −0.5 / −0.5 |
| 20 / 25 | 20.7 / 24.3 | +0.7 / −0.7 |
| 20 / 35 | 21.5 / 32.5 | +1.5 / −2.5 |
| 30 / 40 | 31.2 / 38.2 | +1.2 / −1.8 |
| 35 / 50 | 36.8 / 47.3 | +1.8 / −2.7 |
| 40 / 65 | 43.0 / 60.4 | +3.0 / −4.6 |

Le sens est constant — β1 trop grand, β2 trop petit — et l'écart croît avec le
vrillage : la lecture en restitue 60 à 90 %. Évaluer la fenêtre à son bord par
une droite des moindres carrés plutôt qu'en son milieu gagne un degré au-delà de
50°, mais **en perd deux sous 20°**, là où sont les aubes de pompe : la
correction a été mesurée puis écartée. Le biais est donc chiffré et dit — le
rapport porte une note qui encadre le vrillage réel — plutôt que déplacé.
Au-delà, le corriger demanderait de caler une loi d'aube, et celle des roues de
synthèse n'est pas celle des roues réelles.

### Les réserves d'une lecture écartée portent maintenant leur étiquette

Sur une roue toroïdale, la cambrure est mise de côté au profit des normales,
mais elle a déjà produit ses réserves — dont « β2 = 86° : aubes quasi radiales »,
imprimé à côté d'un tableau qui publie **3.6°**. Un paragraphe d'introduction
n'y suffisait pas : qui parcourt les puces lit les deux chiffres et ne sait
lequel croire. Chaque réserve issue de la lecture écartée porte désormais
`[lecture par la cambrure, écartée]` en tête.

Et la confrontation entre le sens imposé et celui que suggère la géométrie était
faite **avant** le repli : elle nommait donc le sens lu par la cambrure, alors
que le rapport publie celui des normales — les deux pouvant différer. Elle est
maintenant rendue une fois la lecture des angles arrêtée.

### Le domaine des entrées est contrôlé à la porte

`--grille 0 0` sortait en `ZeroDivisionError`, `--secteurs 1` en « inf n'est pas
sérialisable en JSON », et un régime négatif passait pour une « géométrie
dégénérée ». Chaque borne est maintenant celle d'un modèle nommé — domaine de la
corrélation d'Antoine pour la température, troposphère du modèle d'atmosphère
OACI pour l'altitude, finesse sous laquelle la carte d'occupation ne résout plus
rien pour la grille — et le contrôle est porté par `Options.check`, appelé par
`run` : la ligne de commande, la page web et l'appel direct passent par la même
définition.

`python -m impeller_analyzer.cli roue.stl` — faute de frappe naturelle pour
`python -m impeller_analyzer roue.stl` — rendait la main sans rien dire ni rien
écrire, code de sortie zéro. Le module a maintenant son garde `__main__`.

## Stress : ce que la deuxième relecture a trouvé

Une seconde campagne a soumis l'outil à ce qu'il rencontrera hors du banc :
fichiers corrompus ou dans des variantes de format légales, pièces qui ne sont
pas des roues, roues exportées couchées, scannées, loin de l'origine, arguments
absurdes, et le serveur local à des requêtes hostiles. Quatre suites de tests
(`test_stress_*.py`) en gardent la trace ; voici ce qu'elles ont fait corriger.

### Une pièce qui n'est pas une roue ne reçoit plus de performance

Une roue exportée avec un corps parasite à côté recevait un débit de 9,3 m³/h et
une hauteur de 8,5 m, en confiance « faible » — mot qui dit « estimation
incertaine », là où il n'y avait pas d'estimation du tout. Deux roues côte à côte,
pareil. L'analyse s'arrête désormais avant tout calcul hydraulique, avec
l'occupation et la topologie pour comprendre pourquoi, si l'une de deux
propriétés manque :

- **la périodicité** : une roue de N pales se superpose à elle-même après une
  rotation de 2π/N. Mesuré sur toutes les roues valides du banc, y compris une roue
  privée de 5 % de ses triangles, l'écart reste sous 1,8 % du rayon extérieur ; les
  pièces parasites donnent 29 à 30 %. Le seuil de refus est à 10 % (`SYM_REJECT`) ;
- **un anneau de matière** faisant le tour de l'axe — moyeu, alésage, flasque ou
  jante. Deux roues identiques côte à côte sont réellement symétriques d'ordre 2
  autour de l'axe médian : seule l'absence d'anneau les distingue d'un rotor.

`--sans-controle-symetrie` lève ce contrôle pour une roue volontairement
irrégulière, et le rapport le dit.

### L'axe d'une roue à deux pales

L'inertie désigne l'axe des deux moments égaux. Pour deux pales élancées, c'est
**l'envergure** : le maillage était basculé de 90 degrés, en confiance haute.
Quand l'inertie hésite ou contredit Z, chaque axe principal est maintenant mis à
l'épreuve de la périodicité, et l'axe d'ordre le plus élevé l'emporte — une roue
à aubes hélicoïdales est aussi symétrique d'ordre 2 autour d'axes perpendiculaires
au sien. À deux pales, les trois axes sont d'ordre 2 et rien ne les départage :
Z est gardé, et l'axe est publié **non vérifié**, en confiance moyenne.

### Le nombre de pales se vérifie sur la pièce

Le spectre angulaire lisait 9 pales sur une roue de 3 aubes en boucle : un
harmonique de la forme des aubes. Quand la rotation de 2π/N ne superpose pas la
pièce, les diviseurs de N et les pics les plus forts sont essayés, et le plus grand
ordre qui tient est retenu, en confiance moyenne. Un nombre **imposé** par
`--blades` est lui aussi confronté à la pièce : démenti, il perd la confiance
haute qu'il recevait d'office.

### Un scan ne fabrique plus d'aubes en boucle

0,1 mm de bruit sur les sommets — ce que rend un scan 3D — suffisait à faire lire
des aubes en boucle sur une roue centrifuge ordinaire, et β2 tombait de 24,5 à
6,6 degrés. La face bosselée d'un moyeu plein n'est plus pleine qu'à 97 %, bascule
dans le masque des pales, et le dessus et le dessous du moyeu se lisaient comme
les deux brins d'une boucle. Un tronçon d'au plus deux cellules collé à de la
matière pleine est désormais une frange, pas un brin (`LOOP_FRINGE_ROWS`). Jusqu'à
0,4 mm de bruit, β2 reste à moins de trois degrés, et la confiance sur les angles
descend d'elle-même à moyenne.

### Ce qui n'est pas lu n'est pas publié comme lu

Une plaque à deux lobes recevait « 2 pales, roue axiale, moyeu plein » en
confiance haute, des angles « 0,0 / 0,0 » de provenance « mesuré », une vitesse
maximale de « 0 tr/min », et l'assurance que « le débit et le NPSHr restent
fiables » — sans qu'aucun débit n'ait été calculé. Sans aube lue, nombre et type
sont plafonnés à moyenne, les angles sont « non lus », la vitesse maximale « non
calculable ». Zéro tr/min reste une réponse quand le NPSH disponible est négatif.

### Import

Un STL ASCII indenté était lu comme un binaire tronqué ; un binaire SolidWorks à
en-tête « solid » suivi d'octets de bourrage, refusé ; un OFF coloré, décalé ; un
PLY à coordonnées de texture, lu sur la mauvaise liste. Coordonnées infinies et
indices hors bornes (OBJ numéroté à partir de zéro) sortaient en trace de pile :
ils sont refusés à l'import avec leur cause. Un facteur d'unité `nan` rendait tout
le maillage NaN. Une pièce à deux kilomètres de l'origine est signalée : un export
en simple précision n'y résout que 0,1 mm. Les sept formats lisibles rendent
exactement la même analyse de la même roue.

### Entrées et sorties

Hauteur et pertes d'aspiration, rayon d'aspiration : `nan` sortait en trace de
pile, des pertes négatives ajoutaient de l'énergie, un rayon d'aspiration de 10 m
passait sur une roue de 10 cm. 10⁸ secteurs angulaires épuisaient la mémoire. Un
dossier de sortie impossible produisait une trace. Chaque cas donne maintenant une
phrase ; et une exception imprévue est annoncée comme un défaut **de l'outil**,
code de sortie 3, trace jointe.

### Mode composants

Le seuil d'anisotropie des tranches fluides valait 2,0 — exactement la limite
qu'une tranche mince atteint sans jamais la dépasser (I_n = I₁ + I₂) : chaque
analyse avertissait d'une normale « incertaine » sur des tranches parfaites. Il
vaut 1,5, soit une tranche aussi épaisse que son rayon. Une erreur d'import nomme
maintenant l'emplacement fautif.

### Serveur local et page

Le serveur vérifie que la requête lui est adressée (`Host`) et vient de sa propre
page (`Origin`) : une page malveillante ne peut plus lui poster de fichier ni se
faire passer pour lui par rebinding DNS. Un nombre non fini dans le formulaire
donne un 400, plus un 500. La charge JSON de la page est écrite en échappements
unicode pour `<`, `>` et `&` : aucune chaîne ne peut plus fermer l'élément
`<script>`.

### Ce que la campagne a vérifié sans rien trouver

Lois de similitude exactes (Q ∝ n, H ∝ n², P ∝ n³), effet d'échelle exact (taille
×10 à régime /10 : débit ×100, hauteur identique), NPSH disponible décroissant
avec la température et l'altitude, vitesse maximale décroissante avec les pertes,
bilan d'énergie refermé à 10⁻¹⁶. Même roue tournée de 37 degrés, exportée Y en
haut, en mm ou en m, faces permutées, normales inversées, en ASCII : même roue, à
0,2 degré près. Son miroir tourne en sens inverse.

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

### Le fichier d'essai réel

[`examples/Drawing1.stl`](examples/Drawing1.stl) est la seule géométrie non
synthétique du dépôt : une hélice à cinq aubes en boucle, exportée d'AutoCAD en
centimètres, à 29,5 m de l'origine CAO. Elle sert de banc, pas de cas de
référence : aucune performance mesurée ne l'accompagne.

```bash
python -m impeller_analyzer examples/Drawing1.stl --type-de-roue axiale
```

| Grandeur | Mesurée hors de l'outil | Lue par l'outil |
|---|---|---|
| Triangles | 29 598 | 29 598 |
| Arêtes de bord après réparation | 60 | 60 |
| Aire de surface | 5 839 cm² | 5 839,4 cm² |
| Volume | 640 cm³ | 640,2 cm³ |
| Nombre de pales | 5 | 5, confiance haute |
| Topologie | toroïdale | toroïdale, genre 6 |
| Diamètre | 267,2 mm | 270,5 mm (+1,2 %) |
| Rayon intérieur de matière | 26,7 mm | **37,4 mm — écart non résolu** |

Le dernier écart n'a pas été réglé en ajustant l'outil : aucun sommet du
maillage n'approche l'axe à moins de 37 mm, quel que soit le centrage retenu. La
valeur de référence vaut exactement le dixième du diamètre de référence, ce qui
ressemble à une cote dérivée plutôt que mesurée. Elle reste à vérifier sur la
pièce.

## Calage

Toutes les constantes physiques, tous les coefficients empiriques et tous les
seuils numériques sont dans [`impeller_analyzer/config.py`](impeller_analyzer/config.py),
une par ligne, commentée avec sa source. Aucune valeur numérique significative
n'apparaît ailleurs. Pour recaler l'outil, on ne modifie que ce fichier.

## Tests

```bash
python -m unittest discover -s tests -t tests
```

411 tests : une phase par module, le banc d'audit qui balaie des roues entières
et confronte chaque grandeur relue au dessin, les invariants de l'analyse en
hélice libre, ce que l'outil a le droit d'affirmer, l'import par composants
déclarés avec ses sept contrôles, les sorties visuelles — palette unique,
autonomie de la page, recalcul des deux curseurs — et les trois suites de stress :
fichiers corrompus, pièces hors domaine, entrées hostiles. Ils passent aussi sous
`pytest` si vous l'avez : ce sont des `unittest.TestCase`.

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
    ├── energy.py        # bilan d'Euler : d'où vient la hauteur
    ├── cavitation.py    # NPSHr, NPSHa, vitesse maximale
    └── similarity.py    # lois de similitude, vitesse spécifique
```

## Limites

Modèle 1D ligne moyenne : **hauteur ±18 %, débit ±25 %, NPSHr ±30 %** — et le
cas de référence ci-dessus suggère que l'écart sur la hauteur peut être plus
grand encore sur une pompe réelle tant que le calage n'a pas été fait.

Limites de lecture connues, que l'outil signale quand il les rencontre :

- le nombre de pales se cherche entre 2 et 12 ; au-delà, il faut l'imposer
  (`--blades`, jusqu'à 24) ;
- l'axe d'une roue à **deux** pales hélicoïdales ne se déduit pas de la géométrie :
  exportez-la avec son axe sur Z ;
- une pièce volontairement irrégulière — pas variable, pale cassée — est refusée
  par le contrôle de périodicité, qu'il faut lever sciemment.

**À vérifier par essai sur banc avant toute décision d'achat ou de
dimensionnement.**
