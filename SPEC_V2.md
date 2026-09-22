# SPEC v2 — Import par composants déclarés
Extension de impeller-analyzer. Ne remplace pas le mode d'import actuel : l'ajoute à côté.

---

## 1. Pourquoi
Tous les échecs constatés sur un fichier réel viennent de l'inférence, aucun du calcul.

| Échec observé | Valeur inférée | Ce qu'elle devient |
|---|---|---|
| Nombre de pales | 9 au lieu de 5 | déclaré |
| Type de roue | centrifuge, sur une hélice axiale | déclaré |
| Rayon de moyeu r1h | 0.00 mm, confiance « haute » | mesuré sur le STL du moyeu |
| Section d'entrée A1 | déduite de rayons fragiles | mesurée sur le solide d'entrée |
| Sens de sortie du fluide | supposé radial | donné par le solide de sortie |
| Sens de rotation | indéterminé | déduit sans ambiguïté (§5) |
| Topologie en boucle | verdict suspendu | déclarée |
| Mode hélice libre | refusé par la classification | déclaré |

L'outil cesse d'être un devineur pour devenir un vérificateur. C'est un changement de rôle, pas un ajout de fonctionnalité.

---

## 2. Les déclarations
Trois champs, saisis avant tout import :

| Champ | Valeurs | Effet |
|---|---|---|
| `mode` | `pompe_carenee` \| `helice_libre` | choisit le modèle hydraulique ; supprime la classification automatique |
| `topologie_pale` | `conventionnelle` \| `toroidale` | choisit la logique de coupe (une ou deux traversées par azimut) |
| `nombre_de_pales` | entier 2 à 24 | supprime l'analyse de Fourier |

Ces trois valeurs ne sont jamais inférées en mode composants. Elles sont vérifiées (§6), jamais devinées.

---

## 3. Les emplacements d'import

| Emplacement | Statut | Ce que l'outil en tire |
|---|---|---|
| `coque` | optionnel | rayon de carter, jeu en bout de pale, présence d'un flasque |
| `moyeu` | recommandé | r1h, r_hub(z), longueur de moyeu |
| `entree_fluide` | obligatoire | plan d'entrée, aire A1, normale, position |
| `sortie_fluide` | obligatoire | plan de sortie, aire A2, normale, position |
| `pale` | obligatoire | géométrie de pale unique, isolée |

Chaque emplacement accepte un `.stl`, `.obj`, `.ply`, `.off` ou `.step`. Un emplacement vide n'est pas une erreur : l'outil signale simplement les grandeurs qu'il ne pourra pas produire.

### 3.1 Une seule pale par défaut
Le mode normal importe une pale et reconstruit les N−1 autres par rotation de 2π/N autour de Z. Prévoir un mode `pales_distinctes` acceptant N fichiers, pour les roues à pas ou espacement variables. Ne jamais imposer N imports quand la roue est régulière : chaque import supplémentaire est une occasion d'export raté.

### 3.2 Les solides fluide
Les modéliser en CAO comme des tranches minces, épaisseur 2 à 5 mm, placées exactement dans le plan d'entrée et dans le plan de sortie. Pas des tubes longs : un tube ne dit pas où se trouve la section de référence.

De chaque tranche, l'outil tire :
* aire = volume / epaisseur, mesurée et non calculée depuis des rayons
* normale = vecteur propre isolé du tenseur d'inertie de la tranche
* centroide = position du plan
* sens_debitant = de `entree_fluide.centroide` vers `sortie_fluide.centroide`

---

## 4. Règle d'export imposée à l'utilisateur
**Toutes les pièces sont exportées depuis le même repère CAO, sans recentrage, avec Z pour axe de rotation, en centimètres.**

Conséquences directes :
* la détection d'axe est supprimée en mode composants — Z est l'axe, par contrat
* le recentrage est supprimé — l'origine est celle de la CAO
* les pièces sont déjà assemblées entre elles, aucun recalage n'est nécessaire

C'est aussi le seul risque réel de la méthode : si un exportateur recentre une pièce sur l'origine, l'assemblage est silencieusement détruit. L'outil doit donc le détecter (§6).

Afficher cette règle dans l'interface d'import, au-dessus des emplacements.

---

## 5. Ce que les déclarations débloquent

### 5.1 Sens de rotation, enfin déterminé
Le critère reste `signe(ω) = signe(∂z/∂θ)` mesuré sur la surface moyenne de pale. Ce qui change : la direction de refoulement n'est plus une convention (vers −Z) mais un vecteur mesuré, `entree_fluide → sortie_fluide`. On projette k sur cette direction. L'indétermination disparaît, et `--rotation` devient facultatif.

Pour une roue centrifuge, où le refoulement est radial, la même logique s'applique dans le plan méridien : la pale doit pousser le fluide de la section d'entrée vers la section de sortie.

### 5.2 Cambrure, enfin extractible
Une coupe cylindrique d'une pale isolée donne un seul profil fermé. Plus aucune segmentation à faire.

Pour `topologie_pale = toroidale`, la coupe donne deux profils : la branche aller et la branche retour de la boucle. Les traiter comme deux cascades en série, et sommer leurs déviations. Ne plus tenter de les fusionner.

### 5.3 Débit, enfin indépendant de β1
`Q = cm1 · A1` avec A1 mesurée sur le solide d'entrée. Conserver le calcul par incidence nulle `cm1 = u1 · tan(β1)` comme aujourd'hui, mais A1 cesse d'être suspendue à la détection de r1s et r1h.

---

## 6. L'outil vérifie, il n'accepte pas
Une déclaration est une entrée, jamais une dispense de contrôle. Le système de confiance reste en place, à l'identique. Contrôles obligatoires :

1. **Repère commun.** Si le centroïde d'une pièce est à moins de 1 % du rayon extérieur de l'origine alors que les autres en sont loin, cette pièce a probablement été recentrée à l'export. Avertir et bloquer.
2. **Non-interpénétration.** Les boîtes englobantes de moyeu, pale et coque ne doivent pas se recouvrir de plus de 5 %. Au-delà, avertir.
3. **Position de la pale.** Le centroïde de la pale doit se trouver entre le plan d'entrée et le plan de sortie, le long du sens débitant. Sinon, les deux solides fluide sont probablement inversés : avertir explicitement.
4. **Nombre de pales déclaré.** Reconstruire l'assemblage par N rotations de la pale et vérifier qu'aucune copie n'en recoupe une autre. Un recouvrement signifie que N est trop grand pour cette pale : avertir.
5. **Topologie déclarée.** Si toroidale est coché, la pale doit être de genre topologique 1 — une boucle, donc une anse. Compter les composantes de bord et la caractéristique d'Euler du maillage. Si la pale est simplement connexe, avertir que la déclaration contrait la géométrie.
6. **Étanchéité.** Inchangé : plafonner la confiance des grandeurs volumiques si le maillage n'est pas fermé, et donner le nombre d'arêtes de bord restantes.

Aucun de ces contrôles n'annule une déclaration. Ils l'accompagnent d'un avertissement, et la confiance globale reste le minimum de la table.

---

## 7. Ce qui est supprimé en mode composants
* détection d'axe et recentrage — remplacés par le contrat d'export
* analyse de Fourier du nombre de pales — remplacée par la déclaration
* détection de moyeu par histogramme — remplacée par le STL de moyeu
* classification axiale / mixte / centrifuge — remplacée par mode
* détection du côté d'aspiration — remplacée par les solides fluide

Ces fonctions restent en place pour le mode d'import global, elles ne sont pas retirées du code.

---

## 8. Ce qui n'est pas supprimé
Le mode d'import en vrac reste disponible, avec le comportement actuel et sa confiance plafonnée. Les deux modes coexistent, le rapport indique lequel a servi.

---

## 9. Ordre d'implémentation demandé
Une étape à la fois, chacune validée avant la suivante.
1. Structure d'import multi-emplacements et les trois déclarations, sans aucun calcul nouveau : l'outil charge les pièces, affiche leurs cotes, et applique les contrôles du §6. Rien d'autre.
2. Extraction des plans d'entrée et de sortie depuis les solides fluide : aire, normale, centroïde, sens débitant.
3. Sens de rotation ancré sur le sens débitant (§5.1).
4. Cambrure sur pale isolée (§5.2), topologie conventionnelle d'abord.
5. Cas toroïdal : deux profils par coupe, cascades en série.
6. Branchement de A1 et A2 mesurées dans le modèle hydraulique.
