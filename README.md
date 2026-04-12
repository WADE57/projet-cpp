# Simulation thermique 1D en C++

Projet de simulation thermique avec deux modeles:

- **stationnaire** (regime permanent)
- **instationnaire** (evolution temporelle)

Le code resout un systeme tridiagonal par decomposition LU et exporte les resultats en **CSV** et **VTK**.

## Prérequis

- **Compilateur** : `g++` (suite GNU C++)
- **Norme C++** : C++23
- **Système d'exploitation** : Linux, macOS, Windows (MinGW/MSYS2)
- **Outils optionnels** :
  - [Paraview](https://www.paraview.org/) pour visualiser les fichiers VTK
  - Python (pour traiter les fichiers CSV)

## Description scientifique

### Problème physique

Le projet simule la **conduction thermique 1D** dans un domaine solide avec échange convectif aux frontières.

**Équation de conduction** (loi de Fourier) :
$$\rho C_p \frac{\partial T}{\partial t} = k \frac{\partial^2 T}{\partial x^2} + S$$

Où :

- $\rho$ : densité du matériau (kg/m³)
- $C_p$ : capacité calorifique massique (J/kg·K)
- $k$ : conductivité thermique (W/m·K)
- $S$ : source volumique (W/m³)
- $T$ : température (K)

**Conditions aux limites** :

- **Flux prescrit** : $-k \frac{\partial T}{\partial x} = \Phi$ (ou transfert convectif : $h(T - T_e)$)
- **Température imposée** : $T = T_0$

### Approches numériques

**Cas stationnaire** : Résout $\frac{d^2T}{dx^2} = 0$ (équilibre thermique)

**Cas instationnaire** : Utilise un schéma implicite Euler pour la discrétisation temporelle :
$$\frac{T^{n+1} - T^n}{\Delta t} = \frac{k}{\rho C_p} \frac{\partial^2 T^{n+1}}{\partial x^2}$$

## Structure du projet

```text
CSV/
INCLUDE/
PARAMETRE/
RAPPORT/
SRC/
VTK/
````

- `SRC/main.cpp` : point d'entree, lance les deux modeles.
- `SRC/matrice.cpp` + `SRC/resolution.cpp` : assemblage et resolution du cas stationnaire.
- `SRC/matInst.cpp` + `SRC/resolutionInst.cpp` : assemblage et resolution du cas instationnaire.
- `SRC/parametres.cpp` : lecture de `PARAMETRE/simu.cfg`.
- `PARAMETRE/simu.cfg` : parametres physiques et numeriques.
- `CSV/` : resultats tabulaires (stationnaire.csv, instationnaire.csv).
- `VTK/` : sorties pour visualisation.
- `RAPPORT/rapport.pdf` : rapport du projet.

## Parametres utilises

Le fichier `PARAMETRE/simu.cfg` est lu dans un **ordre fixe**, par exemple:

```properties
Lx 0.04 Ly 0.004 Lz 0.05
M 10000
Phip 125000
hc 200
rho 2700
Cp 940
kapa 164
Te 20
stationary 0
TFinal 300
N 600
Mx 50 My 10 Mz 30
```

## Parametres détaillés

| Paramètre | Type | Unité | Description |
|---|---|---|---|
| `Lx`, `Ly`, `Lz` | float | m | Dimensions du domaine (longueur, largeur, profondeur) |
| `Mx`, `My`, `Mz` | int | - | Nombre de points dans chaque direction |
| `M` | int | - | Nombre total de points de discrétisation |
| `rho` | float | kg/m³ | Densité du matériau |
| `Cp` | float | J/kg·K | Capacité calorifique massique |
| `kapa` | float | W/m·K | Conductivité thermique (kappa) |
| `Phip` | float | W/m² | Flux de chaleur prescrit à la frontière |
| `hc` | float | W/m²·K | Coefficient d'échange convectif |
| `Te` | float | K ou °C | Température extérieure (ambiance) |
| `TFinal` | float | s | Temps final de simulation |
| `N` | int | - | Nombre de pas de temps |
| `stationary` | int | - | **0** : mode instationnaire, **1** : mode stationnaire uniquement |

## Architecture du code

### Modules principaux

- **`parametres.cpp/hpp`** : Lecture et stockage des paramètres depuis `simu.cfg`
- **`matrice.cpp/hpp`** : Construction de la matrice tridiagonale pour le cas stationnaire
- **`resolution.cpp/hpp`** : Factorisation LU et résolution du système stationnaire
- **`matInst.cpp/hpp`** : Construction de la matrice pour le schéma implicite (instationnaire)
- **`resolutionInst.cpp/hpp`** : Boucle temporelle et résolution itérative
- **`main.cpp`** : Orchestration : orchestration des deux cas, I/O CSV/VTK

### Flux d'exécution

1. **Initialisation** : Lecture de `simu.cfg` par `parametres.cpp`
2. **Cas stationnaire** (si `stationary >= 0`) :
   - Assemblage matrice `matrice.cpp`
   - Résolution LU `resolution.cpp`
   - Export CSV et VTK
3. **Cas instationnaire** (si `stationary == 0`) :
   - Boucle temporelle `resolutionInst.cpp`
   - Export aux temps définis : t = 15s, 30s, 60s, 90s, 150s, 210s

## Compilation et execution

### Compilateur

`g++` (C++23)

Depuis la racine du projet:

```bash
g++ -std=c++23 -O2 SRC/main.cpp SRC/parametres.cpp SRC/matrice.cpp SRC/resolution.cpp SRC/matInst.cpp SRC/resolutionInst.cpp -o SRC/main
Version debug(affichage supplementaire):
```

Version debug:

```bash
g++ -std=c++23 -DDEBUG -g SRC/main.cpp SRC/parametres.cpp SRC/matrice.cpp SRC/resolution.cpp SRC/matInst.cpp SRC/resolutionInst.cpp -o SRC/main
```

### Exécution

**Important** : Lancer l'exécutable depuis le dossier `SRC/`, car les chemins sont relatifs (`../PARAMETRE/simu.cfg`, `../CSV/`, `../VTK/`).

```bash
cd SRC
./main ../PARAMETRE/simu.cfg
```

#### Exemple complet (depuis la racine du projet)

```bash
# 1. Compiler
g++ -std=c++23 -O2 SRC/main.cpp SRC/parametres.cpp SRC/matrice.cpp \
    SRC/resolution.cpp SRC/matInst.cpp SRC/resolutionInst.cpp -o SRC/main

# 2. Exécuter
cd SRC
./main ../PARAMETRE/simu.cfg

# 3. Vérifier les résultats
ls -la ../CSV/       # Résultats tabulaires
ls -la ../VTK/       # Résultats pour Paraview
```

**Sortie attendue** :

```
Lecture des parametres OK
Simulation stationnaire lancee...
Simulation instationnaire lancee...
Export CSV et VTK termines.
```

## Sorties générées et visualisation

### Format CSV

#### Cas stationnaire : `CSV/stationnaire.csv`

```
x,numeric,exacte
0.0,58.4544,58.4487
0.001,60.2345,60.2234
...
```

#### Cas instationnaire : `CSV/instationnaire.csv`

```
x,t=15s,t=30s,t=60s,t=90s,t=150s,t=210s
0.0,37.0699,42.1234,48.5678,52.3456,56.1234,58.1065
...
```

### Format VTK pour Paraview

- **Stationnaire** : `VTK/stationnaire.vtk` (structure 3D avec champ de température)
- **Instationnaire** : `VTK/instationnaire/masolution_*.vtk` (un fichier par pas de temps pour animation)

### Ouvrir dans Paraview

1. Télécharger [Paraview](https://www.paraview.org/)
2. Ouvrir le fichier VTK :

   ```bash
   paraview ../VTK/stationnaire.vtk
   ```

3. Appliquer un filtre de visualisation (Colors, Representation)
4. Pour l'instationnaire, charger la série `masolution_*.vtk` pour une animation temporelle

## Méthodes numériques

### Factorisation LU pour systèmes tridiagonaux

Pour un système tridiagonal $AX = F$ :

- **Décomposition** : $A = LU$ (L inférieure, U supérieure, adaptée à la structure tridiagonale)
- **Résolution** en deux étapes :
  - **Descente** : $LY = F$ (substitution avant)
  - **Remontée** : $UX = Y$ (substitution arrière)
- **Avantage** : Complexité $O(n)$ au lieu de $O(n^3)$ pour Gauss standard

### Cas stationnaire

- Assemblage d'une matrice tridiagonale pour l'équation de conduction 1D avec échange convectif.
- Conditions aux limites incorporées directement.

### Cas instationnaire

- Schéma implicite **Euler rétrograde** (stabilité inconditionnelle)
- Boucle temporelle : pour chaque pas $n = 1, 2, ..., N_t$
  - Assembler la matrice avec pas de temps $\Delta t$
  - Résoudre le système linéaire
  - Exporter si temps d'intérêt ($t = 15, 30, 60, 90, 150, 210$ s)

## Validation et tests

### Cas stationnaire

Le code compare la solution numérique avec une **solution analytique** (disponible pour conduction sans source avec conditions limites simples).

**Métrique d'erreur** :
$$L^\infty = \max_i |T_{num,i} - T_{exact,i}|$$
$$L^2 = \sqrt{\sum_i (T_{num,i} - T_{exact,i})^2}$$

**Exemple de résultats** :

- x = 0 : numeric = 58.4544, exacte = 58.4487 → erreur ≈ 0.0057°C

### Cas instationnaire

Validation qualitative par :

- Monotonie de la température (évolution physique correcte)
- Convergence vers la solution stationnaire à temps infini
- Respect du bilan énergétique (conservation de l'énergie)

**Exemple** :

- x = 0 : t=15s → 37.0699°C, t=210s → 58.1065°C (approche vers 58.45°C du stationnaire)

### Recommandations de tests

1. Vérifier que les résultats CSV ne contiennent pas de `NaN` ou `inf`
2. Comparer `CSV/stationnaire.csv` et `CSV/instationnaire.csv` à t=TFinal
3. Tracer les courbes avec Python :

   ```python
   import pandas as pd
   import matplotlib.pyplot as plt
   df = pd.read_csv('../CSV/stationnaire.csv')
   plt.plot(df['x'], df['numeric'], 'b-', label='Numérique')
   plt.plot(df['x'], df['exacte'], 'r--', label='Exacte')
   plt.legend()
   plt.show()
   ```

## Documentation détaillée

Pour une analyse approfondie (équations, résultats détaillés, graphiques, discussions) :

- **Rapport complet** : [`RAPPORT/rapport.pdf`](RAPPORT/rapport.pdf)

## Licence

MIT License - Ce projet est libre d'utilisation, modification et distribution.

## Auteur

- [Wade](https://github.com/WADE57)
- **Dépôt** : <https://github.com/WADE57/projet-cpp>
- **Domaine** : Master CSMI, M1, Formation en simulation numérique thermique
