# Simulation thermique 1D en C++

Projet de simulation thermique avec deux modeles:
- **stationnaire** (regime permanent)
- **instationnaire** (evolution temporelle)

Le code resout un systeme tridiagonal par decomposition LU et exporte les resultats en **CSV** et **VTK**.

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

Le fichier `PARAMETRE/simu.cfg` contient notamment:

- Geometrie: `Lx`, `Ly`, `Lz`
- Discretisation spatiale et temporelle: `M`, `N`
- Maillage de visualisation VTK: `Mx`, `My`, `Mz`
- Proprietes physiques: `rho`, `Cp`, `kappa`, `hc`, `Te`, `Phip`
- Temps final: `TFinal`

Exemple de configuration:

```ini
Lx = 0.04 
Ly = 0.004 
Lz = 0.05
M = 10000
Phip = 125000
hc = 200
rho = 2700
Cp = 940
kappa = 164
Te = 20
stationary = 0
TFinal = 300
N = 600
Mx = 50 
My = 10 
Mz = 30
```
## Compilation et execution

- Compilateur: g++ (C++23)

Depuis la racine du projet:

```bash
g++ -std=c++23 -O2 main.cpp parametres.cpp matrice.cpp resolution.cpp matInst.cpp resolutionInst.cpp -o SRC/main
```
Version debug(affichage supplementaire):

```bash
g++ -std=c++23 -DDEBUG -g main.cpp parametres.cpp matrice.cpp resolution.cpp matInst.cpp resolutionInst.cpp -o SRC/main
```

- Execution:

Inportant: lancer l'executable depuis le dossier `SRC/`, car les chemins sont relatifs ( `../PARAMETRE/simu.cfg` , `../CSV/`, `../VTK/` ).   
```bash
cd SRC
./main ../PARAMETRE/simu.cfg
``` 
## Sorties générées

### Cas stationnaire:

- `CSV/stationnaire.csv` : avec les colonnes 

    - `x`
    - `numeric`
    - `exacte`
- `VTK/stationnaire.vtk` : pour visualisation 3D (ex: Paraview)
### Cas instationnaire:
- `CSV/instationnaire.csv` : avec des profils à plusieurs instants (`t=15s`, `30s`, `60s`, `90s`, `150s`, `210s`).
- `VTK/instationnaire/masolution*.vtk` : (un fichier par pas de temps)

## Methode numerique

- Assemblage d'une matrice tridiagonale pour l'equation de conduction 1D avec echange convectif.
- Factorisation LU adaptee au cas tridiagonal.
- Resolution en deux etapes:
  - descente: `L Y = F`
  - remontee: `U X = Y`
- Cas instationnaire: schema implicite avec boucle en temps.
## Extraits de resultats

Exemple stationnaire (`CSV/stationnaire.csv`):

- x = 0 : numeric = 58.4544, exacte = 58.4487

Exemple instationnaire (`CSV/instationnaire.csv`):

- x = 0 : t=15s -> 37.0699, t=210s -> 58.1065
## Rapport
Le rapport du projet est disponible dans:

- `RAPPORT/rapport.pdf`
## Auteur
- [Wade](https://github.com/WADE57)
- Depot: https://github.com/WADE57/projet-cpp