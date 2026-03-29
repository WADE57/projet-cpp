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

