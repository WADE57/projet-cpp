#ifndef RESOLUTION_HPP
#define RESOLUTION_HPP

#include "matrice.hpp"
#include <vector>
#include "parametres.hpp"

class Resolution
{
private:
    std::vector<double> diagonaleL;    /**< Vecteur représentant la diagonale de la matrice */
    std::vector<double> surDiagonaleU; /**< Vecteur représentant la sur-diagonale de la matrice */

    // Membres pour les méthodes
    double h, hc, p, S, Te, Phi_p, kappa, Lx;
    int M, MxSta, MySta, MzSta;

public:
    // Constructeur principal : passage de Parametre
    Resolution(const Matrice &matrice, int taille, const Parametre &param);

    // getters et setters
    double getDiagonaleL(int i) const { return diagonaleL[i]; }
    double getSurDiagonaleU(int i) const { return surDiagonaleU[i]; }
    void setDiagonaleL(int i, double val) { diagonaleL[i] = val; }
    void setSurDiagonaleU(int i, double val) { surDiagonaleU[i] = val; }

    // Méthodes de décomposition LU
    void decompositionLU(const Matrice &matrice, int taille);

    // Méthodes de résolution LY = F
    std::vector<double> resoudreLYF(const Matrice &matrice, int taille);

    // Méthodes de résolution UX = Y
    std::vector<double> resoudreUXY(const std::vector<double> &Y, int taille);

    // Méthode pour écrire dans un fichier CSV
    void ecrireFichier(const std::vector<double> &vecteur, const std::string &nomFichier, const std::vector<double> &TExacte);

    // Méthode pour calculer la solution exacte
    std::vector<double> solutionExacte();

    // Méthode pour écrire dans un fichier VTK
    void writeVTKFile(std::vector<double> vecteur, std::string fileName);

    // Destructeur
    ~Resolution()
    {
        diagonaleL.clear();
        surDiagonaleU.clear();
    }
};

#endif