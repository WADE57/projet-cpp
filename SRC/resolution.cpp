#include <iostream>
#include "../INCLUDE/resolution.hpp"
#include <vector>
#include "../INCLUDE/matrice.hpp"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

// Affichage pour debug
void Aff(std::vector<std::vector<double>> matrice)
{
    for (int i = 0; i < matrice.size(); i++)
    {
        for (int j = 0; j < matrice.size(); j++)
        {
            std::cout << std::setw(12) << matrice[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\ntaille de la matrice = " << matrice.size() << std::endl;
}

void AffVec(std::vector<double> vecteur)
{
    for (int i = 0; i < vecteur.size(); i++)
    {
        std::cout << vecteur[i] << " ";
    }
    std::cout << "\ntaille du vecteur = " << vecteur.size() << std::endl;
}

// Constructeur Resolution prenant les paramètres nécessaires
Resolution::Resolution(const Matrice &matrice, int taille, const Parametre &param)
{
    // Stockage des paramètres utiles pour les méthodes membres
    this->h = param.Lx / param.M;
    this->M = param.M;
    this->hc = param.hc;
    this->p = 2 * (param.Ly + param.Lz);
    this->S = param.Ly * param.Lz;
    this->Te = param.Te;
    this->Phi_p = param.Phi_p;
    this->kappa = param.kappa;
    this->Lx = param.Lx;
    this->MxSta = param.Mx;
    this->MySta = param.My;
    this->MzSta = param.Mz;

    std::vector<double> b_star;
    std::vector<double> c_star;

    b_star.push_back(matrice.getDiagonale(0));
    c_star.push_back(matrice.getSurDiagonale(0) / b_star[0]);

    for (int k = 1; k < taille - 1; k++)
    {
        b_star.push_back(matrice.getDiagonale(k) - matrice.getSousDiagonale(k - 1) * c_star[k - 1]);
        c_star.push_back(matrice.getSurDiagonale(k) / b_star[k]);
    }
    b_star.push_back(matrice.getDiagonale(taille - 1) - matrice.getSousDiagonale(taille - 2) * c_star[taille - 2]);

#ifdef DEBUG
    std::cout << "\nb_star : ";
    AffVec(b_star);
    std::cout << "\nc_star : ";
    AffVec(c_star);
#endif

    for (int i = 0; i < taille; i++)
    {
        diagonaleL.push_back(b_star[i]);
    }
    for (int i = 0; i < taille - 1; i++)
    {
        surDiagonaleU.push_back(c_star[i]);
    }
#ifdef DEBUG
    std::cout << "\ndiagonaleL : ";
    AffVec(diagonaleL);
    std::cout << "\nsurDiagonaleU : ";
    AffVec(surDiagonaleU);
#endif
}

void Resolution::decompositionLU(const Matrice &matrice, int taille)
{
    std::vector<std::vector<double>> L(taille, std::vector<double>(taille, 0));
    std::vector<std::vector<double>> U(taille, std::vector<double>(taille, 0));

    for (int i = 0; i < taille; i++)
    {
        L[i][i] = diagonaleL[i];
    }
    for (int i = 1; i < taille; i++)
    {
        L[i][i - 1] = matrice.getSousDiagonale(i - 1);
    }
    for (int i = 0; i < taille; i++)
    {
        U[i][i] = 1;
    }
    for (int i = 0; i < taille - 1; i++)
    {
        U[i][i + 1] = surDiagonaleU[i];
    }
#ifdef DEBUG
    std::cout << "\nMatrice L : " << std::endl;
    Aff(L);
    std::cout << "\nMatrice U : " << std::endl;
    Aff(U);
#endif
}

std::vector<double> Resolution::resoudreLYF(const Matrice &matrice, int taille)
{
    std::vector<double> Y;
    Y.push_back(matrice.secondMembre()[0] / diagonaleL[0]);
    for (int k = 1; k < taille; k++)
    {
        Y.push_back((matrice.secondMembre()[k] - matrice.getSousDiagonale(k - 1) * Y[k - 1]) / diagonaleL[k]);
    }
#ifdef DEBUG
    std::cout << "\nVecteur Y : " << std::endl;
    AffVec(Y);
#endif
    return Y;
}

std::vector<double> Resolution::resoudreUXY(const std::vector<double> &Y, int taille)
{
    std::vector<double> X(taille, 0);
    X[taille - 1] = Y[taille - 1];
    for (int k = taille - 2; k >= 0; k--)
    {
        X[k] = Y[k] - surDiagonaleU[k] * X[k + 1];
    }
#ifdef DEBUG
    std::cout << "\nVecteur X : " << std::endl;
    AffVec(X);
#endif
    return X;
}

void Resolution::ecrireFichier(const std::vector<double> &vecteur, const std::string &nomFichier, const std::vector<double> &TExacte)
{
    std::ofstream fichier(nomFichier, std::ios::out | std::ios::trunc);
    if (fichier)
    {
        fichier << "x,numeric,exacte" << std::endl;
        for (int i = 0; i < vecteur.size(); i++)
        {
            fichier << i * h << "," << vecteur[i] << "," << TExacte[i] << std::endl;
        }
        fichier.close();
        std::cout << "Ecriture dans le fichier " << nomFichier << " réussie !" << std::endl;
    }
    else
    {
        std::cout << "Erreur à l'ouverture !" << std::endl;
    }
}

std::vector<double> Resolution::solutionExacte()
{
    double a = (hc * p) / (kappa * S);
    double b = (Phi_p / kappa) * ((cosh(sqrt(a) * Lx)) / (sqrt(a) * sinh(sqrt(a) * Lx)));
    std::vector<double> TExacte;
    for (int i = 0; i < M; i++)
    {
        TExacte.push_back(Te + b * ((cosh(sqrt(a) * (Lx - i * h))) / (cosh(sqrt(a) * Lx))));
    }
    return TExacte;
}

void Resolution::writeVTKFile(std::vector<double> vecteur, std::string fileName)
{
    std::ofstream fichier(fileName, std::ios::out | std::ios::trunc);
    if (fichier)
    {
        fichier << "# vtk DataFile Version 2.0" << std::endl;
        fichier << "vtk output" << std::endl;
        fichier << "ASCII" << std::endl;
        fichier << "DATASET STRUCTURED_GRID" << std::endl;
        fichier << "DIMENSIONS " << MxSta << " " << MySta << " " << MzSta << std::endl;
        fichier << "POINTS " << MxSta * MySta * MzSta << " float" << std::endl;
        for (int k = 0; k < MzSta; k++)
        {
            for (int j = 0; j < MySta; j++)
            {
                for (int i = 0; i < MxSta; i++)
                {
                    fichier << i << " " << j << " " << k << std::endl;
                }
            }
        }
        fichier << "POINT_DATA " << MxSta * MySta * MzSta << std::endl;
        fichier << "FIELD FieldData 1" << std::endl;
        fichier << "sol1 " << 1 << " " << MxSta * MySta * MzSta << " float" << std::endl;
        int k = 0, xk;
        double a, b;
        double hx = Lx / MxSta;
        for (int i = 0; i < MxSta; i++)
        {
            double xi00 = i * hx;
            for (int j = 0; j < M; j++)
            {
                double xk = j * h;
                if (xi00 >= xk && xi00 <= xk + h)
                {
                    k = j;
                    break;
                }
            }
            a = (vecteur[k + 1] - vecteur[k]) / (h);
            b = vecteur[k] - a * k * h;
            vecteur[i] = a * xi00 + b;
        }
        for (int k = 0; k < MzSta; k++)
        {
            for (int j = 0; j < MySta; j++)
            {
                for (int i = 0; i < MxSta; i++)
                {
                    fichier << vecteur[i] << std::endl;
                }
            }
        }
        fichier.close();
        std::cout << "Ecriture dans le fichier " << fileName << " réussie !" << std::endl;
    }
    else
    {
        std::cout << "Erreur à l'ouverture !" << std::endl;
    }
}