#include "../INCLUDE/matInst.hpp"
#include <vector>
#include "../INCLUDE/resolution.hpp"
#include "../INCLUDE/matrice.hpp"

// Nouveau constructeur prenant Parametre en argument
MatInst::MatInst(const Parametre &para, std::vector<double> aInst, std::vector<double> bInst, std::vector<double> cInst)
{
    int MInst = para.M;
    double kappaInst = para.kappa;
    double hInst = para.Lx / para.M;
    double hcInst = para.hc;
    double SInst = para.Ly * para.Lz;
    double pInst = 2 * (para.Ly + para.Lz);
    double TeInst = para.Te;
    double Phi_pInst = para.Phi_p;
    double rhoIsnt = para.rho;
    double CpIsnt = para.Cp;
    double TFinal = para.TFinal;
    int N = para.N;
    double dt = TFinal / N;

    // remplissage des vecteurs a, b et c
    for (int i = 1; i < MInst - 1; i++)
    {
        aInst.push_back(-kappaInst / (hInst * hInst));
    }
    aInst.push_back(-kappaInst / hInst);

    bInst.push_back(kappaInst / hInst);
    for (int i = 1; i < MInst - 1; i++)
    {
        bInst.push_back(2 * kappaInst / (hInst * hInst) + ((hcInst * pInst) / SInst) + ((rhoIsnt * CpIsnt) / dt));
    }
    bInst.push_back(kappaInst / hInst);

    cInst.push_back(-kappaInst / hInst);
    for (int i = 1; i < MInst - 1; i++)
    {
        cInst.push_back(-kappaInst / (hInst * hInst));
    }

    for (int i = 0; i < MInst; i++)
    {
        diagInst.push_back(bInst[i]);
        sousDiagInst.push_back(aInst[i]);
        surDiagInst.push_back(cInst[i]);
    }

    // Stocker les paramètres utiles pour FInst
    this->MInst = MInst;
    this->hcInst = hcInst;
    this->pInst = pInst;
    this->TeInst = TeInst;
    this->SInst = SInst;
    this->Phi_pInst = Phi_pInst;
}

// matriceTridiagInst ne change pas
std::vector<std::vector<double>> MatInst::matriceTridiagInst() const
{
    std::vector<std::vector<double>> matrixInst(MInst, std::vector<double>(MInst, 0));
    for (int i = 0; i < MInst; i++)
    {
        matrixInst[i][i] = diagInst[i];
    }
    for (int i = 0; i < MInst - 1; i++)
    {
        matrixInst[i][i + 1] = surDiagInst[i];
    }
    for (int i = 1; i < MInst; i++)
    {
        matrixInst[i][i - 1] = sousDiagInst[i];
    }
    return matrixInst;
}

// FInst utilise les membres stockés dans l'objet
std::vector<double> MatInst::FInst() const
{
    std::vector<double> FInst;
    FInst.push_back(Phi_pInst);
    for (int i = 1; i < MInst - 1; i++)
    {
        FInst.push_back((hcInst * pInst * TeInst) / SInst);
    }
    FInst.push_back(0);
    return FInst;
}