#include <iostream>
#include "matricebande.h"
#include "vecteur_template.h"
using namespace std;

int main() {
    int N = 20;
    
    // CAS (1.1) : rho = 1  --> Laplacien standard
    cout << "===== CAS (1.1) : rho = 1 =====\n";
    matricebande B;
    matricebande A = B.laplacien(N);   

    vector<double> b = smbr(N);        // second membre
    vector<double> x(A.dim1(), 0.0);

    double tol = 1e-6;
    int it;
    vector<double> res;
    // -------- C = I --------
    A.assemblageT(1);
    x.assign(A.dim1(), 0.0);
    res.clear();

    A.PCG(b, x, tol, it, res);
    cout << "C = I : " << it << endl;
    
    {
    ofstream fout("residus_PCG_CI.txt");
    if (!fout)
        throw runtime_error("Erreur ouverture residus_PCG_CI.txt");

    for (size_t k = 0; k < res.size(); ++k)
        fout << k << " " << res[k] << endl;

    fout.close();
    }    

    // -------- C = diag(A) --------
    A.assemblageT(2);
    x.assign(A.dim1(), 0.0);
    res.clear();

    A.PCG(b, x, tol, it, res);
    cout << "C = diag(A) : " << it << endl;

    {
    ofstream fout("residus_PCG_diag.txt");
    if (!fout)
        throw runtime_error("Erreur ouverture residus_PCG_diag.txt");

    for (size_t k = 0; k < res.size(); ++k)
        fout << k << " " << res[k] << endl;

    fout.close();
    }

    // -------- C = IC(0) --------
    A.assemblageT(3,0.0);
    x.assign(A.dim1(), 0.0);
    res.clear();

    A.PCG(b, x, tol, it, res);
    cout << "C = IC(0) : " << it << endl;

    {
    ofstream fout("residus_PCG_IC0_omega_0.txt");
    if (!fout)
        throw runtime_error("Erreur ouverture residus_PCG_IC0_omega_0.txt");

    for (size_t k = 0; k < res.size(); ++k)
        fout << k << " " << res[k] << endl;

    fout.close();
    }

    
    // CAS (1.2) : rho variable (gâteau)
    cout << "===== CAS (1.2) : rho variable (gateau) =====\n";
    matricebande C;
    matricebande A2 = C.laplacien_rho(N);

    vector<double> b2 = smbr_rho(N);        // second membre
    vector<double> x2(A2.dim1(), 0.0);

    double tol2 = 1e-6;
    int it2;
    vector<double> res2;

    // -------- C = I --------
    A2.assemblageT(1);
    x2.assign(A2.dim1(), 0.0);
    res2.clear();

    A2.PCG(b2, x2, tol2, it2, res2);
    cout << "C = I : " << it2 << endl;

    {
        ofstream fout("residus_PCG_rho_CI.txt");
        if (fout)
            for (size_t k = 0; k < res2.size(); ++k)
                fout << k << " " << res2[k] << endl;
    }


    // -------- C = diag(A) --------
    A2.assemblageT(2);
    x2.assign(A2.dim1(), 0.0);
    res2.clear();

    A2.PCG(b2, x2, tol2, it2, res2);
    cout << "C = diag(A) : " << it2 << endl;

    {
        ofstream fout("residus_PCG_rho_diag.txt");
        if (fout)
            for (size_t k = 0; k < res2.size(); ++k)
                fout << k << " " << res2[k] << endl;
    }


    // -------- C = IC(0) --------
    A2.assemblageT(3,0.0);
    x2.assign(A2.dim1(), 0.0);
    res2.clear();

    A2.PCG(b2, x2, tol2, it2, res2);
    cout << "C = IC(0) : " << it2 << endl;
     
    {
        ofstream fout("residus_PCG_rho_IC0_omega_0.txt");
        if (fout)
            for (size_t k = 0; k < res2.size(); ++k)
                fout << k << " " << res2[k] << endl;
    }

    return 0;
}