#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "matricebande.h"
#include "vecteur_template.h"

using namespace std;

double Topt(double, double)
{
    return 200.0;
}

bool is_in_gateau(int i, int j, int N)
{
    double x = (double)i / (N + 1);
    double y = (double)j / (N + 1);
    double cx = 0.5, cy = 0.5, r = 0.25;
    return ( (x - cx)*(x - cx) + (y - cy)*(y - cy) <= r*r );
}

int main()
{
    int N = 20;
    double tol = 1e-6;

    cout << "===== PROBLEME INVERSE : determination des resistances =====\n";

    matricebande B;

    // 1. T1 : chauffage par le bas seul
    cout << "\n[1] Calcul de T1 (bas seul)\n";

    matricebande A1 = B.laplacien_rho(N);
    vector<double> b1 = smbr_rho_bas_seul(N);
    vector<double> T1(A1.dim1(), 0.0);

    int it;
    vector<double> res;

    A1.assemblageT(3, 0.0);
    A1.PCG(b1, T1, tol, it, res);
    cout << "PCG T1 : " << it << " iterations\n";

    // 2. T2 : chauffage par le haut seul
    cout << "\n[2] Calcul de T2 (haut seul)\n";

    matricebande A2 = B.laplacien_rho(N);
    vector<double> b2 = smbr_rho_haut_seul(N);
    vector<double> T2(A2.dim1(), 0.0);

    res.clear();
    A2.assemblageT(3, 0.0);
    A2.PCG(b2, T2, tol, it, res);
    cout << "PCG T2 : " << it << " iterations\n";

    // 3. Construction du système A alpha = b
    cout << "\n[3] Construction du systeme 2x2\n";

    double A11 = 0.0, A12 = 0.0, A22 = 0.0;
    double b_1 = 0.0, b_2 = 0.0;

    int npts = N + 1;

    for (int j = 1; j <= N; ++j)
    {
        for (int i = 1; i <= N; ++i)
        {
            if (is_in_gateau(i, j, N))
            {
                int k = i + j * npts;
                double T1k = T1[k];
                double T2k = T2[k];
                double Toptk = 200.0;

                A11 += T1k * T1k;
                A12 += T1k * T2k;
                A22 += T2k * T2k;

                b_1 += T1k * Toptk;
                b_2 += T2k * Toptk;
            }
        }
    }

    // régularisation énergie
    A11 += (1.0 / 100.0) * 20.0 * 20.0;
    A22 += (1.0 / 100.0) * 180.0 * 180.0;

    // 4. Résolution du système 2x2
    double det = A11 * A22 - A12 * A12;
    if (fabs(det) < 1e-12)
        throw runtime_error("Matrice A singuliere");

    double alpha1 = ( b_1 * A22 - b_2 * A12 ) / det;
    double alpha2 = ( A11 * b_2 - A12 * b_1 ) / det;

    cout << "\nalpha bas  = " << alpha1 << endl;
    cout << "alpha haut = " << alpha2 << endl;

    // 5. Température finale
    vector<double> T_final(T1.size());
    for (size_t k = 0; k < T1.size(); ++k)
        T_final[k] = alpha1 * T1[k] + alpha2 * T2[k];

    cout << "\nTemperature finale calculee.\n";
    cout << "===== FIN DU PROBLEME INVERSE =====\n";

    ofstream fout("T_final.txt");

    for (int j = 0; j <= N; ++j)
    {
        for (int i = 0; i <= N; ++i)
        {
            int k = i + j * npts;
            fout << T_final[k] << " ";
        }
        fout << endl;
    }
    fout.close();

    return 0;
}
