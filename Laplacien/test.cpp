#include <iostream>
#include "matricebande.h"
#include "vecteur_template.h"
using namespace std;

// ce fichier teste les methodes assemblageT, preconditionne et PCG de la classe matricebande

int main(){
    try
    {
        int N = 3;
        cout << "=== TEST assemblageT ===\n";
        cout << "N = " << N << endl;

        // 1) Construire la matrice de Laplacien

        matricebande B;
        matricebande A = B.laplacien(N);

        cout << "\nMatrice A (format bande) :\n";
        cout << A << endl;

        
        // CHOIX 1 : C = I
        cout << "\n---- Choix 1 : C = I ----\n";
        A.assemblageT(1);

        cout << "TIND : ";
        for (int k : A.get_TIND())
            cout << k << " ";
        cout << endl;

        cout << "TTAB :\n";
        cout << A.get_TTAB() << endl;

        // CHOIX 2 : C = diag(A)
        cout << "\n---- Choix 2 : C = diag(A) ----\n";
        A.assemblageT(2);

        cout << "TIND : ";
        for (int k : A.get_TIND())
            cout << k << " ";
        cout << endl;

        cout << "TTAB :\n";
        cout << A.get_TTAB() << endl;

        // CHOIX 3 : IC(0)
       
        cout << "\n---- Choix 3 : IC(0) ----\n";
        double omega = 0.0;
        A.assemblageT(3, omega);

        cout << "TIND : ";
        for (int k : A.get_TIND())
            cout << k << " ";
        cout << endl;

        cout << "TTAB :\n";
        cout << A.get_TTAB() << endl;

        cout << "\n=== FIN TEST assemblageT ===\n";
    }
    catch (const exception& e)
    {
        cerr << "ERREUR : " << e.what() << endl;
    }


    printf("\n"); 
    
    try
    {
        cout << "=== TEST preconditionne(r) ===\n";

        int N = 3;
        cout << "N = " << N << endl;

        matricebande B;
        matricebande A = B.laplacien(N);

        int n = A.dim1();

        // Vecteur r test (simple)
        vector<double> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = i + 1;   // r = (1,2,3,...)

        cout << "\nVecteur r :\n";
        for (double ri : r) cout << ri << " ";
        cout << endl;

       
        // TEST 1 : Choix 1 (C = I)
        cout << "\n---- Choix 1 : C = I ----\n";
        A.assemblageT(1);

        vector<double> z1 = A.preconditionne(r);

        cout << "z = C^{-1} r (doit etre egal a r) :\n";
        for (double zi : z1) cout << zi << " ";
        cout << endl;

        // TEST 2 : Choix 2 (C = diag(A))
        cout << "\n---- Choix 2 : C = diag(A) ----\n";
        A.assemblageT(2);

        vector<double> z2 = A.preconditionne(r);

        cout << "z = C^{-1} r :\n";
        for (double zi : z2) cout << zi << " ";
        cout << endl;

        // TEST 3 : Choix 3 (IC(0))
        cout << "\n---- Choix 3 : IC(0) ----\n";
        double omega = 0.0;
        A.assemblageT(3, omega);

        vector<double> z3 = A.preconditionne(r);

        cout << "z = C^{-1} r :\n";
        for (double zi : z3) cout << zi << " ";
        cout << endl;

        cout << "\n=== FIN TEST preconditionne ===\n";
    }
    catch (const exception& e)
    {
        cerr << "ERREUR : " << e.what() << endl;
    }

    printf("\n"); 
    
    try
    {
        cout << "=== TEST PCG ===\n";

        int N = 3;
        cout << "N = " << N << endl;

        matricebande B;
        matricebande A = B.laplacien(N);

        int n = A.dim1();

        // Second membre b
        // (exemple simple : b = (1,1,...,1))
        vector<double> b(n, 1.0);

        // Vecteur initial x0 = 0
        vector<double> x(n, 0.0);

        double tol = 1e-6;
        int nb_iter;
        vector<double> res_hist;

        // TEST 1 : PCG avec C = I
        cout << "\n---- PCG : C = I ----\n";
        A.assemblageT(1);

        x.assign(n, 0.0);
        res_hist.clear();

        A.PCG(b, x, tol, nb_iter, res_hist);

        cout << "Iterations : " << nb_iter << endl;
        cout << "Dernier residu relatif : " << res_hist.back() << endl;

        // TEST 2 : PCG avec C = diag(A)
        cout << "\n---- PCG : C = diag(A) ----\n";
        A.assemblageT(2);

        x.assign(n, 0.0);
        res_hist.clear();

        A.PCG(b, x, tol, nb_iter, res_hist);

        cout << "Iterations : " << nb_iter << endl;
        cout << "Dernier residu relatif : " << res_hist.back() << endl;

        // TEST 3 : PCG avec IC(0)
        cout << "\n---- PCG : IC(0) ----\n";
        double omega = 0.0;
        A.assemblageT(3, omega);

        x.assign(n, 0.0);
        res_hist.clear();

        A.PCG(b, x, tol, nb_iter, res_hist);

        cout << "Iterations : " << nb_iter << endl;
        cout << "Dernier residu relatif : " << res_hist.back() << endl;

        // Affichage de l'historique des résidus
        cout << "\nHistorique des residus (IC(0)) :\n";
        for (size_t k = 0; k < res_hist.size(); ++k)
            cout << "k = " << k
                 << "   ||r_k|| / ||b|| = "
                 << res_hist[k] << endl;

        cout << "\n=== FIN TEST PCG ===\n";
    }
    catch (const exception& e)
    {
        cerr << "ERREUR : " << e.what() << endl;
    }
    return 0;

}