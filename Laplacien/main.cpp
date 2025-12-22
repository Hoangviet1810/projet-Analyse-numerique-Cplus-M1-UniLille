#include <iostream>
#include "matricebande.h"
#include "vecteur_template.h"
using namespace std;


int main()
{
    int N = 3;   
    cout << "Test (1.8) : assemblage de A et Ax" << endl;
    cout << "N = " << N << endl;

    matricebande B;
    matricebande A = B.laplacien(N);

    cout << "\nVecteur IND :" << endl;
    vector<int> IND = A.ind();
    for (int k : IND)
        cout << k << " ";
    cout << endl;

    cout << "\nMatrice bande ATAB :" << endl;
    cout << A << endl;

    int n = A.dim1();
    vector<double> x(n, 1.0);

    cout << "\nVecteur x = (1,...,1) :" << endl;
    for (double xi : x)
        cout << xi << " ";
    cout << endl;

    //Produit y = A x  (question 1.7)
    vector<double> y = A.prod(x);

    cout << "\nVecteur y = A x :" << endl;
    for (double yi : y)
        cout << yi << " ";
    cout << endl;


    vector<double> b = smbr(N);       
    cout << "\nVecteur b :" << endl;
    for (double bi : b){
        cout << bi << " ";
    }
    cout << endl;
   

    printf("\n");

    vector<double> x1(A.dim1(), 0.0);   // x0 = 0
    int nb_iter;
    vector<double> res_hist;
    A.steepest_descent(b, x1, 1e-4, nb_iter, res_hist);
    cout << "Nombre d'iterations : " << nb_iter << endl;
    cout << "Dernier residu relatif : " << res_hist.back() << endl;

    printf("\n"); 
    cout << "\nVecteur des residus relatifs ||r_k|| / ||b|| :\n"; 
    for (size_t k = 0; k < res_hist.size(); ++k){ 
        cout << "k = " << k << " ||r_k|| / ||b|| = " << res_hist[k] << endl; 
    } 
    // printf("\n");
    // for (int k = 0; k < res_hist.size(); ++k)
    //     cout << k << " " << res_hist[k] << endl;

    printf("\n");
    vector<double> Ax = A.prod(x1);
    cout << "*) On calcul Norme ||Ax - b|| avec plusieurs types de normes : " << endl;
    cout << " Norme L2 : " << norme_l2(Ax-b) << endl;
    cout << " Norme 2  : " << norme_2(Ax-b) << endl;
    cout << " Norme L2 (euclidienne) : " << norme_eculidienne(Ax-b) << endl;
    cout << " Norme L1 : " << norme_l1(Ax-b) << endl;
    cout << " Norme L∞ : " << norme_linf(Ax-b) << endl;
    cout << "------------------------\n";

    // cout << "\nSolution approchée x_k (vecteur x1) :\n";
    // for (int i = 0; i < x1.size(); ++i)
    //     cout << "x[" << i << "] = " << x1[i] << endl;

    cout << "\nSolution x_k sous forme de maillage 2D :\n";
    int npts = N + 1;
    for (int j = 0; j < npts; ++j)
    {
        for (int i = 0; i < npts; ++i)
        {
            int k = i + j * npts;
            cout << x1[k] << "\t";
        }
        cout << endl;
    }


    // matricebande C;
    // matricebande A2 = C.laplacien_rho(N);
    // cout << "\nVecteur IND :" << endl;
    // vector<int> IND2 = A2.ind();
    // for (int k : IND2)
    //     cout << k << " ";
    // cout << endl;

    // cout << "\nMatrice bande ATAB :" << endl;
    // cout << A2 << endl;

    // int n2 = A2.dim1();
    // vector<double> x2(n2, 1.0);

    // // Produit y = A x  (question 1.7)
    // vector<double> y2 = A2.prod(x2);
    // cout << "\nVecteur y = A2 x2 :" << endl;
    // for (double yi : y2)
    //     cout << yi << " ";
    // cout << endl;

    // vector<double> b2 = smbr_rho(N);
    // cout << "\nVecteur b2 :" << endl;
    // for (double bi : b2){
    //     cout << bi << " ";
    // }
    // cout << endl;

    ofstream f("residus_steepest_descent.txt");
    for (size_t k = 0; k < res_hist.size(); ++k){
        f << k << " " << res_hist[k] << endl;
    }
    f.close();


    double diag_max = 0.0;
    for (int i = 0; i < A.dim1(); ++i){
        diag_max = max(diag_max, A(i,2));   // diagonale principale
    }
    double alpha_opt = 1.0 / diag_max;

    cout << "\nComparaison steepest descent / gradient pas fixe\n";
    cout << "alpha_opt = " << alpha_opt << "\n\n";

    // cout << " j   alpha         nb_iter\n";
    // cout << "-----------------------------\n";

    // for (int j = 1; j <= 6; ++j){
    //     double alpha = j * alpha_opt / 5.0;
    //     vector<double> x_pf(A.dim1(), 0.0);
    //     vector<double> res_pf;
    //     int it_pf;

    //     A.gradient_pas_fixe(b, x_pf, alpha, 1e-4, it_pf, res_pf);

    //     cout << " " << j << "   " << alpha << "   "  << it_pf << endl;
    // }
    
    for (int j = 1; j <= 6; ++j){
        double alpha = j * alpha_opt / 5.0;

        vector<double> x_pf(A.dim1(), 0.0);
        vector<double> res_pf;
        int it_pf;

        A.gradient_pas_fixe(b, x_pf, alpha, 1e-4, it_pf, res_pf);

        cout << " " << j << "   " << alpha << "   " << it_pf << endl;

        string filename = "residus_pas_fixe_j" + to_string(j) + ".txt";
        ofstream fout(filename);

        if (!fout)
             throw runtime_error("Erreur ouverture fichier " + filename);

        for (size_t k = 0; k < res_pf.size(); ++k)
             fout << k << " " << res_pf[k] << endl;

        fout.close();   
    }

    printf("\n"); 
    // Cette section est pour le test des 3 fonctions : assemblageT, preconditionne et PCG
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
