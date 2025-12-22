#include "matricebande.h"
#include "vecteur_template.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

matricebande::matricebande(int n, int d) : matrice(n,d){
    this->indice.resize(d);
}

matricebande::matricebande(const matricebande& A) : matrice(A){
    this->indice=A.ind();}

vector<int> matricebande::ind(vector<int> v){
    return this->indice=v;
}

matricebande matricebande::laplacien(int n){
    
    int N=(n+2)*n;
    double h = 1./(n+1), hh=pow(h,2);

    vector<int> v;
    v.push_back(-n-2);v.push_back(-1); v.push_back(0);
    
    matricebande A(N,3);
    A.ind(v);
        

    int ligne;
 
    // Ici j==0
    for (int i=0; i<n+2; i++){ // CL de Neumann sur les bords verticaux
        ligne=i;
        if (i==0)
            A(ligne,2)= 4./hh;
        else if (i<n+1){
            A(ligne,1)=-1./hh;
            A(ligne,2)= 4./hh;
        }
        else{
            A(ligne,1)=-2./hh;
            A(ligne,2)= 4./hh;
        }
    }
    // CL de Dirichlet sur les bords horizontaux, j==0 déjà traité
    for (int j=1; j<n; j++){
        for (int i=0; i<n+2; i++){ // CL de Neumann sur les bords verticaux
            ligne=(n+2)*j+i;
            if (i==0){
                A(ligne,0)=-1./hh;
                A(ligne,2)= 4./hh;
            }
            else if (i<n+1){
                A(ligne,0)=-1./hh;
                A(ligne,1)=-1./hh;
                A(ligne,2)= 4./hh;
            }
            else{
                A(ligne,0)=-1./hh;
                A(ligne,1)=-2./hh;
                A(ligne,2)= 4./hh;
            }
        }
    }
    //cout << "matrice du laplacien avant symétrisation \n";
    //cout << A << endl;

    // Symetrisation de A
    for (int j=0; j<n; j++){
        ligne=(n+2)*j;
        A[ligne]=A[ligne]*0.5; // modification de toute la ligne
        ligne=(n+2)*j+n+1;
        A[ligne]=A[ligne]*0.5;
    }
    
    return A;
}

double rho(double x, double y)
{
    if (x>=0.2 && x<=0.8 && y>=0.3 && y<=0.4)
        return 100.0;
    return 1.0;
}

// (1.6)
matricebande matricebande::laplacien_rho(int n)
{
    int N = (n + 2) * n;
    double h  = 1.0 / (n + 1);
    double hh = h * h;

    // stockage bande : voisin bas, gauche, diagonale
    vector<int> v;
    v.push_back(-(n + 2)); // voisin bas
    v.push_back(-1);       // voisin gauche
    v.push_back(0);        // diagonale principale

    matricebande A(N, 3);
    A.ind(v);

    int ligne;

    // Boucle sur le maillage
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n + 2; ++i) {

            ligne = i + j * (n + 2);

            double x = i * h;
            double y = (j + 1) * h;

            // rho aux demi-mailles (prolongement symétrique)
            double rho_ip = rho(x + 0.5*h, y);
            double rho_im = (i > 0) ? rho(x - 0.5*h, y) : rho_ip;

            double rho_jp = rho(x, y + 0.5*h);
            double rho_jm = (j > 0) ? rho(x, y - 0.5*h) : rho_jp;

            // diagonale principale
            A(ligne, 2) = (rho_ip + rho_im + rho_jp + rho_jm) / hh;

            // voisin gauche
            if (i > 0)
                A(ligne, 1) = -rho_im / hh;

            // voisin bas
            if (j > 0)
                A(ligne, 0) = -rho_jm / hh;
        }
    }
    //cout << "matrice du laplacien_rho avant symétrisation \n";
    //cout << A << endl;

    // Symétrisation (comme dans laplacien)
    for (int j = 0; j < n; ++j) {
        ligne = j * (n + 2);
        A[ligne] = A[ligne] * 0.5;

        ligne = j * (n + 2) + (n + 1);
        A[ligne] = A[ligne] * 0.5;
    }

    return A;
}


// (1.7)
vector<double> matricebande::prod(const vector<double>& x) const{
    int n = this->dim1();
    int d = this->indice.size();

    if ((int)x.size() != n)
        throw runtime_error("Erreur: taille du vecteur incompatible avec la matrice.");

    vector<double> y(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < d; ++k) {

            int offset = indice[k];
            int j = i + offset;

            if (j >= 0 && j < n) {

                double aij = (*this)(i, k);

                // contribution directe
                y[i] += aij * x[j];

                // contribution symétrique
                if (offset < 0 && j != i)
                    y[j] += aij * x[i];
            }
        }
    }

    return y;
}



//(1.9)
vector<double> smbr(int n){
    int N = (n + 2) * n;
    double h  = 1.0 / (n + 1);
    double hh = h * h;
    vector<double> b(N, 0.0);
    // ligne du bas : j = 0
    for (int i = 0; i < n + 2; ++i) {
        int k = i;
        b[k] = 20.0 / hh;
    }
    // ligne du haut : j = n-1
    for (int i = 0; i < n + 2; ++i) {
        int k = i + (n - 1) * (n + 2);
        b[k] = 180.0 / hh;
    }
    return b;
}


vector<double> smbr_rho(int n){
    int N = (n + 2) * n;
    double h  = 1.0 / (n + 1);
    double hh = h * h;
    vector<double> b(N, 0.0);
    // ligne du bas : j = 0
    for (int i = 0; i < n + 2; ++i) {
        double x = i * h;
        double y = 0.0;
        double rho_jm = rho(x, y + 0.5 * h);
        int k = i;
        b[k] = rho_jm * 20.0 / hh;
    }
    // ligne du haut : j = n-1
    for (int i = 0; i < n + 2; ++i) {
        double x = i * h;
        double y = n * h;
        double rho_jp = rho(x, y - 0.5 * h);
        int k = i + (n - 1) * (n + 2);
        b[k] = rho_jp * 180.0 / hh;
    }
    return b;
}

// void matricebande::steepest_descent(const vector<double>& b,vector<double>& x,double tol,int& nb_iter,vector<double>& res_hist) const {
//     int n = this->dim1();
//     if ((int)b.size() != n || (int)x.size() != n)
//         throw runtime_error("Erreur de dimension dans steepest_descent");
//     // r0 = b - A x0
//     vector<double> r = b - (*this) * x;
//     double norm_r0 = 0.0;
//     for (double v : r) norm_r0 += v * v;
//     norm_r0 = sqrt(norm_r0);
//     double norm_r = norm_r0;
//     res_hist.clear();
//     res_hist.push_back(1.0); // ||r0|| / ||r0||
//     int kmax = 10 * n;
//     nb_iter = 0;

//     while (norm_r > tol * norm_r0 && nb_iter < kmax){
//         // Ar = A r
//         vector<double> Ar = (*this) * r;
//         // calcul alpha_k
//         double num = 0.0, den = 0.0;
//         for (int i = 0; i < n; ++i) {
//             num += r[i] * r[i];
//             den += r[i] * Ar[i];
//         }
//         if (den <= 0.0)
//             throw runtime_error("Matrice non définie positive");
//         double alpha = num / den;
//         // x_{k+1}
//         for (int i = 0; i < n; ++i)
//             x[i] += alpha * r[i];
//         // r_{k+1}
//         for (int i = 0; i < n; ++i)
//             r[i] -= alpha * Ar[i];
//         // norme du résidu
//         norm_r = 0.0;
//         for (double v : r) norm_r += v * v;
//         norm_r = sqrt(norm_r);
//         res_hist.push_back(norm_r / norm_r0);
//         // test overflow
//         if (norm_r >= norm_r0 / tol) {
//             cerr << "OVERFLOW dans steepest_descent" << endl;
//             break;
//         }
//         nb_iter++;
//     }
// }

void matricebande::steepest_descent(const vector<double>& b,vector<double>& x,double tol,int& nb_iter,vector<double>& res_hist) const{
    int n = this->dim1();
    if ((int)b.size() != n || (int)x.size() != n)
        throw runtime_error("Erreur de dimension dans steepest_descent");

    // r0 = b - A x0
    vector<double> Ax = this->prod(x);
    vector<double> r(n);
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - Ax[i];

    double norm_r0 = 0.0;
    for (double v : r) norm_r0 += v * v;
    norm_r0 = sqrt(norm_r0);

    double norm_r = norm_r0;

    res_hist.clear();
    res_hist.push_back(1.0);

    int kmax = 10 * n;
    nb_iter = 0;

    while (norm_r > tol * norm_r0 && nb_iter < kmax)
    {
        // y = A r_k
        vector<double> y = this->prod(r);

        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; ++i) {
            num += r[i] * r[i];
            den += r[i] * y[i];
        }

        if (den <= 0.0)
            throw runtime_error("Matrice non définie positive");

        double alpha = num / den;

        for (int i = 0; i < n; ++i)
            x[i] += alpha * r[i];

        for (int i = 0; i < n; ++i)
            r[i] -= alpha * y[i];

        norm_r = 0.0;
        for (double v : r) norm_r += v * v;
        norm_r = sqrt(norm_r);

        res_hist.push_back(norm_r / norm_r0);

        if (norm_r >= norm_r0 / tol) {
            cerr << "OVERFLOW dans steepest_descent" << endl;
            break;
        }

        nb_iter++;
    }
}


void matricebande::gradient_pas_fixe(const vector<double>& b,vector<double>& x,double alpha,double tol,int& nb_iter,vector<double>& res_hist) const{
    int n = this->dim1();
    if ((int)b.size() != n || (int)x.size() != n)
        throw runtime_error("Erreur de dimension dans gradient_pas_fixe");

    vector<double> r(n);
    vector<double> Ar(n);

    // r0 = b - A x0
    vector<double> Ax = this->prod(x);
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - Ax[i];

    double norm_r0 = 0.0;
    for (double v : r) norm_r0 += v * v;
    norm_r0 = sqrt(norm_r0);

    double norm_r = norm_r0;
    res_hist.clear();
    res_hist.push_back(1.0);

    int kmax = 10 * n;
    nb_iter = 0;

    while (norm_r > tol * norm_r0 && nb_iter < kmax)
    {
        // Ar = A r
        Ar = this->prod(r);

        // x_{k+1} = x_k + alpha r_k
        for (int i = 0; i < n; ++i)
            x[i] += alpha * r[i];

        // r_{k+1} = r_k - alpha A r_k
        for (int i = 0; i < n; ++i)
            r[i] -= alpha * Ar[i];

        norm_r = 0.0;
        for (double v : r) norm_r += v * v;
        norm_r = sqrt(norm_r);

        res_hist.push_back(norm_r / norm_r0);

        if (norm_r >= norm_r0 / tol) {
            cerr << "OVERFLOW dans gradient_pas_fixe\n";
            break;
        }

        nb_iter++;
    }
}


void matricebande::assemblageT(int choix, double omega){

    int n = dim1();
    if (n <= 0)
        throw std::runtime_error("assemblageT : matrice vide");

    // ======================================================
    // CHOIX 1 : C = I  ==>  T = I
    // ======================================================
    if (choix == 1)
    {
        TIND.clear();
        TIND.push_back(0);          // diagonale principale seule

        TTAB.resize(n, 1);
        TTAB.zero();

        for (int i = 0; i < n; ++i)
            TTAB(i, 0) = 1.0;

        return;
    }
    // ======================================================
    // CHOIX 2 : C = diag(A)  ==>  T = sqrt(diag(A))
    // ======================================================
    if (choix == 2)
    {
        TIND.clear();
        TIND.push_back(0);

        TTAB.resize(n, 1);
        TTAB.zero();

        int diagA = dim2() - 1;     // diagonale principale de A

        for (int i = 0; i < n; ++i)
        {
            double Aii = (*this)(i, diagA);
            if (Aii <= 0.0)
                throw std::runtime_error("assemblageT : Aii <= 0");

            TTAB(i, 0) = std::sqrt(Aii);
        }

        return;
    }
    // ======================================================
    // CHOIX 3 : IC(0)
    // ======================================================
    if (choix == 3)
    {
        int d = ind().size();
        if (d == 0)
            throw std::runtime_error("assemblageT : structure bande vide");

        TIND = ind();                 // même structure bande que A

        TTAB.resize(n, d);
        TTAB.zero();

        int diag = d - 1;           // position de la diagonale principale

        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;

            // éléments hors diagonale (j < i)
            for (int k = 0; k < diag; ++k)
            {
                int j = i + TIND[k];
                if (j >= 0 && j < i)
                {
                    double Aij = (*this)(i, k);
                    double Tjj = TTAB(j, diag);

                    if (Tjj == 0.0)
                        throw std::runtime_error("assemblageT IC(0) : Tjj = 0");

                    double Tij = Aij / Tjj;
                    TTAB(i, k) = Tij;
                    sum += Tij * Tij;
                }
            }

            double Aii = (*this)(i, diag);
            double val = Aii * (1.0 + omega) - sum;

            if (val <= 0.0)
                throw std::runtime_error("assemblageT IC(0) : racine negative");

            TTAB(i, diag) = std::sqrt(val);
        }

        return;
    }
    // ======================================================
    // CHOIX INVALIDE
    // ======================================================
    throw std::runtime_error("assemblageT : choix invalide");
}

vector<double> matricebande::preconditionne(const vector<double>& r) const{
    int n = dim1();
    assert((int)r.size() == n);

    vector<double> y(n, 0.0);
    vector<double> z(n, 0.0);

    int d = TIND.size();
    int diag = d - 1;   // position de la diagonale principale

    // ==================================================
    // 1) Descente : T y = r
    // ==================================================
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;

        for (int k = 0; k < diag; ++k)
        {
            int j = i + TIND[k];
            if (j >= 0 && j < i)
                sum += TTAB(i, k) * y[j];
        }

        double Tii = TTAB(i, diag);
        assert(Tii != 0.0);

        y[i] = (r[i] - sum) / Tii;
    }

    // ==================================================
    // 2) Remontée : Tᵗ z = y
    // ==================================================
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0.0;

        for (int k = 0; k < diag; ++k)
        {
            int j = i - TIND[k];   // car TIND[k] < 0
            if (j > i && j < n)
                sum += TTAB(j, k) * z[j];
        }

        double Tii = TTAB(i, diag);
        assert(Tii != 0.0);

        z[i] = (y[i] - sum) / Tii;
    }

    return z;
}

void matricebande::PCG(const vector<double>& b,vector<double>& x,double tol,int& nb_iter,vector<double>& res_hist) const{
    int n = dim1();
    assert((int)b.size() == n);
    assert((int)x.size() == n);

    // r = b - A x
    vector<double> r = b - this->prod(x);

    // norme de b
    double normb = 0.0;
    for (double bi : b) normb += bi * bi;
    normb = std::sqrt(normb);

    // z = C^{-1} r
    vector<double> z = preconditionne(r);

    // p = z
    vector<double> p = z;

    // gamma = z^T r
    double gamma = 0.0;
    for (int i = 0; i < n; ++i)
        gamma += z[i] * r[i];

    double gamma0 = gamma;

    res_hist.clear();
    res_hist.push_back(std::sqrt(gamma) / normb);

    nb_iter = 0;

    // ==============================
    // Boucle PCG
    // ==============================
    while (gamma > tol * tol * gamma0)
    {
        // q = A p
        vector<double> q = this->prod(p);

        // alpha = gamma / (p^T q)
        double pq = 0.0;
        for (int i = 0; i < n; ++i)
            pq += p[i] * q[i];

        assert(pq != 0.0);

        double alpha = gamma / pq;

        // x = x + alpha p
        for (int i = 0; i < n; ++i)
            x[i] += alpha * p[i];

        // r = r - alpha q
        for (int i = 0; i < n; ++i)
            r[i] -= alpha * q[i];

        // z = C^{-1} r
        z = preconditionne(r);

        // gamma_new = z^T r
        double gamma_new = 0.0;
        for (int i = 0; i < n; ++i)
            gamma_new += z[i] * r[i];

        // historique des résidus
        res_hist.push_back(std::sqrt(gamma_new) / normb);

        // beta = gamma_new / gamma
        double beta = gamma_new / gamma;

        // p = z + beta p
        for (int i = 0; i < n; ++i)
            p[i] = z[i] + beta * p[i];

        gamma = gamma_new;
        nb_iter++;
    }
}

vector<double> smbr_rho_bas_seul(int n)
{
    int N = (n + 2) * n;
    double h  = 1.0 / (n + 1);
    double hh = h * h;
    vector<double> b(N, 0.0);

    for (int i = 0; i < n + 2; ++i)
    {
        double x = i * h;
        double y = 0.0;
        double rho_jm = rho(x, y + 0.5 * h);
        int k = i;
        b[k] = rho_jm * 20.0 / hh;   // bas seul
    }
    return b;
}

vector<double> smbr_rho_haut_seul(int n)
{
    int N = (n + 2) * n;
    double h  = 1.0 / (n + 1);
    double hh = h * h;
    vector<double> b(N, 0.0);

    for (int i = 0; i < n + 2; ++i)
    {
        double x = i * h;
        double y = n * h;
        double rho_jp = rho(x, y - 0.5 * h);
        int k = i + (n - 1) * (n + 2);
        b[k] = rho_jp * 180.0 / hh;   // haut seul
    }
    return b;
}
