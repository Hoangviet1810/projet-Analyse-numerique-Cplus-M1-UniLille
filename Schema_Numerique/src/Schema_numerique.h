#ifndef SCHEMA_NUMERIQUE_H
#define SCHEMA_NUMERIQUE_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "vecteur_template.h" // phải nằm cùng thư mục

using namespace std;

// ======= Struct =======
struct mat_bande {
    int n;
    vector<int> IND;                 // indices diagonale (-1,0,1)
    vector<vector<double>> TAB;      // tableau des coefficients
};

// ======= Fonctions =======

// Construction matrice bande
mat_bande matrice_1D_bande(int N, double c0, double c1);

// Produit matrice bande * vecteur
template <typename T>
vector<T> prod_mat_vect(const mat_bande& A, const vector<T>& x){
    int nd = A.IND.size();
    int n = A.TAB.size();
    assert(x.size() == n);
    vector<T> y(n, 0.0); // initialiser le vecteur resultat
    
    for ( int i=0; i<n; i++){
        int j;
        j = i + A.IND[0];
        for(int k = max(0, -j);k<= min(nd-1,n-1-j); k++){
            y[i] += A.TAB[i][k] * x[i+A.IND[k]];
        }
    }
    return y;
}
// Factorisation LU tridiag
void factorisation_LU_tridiag(mat_bande& A);

// Résolution système tridiagonal LU
template <typename T>
vector<T> resoudre_LU_tridiag(const mat_bande& A, const vector<T>& b){
    int n = A.TAB.size();
    vector<T> x = b;
    for(int i=1; i<n; i++){
        x[i] -= A.TAB[i][0] * x[i-1];
    }
    x[n-1] /= A.TAB[n-1][1];
    for(int i=n-2; i>=0; i--){
        x[i] = (x[i] - A.TAB[i][2] * x[i+1]) / A.TAB[i][1];
    }
    return x;
}
// Fonctions exactes / conditions aux bords / source
double v(double x);
double uex(double t, double x);
double u0(double x);
double f(double t, double x);
double g0(double t);
double g1(double t);

// Schéma implicite
template <typename T>
vector< T> chaleur_implicite(const mat_bande&A, int N, double h, const vector<T>& nu, double Tfinal, const vector<T>& u0) {
    double k = h;                        
    int K = static_cast<int>(Tfinal / k); 

    // Initialisation
    vector<T> U = u0;

    // Construction de la matrice bande M = I + νkA
    mat_bande M = A; // copie pour modification
    for (int i = 0; i < N; ++i) {
        double alpha = nu[i] * k / (h * h);
        M.TAB[i][1] = 1.0 + 2.0 * alpha;    // diagonale principale
        if (i > 0)   M.TAB[i][0] = -alpha;  // sous-diagonale
        if (i < N-1) M.TAB[i][2] = -alpha;  // sur-diagonale
    }
   
    M.TAB[0][0] = 0.0; M.TAB[0][1] = 1.0; M.TAB[0][2] = 0.0;
    M.TAB[N-1][0] = 0.0; M.TAB[N-1][1] = 1.0; M.TAB[N-1][2] = 0.0;

    mat_bande M_copy = M; // copie pour factorisation
    factorisation_LU_tridiag(M_copy);

    for (int n = 0; n < K; ++n) {
      vector<T> B(N);
      double t_next = (n+1) * k;   
      for (int i = 0; i < N; ++i) {
          double xi = (i + 1) * h;
          B[i] = U[i] + k * f(t_next, xi);   // f(t^{n+1}, x_i)
       }
      B[0]   = g0(t_next);       // Conditions aux bords
      B[N-1] = g1(t_next);
      U = resoudre_LU_tridiag(M_copy, B);
      U[0]   = g0(t_next);      // impose de nouveau les conditions aux bords
      U[N-1] = g1(t_next);
    }

    return U; // approximation au temps final T = K·k
}

// Schéma explicite
template <typename T>
vector<T> chaleur_explicite(const mat_bande& A, int N, double h, const vector<T>& nu, double Tfinal, const vector<T>& u0) {
    double k= (h*h)/(2.0*(*max_element(nu.begin(), nu.end()))); // pas de temps stable
    //double k =h;
    int K = static_cast<int>(Tfinal / k);
    vector<T> U = u0;

    for(int n=0; n<K; n++){
        double t_n = n*k;
        vector<T> Fn(N);
        for(int i=0;i<N;i++) Fn[i] = f(t_n, (i+1)*h);

        // calcul (nu * k * A) * U^n
        mat_bande L = A;
        for(int i=0;i<N;i++){
            L.TAB[i][1] = 2.0*nu[i]*k/(h*h);
            if(i>0) L.TAB[i][0] = -nu[i]*k/(h*h);
            if(i<N-1) L.TAB[i][2] = -nu[i]*k/(h*h);
        }
        vector<T> LU = prod_mat_vect(L, U);
        // U^{n+1} = U^n - (nu k A) U^n + k F^n
        for(int i=0;i<N;i++) U[i] = U[i] - LU[i] + k*Fn[i];

        // conditions aux bords
        U[0] = g0(t_n+k);
        U[N-1] = g1(t_n+k);
    }

    return U;
}

// Vérification schéma implicite / explicite
void verif_chaleur_implicite(double Tfinal, int N, const string& file);
void verif_chaleur_explicite(double Tfinal, int N, const string& file);

#endif // SCHEMA_NUMERIQUE_H
