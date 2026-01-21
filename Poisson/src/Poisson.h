#ifndef POISSON_H
#define POISSON_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "vecteur_template.h"

using namespace std;

// ==================== Struct ====================
struct mat_profil {
    vector<int> INDIAG;     // indices diagonale principale
    vector<double> ATAB;    // tableau des coefficients
};

// ==================== Fonctions ====================

// Construction matrice profilée tridiagonale
mat_profil matrice_1D_profil(int N, double c0, double c1);

// Factorisation LDLt
void factorisation_LDLt(mat_profil& A);

// Résolution système avec LDLt (template)
// -------------------- template functions --------------------
template <typename T>
vector<T> resoudre_LDLt(const mat_profil& A, const vector<T>& b) {
    int n = A.INDIAG.size();
    assert((int)b.size() == n);

    vector<T> x = b; 
    for (int i = 1; i < n; i++) {
        int ji = i - A.INDIAG[i] + A.INDIAG[i-1] + 1;
        for (int j = ji; j <= i-1; j++) {
            x[i] -= A.ATAB[A.INDIAG[i] - i + j] * x[j];
        }
    }

    for (int i = 0; i < n; i++) {
        x[i] /= A.ATAB[A.INDIAG[i]];
    }

    for (int i = n-1; i >= 1; i--) {
        int ji = i - A.INDIAG[i] + A.INDIAG[i-1] + 1;
        for (int j = ji; j <= i-1; j++) {
            x[j] -= A.ATAB[A.INDIAG[i] - i + j] * x[i];
        }
    }
    return x;
}

// Vérification solution de Poisson (template)
template <typename T>
void verif_poisson(const mat_profil& A,const vector<T>& b,int N,double h, double alpha, const string& file) {
    mat_profil Acopy = A; 
    factorisation_LDLt(Acopy);
    vector<double> x = resoudre_LDLt(Acopy, b);
    
    double erreur_max = 0.0;
    for (int i = 0; i < N; i++) {
        double x_exact = sin(2*M_PI*i*h);
        double err = fabs(x[i] - x_exact);
        erreur_max = max(erreur_max, err);
    }
    
    cout << "Erreur maximale = " << erreur_max << endl;

    ofstream fout(file, ios::app);
    if (fout.is_open()) {
        fout << N << "\t" << erreur_max << endl;
        fout.close();
    } else {
        cerr << "Impossible d'ouvrir le fichier " << file << endl;
    }
}
#endif // POISSON_H
