#include<iostream> 
#include <fstream>
#include <assert.h>
#include <cmath>
#include <string>
#include <vector>
#include "vecteur_template.h" 

using namespace std;

struct mat_profil {
    vector<int> INDIAG;     
    vector<double> ATAB;    
};

mat_profil matrice_1D_profil(int N, double c0, double c1) {
    mat_profil A;
    A.INDIAG.resize(N);
    A.ATAB.resize(3*N); // taille provisoire

    int k = 0;

    // première ligne
    A.ATAB[k++] = c0; // diagonale principale
    A.INDIAG[0] = 0;

    // lignes suivantes
    for (int i = 1; i < N-1; i++) {
        A.ATAB[k++] = c1; // sous-diagonal
        A.ATAB[k++] = c0; // diagonale principale
        A.INDIAG[i] = k-1; // index de la diagonale principale pour cette ligne
    }
    // dernière ligne de numero N-1
    A.ATAB[k++] = c1; // sous-diagonal
    for (int i = 1; i < N-2; i++) {
        A.ATAB[k++] = 0; // diagonale principale
    }   
    A.ATAB[k++] = c1;
    A.ATAB[k++] = c0;
    A.INDIAG[N-1] = A.INDIAG[N-2]+ N; // dernier élément sous-
    A.ATAB.resize(k); // ajuster la taille
    return A;
}


void factorisation_LDLt(mat_profil& A){
    
    int n = A.INDIAG.size();
    for (int i = 0; i < n; i++) {
        int ji;
        if (i == 0) ji = 0;
        else ji = i - A.INDIAG[i] + A.INDIAG[i-1] + 1;

        for (int j = ji; j <= i - 1; j++) {
            int jj;
            if (j == 0) jj = 0;
            else jj = j - A.INDIAG[j] + A.INDIAG[j-1] + 1;

            // etap 5–7
            for (int k = max(ji, jj); k <= j - 1; k++) {
                A.ATAB[A.INDIAG[i] - i + j] -= A.ATAB[A.INDIAG[i] - i + k] * A.ATAB[A.INDIAG[j] - j + k];
            }
            // etap 9–13
            double S = A.ATAB[A.INDIAG[i] - i + j] / A.ATAB[A.INDIAG[j]];
            A.ATAB[A.INDIAG[i]] -= S * A.ATAB[A.INDIAG[i] - i + j];
            A.ATAB[A.INDIAG[i] - i + j] = S;
        }
    } // A  maintenant contient L et D (LDLt)
}

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

template <typename T>
void verif_poisson(const mat_profil& A,const vector<T>& b,int N,double h, double alpha, const string& file) {
    
    mat_profil Acopy = A; 
    factorisation_LDLt(Acopy);
    vector<double> x = resoudre_LDLt(Acopy, b);
    
    // cout << "i  \t  x_i  \t  Approx  \t  exact  \t  erreur" << endl;
    double erreur_max = 0.0;
    for (int i = 0; i < N; i++) {
        double x_exact = sin(2*M_PI*i*h);
        double err = fabs(x[i] - x_exact);
        erreur_max = max(erreur_max, err);
        // cout << i << "\t" << i*h << "  \t  " << x[i] << "  \t  " << x_exact << "   \t   " << err << endl;
    }
    cout << "Erreur maximale = " << erreur_max << endl;
    
    // Écriture dans un fichier
    ofstream fout(file, ios::app); // append
    if (fout.is_open()) {
        fout << N << "\t" << erreur_max << endl;
        fout.close();
    } else {
        cerr << "Erreur: impossible d'ouvrir le fichier " << file << endl;
    }
}

// int main() {
//     double alpha = 0.0;              
//     string fichier = "erreur.txt";
//     ofstream fout(fichier);
//     fout << "N \tErreur_max s\n";
//     fout.close();
//     // Tests pour différentes tailles de matrice
//     vector<int> Ns = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
//     for (int N : Ns) {
//         double h = 1.0 / N;
//         double c0 = 2.0/(h*h) + alpha;
//         double c1 = -1.0/(h*h);
        
//         mat_profil A = matrice_1D_profil(N, c0, c1);
//         vector<double> b(N);
//         for (int i = 0; i < N; i++) {
//             double x_i = i*h;
//             b[i] = (4*M_PI*M_PI + alpha) * sin(2*M_PI*x_i);
//         }
//         cout << "Résolution pour N = " << N << endl;
//         verif_poisson(A,b, N,h, alpha, fichier);
//     }
//     cout << "\nRésultats sauvegardés dans : " << fichier << endl;
//     return 0; 
// }