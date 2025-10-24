#include<iostream> 
#include <fstream>
#include <assert.h>
#include <cmath>
#include <string>
#include <vector>
#include "vecteur_template.h" // doit etre dans le meme repertoire

using namespace std;


struct mat_bande{
    int n;
    vector<int> IND;
    vector<vector<double>> TAB;
};

mat_bande matrice_1D_bande(int N, double c0, double c1) {
    mat_bande A;
    A.n = N;
    int d = 3;  // tridiagonal
    A.IND.resize(d);
    A.TAB.resize(N, vector<double>(d,0.0));

    for(int k=0;k<d;k++)
        A.IND[k] = k - (d-1)/2;  // -1,0,1

    // Tab
    for(int i=0;i<N;i++){
        for(int k=0;k<d;k++){
            int j = i + A.IND[k]; 
            if(j>=0 && j<N){
                if(k==0) A.TAB[i][k] = c1; // sub-diagonal
                else if(k==1) A.TAB[i][k] = c0; // diagonal
                else if(k==2) A.TAB[i][k] = c1; // super-diagonal
            }
            else A.TAB[i][k] = 0; 
        }
    }

    return A;
}

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

void factorisation_LU_tridiag(mat_bande& A){
    int n = A.TAB.size();
    
    A.TAB[0][1] = A.TAB[0][1]; // beta_1
    A.TAB[0][2] = A.TAB[0][2]; // gamma_1
    if(n>1) A.TAB[1][0] = A.TAB[1][0] / A.TAB[0][1]; // alpha_2

    for(int i=1; i<n; i++){
        double gamma_prev = A.TAB[i-1][2];                   // gamma_{i-1} = c_{i-1}
        double alpha_i = A.TAB[i][0];                        // alpha_i = A.TAB[i][0] 
        A.TAB[i][1] = A.TAB[i][1] - alpha_i * gamma_prev;    // beta_i = b_i - alpha_i * gamma_{i-1}
        if(i < n-1)                                          // alpha_{i+1} = a_{i+1} / beta_i 
            A.TAB[i+1][0] = A.TAB[i+1][0] / A.TAB[i][1];
        A.TAB[i][2] = A.TAB[i][2];                           // gamma_i = c_i 
    }
}

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

double v(double x){
    if(x < 0.5){
         return 1.0;
    }
    else return 10.0;
}

double uex(double t, double x) {
    return cos(t) * cos(M_PI*x);  
}

double u0(double x) {
    return uex(0.0, x);
}

double f(double t, double x) {
    double v_val = v(x);
    return -sin(t) * cos(M_PI * x)
           + M_PI * M_PI * v_val * cos(t) * cos(M_PI * x);
}

double g0(double t) {
    return cos(t) * cos(0.0);
}

double g1(double t) {
    return cos(t) * cos(M_PI*1.0);
}

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


void verif_chaleur_implicite(double Tfinal, int N, const string& file) {
    double h = 1.0 / (N + 1);
    double c0 = 2.0/(h*h);
    double c1 = -1.0/(h*h);
    
    vector<double> u0_vec(N);
    vector<double> nu_vec(N);
    for(int i=0;i<N;i++){
        double xi = i*h;
        u0_vec[i] = u0(xi);
        nu_vec[i] = v(xi);
    }

    mat_bande A = matrice_1D_bande(N, c0, c1);
    vector<double> U = chaleur_implicite(A, N, h, nu_vec, Tfinal, u0_vec);  // U(T)

    vector<double> Uex(N);
    double errmax = 0.0;
    for (int j = 0; j < N; ++j) {
        double xj = (j + 1) * h;
        Uex[j] = uex(Tfinal, xj);
        errmax = max(errmax, fabs(U[j] - Uex[j]));
    }
    cout << " Erreur maximale = " << errmax << endl;

    // Écriture dans un fichier
    ofstream fout(file, ios::app); // append
    if (fout.is_open()) {
        fout << N << "    \t  " << errmax << endl;
        fout.close();
    } else {
        cerr << "Erreur: impossible d'ouvrir le fichier " << file << endl;
    }
}

// Euler explicite (Sexpl) utilisant prod_mat_vect
// 
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



// template <typename T>
// vector<T> chaleur_explicite(const mat_bande& A, int N, double h, const vector<T>& nu, double Tfinal, const vector<T>& u0) {
//     double k= (h*h)/(2.0*(*max_element(nu.begin(), nu.end()))); // pas de temps stable
//     int K = static_cast<int>(Tfinal / k);
//     vector<T> U = u0;

//     for(int n=0; n<K; n++){
//         double t_n = n*k;
//         vector<T> Fn(N);
//         for(int i=0;i<N;i++) Fn[i] = f(t_n, (i+1)*h);

//         // calcul (nu * k * A) * U^n
//         mat_bande L = A;
//         for(int i=0;i<N;i++){
//             L.TAB[i][1] = 2.0*nu[i]*k/(h*h);
//             if(i>0) L.TAB[i][0] = -nu[i]*k/(h*h);
//             if(i<N-1) L.TAB[i][2] = -nu[i]*k/(h*h);
//         }

//         vector<T> LU = prod_mat_vect(L, U);

//         // U^{n+1} = U^n - (nu k A) U^n + k F^n
//         for(int i=0;i<N;i++) U[i] = U[i] - LU[i] + k*Fn[i];

//         // conditions aux bords
//         U[0] = g0(t_n+k);
//         U[N-1] = g1(t_n+k);
//     }

//     return U;
// }


void verif_chaleur_explicite(double Tfinal, int N, const string& file) {
    double h = 1.0 / (N + 1);
    double c0 = 2.0/(h*h);
    double c1 = -1.0/(h*h);
    
    vector<double> u0_vec(N);
    vector<double> nu_vec(N);
    for(int i=0;i<N;i++){
        double xi = i*h;
        u0_vec[i] = u0(xi);
        nu_vec[i] = v(xi);
    }

    mat_bande A = matrice_1D_bande(N, c0, c1);
    vector<double> U = chaleur_explicite(A, N, h, nu_vec, Tfinal, u0_vec);  // U(T)

    vector<double> Uex(N);
    double errmax = 0.0;
    for (int j = 0; j < N; ++j) {
        double xj = (j + 1) * h;
        Uex[j] = uex(Tfinal, xj);
        errmax = max(errmax, fabs(U[j] - Uex[j]));
    }
    cout << " Erreur maximale = " << errmax << endl;

    // Écriture dans un fichier
    ofstream fout(file, ios::app); // append
    if (fout.is_open()) {
        fout << N << "    \t  " << errmax << endl;
        fout.close();
    } else {
        cerr << "Erreur: impossible d'ouvrir le fichier " << file << endl;
    }
}

 

int main() {
   
    double Tfinal = 1.0; 
    string fichier = "erreur2.txt";

    ofstream fout(fichier);
    fout << "N \t Erreur_max_implicite \n";
    fout.close();

    vector<int> Ns = {8, 16, 32, 64, 128, 256};
    
    for (int N : Ns) {
        cout << "Résolution pour N = " << N  << endl;
        //verif_chaleur_explicite(Tfinal, N, fichier);
        verif_chaleur_implicite(Tfinal, N, fichier);
    }

    cout << "\nRésultats sauvegardés dans : " << fichier << endl;


    return 0;
}