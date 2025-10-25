#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "vecteur_template.h" 

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

// on cree vecteur X de taille M=N+1 car on inclut les deux bornes 0 et L 
vector<double> generer_X(double L, int N) {
    int M = N+1;
    vector<double> X(M);
    double h = L / double(N); // note: h = L/N
    for(int i=0;i<M;i++) X[i] = i * h;
    return X;
}

// strike K(x)
double K_of_x(double x) { return 0.95 * x; }

// payoff initial h(x) = max(x - K(x), 0)
vector<double> u0_call(const vector<double>& X) {
    int M = X.size();
    vector<double> u0(M);
    for(int i=0;i<M;i++){
        u0[i] = max(X[i] - K_of_x(X[i]), 0.0);
    }
    return u0;
}

// Matrice bande pour le schéma implicite Black-Scholes
// M = size M (number of points), h space step, k time step, sigma vol, X positions
mat_bande matrice_BS_bande(const vector<double>& X, double h, double k, double sigma){
    int M = X.size();
    mat_bande A;
    A.n = M;
    A.IND = { -1, 0, 1 };
    A.TAB.assign(M, vector<double>(3, 0.0));

    
    for(int i=0;i<M;i++){
        double nu = 0.5 * sigma * sigma * X[i] * X[i]; // nu_i
        double factor = k * nu / (h*h);
        if(i==0){
            // Dirichlet sur x=0: u=0 
            A.TAB[i][0] = 0.0;
            A.TAB[i][1] = 1.0;
            A.TAB[i][2] = 0.0;
        } else if(i == M-1){
            // Neumann sur x=L: du/dx=0  => u_{M-1} = u_{M-2}
            // on utilise une condition fantome pour eliminer u_{M}
            A.TAB[i][0] = -2.0 * factor;      // sub-diagonal (connect to M-2)
            A.TAB[i][1] = 1.0 + 2.0 * factor; // diagonal
            A.TAB[i][2] = 0.0;                // no super-diagonal
        } else {
            // standard interior
            A.TAB[i][0] = -factor;            // a_i
            A.TAB[i][1] = 1.0 + 2.0*factor;   // b_i
            A.TAB[i][2] = -factor;            // c_i
        }
    }
    return A;
}

vector<double> BS_implicite(const vector<double>& X, double L, double T, double sigma, int N, double h, double k, const vector<double>& u0, const string& file_out){
    int M = X.size();
    assert(M == N+1);
    vector<double> U = u0;

    mat_bande Mmat = matrice_BS_bande(X, h, k, sigma);

    mat_bande Mcopy = Mmat;
    factorisation_LU_tridiag(Mcopy);

    int Ktime = static_cast<int>(round(T / k)); 
    ofstream fout(file_out, ios::out | ios::trunc);
    if(!fout.is_open()){
        cerr << "Impossible d'ouvrir " << file_out << endl;
        return U;
    }
    fout << "# t ";
    for(int i=0;i<M;i++) fout << " x_" << i;
    fout << "\n";
    for(int n=0;n<=Ktime;n++){
        double t = n * k;
        fout << t;
        for(int i=0;i<M;i++) fout << " " << U[i];
        fout << "\n";
        if(n < Ktime){
            vector<double> B = U;
            B[0] = 0.0;
            U = resoudre_LU_tridiag(Mcopy, B);
            U[0] = 0.0;
            U[M-1] = U[M-2];
        }
    }
    fout.close();
    return U;
}

int main(){
    double L = 10000.0;   
    double T = 365.0;     
    double sigma = 0.2;   
    int N = 200;          
    double h = L / double(N+1); 
    double k = h;        

    vector<double> X = generer_X(L, N);   
    vector<double> u0 = u0_call(X);

    string file_out = "BS_implicit.txt";
    vector<double> Ufinal = BS_implicite(X, L, T, sigma, N, h, k, u0, file_out);

    cout << " Solution calculée au temps final T = " << T << " (Resultats sauvegardes dans " << file_out << ")\n";
    cout << "x   U(T,x)\n";
    for(size_t i=0;i<X.size();i++){
        cout << X[i] << "   " << Ufinal[i] << "\n";
    }

    return 0;
}
