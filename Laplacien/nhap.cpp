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
    
    int N=(n+1)*(n+1);
    double h = 1./(n+1);
    double hh = h * h;

    vector<int> v;
    v.push_back(-(n+1));
    v.push_back(-1);
    v.push_back(0);
    
    matricebande A(N,3);
    A.ind(v);
        
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i) {

            int ligne = i + j * (n + 1);

            // Diagonale principale
            A(ligne, 2) = 4.0 / hh;

            // Voisin gauche : T_{i-1,j}
            if (i > 0)
                A(ligne, 1) = -1.0 / hh;

            // Voisin bas : T_{i,j-1}
            if (j > 0)
                A(ligne, 0) = -1.0 / hh;
        }
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
    // n = N dans l'énoncé
    int N = (n + 1) * (n + 1);

    double h  = 1.0 / (n + 1);
    double hh = h * h;

    matricebande A(N, 3);

    // diagonales de la partie triangulaire inférieure
    vector<int> IND;
    IND.push_back(-(n + 1)); // voisin bas
    IND.push_back(-1);       // voisin gauche
    IND.push_back(0);        // diagonale principale
    A.ind(IND);

    // Boucle sur tous les points (i,j)
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i) {

            int k = i + j * (n + 1);

            double x = i * h;
            double y = j * h;

            // rho aux demi-mailles
            double rho_ip = rho(x + 0.5*h, y);
            double rho_im = (i > 0) ? rho(x - 0.5*h, y) : rho_ip;

            double rho_jp = rho(x, y + 0.5*h);
            double rho_jm = (j > 0) ? rho(x, y - 0.5*h) : rho_jp;

            // diagonale principale
            A(k,2) = (rho_ip + rho_im + rho_jp + rho_jm) / hh;

            // voisin gauche
            if (i > 0)
                A(k,1) = -rho_im / hh;

            // voisin bas
            if (j > 0)
                A(k,0) = -rho_jm / hh;
        }
    }

    return A;
}

// (1.7)
vector<double> matricebande::prod(const vector<double>& x) const{
    int n = this->dim1();

    // Vérification de la taille
    if (x.size() != n)
        throw std::runtime_error("Erreur: taille incompatible dans A * x");

    vector<double> y(n);
    fill(y.begin(), y.end(), 0.0);

    int d = indice.size();  // nombre de diagonales stockées

    for (int i = 0; i < n; ++i)
    {
        for (int k = 0; k < d; ++k)
        {
            int j = i + indice[k];

            if (j < 0 || j >= n)
                continue;

            double aij = (*this)(i, k);  // coefficient stocké

            // Contribution à y_i
            y[i] += aij * x[j];

            // Contribution symétrique à y_j
            if (j != i)
                y[j] += aij * x[i];
        }
    }

    return y;
}


//(1.9)
vector<double> smbr(int N)
{
    int n = (N + 1) * (N + 1);
    vector<double> b(n, 0.0);

    // Dirichlet bas : T = 20
    for (int i = 0; i <= N; ++i)
    {
        int k = i;           // j = 0
        b[k] = 20.0;
    }

    // Dirichlet haut : T = 180
    for (int i = 0; i <= N; ++i)
    {
        int k = i + N * (N + 1);   // j = N
        b[k] = 180.0;
    }

    return b;
}

vector<double> smbr_rho(int N)
{
    int n = (N + 1) * (N + 1);
    vector<double> b(n, 0.0);

    double h  = 1.0 / (N + 1);
    double hh = h * h;

    // Dirichlet bas : T = 20
    for (int i = 0; i <= N; ++i)
    {
        double x = i * h;
        double y = 0.0;

        double rho_jm = rho(x, y - 0.5 * h);

        int k = i;
        b[k] = rho_jm * 20.0 / hh;
    }

    // Dirichlet haut : T = 180
    for (int i = 0; i <= N; ++i)
    {
        double x = i * h;
        double y = N * h;

        double rho_jp = rho(x, y + 0.5 * h);

        int k = i + N * (N + 1);
        b[k] = rho_jp * 180.0 / hh;
    }

    return b;
}