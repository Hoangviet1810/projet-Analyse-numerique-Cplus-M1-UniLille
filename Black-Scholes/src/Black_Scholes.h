#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H

#include <vector>
#include <string>

// ===== structure =====
struct mat_bande {
    int n;
    std::vector<int> IND;
    std::vector<std::vector<double>> TAB;
};

// ===== fonctions matrices =====
mat_bande matrice_1D_bande(int N, double c0, double c1);
void factorisation_LU_tridiag(mat_bande& A);

template <typename T>
std::vector<T> resoudre_LU_tridiag(const mat_bande& A,
                                  const std::vector<T>& b);

// ===== utilitaires =====
std::vector<double> generer_X(double L, int N);
double K_of_x(double x);
std::vector<double> u0_call(const std::vector<double>& X);

// ===== Black–Scholes =====
mat_bande matrice_BS_bande(const std::vector<double>& X,double h,double k,double sigma);

std::vector<double> BS_implicite(const std::vector<double>& X,double L,double T,double sigma,int N,double h,double k,const std::vector<double>& u0,const std::string& file_out);

#endif // BLACK_SCHOLES_H
