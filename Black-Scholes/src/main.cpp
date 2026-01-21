#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "vecteur_template.h" 
#include "Black_Scholes.h"

using namespace std;

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