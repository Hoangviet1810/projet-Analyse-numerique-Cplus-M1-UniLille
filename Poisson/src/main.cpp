#include<iostream> 
#include <fstream>
#include <assert.h>
#include <cmath>
#include <string>
#include <vector>
#include "vecteur_template.h" 
#include "Poisson.h"

using namespace std;

int main() {
    // ----------------------- alpha = 0 -----------------------
    double alpha1 = 0.0;
    string fichier1 = "erreur.txt";
    
    // Création du fichier avec header
    ofstream fout1(fichier1);
    fout1 << "N \tErreur_max\n";
    fout1.close();

    vector<int> Ns1 = {4, 8, 16, 32, 64, 128, 256, 512}; // etc
    for(int N : Ns1){
        double h = 1.0 / N;
        double c0 = 2.0/(h*h) + alpha1;
        double c1 = -1.0/(h*h);

        mat_profil A = matrice_1D_profil(N, c0, c1);
        vector<double> b(N);
        for(int i=0;i<N;i++) b[i] = (4*M_PI*M_PI + alpha1)*sin(2*M_PI*i*h);

        verif_poisson(A, b, N, h, alpha1, fichier1); // passe le fichier spécifique
    }

    // ----------------------- alpha = 10 -----------------------
    double alpha2 = 10.0;
    string fichier2 = "erreur_alpha10.txt";

    ofstream fout2(fichier2);
    fout2 << "N \tErreur_max\n";
    fout2.close();

    vector<int> Ns2 = {4, 8, 16, 32, 64, 128, 256, 512}; // etc
    for(int N : Ns2){
        double h = 1.0 / N;
        double c0 = 2.0/(h*h) + alpha2;
        double c1 = -1.0/(h*h);

        mat_profil A = matrice_1D_profil(N, c0, c1);
        vector<double> b(N);
        for(int i=0;i<N;i++) b[i] = (4*M_PI*M_PI + alpha2)*sin(2*M_PI*i*h);

        verif_poisson(A, b, N, h, alpha2, fichier2);
    }

    cout << "\nRésultats sauvegardés dans : " << fichier1 << " et " << fichier2 << endl;
    return 0;
}
