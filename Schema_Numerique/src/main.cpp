#include<iostream> 
#include <fstream>
#include <assert.h>
#include <cmath>
#include <string>
#include <vector>
#include "vecteur_template.h" // doit etre dans le meme repertoire
#include "Schema_numerique.h"
using namespace std;

int main() {
    double Tfinal = 1.0;

    // ------------------- Euler explicite -------------------
    string fichier_exp = "erreur_explicite.txt";
    ofstream fout1(fichier_exp);
    fout1 << "N \tErreur_max\n";
    fout1.close();

    vector<int> Ns = {4, 8, 16, 32, 64};
    for(int N : Ns){
        cout << "[Explicite] Résolution pour N = " << N << endl;
        verif_chaleur_explicite(Tfinal, N, fichier_exp);
    }

    // ------------------- Euler implicite -------------------
    string fichier_imp = "erreur_implicite.txt";
    ofstream fout2(fichier_imp);
    fout2 << "N \tErreur_max\n";
    fout2.close();

    for(int N : Ns){
        cout << "[Implicite] Résolution pour N = " << N << endl;
        verif_chaleur_implicite(Tfinal, N, fichier_imp);
    }

    cout << "\nRésultats sauvegardés dans : " << fichier_exp << " et " << fichier_imp << endl;
    return 0;
}
