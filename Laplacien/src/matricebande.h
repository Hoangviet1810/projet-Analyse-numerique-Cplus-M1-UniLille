#ifndef MATRICEBANDE_H
#define MATRICEBANDE_H

#include "matrice.h"
#include <vector>
#include <iostream>
#include <stdexcept>    

class matricebande : public matrice{
private:
    vector<int> indice;
    vector<int> TIND;
    matrice TTAB;
public:
    matricebande() {};                         // constructeur par defaut
    matricebande(int,int);                     // constructeur en donnant 2 dimensions
    matricebande(const matricebande&);         // constructeur par recopie
    ~matricebande() {};                        // destructeur par defaut
    
    vector<int> ind() const {return this->indice;}    // retourne le champ privé indice
    vector<int> ind(vector<int>);                     // remplit le champ privé indice

    matricebande laplacien(int);                      // assemblage du laplacien TP3, question (1.1)
    matricebande laplacien_rho(int);
    vector<double> prod(const vector<double>& x) const;


    friend vector<double> smbr(int N);
    friend vector<double> smbr_rho(int N);

    void steepest_descent(const vector<double>& b, vector<double>& x, double tol, int& nb_iter, vector<double>& res_hist) const;

    void gradient_pas_fixe(const vector<double>& b,vector<double>& x,double alpha,double tol,int& nb_iter,vector<double>& res_hist) const;

    void assemblageT(int choix, double omega = 0.0);
    const vector<int>& get_TIND() const { return TIND; } // co the xoa luc sau
    const matrice& get_TTAB() const { return TTAB; }

    vector<double> preconditionne(const vector<double>& r) const;


    void PCG(const vector<double>& b,vector<double>& x,double tol,int& nb_iter,vector<double>& res_hist) const;

};

// Déclarations des seconds membres (1.9)
vector<double> smbr(int N);
vector<double> smbr_rho(int N);

vector<double> smbr_rho_bas_seul(int n);
vector<double> smbr_rho_haut_seul(int n);

#endif // MATRICEBANDE_H
