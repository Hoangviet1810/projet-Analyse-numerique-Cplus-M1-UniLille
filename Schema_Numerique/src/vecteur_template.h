/*
*
*  Created by caterina calgaro on 08/09/2020.
*  Copyright 2020 INRIA LILLE NORD EUROPE. All rights reserved.
*
*/

#if !defined (__IOSTREAM_H)
#include <iostream>
#endif
#if !defined (__FSTREAM_H)
#include <fstream>
#endif
#if !defined (__ASSERT_H)
#include <assert.h>
#endif
#if !defined(__CMATH_H)
#include <cmath>
#endif
#if !defined (__VECTOR_H)
#include <vector>
#endif

#if !defined(VECTEUR_TEMPLATE_H)
#define VECTEUR_TEMPLATE_H

using namespace std;

// Surcharge operator >> (vecteur vide au départ)
template <class T>
istream& operator >> (istream& s, vector<T> & vect){
    assert(vect.size()>0);
    cout << "Entrer au clavier un vecteur de taille " << vect.size() << endl;
    for (int i=0; i<vect.size(); i++)
        s >> vect[i];
    return s;
}

// Surcharge operator <<
template <class T>
ostream& operator << (ostream& s, const vector<T> & vect){
    assert(vect.size()>0);
    //cout << "Affiche le vecteur de taille " << vect.size() << endl;
    for (int i=0; i<vect.size(); i++)
        s << vect[i] << " ";
    s << endl;
    return s;
}

// Lire un vecteur depuis un fichier
template <class T>
void lire_vect(string nom_fichier, vector<T> & vect){
    ifstream file;
    T val;
    file.open(nom_fichier,ios::in); // ouverture du fichier
    assert(!file.fail()); // verifier ouverture correcte du fichier
    while (!file.eof()){ // on lit tant qu'on est pas à la fin du fichier
        file >> val;
        vect.push_back(val);
    }
    file.close();
    //return vect;
    }

// Ecrire un vecteur dans un fichier
template <class T>
void ecrire_vect(string nom_fichier, vector<T> & vect){
    ofstream file;
    file.open(nom_fichier,ios::out); // ouverture du fichier
    assert(!file.fail()); // verifier ouverture correcte du fichier
    for (int i=0; i<vect.size(); i++)
        file << " " << vect[i];
    file.close();
}

// Surcharge operateur + : somme de 2 vecteurs
template <class T>
vector<T> operator + (const vector<T>& v1, const vector<T>& v2){
    assert((v1.size()>0)&&(v1.size()==v2.size()));
    vector<T> w(v1.size(),0);
    for (int i=0; i<v1.size(); i++)
         w[i]= v1[i]+v2[i];
    return w;
}

// Surcharge operateur - : difference de 2 vecteurs
template <class T>
vector<T> operator - (const vector<T>& v1, const vector<T>& v2){
    assert((v1.size()>0)&&(v1.size()==v2.size()));
    vector<T> w(v1.size(),0);
    for (int i=0; i<v1.size(); i++)
         w[i]= v1[i]-v2[i];
    return w;
}

// Surcharge operateur * : produit scalaire
template <class T>
double operator * (const vector<T>& v1, const vector<T>& v2){
    assert((v1.size()>0)&&(v1.size()==v2.size()));
    double ps=0.0;
    for (int i=0; i<v1.size(); i++)
         ps += v1[i]*v2[i];
    return ps;
}

// Surcharge operateur * : produit vecteur x cte
template <class T>
vector<T> operator * (const vector<T>& v, const T& cte){
    assert(v.size()>0);
    vector<T> w(v.size(),0);
    for (int i=0; i<v.size(); i++)
         w[i]= v[i]*cte;
    return w;
}

// le produit scalaire de deux vecteurs
template <typename T>
T operator^ (const vector<T>& v1, const vector<T>& v2) {
    assert(v1.size() == v2.size());
    T res = 0;
    for (size_t i = 0; i < v1.size(); ++i)
        res += v1[i] * v2[i];
    return res;
}

// Norme_inf
template <class T>
T norme_inf(const vector<T> & vect){
    T nrm=-1;
    assert(vect.size()>0); // verifier si vecteur non vide
    for (int i=0; i<vect.size(); i++)
        if (abs(vect[i])>nrm) nrm = abs(vect[i]);
    return nrm;
}

// norme l-infinin de un vecteur cái này của mình 
template <typename T>
T norme_linf(const vector<T>& v){
    T max_val = fabs(v[0]);
    for (size_t i = 1; i < v.size(); ++i) {
        if (fabs(v[i]) > max_val) {
            max_val = fabs(v[i]);
        }
    }
    return max_val; 
}


// Norme_1
template <class T>
T norme_1(const vector<T> & vect){
    T nrm=0.0;
    assert(vect.size()>0);
    for (int i=0; i<vect.size(); i++)
        nrm += abs(vect[i]);
    return nrm;
}

// norme l1 de un vecteur
template <typename T>
T norme_l1(const vector<T>& v) {
    T sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += fabs(v[i]);
    }
    return sum;
}

// Norme_2
template <class T>
double norme_2(const vector<T> & vect){
    assert(vect.size()>0);
    return sqrt(vect*vect);
}


// norme euclidienne l2 de un vecteur
template <typename T>
T norme_eculidienne(const vector<T>& v) {
    T sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += fabs(v[i]) * fabs(v[i]);
    }
    return sqrt(sum);
}
// Norme_l2
template <class T>
double norme_l2(const vector<T> & vect){
    assert(vect.size()>0);
    return sqrt(vect*vect/vect.size());
}

// norme l2 d'un vecteur
template <typename T>
T norm_l2(const vector<T>& v) {
    T sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += (fabs(v[i]) * fabs(v[i])) * 1 / v.size();
    }
    return sqrt (sum);
}

#endif // !defined VECTEUR_TEMPLATE_H
