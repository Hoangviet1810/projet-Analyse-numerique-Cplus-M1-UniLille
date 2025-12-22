#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

# ifndef matrice_h
# define matrice_h

using namespace std;

class matrice
{
private:
    int size1, size2;
    vector<double>* matrix;           
public:
    matrice() {this->size1=0; this->size2=0; this->matrix=nullptr;};  // constructeur par défaut
    matrice(int, int);          // constructeur en donnant les dimensions
    matrice(const matrice &);   // constructeur par copie
    ~matrice();                 // destructeur

    int dim1() const { return this->size1; }
    int dim2() const { return this->size2; }
    vector<double>& operator[](int i);
    const vector<double>& operator[](int i) const;
    double& operator()(int, int );
    const double& operator()(int, int ) const;

    matrice& operator=(const matrice &);

    matrice transpose();
    vector<double> colonne(int) const;
    vector<double> ligne(int) const;

    void zero();
    void set_constant(double);

    void resize(int n, int m);

    matrice operator*(const matrice &) const;
    matrice operator+(const matrice &) const;
    matrice operator-(const matrice &) const;

    matrice operator*(const double &);
    vector<double> operator*(const vector<double> &) const;
    friend vector<double> operator*(const vector<double> &, const matrice &);

    friend ostream &operator<<(ostream &o, const matrice &m);
    friend istream &operator>>(istream &i, matrice &m);

    friend void save_matrice(const string &filename, const matrice &);
    friend matrice read_matrice(const char* filename);
    
};
matrice read_matrice(const char* filename);
void save_matrice(const string &filename, const matrice &m);
#endif