#include "matrice.h"
#include <vector>
//#include "vecteur_template.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cmath>
using namespace std;

// constructeur en donnant les dimensions
matrice::matrice(int n, int m) {
    assert(n > 0 && m > 0);
    this->size1 = n; 
    this->size2 = m;
    this->matrix = new vector<double>[n];
    vector<double> v(m,0.0);  
    for (int i = 0; i < n; i++) {
        this->matrix[i] = v;  
    }
}  

// constructeur par copie
matrice::matrice(matrice const & mat) {   
    assert((mat.size1 > 0) && (mat.size2 > 0));
    this->size1 = mat.size1;
    this->size2 = mat.size2;
    this->matrix = new vector<double>[size1];
    for (int i = 0; i < size1; i++) {
        this->matrix[i] = mat.matrix[i];
    }
}

// destructeur
matrice::~matrice(){
    if (this->size1||this->size2) delete[] this->matrix;
    this->size1 = 0;
    this->size2 = 0;
}

// Retour la ligne i
vector<double>& matrice::operator[](int i){
    assert(i >= 0 && i < this->size1);
    return this->matrix[i];
}


const vector<double>& matrice::operator[](int i) const {
    assert(i >= 0 && i < this->size1);
    return matrix[i];
}

double& matrice::operator()(int i, int j) {
    assert(i >= 0 && i < size1 && j >= 0 && j < size2);
    return matrix[i][j]; 
}


const double& matrice::operator()(int i, int j) const {
    assert(i >= 0 && i < size1 && j >= 0 && j < size2);
    return matrix[i][j];  // trả về giá trị (copy) tại (i,j)
}

matrice& matrice::operator=(const matrice  &m){
        if (this != &m){
            delete[] matrix;
            size1 = m.dim1();
            size2 = m.dim2();
            matrix = new vector<double>[size1];
            for (int i = 0; i < size1; i++)
            {
                matrix[i] = vector<double>(size2, 0.0);
                for (int j = 0; j < size2; j++)
                {
                    matrix[i][j] = m[i][j];
                }
            }
        }
        return *this;
}


matrice matrice::transpose(){
    matrice T(this->size2, this->size1);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            T[j][i] = this->matrix[i][j]; 
        }    
    }
    return T;        
}


vector<double> matrice::colonne(int j) const{ 
    assert(j >= 0 && j < this->size2);
    vector<double> col(this->size1);  
    for (int i = 0; i < this->size1; i++){
        col[i] = this->matrix[i][j];
    }
    return col;
}

vector<double> matrice::ligne(int i) const{
    assert(i >= 0 && i < this->size1);
    vector<double> row(this->size2);   
    for (int j = 0; j < this->size2; j++){
        row[j] = this->matrix[i][j];
    }
    return row;
}

void matrice::zero(){
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            this->matrix[i][j] = 0.0;
        }
    }
}
void matrice::set_constant(double d){
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            this->matrix[i][j] = d;
        }
    }
}

matrice matrice::operator*(const matrice &m ) const{
    assert(this->size2 == m.size1);
    matrice result(this->size1, m.size2);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < m.size2; j++){
            for (int k = 0; k < this->size2; k++){
                result[i][j] += this->matrix[i][k] * m.matrix[k][j];
            }               
        }    
    }
    return result;
}

matrice matrice::operator+(const matrice &m ) const{
    assert(this->size1 == m.size1 && this->size2 == m.size2);
    matrice result(this->size1, this->size2);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            result[i][j] = this->matrix[i][j] + m.matrix[i][j];
        }    
    }
    return result;
}

matrice matrice::operator-(const matrice &m ) const{
    assert(this->size1 == m.size1 && this->size2 == m.size2);
    matrice result(this->size1, this->size2);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            result[i][j] = this->matrix[i][j] - m.matrix[i][j];
        }    
    }
    return result;
}

matrice matrice::operator*(const double &valeur){
    matrice result(this->size1, this ->size2);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){ 
            result[i][j] = this->matrix[i][j] * valeur;
        }
    }
    return result;      
}

vector<double> matrice::operator*(const vector<double> &v) const{
    assert(this->size2 == v.size()); 
    vector<double> result(this->size1, 0.0);
    for (int i = 0; i < this->size1; i++){
        for (int j = 0; j < this->size2; j++){
            result[i] += this->matrix[i][j] * v[j];
        }
    }
    return result;
}

vector<double> operator*(const vector<double> &v, const matrice &mat){
    assert (v.size() == mat.dim1());
    vector<double> result(mat.dim2(), 0.0);
    for (int j = 0; j < mat.dim2(); j++){
        for (int i = 0; i < mat.dim1(); i++){
            result[j] += v[i] * mat.matrix[i][j];   
        }
    }
    return result;
}

ostream &operator<<(ostream &o, matrice const &mat){
    assert (mat.dim1() > 0 && mat.dim2() > 0);    
    for (int i = 0; i < mat.dim1(); i++){
        for (int j = 0; j < mat.dim2(); j++){
            o << mat[i][j] << " ";
        }
        o << endl;
    }
    return o;
}

istream &operator>>(istream &i, matrice &mat){
    assert(mat.dim1() > 0 && mat.dim2() > 0);
    for (int r = 0; r < mat.dim1(); r++){
        for (int j = 0; j < mat.dim2(); j++){
            i >> mat[r][j] ;
        }
    }
    return i;
}

void save_matrice(const string &filename, const matrice &mat){
    ofstream file;   
    file.open(filename, ios::out);
    assert(!file.fail());
    
    file << mat.dim1() << " " << mat.dim2() << endl;

    for (int i = 0; i < mat.dim1(); i++){
        for (int j = 0; j < mat.dim2(); j++){
            file << mat[i][j] << " ";  
        }
        file << endl;                   
    }

    file.close();   
}

matrice read_matrice(const char* filename){
    ifstream file(filename);   
    assert(file.is_open());

    int rows, cols;
    file >> rows >> cols;     

    matrice mat(rows, cols);  

   
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            file >> mat[i][j];   
        }
    }

    file.close();   
    return mat;    
}

void matrice::resize(int n, int m){
    if (matrix != nullptr)
        delete[] matrix;

    size1 = n;
    size2 = m;

    matrix = new vector<double>[n];

    for (int i = 0; i < n; ++i)
        matrix[i].assign(m, 0.0);   
}

