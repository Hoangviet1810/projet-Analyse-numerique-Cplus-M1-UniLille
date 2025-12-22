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
    vector<double> v(m,0.0);  // tạo vector tạm với kích thước m và giá trị 0.0
    for (int i = 0; i < n; i++) {
        this->matrix[i] = v;  // gán vector tạm cho mỗi hàng của ma trận
    }
}  

// constructeur par copie Lille
matrice::matrice(matrice const & mat) {   // đây là copy trực tiếp từ vecteur 
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

// matrice::~matrice()
// {
//     cout << "Destruct matrice " << this << " matrix=" << matrix << endl;
//     delete[] matrix;   // delete nullptr là an toàn
//     matrix = nullptr;
//     size1 = 0;
//     size2 = 0;
// }

// Retour la ligne i
vector<double>& matrice::operator[](int i){
    assert(i >= 0 && i < this->size1);
    return this->matrix[i];
}

// Phiên bản chỉ đọc
const vector<double>& matrice::operator[](int i) const {
    assert(i >= 0 && i < this->size1);
    return matrix[i];  // trả về vector dòng const
}

// Phiên bản đọc/ghi
double& matrice::operator()(int i, int j) {
    assert(i >= 0 && i < size1 && j >= 0 && j < size2);
    return matrix[i][j];  // trả về phần tử tại (i,j)
}

// Phiên bản chỉ đọc
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


vector<double> matrice::colonne(int j) const{ // ở đây int j là chỉ số cột
    assert(j >= 0 && j < this->size2);
    vector<double> col(this->size1);   // lấy size1 do nó phù hơp với số dòng
    for (int i = 0; i < this->size1; i++){
        col[i] = this->matrix[i][j];
    }
    return col;
}

vector<double> matrice::ligne(int i) const{
    assert(i >= 0 && i < this->size1);
    vector<double> row(this->size2);   // lấy size2 do nó phù hơp với số cột
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
    ofstream file;   // mở file để ghi
    file.open(filename, ios::out);
    assert(!file.fail());
    // ghi kích thước ma trận trước, ví dụ: "2 3"
    file << mat.dim1() << " " << mat.dim2() << endl;

    // ghi từng phần tử ma trận
    for (int i = 0; i < mat.dim1(); i++){
        for (int j = 0; j < mat.dim2(); j++){
            file << mat[i][j] << " ";   // ghi các phần tử của hàng
        }
        file << endl;                   // xuống dòng sau mỗi hàng
    }

    file.close();   // đóng file sau khi ghi xong
}

matrice read_matrice(const char* filename){
    ifstream file(filename);   // mở file để đọc
    assert(file.is_open());

    int rows, cols;
    file >> rows >> cols;     // đọc kích thước ma trận

    matrice mat(rows, cols);  // tạo ma trận với kích thước đã đọc

    // đọc từng phần tử ma trận
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            file >> mat[i][j];   // đọc các phần tử của ma trận
        }
    }

    file.close();   // đóng file sau khi đọc xong
    return mat;     // trả về ma trận đã đọc
}

// matrice read_matrice(const char* filename){
//     matrice A;
//     ifstream fh(filename);
//     fh >> A;
//     return A;
// }

void matrice::resize(int n, int m)
{
    // giải phóng đúng cách
    if (matrix != nullptr)
        delete[] matrix;

    size1 = n;
    size2 = m;

    matrix = new vector<double>[n];

    for (int i = 0; i < n; ++i)
        matrix[i].assign(m, 0.0);   
}

// matrice& matrice::operator=(const matrice& other)
// {
//     if (this == &other)
//         return *this;

//     delete[] matrix;

//     size1 = other.size1;
//     size2 = other.size2;

//     if (other.matrix)
//     {
//         matrix = new std::vector<double>[size1];
//         for (int i = 0; i < size1; ++i)
//             matrix[i] = other.matrix[i];
//     }
//     else
//         matrix = nullptr;

//     return *this;
// }

// matrice::matrice(const matrice& other)
//     : size1(other.size1), size2(other.size2)
// {
//     if (other.matrix)
//     {
//         matrix = new std::vector<double>[size1];
//         for (int i = 0; i < size1; ++i)
//             matrix[i] = other.matrix[i];
//     }
//     else
//         matrix = nullptr;
// }
