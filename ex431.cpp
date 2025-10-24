#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

// Tridiagonal matrix stored as band
struct mat_bande{
    int n;
    vector<int> IND;              // offsets [-1,0,1]
    vector<vector<double>> TAB;   // n x 3
};

// Factorisation LU for tridiagonal
void factorisation_LU_tridiag(mat_bande& A){
    int n = A.n;
    if(n < 2) return;
    A.TAB[1][0] /= A.TAB[0][1];
    A.TAB[1][1] -= A.TAB[1][0]*A.TAB[0][2];
    for(int i=2;i<n;i++){
        A.TAB[i][0] /= A.TAB[i-1][1];
        A.TAB[i][1] -= A.TAB[i][0]*A.TAB[i-1][2];
    }
}

// Solve LU tridiagonal
vector<double> resoudre_LU_tridiag(const mat_bande& A, const vector<double>& b){
    int n = A.n;
    vector<double> x = b;
    // forward substitution
    for(int i=1;i<n;i++) x[i] -= A.TAB[i][0]*x[i-1];
    // backward substitution
    x[n-1] /= A.TAB[n-1][1];
    for(int i=n-2;i>=0;i--) x[i] = (x[i] - A.TAB[i][2]*x[i+1])/A.TAB[i][1];
    return x;
}

// Generate X vector
vector<double> generer_X(double L, int N){
    vector<double> X(N+1);
    double h = L/double(N);
    for(int i=0;i<=N;i++) X[i] = i*h;
    return X;
}

// Strike function
double K_of_x(double x){ return 0.95*x; }

// Payoff
vector<double> u0_call(const vector<double>& X){
    vector<double> u0(X.size());
    for(size_t i=0;i<X.size();i++) u0[i] = max(X[i]-K_of_x(X[i]),0.0);
    return u0;
}

// Build matrix for implicit BS
mat_bande matrice_BS_bande(const vector<double>& X, double h, double k, double sigma){
    int M = X.size();
    mat_bande A;
    A.n = M;
    A.IND = {-1,0,1};
    A.TAB.assign(M, vector<double>(3,0.0));

    for(int i=0;i<M;i++){
        double nu = 0.5*sigma*sigma*X[i]*X[i];
        double factor = k*nu/(h*h);
        if(i==0){
            // Dirichlet u=0
            A.TAB[i][0] = 0.0; A.TAB[i][1] = 1.0; A.TAB[i][2] = 0.0;
        } else if(i==M-1){
            // Neumann du/dx=0
            A.TAB[i][0] = -factor;
            A.TAB[i][1] = 1.0 + factor;
            A.TAB[i][2] = 0.0;
        } else {
            A.TAB[i][0] = -factor;
            A.TAB[i][1] = 1.0 + 2.0*factor;
            A.TAB[i][2] = -factor;
        }
    }
    return A;
}

// Implicit scheme
vector<double> BS_implicite(const vector<double>& X, double L, double T, double sigma, int N, double h, double k, const vector<double>& u0, const string& file_out){
    int M = X.size();
    vector<double> U = u0;

    mat_bande Mmat = matrice_BS_bande(X,h,k,sigma);
    mat_bande Mcopy = Mmat;
    factorisation_LU_tridiag(Mcopy);

    int Ktime = int(T/k);
    ofstream fout(file_out);
    fout << "# t ";
    for(int i=0;i<M;i++) fout << " x_"<<i;
    fout<<"\n";

    for(int n=0;n<=Ktime;n++){
        double t = n*k;
        fout<<t;
        for(int i=0;i<M;i++) fout<<" "<<U[i];
        fout<<"\n";

        if(n<Ktime){
            vector<double> B = U;
            B[0] = 0.0; // Dirichlet
            U = resoudre_LU_tridiag(Mcopy,B);
        }
    }

    fout.close();
    return U;
}

int main(){
    double L = 10000.0;
    double T = 365.0;
    double sigma = 0.2;
    int N = 200;
    double h = L/double(N);
    double k = 1.0; // step 1 day

    vector<double> X = generer_X(L,N);
    vector<double> u0 = u0_call(X);

    string file_out = "BS_implicit.txt";
    vector<double> Ufinal = BS_implicite(X,L,T,sigma,N,h,k,u0,file_out);

    cout<<"Computed solution at final time T = "<<T<<" (written in "<<file_out<<")\n";
    cout<<"x   U(T,x)\n";
    for(size_t i=0;i<X.size();i++) cout<<X[i]<<"   "<<Ufinal[i]<<"\n";

    return 0;
}
