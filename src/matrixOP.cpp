#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/matrixOP.h"
#include "../header/global.h"


void generate_Pauli(){
    for(int i=0; i<2; i++){
        for(int j = 0; j<2; j++){
            if(i==j){
            identity2(i,j)=std::complex<double> {1.0,0.0};
            Pauli1(i,j)=std::complex<double> {0.0,0.0};
            Pauli2(i,j)=std::complex<double> {0.0,0.0};
            }
            else{
                identity2(i,j)=std::complex<double> {0.0,0.0};
                Pauli1(i,j)=std::complex<double> {1.0,0.0};
                Pauli3(i,j)=std::complex<double> {0.0,0.0};

            }
        }
    }

    Pauli2(0,1)=std::complex<double> {0.0,-1.0};
    Pauli2(1,0)=std::complex<double> {0.0,1.0};
    Pauli3(0,0)=std::complex<double> {1.0,0.0};
    Pauli3(1,1)=std::complex<double> {-1.0,0.0};

    pauliMatrices.push_back(identity2);
    pauliMatrices.push_back(Pauli1);
    pauliMatrices.push_back(Pauli2);
    pauliMatrices.push_back(Pauli3);
    
}

Matrix<rSU,cSU> matrix_multiplication(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            std::complex<double> sum;
            for(int k = 0; k<rSU; k++){
                sum += A(i,k)*B(k,j);

                }
                C(i,j)= sum;
            }

        }
    return C;
}


Matrix<rSU,cSU> matrix_addition(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){

            C(i,j)= A(i,j)+B(i,j);

            }

        }
    return C;
}

Matrix<rSU,cSU> matrix_subtraction(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){

            C(i,j)= A(i,j)-B(i,j);

            }

        }
    return C;
}

Matrix<rSU,cSU> matrix_hermitean_conjugate(const Matrix<rSU,rSU>& A){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){

            C(i,j)= std::conj(A(j,i));

            }

        }
    return C;
}

std::complex<double> matrix_trace(const Matrix<rSU,rSU>& A){
    std::complex<double> trace = {0.0,0.0}; 
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i==j){
                trace += A(i,j);
            }

            }

        }
    return trace;
}



//generates a rSU x cSU identity matrix
void generate_identity(){
    for(int i = 0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i==j){
                std::complex<double> z = {1.0, 0.0};
                identityMatrix(i,j)= z;
            }
            else{

                std::complex<double> z1 = {0.0, 0.0};
                identityMatrix(i,j)=z1;
            }
        }
    }
}