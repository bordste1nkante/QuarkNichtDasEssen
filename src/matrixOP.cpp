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


            }
        }
    }
    // this should be 0 & 1 as indices?
            Pauli2(1,2)=std::complex<double> {0.0,-1.0};
            Pauli2(2,1)=std::complex<double> {0.0,1.0};
            Pauli3(1,1)=std::complex<double> {1.0,0.0};
            Pauli3(2,2)=std::complex<double> {-1.0,0.0};
            // was ist mit den off diagonal elementen für sigma_3? -> diagonalelemente für sigma_2?
}

Matrix<rSU,cSU> matrix_multiplication(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            std::complex<double> sum;// ich check die funktion von diesem sum nicht
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

Matrix<rSU,cSU> matrix_conjugate(const Matrix<rSU,rSU>& A){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){

            C(i,j)= std::conj(A(j,i));//transponiert nicht?

            }

        }
    return C;
}

std::complex<double> matrix_trace(const Matrix<rSU,rSU>& A){
    std::complex<double> trace; //{0.0,0.0}; //-> muss ich dem Lachs nicht erst nen Wert zuweisen?
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i==j){//don't assign i to j
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