#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/matrixOP.h"
#include "../header/global.h"

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

Matrix<rSU,cSU> matrix_conjugate(const Matrix<rSU,rSU>& A){
    Matrix<rSU,cSU> C;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){

            C(i,j)= std::conj(A(j,i));

            }

        }
    return C;
}

std::complex<double> matrix_trace(const Matrix<rSU,rSU>& A){
    std::complex<double> trace;
    for(int i=0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i=j){
                trace += A(i,j);
            }

            }

        }
    return trace;
}