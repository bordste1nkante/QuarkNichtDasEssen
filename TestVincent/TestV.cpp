#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"
#include "../header/observables.h"


template <size_t R, size_t C>
void printArray(Matrix<rSU,cSU> A) {
    for (size_t i = 0; i < rSU; ++i) {
        for (size_t j = 0; j < cSU; ++j) {
            std::cout << A(i,j)<<'\n';
        }
        
    }
    std::cout << "done" << std::endl;
}

template <size_t R, size_t C>
void printArray2(Matrix<2,2> A) {
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            std::cout << A(i,j)<<'\n';
        }
        
    }
    std::cout << "done" << std::endl;
}


int main(){


    Matrix <rSU,cSU> A, B, C;

    Matrix <rSU,cSU> A2, B2, C2;

    Matrix <rSU,cSU> hA, hB, hC;

    Matrix <rSU,cSU> ApB, BpC;

    Matrix <rSU,cSU> BmC;

    Matrix <rSU,cSU> AB,BA, BC, CC;



    for(int i= 0; i<rSU; i++){
        for(int j=0; j<rSU; j++){

            if(j==i){
                A(i,j)= std::complex<double>{1.0,0};

            }
            else{
                A(i,j)= std::complex<double>{0,0};
                B(i,j)= std::complex<double>{0,0};
                C(i,j)= std::complex<double>{0,0};
            }
        }
    }
    B(2,2)= std::complex<double>{1.0,0};
    B(0,0)= std::complex<double>{-0.5,std::sqrt(3)/2};
    B(1,1)= std::complex<double>{-0.5,-std::sqrt(3)/2};

    C(2,2)= std::complex<double>{1.0,0};
    C(1,1)= std::complex<double>{-0.5,std::sqrt(3)/2};
    C(0,0)= std::complex<double>{-0.5,-std::sqrt(3)/2};


    generate_identity();
    generate_Pauli();
    //std::complex<double> detA,detB,detC;
    //std::complex<double> traceA,traceB,traceC;

    //detA = det_A(A);
    //detB=det_A(B);
    //detC= det_A(C);
//
    ////std::cout << detA << std::endl;
//
    ////std::cout << detB << std::endl;
//
    ////std::cout << detC << std::endl;
//
    //traceA = matrix_trace(A);
    //traceB=matrix_trace(B);
    //traceC= matrix_trace(C);
//
//
    ////std::cout << traceA << std::endl;
    ////std::cout << traceB << std::endl;
    ////std::cout << traceC << std::endl;
//
    //ApB = matrix_addition(A,B);
    //BpC = matrix_addition(B,C);
//
    //printArray<rSU,cSU>(ApB);
//
    ////std::cout << "done" << std::endl;
//
    //printArray<rSU,cSU>(BpC);
    ////std::cout << "done" << std::endl;
//
    //BmC = matrix_subtraction(B,C);
    ////std::cout << "done" << std::endl;
//
    //printArray<rSU,cSU>(BmC);
//
    //AB = matrix_multiplication(A,B);
    //BA = matrix_multiplication(B,A);
    //BC= matrix_multiplication(B,C);
    //CC= matrix_multiplication(C,C);
//
//
    //printArray<rSU,cSU>(AB);
    ////std::cout << "done" << std::endl;
    //printArray<rSU,cSU>(BA);
    ////std::cout << "done" << std::endl;
    //printArray<rSU,cSU>(BC);
    ////std::cout << "done" << std::endl;
    //printArray<rSU,cSU>(CC);
    ////std::cout << "done" << std::endl;
//
    
    //hA= matrix_hermitean_conjugate(A);
    //hB = matrix_hermitean_conjugate(B);
    //hC = matrix_hermitean_conjugate(C);
//
    //printArray<rSU,cSU>(hA);
    //printArray<rSU,cSU>(A);
    //printArray<rSU,cSU>(hB);
    //printArray<rSU,cSU>(B);
    //printArray<rSU,cSU>(hC);
    //printArray<rSU,cSU>(C);

    //A2= inverse_A(A);
    //B2= inverse_A(B);
    //C2= inverse_A(C);
//
    //printArray<rSU,cSU>(A2);
    //printArray<rSU,cSU>(A);
    //printArray<rSU,cSU>(B2);
    //printArray<rSU,cSU>(B);
    //printArray<rSU,cSU>(C2);
    //printArray<rSU,cSU>(C);

    //printArray2<2,2>(pauliMatrices[0]);
    //printArray2<2,2>(pauliMatrices[1]);
    //printArray2<2,2>(pauliMatrices[2]);
    //printArray2<2,2>(pauliMatrices[3]);
//
    //printArray<rSU,cSU>(identityMatrix);



    return 0;
}