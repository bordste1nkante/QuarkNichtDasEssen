#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <Eigen/Dense> //requires installation of this specific library on Linux:
#include <H5Cpp.h>  //requires installation of this specific library on Linux: sudo apt install hdf5-tools libhdf5-dev
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"
#include "../header/observables.h"
#include "../header/utils.h"


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


    Matrix <rSU,cSU> A, B, C, X, XInv, I, hX, U, UPrime;


    std::complex<double> det, trace;
    




    generate_identity();
    generate_Pauli();
    X_updateSU3();


    xAxis = 3;
    yAxis= 3;
    zAxis=3;
    tAxis=3;
    std::vector<Matrix<rSU,cSU>> lattice(4*(xAxis*yAxis*zAxis*tAxis));
    //printArray<rSU,cSU>(identityMatrix);

    std::mt19937_64 mega(123);
    std::mt19937_64 mega2(123);
    std::mt19937_64 mega3(111);
    std::uniform_real_distribution<double> distro(0,1);

    std::uniform_int_distribution<int> intro(1,3);


    int k = idx(2,2,0,0,1);


    hot_start(lattice);

    A = determineA(lattice,2,2,0,0,1 );

    U = lattice[k];


    std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(k);
    size_t mu,x,y,z,t;
    mu= std::get<0>(temp);
    x= std::get<1>(temp);
    y= std::get<2>(temp);
    z= std::get<3>(temp);
    t= std::get<4>(temp);


    std::cout << x<< std::endl;
    std::cout << y<< std::endl;
    std::cout << z<< std::endl;
    std::cout << t<< std::endl;
    std::cout << mu<< std::endl;


    //X= XSet[21];
    UPrime=overrelaxation(A,U,intro, mega3);


    //UPrime = matrix_multiplication(X,U);


    bool accept;
    double r;
    double probability;
    
    double SActionDif;
    SActionDif = -beta/(rSU)*(matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real();

   // std::cout << SActionDif << std::endl;
    probability = std::min(1.0, exp(- SActionDif)); // according to my notes -> check in doubt

    //std::cout << probability << std::endl;

    r= distro(mega2);
    if(r<=probability){
        accept = true;
    }
    else{
        accept = false;
    }




    bool accept1 = latticeAction(lattice,U,UPrime,A,mega,distro);

    //std::cout << accept << std::endl;
    //std::cout << accept1 << std::endl;
//














    return 0;
}