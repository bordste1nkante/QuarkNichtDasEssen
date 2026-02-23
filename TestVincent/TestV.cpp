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
    X_updateSU3(7);


    xAxis = 3;
    yAxis= 3;
    zAxis=3;
    tAxis=3;
    std::vector<Matrix<rSU,cSU>> lattice(4*(xAxis*yAxis*zAxis*tAxis));
    //printArray<rSU,cSU>(identityMatrix);

    std::mt19937_64 mega(12003);
    std::mt19937_64 mega2(1013);
    std::mt19937_64 mega3(5479);
    std::uniform_real_distribution<double> distro(0,1);

    std::uniform_int_distribution<int> intro(1,3);

    std::cout << "beta: ";
    std::cin >> beta;
    std::cout << "Epsilon: ";
    std::cin >> epsilon;

    hot_start(lattice);

    epsilonTune(lattice,0,2,1,3);
    ThermalTune(lattice, 0,2,1,3);
    std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
    plaquette(lattice, plaquettes);
    double PTest=average(plaquettes);



    std::cout << PTest << std::endl;





    return 0;
}
//














