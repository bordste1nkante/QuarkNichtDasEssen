#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
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


    Matrix <rSU,cSU> A, B, C, X, XInv, I, hX;





    //for(int i= 0; i<rSU; i++){
    //    for(int j=0; j<rSU; j++){
//
    //        if(j==i){
    //            A(i,j)= std::complex<double>{1.0,0};
//
    //        }
    //        else{
    //            A(i,j)= std::complex<double>{0,0};
    //            B(i,j)= std::complex<double>{0,0};
    //            C(i,j)= std::complex<double>{0,0};
    //        }
    //    }
    //}
    //B(2,2)= std::complex<double>{1.0,0};
    //B(0,0)= std::complex<double>{-0.5,std::sqrt(3)/2};
    //B(1,1)= std::complex<double>{-0.5,-std::sqrt(3)/2};
//
    //C(2,2)= std::complex<double>{1.0,0};
    //C(1,1)= std::complex<double>{-0.5,std::sqrt(3)/2};
    //C(0,0)= std::complex<double>{-0.5,-std::sqrt(3)/2};


    generate_identity();
    generate_Pauli();
    X_updateSU3();


    xAxis = 1;
    yAxis= 1;
    zAxis=1;
    tAxis=1;
    std::vector<Matrix<rSU,cSU>> lattice(4*(xAxis*yAxis*zAxis*tAxis));
    //printArray<rSU,cSU>(identityMatrix);

    hot_start(lattice);
    //for(int k = 0; k<100; k++){
    //for(int mu= 0; mu<1; mu++){
    //    //printArray<rSU,cSU>(lattice[idx(0,0,0,0,mu)]);
//
//
    //    std::complex<double> det = det_A(lattice[idx(0,0,0,0,mu)]);
    //    std::complex<double> trace = matrix_trace(matrix_multiplication( matrix_hermitean_conjugate(lattice[idx(0,0,0,0,mu)]),lattice[idx(0,0,0,0,mu)]));
//
    //    std::cout << det<< std::endl;
    //    std::cout << trace << std::endl;
    //    std::cout << "done" << std::endl;
    //    X = XSet[indexDist(indexing)];
    //    
    //    lattice[idx(0,0,0,0,mu)]= matrix_multiplication(X,lattice[idx(0,0,0,0,mu)]);
    //    normalizeSU3Matrix(lattice[idx(0,0,0,0,mu)]);
    //}
    //}

    //printArray<rSU, cSU>(XSet[2*1+1]);

    //for(int i=0; i<NSetXMatrices; i++){
//
    //    X=XSet[2*i];
    //    XInv =XSet[2*i+1];
//
    //    I= matrix_multiplication(XInv,matrix_hermitean_conjugate(XInv));
//
//
//
    //    std::complex<double> detI, detX, detXInv;
//
//
    //    detI = det_A(I);
    //    detXInv = det_A(XInv);
    //    detX = det_A(X);
//
    //    std::cout << detI << std::endl;
    //    std::cout << detX << std::endl;
    //    std::cout << detXInv << std::endl;
    //    std::cout << "done" << std::endl;
//
//
    //}



    try{
    H5::H5File file("../testV.h5",H5F_ACC_TRUNC);

    std::vector<double> test = {0.1,0.2,0.3,0.4,0.5};

    hsize_t dimsA[1]= {5};
    H5::DataSpace spaceA(1,dimsA);

    H5::Group group =file.createGroup("/lattice");

    H5::DataSet dSetA = file.createDataSet(
        "/lattice/r",
        H5::PredType::NATIVE_DOUBLE,
        spaceA
    );

    dSetA.write(test.data(), H5::PredType::NATIVE_DOUBLE);
    file.close();


    
}
catch(H5::Exception& e){
    e.printErrorStack();
}


    return 0;
}