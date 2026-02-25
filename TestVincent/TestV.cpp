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


    xAxis=4;
    yAxis=4;
    zAxis=4;
    tAxis=4;
    std::vector<Matrix<rSU,cSU>> lattice(linksPerSite*(xAxis*yAxis*zAxis*tAxis));
    //printArray<rSU,cSU>(identityMatrix);

    std::mt19937_64 mega(14703);
    std::mt19937_64 mega2(1613);
    std::mt19937_64 mega3(3479);
    std::uniform_real_distribution<double> distro(0,1);

    std::uniform_int_distribution<int> intro(1,3);

    std::cout << "beta: ";
    std::cin >> beta;
    std::cout << "Epsilon: ";
    std::cin >> epsilon;
    hotEpsilon /= beta;
    //hotEpsilon*=3.5;
    if(hotEpsilon > 0.7){
        hotEpsilon =0.7;
    }
    hot_start(lattice);


std::atomic<int> acceptanceRate=0;
    std::atomic<int> updates=0;

    size_t observableCollected= 0;
    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);

    size_t SweepFactor = 2;



    //thermalization
    for(int p=0; p<400/5; p++){
        //parallelization
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*17007+ (p+2)* 10111);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);

            std::uniform_real_distribution<double> threadAcceptReject(0,1);
        
            std::uniform_int_distribution<int> threadrefelctDist(1,3);


            Matrix<rSU, cSU> UPrime, X, A, U; 

            U= lattice[i];

            std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(i);
            size_t mu,x,y,z,t;
            mu= std::get<0>(temp);
            x= std::get<1>(temp);
            y= std::get<2>(temp);
            z= std::get<3>(temp);
            t= std::get<4>(temp);
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<1; j++){


                //in theory automatically accepted
                if(j%2==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    //std::cout <<  (matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real() << std::endl;
                  //Matrix<rSU,cSU> sub = matrix_subtraction(U, UPrime);
                    //std::cout << "sub: "<< sub(0,0) <<std::endl;
                    

                    //Matrix<rSU,cSU> mult = matrix_multiplication(sub,A);
                    //double trace = matrix_trace(mult).real();     
                    //bool acceptance = latticeAction(lattice, lattice[i], U, A, Threadindexing, threadAcceptReject);
                    normalizeSU3Matrix(U);



                   // if(acceptance==true){
                   //     std::cout << "accepted" << std::endl; 
                   //     }
                   // else{
                   //     std::cout << "false" << std::endl;
                   // }

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,lattice[i]);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, lattice[i], UPrime, A, Threadindexing, threadAcceptReject);


                    if(acceptance==true){
                        U=UPrime;
                        acceptanceRate.fetch_add(1);
                        }
                    updates.fetch_add(1);
                    }
//

                }
            bufferLattice[i]=U;

        });


        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);
            //Update X matrices
        if(p%2 ==0 && p!=0){
            double rate = double(acceptanceRate.load())/double(updates.load());
            epsilon *= 1 + alpha * (rate - target_rate);
            if(epsilon>0.7){
                epsilon=1;
            }
            
            acceptanceRate.store(0);
            updates.store(0);
            X_updateSU3(p*400+ 270);
                        }
        

    }

    epsilonTune(lattice,0,2,1,3);
    ThermalTune(lattice, 0,2,1,3);
    SweepFactor=AutoCorrelationTune(lattice,0,2,1,3);

    std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
    plaquette(lattice, plaquettes);
    //for(int i =0; i< plaquettes.size(); i++){
    //    std::cout << plaquettes[i] << std::endl;
    //}
    double PTest=average(plaquettes);



    std::cout << PTest << std::endl;

    std::cout << SweepFactor << std::endl;





    return 0;
}















