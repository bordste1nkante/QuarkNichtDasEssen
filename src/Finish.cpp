#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <execution>
#include <thread>
#include <complex>
#include <random>
#include <Eigen/Dense>
#include "global.h"
#include "latticeOP.h"
#include "matrixOP.h"
#include "utils.h"

void Simulation(   std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable,
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep ){


    size_t observableCollected= 0;
    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);



    //thermalization
    for(int p=0; p<numberOfThermalSweeps/numberOfMultiHit; p++){
        //parallelization
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            std::mt19937_64 Threadindexing(dindexing*i+p);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);

            std::mt19937_64 Threadaccept(dacceptReject*i+p);
            std::uniform_real_distribution<double> threadAcceptReject(0,1);
        
            std::mt19937_64 Threadreflect(dreflection*i+p);
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

            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadreflect);
                    normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,lattice[i]);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, lattice[i], UPrime, A, Threadaccept, threadAcceptReject);

                    if(acceptance==true){
                        U=UPrime;
                        }

                    }


                }
            bufferLattice[i]=U;

        });

        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);
            //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3(p+numberOfThermalSweeps);
                        }

    }

    size_t k = 0;
    for(int p=0; p<NConfigs*SweepFactor/numberOfMultiHit; p++){

        //parallelization
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            std::mt19937_64 Threadindexing(dindexing*i+pow(p,1.5));
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);

            std::mt19937_64 Threadaccept(dacceptReject*i+pow(p,1.5));
            std::uniform_real_distribution<double> threadAcceptReject(0,1);

            std::mt19937_64 Threadreflect(dreflection*i+pow(p,1.5));
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

            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadreflect);
                    normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,lattice[i]);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, lattice[i], UPrime, A, Threadaccept ,threadAcceptReject);
                    if(acceptance==true){
                        U=UPrime;
                        }

                    }


                }
            bufferLattice[i]=U;

        });

        k+=numberOfMultiHit;

        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);
        
        //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3(p*p+SweepFactor);
                        }
        if(k % SweepFactor==0) {

            std::vector<std::complex<double>> loops;
            std::vector<double> r;
            if(observable == 1){
                wilsonLoop( lattice, loops, r, startingPoint,endPoint );
                std::string pathR = "/Configuration/"+ std::to_string(observableCollected)+ "/r";
                std::string pathL = "/Configuration/"+ std::to_string(observableCollected)+ "/loops";
                saveArrayH5(r,pathR);
                saveArrayH5complex(loops,pathL);
            }
            if(observable == 2){
                polyakovLoop( lattice, loops, r, startingPoint,endPoint );
                std::string pathR = "/Configuration/"+ std::to_string(observableCollected)+ "/r";
                std::string pathL = "/Configuration/"+ std::to_string(observableCollected)+ "/loops";
                saveArrayH5(r,pathR);
                saveArrayH5complex(loops,pathL);
            }
            observableCollected +=1;
            std::cout << "Collected:" <<observableCollected <<"/" <<NConfigs << std::endl;


            
            }
    }
    }