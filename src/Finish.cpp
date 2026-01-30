#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
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

    //thermalization
    for(int p=0; p<numberOfThermalSweeps/numberOfMultiHit; p++){
    for(size_t i = 0; i<lattice.size(); i++ ){

        //this all belongs to the multi-thread part
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


            //in theory atomatically accepted
            if(j%overrelaxationStep && j!=0){
                U = overrelaxation(A, U);
                normalizeSU3Matrix(U);

            }
            else{
                size_t index = indexDist(indexing);
                X=XSet[index];

                UPrime = matrix_multiplication(X,lattice[idx(x,y,z,t, mu)]);
                normalizeSU3Matrix(UPrime);
                bool acceptance = latticeAction(lattice, lattice[idx(x,y,z,t, mu)], UPrime,x,y,z,t,mu, A);

                if(acceptance==true){
                    U=UPrime;
                }

                }
        }

        //multi thread stops

        //ensure that all U are inserted at their location


    }
            //Update X matrices
            if(p%XUpdate ==0 && p!=0){
                X_updateSU3();
                        }

    }

    for(int p=0; p<NConfigs*SweepFactor/numberOfMultiHit; p++){
    for(size_t i = 0; i<lattice.size(); i++ ){

        //this all belongs to the multi-thread part
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


            //in theory atomatically accepted
            if(j%overrelaxationStep && j!=0){
                U = overrelaxation(A, U);
                normalizeSU3Matrix(U);

            }
            else{
                size_t index = indexDist(indexing);
                X=XSet[index];

                UPrime = matrix_multiplication(X,lattice[idx(x,y,z,t, mu)]);
                normalizeSU3Matrix(UPrime);
                bool acceptance = latticeAction(lattice, lattice[idx(x,y,z,t, mu)], UPrime,x,y,z,t,mu, A);

                if(acceptance==true){
                    U=UPrime;
                }

                }
        }

        //multi thread stops

        //ensure that all U are inserted at their location


    }
        //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3();
                        }
        if(p % SweepFactor==0) {

            std::vector<std::complex<double>> loops;
            std::vector<double> r;
            if(observable == 1){
                wilsonLoop( lattice, loops, r, startingPoint,endPoint );
            }
            if(observable == 2){
                polyakovLoop( lattice, loops, r, startingPoint,endPoint );
            }
            observableCollected +=1;
            //saveArrayH5(r,);
            //saveArrayH5(loops,);
            
            }
    }
    }