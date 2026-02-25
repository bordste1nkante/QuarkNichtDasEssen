#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <execution>
#include <thread>
#include <complex>
#include <atomic>
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
    size_t SweepFactor, 
    const size_t observable,
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep ){


    std::atomic<int> acceptanceRate=0;
    std::atomic<int> updates=0;

    size_t observableCollected= 0;
    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);



    //thermalization
    for(int p=0; p<numberOfThermalSweeps/numberOfMultiHit; p++){
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


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    normalizeSU3Matrix(U);

                    }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,U);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, U, UPrime, A, Threadindexing, threadAcceptReject);


                    if(acceptance==true){
                        U=UPrime;
                        acceptanceRate.fetch_add(1);
                        }
                    updates.fetch_add(1);
                    }


                }
            bufferLattice[i]=U;

        });


        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);
            //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            double rate = double(acceptanceRate.load())/double(updates.load());
            std::cout << rate << std::endl;
           // epsilon *= 1 + alpha * (rate - target_rate);
            //if(rate < 0.45){
            //    if(rate < 0.35){
            //        epsilon *=0.9;
            //    }
            //    else{
            //        epsilon *=0.95;
            //    }
//
            //}
            //if(rate>0.6){
            //    if(rate < 0.7){
            //        epsilon*=1.1;
            //    }
            //    else{
            //        epsilon*=1.05;
            //    }
//
            //}
            if(epsilon>0.7){
                epsilon=0.7;
            }
            acceptanceRate.store(0);
            updates.store(0);
            X_updateSU3(p*numberOfThermalSweeps+270);
                        }
        

    }
    //std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
    //plaquette(lattice, plaquettes);
    //for(int i =0; i< plaquettes.size(); i++){
    //    std::cout << plaquettes[i] << std::endl;
    //}
    //double PTest=average(plaquettes);



    //std::cout <<"P: " << PTest << std::endl;
    //double rate = double(acceptanceRate.load())/double(updates.load());
    //std::cout << rate << std::endl;


    //tune params
    //epsilonTune(lattice, numberOfThermalSweeps, XUpdate, numberOfMultiHit,overrelaxationStep);
    //ThermalTune(lattice,numberOfThermalSweeps, XUpdate,numberOfMultiHit,overrelaxationStep);
    //SweepFactor = AutoCorrelationTune(lattice, numberOfThermalSweeps, XUpdate, numberOfMultiHit,overrelaxationStep);
    //SaveTune_H5(const double epsilon, const double avgPlaq, const size_t SweepFactor);


    for(int p=0; p<NConfigs*SweepFactor; p++){

        //parallelization
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*19007+ (p+2)* 12111+ 3011);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);


            std::uniform_real_distribution<double> threadAcceptReject(0,1);
        

            std::uniform_int_distribution<int> threadrefelctDist(1,3);


            Matrix<rSU, cSU> UPrime, X, A, U, UAccepted; 

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
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    //std::cout << matrix_trace(matrix_multiplication(U,matrix_hermitean_conjugate(U))) << std::endl;
                    
                    //bool acceptance = latticeAction(lattice, lattice[i], U, A, Threadindexing, threadAcceptReject);
                    //normalizeSU3Matrix(U);


//
                    //if(acceptance==true){
                    //    std::cout << "acceptedOver" << std::endl; 
                    //    }
                    //else{
                    //    std::cout << "falseOver" << std::endl;
                    //}
                    normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,U);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, U, UPrime, A, Threadindexing, threadAcceptReject);
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
            X_updateSU3(p*p+SweepFactor);
                        }
        if(p % SweepFactor==0) {

            std::vector<std::complex<double>> loops;
            std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
            std::vector<double> r;

            if(observable==0){
                plaquette(lattice,plaquettes);
                std::string pathL = "/Configuration/"+ std::to_string(observableCollected)+ "/plaquette";
                saveArrayH5(plaquettes,pathL);
            }
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