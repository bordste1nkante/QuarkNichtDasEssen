#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "header/global.h"
#include "header/latticeOP.h"
#include "header/matrixOP.h"




int main(){

    //obtain .yaml name
    std::string input;
    std::cout << "Enter name of yaml file: ";
    std::getline(std::cin, input);
    input += ".yaml";

    YAML::Node information = YAML::LoadFile(input);

    // assigned them the correct value by reading it from yaml file
    xAxis = information["lattice"]["x"].as<size_t>(); 
    yAxis = information["lattice"]["y"].as<size_t>();
    zAxis = information["lattice"]["z"].as<size_t>();
    tAxis = information["lattice"]["t"].as<size_t>();


    beta = information["lattice"]["beta"].as<double>();
    a= information["lattice"]["lattice spacing"].as<double>();

    // if true we have a cold start else a hot start, read from yaml
    bool coldOrHot = information["startConfig"].as<bool>();

    size_t numberOfThermalSweeps = information["updates"]["NSweepsThermal"].as<size_t>(); // number of complete lattice updates till you start the data run
    size_t NConfigs = information["updates"]["NSweepsThermal"].as<size_t>(); //Number of configs for analysis
    size_t SweepFactor = information["updates"]["Sweep"].as<size_t>();//autocorrelation needs to be overcome, wait some iterations before collecting the next config.
    size_t roundingFactor = information["updates"]["Rounding"].as<size_t>(); //rounding errors need to be corrected
    size_t XUpdate = information["updates"]["XUpdate"].as<size_t>(); //how often to generate new X

    // our lattice as 1D array of matrices (3x3), factor 4 because every lattice site has 4 link variable (technically 8, but hermitean conjugate reduces it to 4 indepent ones)
    std::vector<Matrix<rSU,rSU>> lattice(4*xAxis*yAxis*zAxis*tAxis);

    Matrix<rSU, cSU> UPrime;
    Matrix<rSU, cSU> X;    

    //generate the pauli matrices 0-3
    generate_Pauli();
    //create the first set of matrices X
    X_updateSU3();

    if(coldOrHot == true){
        cold_start(lattice);
    }
    else{
        hot_start(lattice);
    }

    //thermal sweeps, update the lattice a couple times so that equilibrium distribution can manifest
    for(int p=0; p<numberOfThermalSweeps; p++)
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for (int mu = 0; mu<linksPerSite; mu++){
                        size_t index = indexDist(indexing);
                        X=XSet[index];

                        //matrix multiplication X*U = U'
                        for(int i=0; i<rSU; i++){
                            for(int j=0; j<cSU; j++){
                                std::complex<double> sum;
                                for(int k = 0; k<rSU; k++){
                                    sum += X(i,k)*lattice[idx(l,k,j,i, mu)](k,j);
                                
                                }
                                UPrime(i,j)= sum;
                            }
                        
                        }
                        // calculate lattice action change and change or not change the lattice
                        bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu );
                        if (acceptance == true)
                        {
                            lattice[idx(l,k,j,i, mu)] = UPrime;
                        }
                        

                        //Update X matrices
                        if(p%XUpdate ==0 && p!=0){
                            X_updateSU3();
                        }

                        // from time to time our matrices have to be projected to det=1, rounding errors cause trouble and like X is also not neccesarily det 1, right?
                        if (p% roundingFactor == 0 && p!=0){
                            normalizeSU3(lattice);
                        }
                    }

                }
            }
        }
    }
    //actual data generation save a config after sufficient update to avoid autocorrelation 
    for(int p=0; p<NConfigs*SweepFactor; p++)
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu=0; mu<linksPerSite; mu++){
                        size_t index = indexDist(indexing);
                        X=XSet[index];

                        //matrix multiplication X*U = U'
                        for(int i=0; i<rSU; i++){
                            for(int j=0; j<cSU; j++){
                                std::complex<double> sum;
                                for(int k = 0; k<rSU; k++){
                                    sum += X(i,k)*lattice[idx(l,k,j,i, mu)](k,j);
                                
                                }
                                UPrime(i,j)= sum;
                            }
                        
                        }
                        // calculate lattice action change and change or not change the lattice
                        bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu );
                        if (acceptance == true)
                        {
                            lattice[idx(l,k,j,i, mu)] = UPrime;
                        }
                        

                        //Update X matrices
                        if(p%XUpdate ==0 && p!=0){
                            X_updateSU3();
                        }
                        // from time to time our matrices have to be projected to det=1, rounding errors cause trouble and like X is also not neccesarily det 1, right?
                        if (p% roundingFactor == 0 && p!=0){
                            normalizeSU3(lattice);
                        }
                        //save these configs, the rest is updates and after sufficient sweeps autocorrelation becomes negliable 
                        if(p % SweepFactor==0) {

                        }
                    }

                    }
            }
        }
    }

}