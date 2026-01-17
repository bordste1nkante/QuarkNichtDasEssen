#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"


// run the simulation
void latticeSimulationPureMetropolis(std::vector<Matrix<rSU,cSU>>& lattice, const std::vector<size_t>& startingPoint, const std::vector<size_t>& endPoint, const size_t numberOfThermalSweeps, const size_t roundingFactor, const size_t XUpdate, const size_t NConfigs, const size_t SweepFactor, const size_t observable ){
    Matrix<rSU, cSU> UPrime;
    Matrix<rSU, cSU> X; 
    
    
    //function to establish temporal gauge
    //temporalGauge(lattice);

    //run only over other variables

    //thermal sweeps, update the lattice a couple times so that equilibrium distribution can manifest
    for(int p=0; p<numberOfThermalSweeps; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for (int mu = 0; mu<linksPerSite; mu++){
                        size_t index = indexDist(indexing);
                        X=XSet[index];

                        UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);

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
    }




    //actual data generation save a config after sufficient update to avoid autocorrelation 
    for(int p=0; p<NConfigs*SweepFactor; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu=0; mu<linksPerSite; mu++){
                        size_t index = indexDist(indexing);
                        X=XSet[index];

                        UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);

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


}