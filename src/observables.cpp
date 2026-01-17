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




    //std::vector<int> diff;
//
    //for(int i= 0; i<4; i++){
    //    diff.push_back(endPoint[i]-startingPoint[i]);
    //}


    
    
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
                        if(p% roundingFactor == 0 && p!=0){
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
                            Matrix<rSU, cSU> Path1;
                            Matrix<rSU, cSU> Path2;
                            std::vector<std::complex<double>> loops;
                            std::vector<double> r;

                            //wilson loop
                            if(observable == 1){
                                int dx=(endPoint[0]>startingPoint[0]) ? +1:-1;
                                int dy=(endPoint[1]>startingPoint[1]) ? +1:-1;                              
                                int dz=(endPoint[2]>startingPoint[2]) ? +1:-1;

                                int x= startingPoint[0];
                                int y= startingPoint[1];
                                int z= startingPoint[2];

                                int xDistance=0;
                                int yDistance=0;
                                int zDistance=0;



                                while (x!=endPoint[0])
                                {

                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,0)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,0)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_conjugate(Path1))));
                                    x+=dx;
                                    xDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));

                                }
                                while (y!=endPoint[0])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,1)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,1)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_conjugate(Path1))));
                                    y+=dy;
                                    yDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));                                

                                }
                                while (z!=endPoint[0])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,2)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,2)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_conjugate(Path1))));
                                    z+=dz;
                                    zDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));  
                                }

                                //now save loops and r


                            }
                            //polyakov loop
                            if(observable == 2){

                                int dt=(endPoint[3]>startingPoint[3]) ? +1:-1;
                                int x= startingPoint[0];
                                int y= startingPoint[1];
                                int z= startingPoint[2];
                                int t= startingPoint[3];

                                int tDistance= 0;

                                while(t!=endPoint[3]){
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,t,3)]);
                                    Path2=matrix_multiplication(Path1, lattice[idx(x,y,z,endPoint[3],3)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_conjugate(Path1))));
                                    r.push_back(std::sqrt(tDistance*tDistance));
                                    tDistance +=1;
                                    t+=dt;

                                }

                                //now save loops and r

                            }

                        }
                    }

                    }
            }
        }
    }
    }


}