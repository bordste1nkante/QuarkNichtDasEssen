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
void latticeSimulationPureMetropolis(
    std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable ){
    Matrix<rSU, cSU> UPrime;
    Matrix<rSU, cSU> X; 

    Matrix<rSU,cSU> A;

    size_t observableCollected= 0;


    //thermal sweeps, update the lattice a couple times so that equilibrium distribution can manifest
    for(int p=0; p<numberOfThermalSweeps; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu = 0; mu<linksPerSite; mu++){
                            size_t index = indexDist(indexing);
                            X=XSet[index];

                            UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);
                            normalizeSU3Matrix(UPrime);


                            A=determineA(lattice,l,k,j,i,mu);

                            // calculate lattice action change and change or not change the lattice
                            bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu, A);
                            if (acceptance == true)
                            {
                                lattice[idx(l,k,j,i, mu)] = UPrime;
                            }
                            


                        }

                    }
                }
            }
        }
        //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3();
        }
    }

    std::cout << "Thermalization complete "<<std::endl;



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
                        normalizeSU3Matrix(UPrime);


                        A=determineA(lattice,l,k,j,i,mu);
                        // calculate lattice action change and change or not change the lattice
                        bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu,A );
                        if (acceptance == true)
                        {
                            lattice[idx(l,k,j,i, mu)] = UPrime;
                            }

                        }
                    }
                }
            }
        }
                        //Update X matrices
                        if(p%XUpdate ==0 && p!=0){
                            X_updateSU3();
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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    x+=dx;
                                    xDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));

                                }
                                while (y!=endPoint[1])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,1)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,1)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    y+=dy;
                                    yDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));                                

                                }
                                while (z!=endPoint[2])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,2)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,2)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    z+=dz;
                                    zDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));  
                                }
                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;



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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));

                                    //save t-distance in r array (for simplicity)
                                    r.push_back(std::sqrt(tDistance*tDistance));
                                    tDistance +=1;
                                    t+=dt;

                                }

                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;

                            }

                        }
                    }

                    }





// run the simulation
void latticeSimulationPureMetropolisMultiHit(
    std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable,
    const size_t numberOfMultiHit ){


    Matrix<rSU, cSU> UPrime;
    Matrix<rSU, cSU> X; 

    Matrix<rSU,cSU> A;

    size_t observableCollected= 0;


    //thermal sweeps, update the lattice a couple times so that equilibrium distribution can manifest
    for(int p=0; p<numberOfThermalSweeps/numberOfMultiHit; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu = 0; mu<linksPerSite; mu++){
                            A=determineA(lattice,l,k,j,i,mu);
                            for(int multi = 0; multi < numberOfMultiHit; multi++){

                                size_t index = indexDist(indexing);
                                X=XSet[index];

                                UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);
                                normalizeSU3Matrix(UPrime);




                                // calculate lattice action change and change or not change the lattice
                                bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu, A);
                                if (acceptance == true)
                                {
                                    lattice[idx(l,k,j,i, mu)] = UPrime;
                                }

                            }

                            

                        }

                    }
                }
            }
        }

        //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3();
        }
    }

    std::cout << "Thermalization complete "<<std::endl;



    //actual data generation save a config after sufficient update to avoid autocorrelation 
    for(int p=0; p<NConfigs*SweepFactor/numberOfMultiHit; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu=0; mu<linksPerSite; mu++){
                        A=determineA(lattice,l,k,j,i,mu);

                            for(int multi = 0; multi < numberOfMultiHit; multi++){

                                size_t index = indexDist(indexing);
                                X=XSet[index];

                                UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);
                                normalizeSU3Matrix(UPrime);




                                // calculate lattice action change and change or not change the lattice
                                bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu, A);
                                if (acceptance == true)
                                {
                                    lattice[idx(l,k,j,i, mu)] = UPrime;
                                }

                            }

                        }
                    }
                }
            }
        }
                        //Update X matrices
                        if(p%XUpdate ==0 && p!=0){
                            X_updateSU3();
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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    x+=dx;
                                    xDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));

                                }
                                while (y!=endPoint[1])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,1)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,1)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    y+=dy;
                                    yDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));                                

                                }
                                while (z!=endPoint[2])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,2)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,2)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    z+=dz;
                                    zDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));  
                                }
                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;



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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));

                                    //save t-distance in r array (for simplicity)
                                    r.push_back(std::sqrt(tDistance*tDistance));
                                    tDistance +=1;
                                    t+=dt;

                                }

                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;

                            }

                        }
                    }

                    }















// run the simulation
void latticeSimulationMixedMetropolisMultiHit(
    std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable, const size_t numberOfMultiHit, const size_t overrelaxationStep ){


    Matrix<rSU, cSU> UPrime;
    Matrix<rSU, cSU> X; 

    Matrix<rSU,cSU> A;

    size_t observableCollected= 0;


    //thermal sweeps, update the lattice a couple times so that equilibrium distribution can manifest
    for(int p=0; p<numberOfThermalSweeps/numberOfMultiHit; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu = 0; mu<linksPerSite; mu++){
                            A=determineA(lattice,l,k,j,i,mu);
                            for(int multi = 0; multi < numberOfMultiHit; multi++){

                                if(multi % overrelaxationStep ==0 && multi !=0){
                                    //UPrime = overrelaxation(lattice[idx(l,k,j,i, mu)]);
                                }
                                else{
                                    size_t index = indexDist(indexing);
                                    X=XSet[index];

                                    UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);
                                    normalizeSU3Matrix(UPrime);

                                }






                                // calculate lattice action change and change or not change the lattice
                                bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu, A);
                                if (acceptance == true)
                                {
                                    lattice[idx(l,k,j,i, mu)] = UPrime;
                                }

                            }

                            

                        }

                    }
                }
            }
        }

        //Update X matrices
        if(p%XUpdate ==0 && p!=0){
            X_updateSU3();
        }
    }

    std::cout << "Thermalization complete "<<std::endl;



    //actual data generation save a config after sufficient update to avoid autocorrelation 
    for(int p=0; p<NConfigs*SweepFactor/numberOfMultiHit; p++){
        for(int i = 0; i<tAxis; i++){
            for(int j= 0; j<zAxis; j++){
                for(int k=0; k<yAxis; k++){
                    for(int l=0; l<xAxis; l++){
                        for(int mu=0; mu<linksPerSite; mu++){
                        A=determineA(lattice,l,k,j,i,mu);

                            for(int multi = 0; multi < numberOfMultiHit; multi++){

                                if(multi % overrelaxationStep ==0 && multi !=0){
                                   // UPrime = overrelaxation(lattice[idx(l,k,j,i, mu)]);
                                }
                                else{
                                    size_t index = indexDist(indexing);
                                    X=XSet[index];

                                    UPrime = matrix_multiplication(X,lattice[idx(l,k,j,i, mu)]);
                                    normalizeSU3Matrix(UPrime);

                                }




                                // calculate lattice action change and change or not change the lattice
                                bool acceptance = latticeAction(lattice, lattice[idx(l,k,j,i, mu)], UPrime,l,k,j,i,mu, A);
                                if (acceptance == true)
                                {
                                    lattice[idx(l,k,j,i, mu)] = UPrime;
                                }

                            }

                        }
                    }
                }
            }
        }
                        //Update X matrices
                        if(p%XUpdate ==0 && p!=0){
                            X_updateSU3();
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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    x+=dx;
                                    xDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));

                                }
                                while (y!=endPoint[1])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,1)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,1)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    y+=dy;
                                    yDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));                                

                                }
                                while (z!=endPoint[2])
                                {
                                    Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,2)]);
                                    Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,2)]);
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
                                    z+=dz;
                                    zDistance+=1;
                                    r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));  
                                }
                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;



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
                                    loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));

                                    //save t-distance in r array (for simplicity)
                                    r.push_back(std::sqrt(tDistance*tDistance));
                                    tDistance +=1;
                                    t+=dt;

                                }

                                observableCollected += 1;
                                //now save loops and r
                                std::cout << "Extracted observable: " << observableCollected << "/" << NConfigs <<std::endl;

                            }

                        }
                    }

                    }

