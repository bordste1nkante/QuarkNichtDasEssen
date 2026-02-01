#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <H5Cpp.h> //requires installation of this specific library on Linux:
#include <Eigen/Dense> //requires installation of this specific library on Linux:
#include "header/global.h"
#include "header/latticeOP.h"
#include "header/matrixOP.h"
#include "header/observables.h"
#include "header/utils.h"
#include "header/Finish.h"




int main(){

    //obtain .yaml name
    std::string input;
    std::cout << "Enter name of yaml file: ";
    std::getline(std::cin, input);
    input += ".yaml";
    //std::cout << "input read" << std::endl;

    YAML::Node information = YAML::LoadFile(input);

    // assigned them the correct value by reading it from yaml file
    xAxis = information["lattice"]["x"].as<size_t>(); 
    yAxis = information["lattice"]["y"].as<size_t>();
    zAxis = information["lattice"]["z"].as<size_t>();
    tAxis = information["lattice"]["t"].as<size_t>();


    beta = information["lattice"]["beta"].as<double>();
    


    // if true we have a cold start else a hot start, read from yaml
    bool coldOrHot = information["startConfig"].as<bool>();

    size_t numberOfThermalSweeps = information["updates"]["NSweepsThermal"].as<size_t>(); // number of complete lattice updates till you start the data run
    size_t NConfigs = information["updates"]["NConfig"].as<size_t>(); //Number of configs for analysis
    size_t SweepFactor = information["updates"]["Sweep"].as<size_t>();//autocorrelation needs to be overcome, wait some iterations before collecting the next config.
    size_t XUpdate = information["updates"]["XUpdate"].as<size_t>(); //how often to generate new X
    size_t numberOfMultiHit = information["updates"]["Multi"].as<size_t>(); // how many hits per link in a row
    size_t overrelaxationStep = information["updates"]["overrelaxation"].as<size_t>(); // 



    epsilon = information["epsilon"].as<double>();
    
    //place seeds
    randNumb.seed(information["seeds"]["distEpsilon"].as<double>());
    rng1.seed(information["seeds"]["dist"]["rng1"].as<double>());
    rng2.seed(information["seeds"]["dist"]["rng2"].as<double>());
    rng3.seed(information["seeds"]["dist"]["rng3"].as<double>());
    hotNumb1.seed(information["seeds"]["dist"]["hotNumb1"].as<double>());
    hotNumb2.seed(information["seeds"]["dist"]["hotNumb2"].as<double>());
    hotNumb3.seed(information["seeds"]["dist"]["hotNumb3"].as<double>());
    hotNumb1Extra.seed(information["seeds"]["dist"]["hotNumb1Extra"].as<double>());
    hotNumb2Extra.seed(information["seeds"]["dist"]["hotNumb2Extra"].as<double>());
    hotNumb3Extra.seed(information["seeds"]["dist"]["hotNumb3Extra"].as<double>());
    hotNumbS1.seed(information["seeds"]["dist"]["hotNumbS1"].as<double>());
    hotNumbS2.seed(information["seeds"]["dist"]["hotNumbS2"].as<double>());
    hotNumbS3.seed(information["seeds"]["dist"]["hotNumbS3"].as<double>());
    hotNumbT1.seed(information["seeds"]["dist"]["hotNumbT1"].as<double>());
    hotNumbT2.seed(information["seeds"]["dist"]["hotNumbT2"].as<double>());
    hotNumbT3.seed(information["seeds"]["dist"]["hotNumbT3"].as<double>());
    hotNumbR1.seed(information["seeds"]["dist"]["hotNumbR1"].as<double>());
    hotNumbR2.seed(information["seeds"]["dist"]["hotNumbR2"].as<double>());
    hotNumbR3.seed(information["seeds"]["dist"]["hotNumbR3"].as<double>());
    indexing.seed(information["seeds"]["indexDist"].as<double>());
    acceptReject.seed(information["seeds"]["uniformAcceptReject"].as<double>());
    reflection.seed(information["seeds"]["overrelaxation"].as<double>());


    drandNumb = information["seeds"]["distEpsilon"].as<double>();
    drng1=information["seeds"]["dist"]["rng1"].as<double>();
    drng2=information["seeds"]["dist"]["rng2"].as<double>();
    drng3=information["seeds"]["dist"]["rng3"].as<double>();
    dhotNumb1=information["seeds"]["dist"]["hotNumb1"].as<double>();
    dhotNumb2=information["seeds"]["dist"]["hotNumb2"].as<double>();
    dhotNumb3=information["seeds"]["dist"]["hotNumb3"].as<double>();
    dhotNumb1Extra = information["seeds"]["dist"]["hotNumb1Extra"].as<double>();
    dhotNumb2Extra = information["seeds"]["dist"]["hotNumb2Extra"].as<double>();
    dhotNumb3Extra = information["seeds"]["dist"]["hotNumb3Extra"].as<double>();
    dhotNumbS1= information["seeds"]["dist"]["hotNumbS1"].as<double>();
    dhotNumbS2= information["seeds"]["dist"]["hotNumbS2"].as<double>();
    dhotNumbS3= information["seeds"]["dist"]["hotNumbS3"].as<double>();
    dhotNumbT1= information["seeds"]["dist"]["hotNumbT1"].as<double>();
    dhotNumbT2= information["seeds"]["dist"]["hotNumbT2"].as<double>();
    dhotNumbT3= information["seeds"]["dist"]["hotNumbT3"].as<double>();
    dhotNumbR1= information["seeds"]["dist"]["hotNumbR1"].as<double>();
    dhotNumbR2= information["seeds"]["dist"]["hotNumbR2"].as<double>();
    dhotNumbR3= information["seeds"]["dist"]["hotNumbR3"].as<double>();
    dindexing=information["seeds"]["indexDist"].as<double>();
    dacceptReject=information["seeds"]["uniformAcceptReject"].as<double>();
    dreflection=information["seeds"]["overrelaxation"].as<double>();


    std::vector<size_t> start = information["positions"]["startPoint"].as<std::vector<size_t>>();
    std::vector<size_t> end = information["positions"]["endPoint"].as<std::vector<size_t>>();
    size_t observable = information["observable"].as<size_t>();
    filenameh5 = information["h5Name"].as<std::string>();


 

    bool stop = Setup_H5(coldOrHot, start, end, numberOfThermalSweeps,XUpdate,NConfigs,SweepFactor,observable);

    if(stop == false){
        std::cout << "Something went wrong with the h5 file" << std::endl;
        return 0; 
    }




    // our lattice as 1D array of matrices (3x3), factor 4 because every lattice site has 4 link variable (technically 8, but hermitean conjugate reduces it to 4 indepent ones)
    std::vector<Matrix<rSU,rSU>> lattice(4*xAxis*yAxis*zAxis*tAxis);//why 4* ?


    //std::cout << lattice.size() << std::endl;
    //std::vector<size_t> indices(lattice.size());
    //std::iota(indices.begin(), indices.end(),0);
//
    //for(int i= 0; i<indices.size();i++){
    //    //std::cout << indices[i] << std::endl;
    //    std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(i);
    //    size_t mu,x,y,z,t;
    //    mu= std::get<0>(temp);
    //    x= std::get<1>(temp);
    //    y= std::get<2>(temp);
    //    z= std::get<3>(temp);
    //    t= std::get<4>(temp);
//
    //    if(t>5){
    //        std::cout << "t problem" << std::endl;
    //    }
    //    if(x>5){
    //        std::cout << "x problem" << std::endl;
    //    }
    //    if(y>5){
    //        std::cout << "y problem" << std::endl;
    //    }
    //    if(z>5){
    //        std::cout << "z problem" << std::endl;
    //    }
    //    if(mu>4){
    //        std::cout << "mu problem" << std::endl;
    //    }
//
    //}
//
    //std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx();
    //size_t mu,x,y,z,t;
    //mu= std::get<0>(temp);
    //x= std::get<1>(temp);
    //y= std::get<2>(temp);
    //z= std::get<3>(temp);
    //t= std::get<4>(temp);
//
    //std::cout << mu << << std::endl;
    //std::cout << "lattice vector done" << std::endl;

    //generate the pauli matrices 0-3
    generate_Pauli();

    //std::cout << "Pauli" << std::endl;
    //generate rSU x cSU identity
    generate_identity();
    //std::cout << "identity" << std::endl;
    //create the first set of matrices X
    X_updateSU3();


    //std::cout << XSet.size() << std::endl;

    if(coldOrHot == true){
        cold_start(lattice);
    }
    else{
        hot_start(lattice);
    }
    std::cout << "starting configuration created" << std::endl;


    //latticeSimulationPureMetropolis(lattice,start, end, numberOfThermalSweeps,XUpdate,NConfigs,SweepFactor,observable);
    Simulation(lattice,start,end,numberOfThermalSweeps,XUpdate,NConfigs,SweepFactor,observable,numberOfMultiHit,overrelaxationStep);
    std::cout << "simulation done" << std::endl;



    return 0;
}