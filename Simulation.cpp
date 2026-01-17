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
#include "header/observables.h"




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

    //place seeds
    randNumb.seed(information["seeds"]["distEpsilon"].as<double>());
    rng.seed(information["seeds"]["dist"].as<double>());
    hotNumb.seed(information["seeds"]["hotDist"].as<double>());
    indexing.seed(information["seeds"]["indexDist"].as<double>());
    acceptReject.seed(information["seeds"]["uniformAcceptReject"].as<double>());


    std::vector<size_t> start = information["positions"]["startPoint"].as<std::vector<size_t>>();
    std::vector<size_t> end = information["positions"]["endPoint"].as<std::vector<size_t>>();
    size_t observable = information["observable"].as<size_t>();



    // our lattice as 1D array of matrices (3x3), factor 4 because every lattice site has 4 link variable (technically 8, but hermitean conjugate reduces it to 4 indepent ones)
    std::vector<Matrix<rSU,rSU>> lattice(4*xAxis*yAxis*zAxis*tAxis);

   

    //generate the pauli matrices 0-3
    generate_Pauli();
    //generate rSU x cSU identity
    generate_identity();
    //create the first set of matrices X
    X_updateSU3();

    if(coldOrHot == true){
        cold_start(lattice);
    }
    else{
        hot_start(lattice);
    }
    latticeSimulationPureMetropolis(lattice,start, end, numberOfThermalSweeps,roundingFactor,XUpdate,NConfigs,SweepFactor,observable);

}