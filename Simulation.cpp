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
    //std::cout << "input read" << std::endl;

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
    size_t NConfigs = information["updates"]["NConfig"].as<size_t>(); //Number of configs for analysis
    size_t SweepFactor = information["updates"]["Sweep"].as<size_t>();//autocorrelation needs to be overcome, wait some iterations before collecting the next config.
    size_t roundingFactor = information["updates"]["Rounding"].as<size_t>(); //rounding errors need to be corrected
    size_t XUpdate = information["updates"]["XUpdate"].as<size_t>(); //how often to generate new X

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


    std::vector<size_t> start = information["positions"]["startPoint"].as<std::vector<size_t>>();
    std::vector<size_t> end = information["positions"]["endPoint"].as<std::vector<size_t>>();
    size_t observable = information["observable"].as<size_t>();

    //std::cout << "yaml read" << std::endl;

    // our lattice as 1D array of matrices (3x3), factor 4 because every lattice site has 4 link variable (technically 8, but hermitean conjugate reduces it to 4 indepent ones)
    std::vector<Matrix<rSU,rSU>> lattice(4*xAxis*yAxis*zAxis*tAxis);//why 4* ?

    //std::cout << "lattice vector done" << std::endl;

    //generate the pauli matrices 0-3
    generate_Pauli();

    //std::cout << "Pauli" << std::endl;
    //generate rSU x cSU identity
    generate_identity();
    //std::cout << "identity" << std::endl;
    //create the first set of matrices X
    X_updateSU3();
    //std::cout << "X" << std::endl;

    if(coldOrHot == true){
        cold_start(lattice);
    }
    else{
        hot_start(lattice);
    }
    std::cout << "starting configuration created" << std::endl;
    latticeSimulationPureMetropolis(lattice,start, end, numberOfThermalSweeps,roundingFactor,XUpdate,NConfigs,SweepFactor,observable);

    std::cout << "simulation done" << std::endl;
    //add heatbath with and without overrelacation
    //vielleicht als variable übergeben welcher Algorithmus gewählt werden soll?


    return 0;
}