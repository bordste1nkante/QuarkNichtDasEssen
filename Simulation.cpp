#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <atomic>
#include <H5Cpp.h> //requires installation of this specific library on Linux:
#include <Eigen/Dense> //requires installation of this specific library on Linux:
#include "header/global.h"
#include "header/latticeOP.h"
#include "header/matrixOP.h"
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


    //epsilon spread depends on beta. hopefully speeds up thermalization
    hotEpsilon /= beta;
    

    drng1=information["seeds"]["dist"]["rng1"].as<double>();
    dhotNumb1=information["seeds"]["dist"]["hotNumb1"].as<double>();
    dhotNumb2=information["seeds"]["dist"]["hotNumb2"].as<double>();
    dhotNumb3=information["seeds"]["dist"]["hotNumb3"].as<double>();
    dhotNumb1Extra = information["seeds"]["dist"]["hotNumb1Extra"].as<double>();
    dhotNumb2Extra = information["seeds"]["dist"]["hotNumb2Extra"].as<double>();
    dhotNumb3Extra = information["seeds"]["dist"]["hotNumb3Extra"].as<double>();
    dhotNumbS1= information["seeds"]["dist"]["hotNumbS1"].as<double>();
    //dhotNumbS2= information["seeds"]["dist"]["hotNumbS2"].as<double>();
    //dhotNumbS3= information["seeds"]["dist"]["hotNumbS3"].as<double>();
    dhotNumbT1= information["seeds"]["dist"]["hotNumbT1"].as<double>();
    //dhotNumbT2= information["seeds"]["dist"]["hotNumbT2"].as<double>();
    //dhotNumbT3= information["seeds"]["dist"]["hotNumbT3"].as<double>();
    dhotNumbR1= information["seeds"]["dist"]["hotNumbR1"].as<double>();
    //dhotNumbR2= information["seeds"]["dist"]["hotNumbR2"].as<double>();
    //dhotNumbR3= information["seeds"]["dist"]["hotNumbR3"].as<double>();
    dindexing=information["seeds"]["indexDist"].as<double>();



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

    generate_zero();


    //generate the pauli matrices 0-3
    generate_Pauli();


    //generate rSU x cSU identity
    generate_identity();

    //create the first set of matrices X

    //get a random number
    X_updateSU3(100);



    if(coldOrHot == true){
        cold_start(lattice);
    }
    else{
        hot_start(lattice);
    }
    std::cout << "starting configuration created" << std::endl;


    Simulation(lattice,start,end,numberOfThermalSweeps,XUpdate,NConfigs,SweepFactor,observable,numberOfMultiHit,overrelaxationStep);
    std::cout << "simulation done" << std::endl;



    return 0;
}