#pragma once

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

void ensureGroup(H5::H5File& file, const std::string& path);


//stores essential meta data
bool Setup_H5(bool ColdOrHot, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps,  
    const size_t XUpdate, 
    const size_t NConfigs,  
    const size_t observable);


//stores tuning results
bool SaveTune_H5(const double epsilon, const double avgPlaq, const size_t SweepFactor);

//stores an array in h5
bool saveArrayH5(const std::vector<double>& array, std::string dataSetPath);

//stores a complex array
bool saveArrayH5complex(const std::vector<std::complex<double>>& array, std::string dataSetPath);


//calculate the average of a vector
double average(const std::vector<double>& vec);

//calculate the correlation between to vectors
double correlationFunc(const std::vector<double>& Plaqs,const std::vector<double>&OGPlaqs );


//test for equilibrium via a mean stabilization
void ThermalTune(std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep);

    
//tune epsilon to obtain a wanted acceptance rate
void epsilonTune( std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep);


//determine autocorrelation in computer time and return it for use in data aquisition
size_t AutoCorrelationTune( std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep);


//change matrix style to Eigen
Eigen::Matrix3cd translateMatrices(const Matrix<rSU,cSU>& A);


//translate eigen to our matrices, we use them :(
Matrix<rSU,cSU> retranslateMatrices(const Eigen::Matrix3cd& A);