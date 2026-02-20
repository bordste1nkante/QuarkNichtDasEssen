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

bool Setup_H5(bool ColdOrHot, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps,  
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable);


bool saveArrayH5(const std::vector<double>& array, std::string dataSetPath);

bool saveArrayH5complex(const std::vector<std::complex<double>>& array, std::string dataSetPath);


void ThermalAndCorrelationTune(std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep);

    
void epsilonTune( std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep);

Eigen::Matrix3cd translateMatrices(const Matrix<rSU,cSU>& A);


//translate eigen to our matrices, we use them :(
Matrix<rSU,cSU> retranslateMatrices(const Eigen::Matrix3cd& A);