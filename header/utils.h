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


bool Setup_H5(bool ColdOrHot, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps,  
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable);


bool saveArrayH5(const std::vector<size_t>& array, std::string dataSetPath);


Eigen::Matrix3d translateMatrices(const Matrix<rSU,cSU>& A);


//translate eigen to our matrices, we use them :(
Matrix<rSU,cSU> retranslateMatrices(const Eigen::Matrix3cd& A);