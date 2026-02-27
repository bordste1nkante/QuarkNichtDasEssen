#pragma once

#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <execution>
#include <thread>
#include <complex>
#include <random>
#include <Eigen/Dense>
#include "global.h"



//runs the simulation, including tuning and data aquisition (does not generate starting config)
void Simulation(   std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    size_t SweepFactor, 
    const size_t observable,
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep );

