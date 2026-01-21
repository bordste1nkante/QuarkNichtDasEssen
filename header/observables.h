#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "global.h"
#include "latticeOP.h"
#include "matrixOP.h"


// run the simulation for wilson loop in temporal gauge (all temporal links are identity), also store data
 void latticeSimulationPureMetropolis(
    std::vector<Matrix<rSU,cSU>>& lattice, 
    const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps, 
    const size_t roundingFactor, 
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
    const size_t observable);

//might add stuff for heatbath etc.