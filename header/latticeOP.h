#pragma once
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"




// this mimics the behaviour of a 4D lattice from our 1D array
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu);


//updates the whole set of X matrices
void X_updateSU3();


//rounding errors causes matrices to potentially digress from det = 1, we correct that from time to time.
void normalizeSU3(std::vector<Matrix<rSU,rSU>>& lattice);


//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice);

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice);

//calculates the lattice action change and determines whether new configuration is accepted, U -> U'
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime ,size_t x, size_t y, size_t z, size_t t, size_t mu);


//set all temporal link variables to identity 
void temporalGauge(std::vector<Matrix<rSU,rSU>>& lattice);

//set all spatial link variables to identity
void spacialGauge(std::vector<Matrix<rSU,rSU>>& lattice);
