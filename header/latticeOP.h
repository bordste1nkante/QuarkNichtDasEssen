#pragma once
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <Eigen/Dense>
#include "../header/global.h"




// this mimics the behaviour of a 4D lattice from our 1D array
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu);

std::tuple <size_t, size_t, size_t, size_t, size_t> ReIdx(size_t idx);


//updates the whole set of X matrices
void X_updateSU3();


//rounding errors causes matrices to potentially digress from det = 1, we correct that from time to time, also we project to unitarity.
void normalizeSU3(std::vector<Matrix<rSU,rSU>>& lattice);

void normalizeSU3Matrix(Matrix<rSU,rSU>& U);


//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice);

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice);

//calculates the lattice action change and determines whether new configuration is accepted, U -> U'
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime ,size_t x, size_t y, size_t z, size_t t, size_t mu, const Matrix<rSU,cSU>& A);

//calculates the Plaquette values to determine lattice change
Matrix<rSU,cSU> determineA(const std::vector<Matrix<rSU,rSU>>& lattice ,size_t x, size_t y, size_t z, size_t t, size_t mu);

//determine U' by overrelaxation
Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U);


void wilsonLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint );


void polyakovLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint );



//set all temporal link variables to identity 
void temporalGauge(std::vector<Matrix<rSU,rSU>>& lattice);

//set all spatial link variables to identity
void spacialGauge(std::vector<Matrix<rSU,rSU>>& lattice);



