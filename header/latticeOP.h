#pragma once
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <execution>
#include <thread>
#include <Eigen/Dense>
#include "../header/global.h"




//from 4d lattice (+1 link variable) point obtain index of 1d storage.
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu);

//from storage index obtain link variable
std::tuple <size_t, size_t, size_t, size_t, size_t> ReIdx(size_t idx);


//Updates the set of X matrices from metropolis, XU=U'. We generate three SU2 matrices and combine them
void X_updateSU3(size_t input);


//fix 2x2 matrix to SU2
void normalizeSU2Matrix(Matrix<2,2>& U);
//fix 3x3 matrix to SU3
void normalizeSU3Matrix(Matrix<rSU,rSU>& U);


//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice);

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice);


//calculate lattice action change and determine wheter proposal is accepted
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime , const Matrix<rSU,cSU>& A, std::mt19937_64& ActionAccept, std::uniform_real_distribution<double>& Distribution);

//calculates the stapel
Matrix<rSU,cSU> determineA(const std::vector<Matrix<rSU,rSU>>& lattice ,size_t x, size_t y, size_t z, size_t t, size_t mu);

//perform an overrelaxation step
Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U,std::uniform_int_distribution<int>& distribution,std::mt19937_64& random );


//calculate WilsonLoop
void wilsonLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, const double T );

//calculate polyakovLoop
void polyakovLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint );

//calculate average plaquette for all lattice points
void plaquette(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<double>& plaquette);




