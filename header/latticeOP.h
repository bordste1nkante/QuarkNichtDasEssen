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




// this mimics the behaviour of a 4D lattice from our 1D array
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu);

std::tuple <size_t, size_t, size_t, size_t, size_t> ReIdx(size_t idx);


//parallel
void X_updateSU3(size_t input);



void normalizeSU2Matrix(Matrix<2,2>& U);

void normalizeSU3Matrix(Matrix<rSU,rSU>& U);


//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice);

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice);


//parallel version
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime , const Matrix<rSU,cSU>& A, std::mt19937_64& ActionAccept, std::uniform_real_distribution<double>& Distribution);

//calculates the Plaquette values to determine lattice change
Matrix<rSU,cSU> determineA(const std::vector<Matrix<rSU,rSU>>& lattice ,size_t x, size_t y, size_t z, size_t t, size_t mu);

//parallization
Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U,std::uniform_int_distribution<int>& distribution,std::mt19937_64& random );

void wilsonLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint );


void polyakovLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint );


void plaquette(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<double>& plaquette);




