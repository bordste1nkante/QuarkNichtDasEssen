#pragma once
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>




// storage for a SU(N) matrix
template <size_t row, size_t column>
struct Matrix{
    std::complex<double> data;


    constexpr std::complex<double>& operator()(size_t r, size_t c) noexcept;


    const std::complex<double>& operator()(size_t r, size_t c) const;

};


//axes of lattice and number of independent link variable per site
extern size_t xAxis; 
extern size_t yAxis;
extern size_t zAxis;
extern size_t tAxis;
const size_t linksPerSite = 4;

//coupling strength
extern double beta;

//lattice constant
extern double a;


// row and column of SU(N)
const size_t rSU= 3;
const size_t cSU=3;



// list of N matrices and their inverse to generate new link variables
const int NSetXMatrices = 50;
std::vector<Matrix<rSU,rSU>> XSet;

//required for generation of X. Epsilon affects acceptance rate.
const double epsilon = 0.4;

extern std::mt19937_64 randNumb;
extern std::uniform_real_distribution<double> distEpsilon;

extern std::mt19937_64 rng;
extern std::uniform_real_distribution<double> dist;


//required for hot start, to create a set of not identity matrices, like generating X, but with bigger spread 
const double hotEpsilon = 0.7;


//uniform dist for generation of random variables in hot_start()
extern std::mt19937_64 hotNumb;
extern std::uniform_real_distribution<double> hotDist;



//uniform dist for selecting indece of XSet
extern std::mt19937_64 indexing;
extern std::uniform_int_distribution<int> indexDist;

//uniform for accept/reject of new config
extern std::mt19937_64 acceptReject;
extern std::uniform_real_distribution<double> uniformAcceptReject;



//Pauli matrices and 2x2 identity
std::vector<Matrix<2,2>> pauliMatrices;
Matrix<2,2> identity2, Pauli1,Pauli2, Pauli3;

//rsU x cSU identity
Matrix<rSU,cSU> identityMatrix;

