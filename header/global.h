#pragma once
#include <stdlib.h>
#include <iostream> //std in/out etc.
#include <cmath>
#include <vector>
#include <complex>
#include <random>

// storage for a SU(N) matrix
template <size_t row, size_t column>
struct Matrix{
    std::complex<double> data[row][column];


    constexpr std::complex<double>& operator()(size_t r, size_t c) noexcept{
        return data[r][c];
    }


    const std::complex<double>& operator()(size_t r, size_t c) const{
        return data[r][c];
    }
};

//axes of lattice and number of independent link variable per site
extern size_t xAxis; 
extern size_t yAxis;
extern size_t zAxis;
extern size_t tAxis;
const size_t linksPerSite = 4;

//coupling strength
extern double beta;

//lattice constant -> lattice spacing
extern double a;


// row and column of SU(N)
const size_t rSU= 3;
const size_t cSU=3;



// list of N matrices and their inverse to generate new link variables -> improves efficiency -> avoids regeneration for every update
const int NSetXMatrices = 100;
extern std::vector<Matrix<rSU,rSU>> XSet;

//required for generation of X it controls spread around unity. Epsilon affects acceptance rate.
extern double epsilon;


//relevance?
extern std::uniform_real_distribution<double> dist;


//required for hot start, to create a set of not identity matrices, like generating X
extern double hotEpsilon;


//relevant?
extern std::uniform_real_distribution<double> hotDist;



//Pauli matrices and 2x2 identity
extern std::vector<Matrix<2,2>> pauliMatrices;
extern Matrix<2,2> identity2, Pauli1, Pauli2, Pauli3;

//rsU x cSU identity
extern Matrix<rSU,cSU> identityMatrix;

//zeroMatrices
extern Matrix<rSU,cSU> zeroMatrix;
extern Matrix<2,2> zeroMatrix2;


//h5-filename
extern std::string filenameh5;



//relevant factors for seeding
extern double drng1;
extern double dhotNumbS1;
extern double dhotNumbT1;
extern double dhotNumbR1;
extern double dindexing;



//some condition for equilibrium:

//controls size of checked set, we compare the change of means between sets
extern double PTestSizeGlobal;

//control allowed spread of means between sets
extern double changeRateLow;
extern double changeRateHigh;

//epsilon tuning:

//controls changefactor of epsilon
extern double alpha;
//what would be the ideal rate
extern double target_rate;

//target rate +- rateInterval controls accepted epsilon
extern double rateInterval;
