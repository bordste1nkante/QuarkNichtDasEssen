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
const int NSetXMatrices = 50;
extern std::vector<Matrix<rSU,rSU>> XSet;

//required for generation of X. Epsilon affects acceptance rate.
extern double epsilon;

//generate random values scaled with epsilon
extern std::mt19937_64 randNumb;
extern std::uniform_real_distribution<double> distEpsilon;

//those are the random numbers for updates?
extern std::mt19937_64 rng1;
extern std::mt19937_64 rng2;
extern std::mt19937_64 rng3;
extern std::uniform_real_distribution<double> dist;


//required for hot start, to create a set of not identity matrices, like generating X, but with bigger spread 
extern double hotEpsilon;


//uniform dist for generation of random variables in hot_start()
extern std::mt19937_64 hotNumb1;
extern std::mt19937_64 hotNumb2;
extern std::mt19937_64 hotNumb3;

extern std::mt19937_64 hotNumb1Extra;
extern std::mt19937_64 hotNumb2Extra;
extern std::mt19937_64 hotNumb3Extra;


extern std::mt19937_64 hotNumbT1;
extern std::mt19937_64 hotNumbT2;
extern std::mt19937_64 hotNumbT3;


extern std::mt19937_64 hotNumbS1;
extern std::mt19937_64 hotNumbS2;
extern std::mt19937_64 hotNumbS3;


extern std::mt19937_64 hotNumbR1;
extern std::mt19937_64 hotNumbR2;
extern std::mt19937_64 hotNumbR3;

extern std::uniform_real_distribution<double> hotDist;



//uniform dist for selecting indece of XSet
extern std::mt19937_64 indexing;
extern std::uniform_int_distribution<int> indexDist;

//uniform for accept/reject of new config
extern std::mt19937_64 acceptReject;
extern std::uniform_real_distribution<double> uniformAcceptReject;

//uniform dist for selecting which reflection to use in overrelaxation
extern std::mt19937_64 reflection;
extern std::uniform_int_distribution<int> reflectDist;


//Pauli matrices and 2x2 identity
extern std::vector<Matrix<2,2>> pauliMatrices;
extern Matrix<2,2> identity2, Pauli1, Pauli2, Pauli3;

//rsU x cSU identity
extern Matrix<rSU,cSU> identityMatrix;


//h5-filename
extern std::string filenameh5;




