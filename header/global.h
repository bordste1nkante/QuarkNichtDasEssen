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

//required for generation of X. Epsilon affects acceptance rate.
extern double epsilon;

//generate random values scaled with epsilon

//those are the random numbers for updates?

extern std::uniform_real_distribution<double> dist;


//required for hot start, to create a set of not identity matrices, like generating X, but with bigger spread 
extern double hotEpsilon;


extern std::uniform_real_distribution<double> hotDist;



//Pauli matrices and 2x2 identity
extern std::vector<Matrix<2,2>> pauliMatrices;
extern Matrix<2,2> identity2, Pauli1, Pauli2, Pauli3;

//rsU x cSU identity
extern Matrix<rSU,cSU> identityMatrix;

extern Matrix<rSU,cSU> zeroMatrix;

extern Matrix<2,2> zeroMatrix2;


//h5-filename
extern std::string filenameh5;



extern double drng1;
extern double dhotNumb1;
extern double dhotNumb2;
extern double dhotNumb3;
extern double dhotNumb1Extra;
extern double dhotNumb2Extra;
extern double dhotNumb3Extra;
extern double dhotNumbS1;
extern double dhotNumbT1;
extern double dhotNumbR1;
extern double dindexing;


