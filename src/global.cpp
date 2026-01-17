#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"


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


std::vector<Matrix<rSU,rSU>> XSet(2*NSetXMatrices);



//std::mt19937_64 randNumb(69);
std::uniform_real_distribution<double> distEpsilon(-epsilon,epsilon);

//std::mt19937_64 rng(42);
std::uniform_real_distribution<double> dist(-0.5,0.5);



//uniform dist for generation of random variables in hot_start()
//std::mt19937_64 hotNumb(67);
std::uniform_real_distribution<double> hotDist(-hotEpsilon, hotEpsilon);



//uniform dist for selecting indece of XSet
//std::mt19937_64 indexing(666);
std::uniform_int_distribution<int> indexDist (0, 2*NSetXMatrices);

//uniform for accept/reject of new config
//std::mt19937_64 acceptReject(1949);
std::uniform_real_distribution<double> uniformAcceptReject(0, 1);



//Pauli matrices and 2x2 identity
std::vector<Matrix<2,2>> pauliMatrices;
Matrix<2,2> identity2, Pauli1,Pauli2, Pauli3;