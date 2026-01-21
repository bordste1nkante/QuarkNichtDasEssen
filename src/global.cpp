#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"


size_t xAxis=0; //unsigned integer
size_t yAxis=0;
size_t zAxis=0;
size_t tAxis=0;

//coupling strength
double beta=1;

//lattice spacing -> set to 1 in lattice units?
double a=1;


//100 matrices in total -> half X, half inverse of X
std::vector<Matrix<rSU,rSU>> XSet(2*NSetXMatrices);

double epsilon = 0.4;

std::mt19937_64 randNumb(69);
std::uniform_real_distribution<double> distEpsilon(-epsilon,epsilon);

//RNG in [-0.5, 0.5]
std::mt19937_64 rng(42);
std::uniform_real_distribution<double> dist(-0.5,0.5);

double hotEpsilon=0.7;

//uniform dist for generation of random variables in hot_start()
std::mt19937_64 hotNumb(67);
std::uniform_real_distribution<double> hotDist(-hotEpsilon, hotEpsilon);



//uniform dist for selecting indece of XSet
std::mt19937_64 indexing(666);
std::uniform_int_distribution<int> indexDist (0, 2*NSetXMatrices - 1);

//uniform for accept/reject of new config
std::mt19937_64 acceptReject(1949);
std::uniform_real_distribution<double> uniformAcceptReject(0, 1);

//Pauli matrices and 2x2 identity

Matrix<2,2> identity2, Pauli1,Pauli2, Pauli3;
std::vector<Matrix<2,2>> pauliMatrices;
//rsU x cSU identity
Matrix<rSU,cSU> identityMatrix;


