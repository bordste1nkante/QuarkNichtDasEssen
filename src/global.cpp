#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"


size_t xAxis=0;
size_t yAxis=0;
size_t zAxis=0;
size_t tAxis=0;

//coupling strength
double beta=1;



//100 matrices in total -> half X, half inverse of X
std::vector<Matrix<rSU,rSU>> XSet(2*NSetXMatrices);

double epsilon = 0.4;

//RNG in [-0.5, 0.5]
std::uniform_real_distribution<double> dist(-0.5,0.5);

double hotEpsilon=0.5;

//uniform dist for generation of random variables in hot_start()
std::uniform_real_distribution<double> hotDist(-hotEpsilon, hotEpsilon);



//Pauli matrices and 2x2 identity

Matrix<2,2> identity2, Pauli1,Pauli2, Pauli3;
std::vector<Matrix<2,2>> pauliMatrices;
//rsU x cSU identity
Matrix<rSU,cSU> identityMatrix;

Matrix<rSU,cSU> zeroMatrix;

Matrix<2,2> zeroMatrix2;

double drng1=42;
double dhotNumbS1=20;
double dhotNumbT1=10;
double dhotNumbR1=30;
double dindexing=666;

//equilibrium test
double PTestSizeGlobal = 100;


double changeRateLow = -0.035;

double changeRateHigh = 0.035;

double alpha=0.07;
double target_rate=0.5;
double rateInterval=0.03;


std::string filenameh5 = "mistake";

