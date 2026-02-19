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



//100 matrices in total -> half X, half inverse of X
std::vector<Matrix<rSU,rSU>> XSet(2*NSetXMatrices);

double epsilon = 0.4;

//RNG in [-0.5, 0.5]
std::uniform_real_distribution<double> dist(-0.5,0.5);

double hotEpsilon=0.5;

//uniform dist for generation of random variables in hot_start()
std::uniform_real_distribution<double> hotDist(-hotEpsilon, hotEpsilon);



//uniform dist for selecting indece of XSet

//std::uniform_int_distribution<int> indexDist (0, 2*NSetXMatrices - 1);

//uniform for accept/reject of new config

//std::uniform_real_distribution<double> uniformAcceptReject(0, 1);

//uniform for refelct overrelaxation
//std::uniform_int_distribution<int> reflectDist (1,3);


//Pauli matrices and 2x2 identity

Matrix<2,2> identity2, Pauli1,Pauli2, Pauli3;
std::vector<Matrix<2,2>> pauliMatrices;
//rsU x cSU identity
Matrix<rSU,cSU> identityMatrix;

Matrix<rSU,cSU> zeroMatrix;

Matrix<2,2> zeroMatrix2;

double drng1=42;
double dhotNumb1=67;
double dhotNumb2=1;
double dhotNumb3=2;
double dhotNumb1Extra=3;
double dhotNumb2Extra=4;
double dhotNumb3Extra=5;
double dhotNumbS1=20;
//double dhotNumbS2=21;
//double dhotNumbS3=22;
double dhotNumbT1=10;
//double dhotNumbT2=11;
//double dhotNumbT3=12;
double dhotNumbR1=30;
//double dhotNumbR2=31;
//double dhotNumbR3=32;
double dindexing=666;



std::string filenameh5 = "mistake";

