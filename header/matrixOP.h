#pragma once
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"

//generates identity and pauli matrices (base of complex 2x2 matrices)
void generate_Pauli();

//A*B A_ik B_kj = C_ij
Matrix<rSU,cSU> matrix_multiplication(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B);

//A+B
Matrix<rSU,cSU> matrix_addition(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B);

// A-B
Matrix<rSU,cSU> matrix_subtraction(const Matrix<rSU,rSU>& A, const Matrix<rSU,rSU>& B);


// returns the complex conjugate of the matrix -> but not the hermitian conjugate, i.e. conjugate & transpose
Matrix<rSU,cSU> matrix_conjugate(const Matrix<rSU,rSU>& A);


// returns the trace of the matrix
std::complex<double> matrix_trace(const Matrix<rSU,rSU>& A);


//generates a rSU x cSU identity matrix
void generate_identity();