#include <stdlib.h>
#include <vector>

size_t xAxis; 
size_t yAxis;
size_t zAxis;
size_t tAxis;


// row and column of SU(N)
const int rSU= 3;
const int cSU=3;




// storage for a SU(3) matrix
template <size_t row, size_t column>
struct Matrix{
    double data[row][column];


    constexpr double& operator()(size_t r, size_t c) noexcept{
        return data[r][c];
    }
    // maybe add easier access and altering options
};





// this mimics the behaviour of a 4D lattice from our 1D array
double idx(size_t x, size_t y, size_t z, size_t t){
    return x+ xAxis*(y+ yAxis*(z+zAxis*t));
}


//sets all matrices to identity
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice){

    Matrix<rSU,cSU> identityMatrix;

    for(int i = 0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i==j){
                identityMatrix(i,j)=1.0;
            }
            else{
                identityMatrix(i,j)=0.0;
            }
        }
    }


    for(int i = 0; i<tAxis; i++){
        for(int j= 0; j<zAxis; j++){
            for(int k=0; k<yAxis; k++){
                for(int l=0; l<xAxis; l++){
                    lattice[idx(l,k,j,i)]=identityMatrix;


                }


            }
        }
    }

    return;
}







int main(){


    // assigned them the correct value by reading it from yaml file
    xAxis = 5; 
    yAxis = 5;
    zAxis = 5;
    tAxis = 5;

    // if true we have a cold start else a hot start, read from yaml
    bool coldOrHot = true;

    // our lattice as 1D array of matrices (3x3)
    std::vector<Matrix<rSU,rSU>> lattice(xAxis*yAxis*zAxis*tAxis);

    if(coldOrHot == true){
        cold_start(lattice);
    }










    


}