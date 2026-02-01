#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <execution>
#include <thread>
#include <algorithm>
#include <Eigen/Dense>
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"
#include "../header/utils.h"

// this mimics the behaviour of a 4D lattice from our 1D array
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu){
    return 4*(x+ xAxis*(y+ yAxis*(z+zAxis*t)))+mu;
}
std::tuple <size_t, size_t, size_t, size_t, size_t> ReIdx(size_t idx){

    size_t mu = idx%linksPerSite;
    idx/=linksPerSite;
    size_t x = idx%xAxis;
    idx/=xAxis;
    size_t y = idx%yAxis;
    idx/=yAxis;
    size_t z= idx%zAxis;
    idx /= zAxis;
    size_t t = idx%tAxis;

    return {mu,x,y,z,t};
}
void X_updateSU3(){
    for(int p=0; p<NSetXMatrices; p++){
        //final matrices
        Matrix<rSU,cSU> X;
        Matrix<rSU,cSU> invX;


        //X=RST
        Matrix<rSU,cSU> R;
        Matrix<rSU,cSU> S;
        Matrix<rSU,cSU> T;

        Matrix<2,2> Rsmall;
        Matrix<2,2> Ssmall;
        Matrix<2,2> Tsmall;

        Matrix<rSU,cSU> ST;

        //std::cout << "Matrix declared" << std::endl;



        //number required to generate 3 SU(2) matrices, from these we form a SU(3)

        double r[3];
        double r0 = dist(rng1); 
        double s[3];
        double s0 = dist(rng2); 
        double t[3];
        double t0 = dist(rng3);

        r[0]= dist(hotNumbR1);
        r[1]= dist(hotNumbR2);
        r[2]= dist(hotNumbR3);  

        double er = dist(hotNumb1Extra);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);


        s[0]= dist(hotNumbS1);
        s[1]= dist(hotNumbS2);
        s[2]= dist(hotNumbS3);                     

        double es = dist(hotNumb2Extra);    

        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);



        t[0]= dist(hotNumbT1);
        t[1]= dist(hotNumbT2);
        t[2]= dist(hotNumbT3);                    

        double et = dist(hotNumb3Extra);
        double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        t0=t0/std::sqrt(t0*t0)*std::sqrt(1-et*et);
  
        for(int i=0; i<3; i++){
            s[i]=es*s[i]/sLength;
            t[i]=et*t[i]/tLength;
            r[i]=er*r[i]/rLength;

        }
        //std::cout << "vectors done" << std::endl;

        for(int i=0; i<4; i++){

            if(i==0){
                for(int j=0; j<2; j++){
                    for(int k=0; k<2; k++){
                        

                        Ssmall(j,k)= pauliMatrices[i](j,k)*s0;
                        Rsmall(j,k)= pauliMatrices[i](j,k)*r0;
                        Tsmall(j,k)= pauliMatrices[i](j,k)*t0;

                    }
                } 
            }
            else{
            for(int j=0; j<2; j++){
                for(int k=0; k<2; k++){
                    Ssmall(j,k)+= pauliMatrices[i](j,k)*s[i-1];
                    Rsmall(j,k)+= pauliMatrices[i](j,k)*r[i-1];
                    Tsmall(j,k)+= pauliMatrices[i](j,k)*t[i-1];
                    }
                }
            }
        }       

        //std::cout << "2D - matrix done" << std::endl;
        //fill R,S,T
        std::complex<double> uno= {1.0,0.0};
        std::complex<double> zero = {0.0,0.0};

        R(0,0)=Rsmall(0,0);
        R(1,0)=Rsmall(1,0);
        R(2,0)=zero;
        R(0,1)=Rsmall(0,1);
        R(1,1)=Rsmall(1,1);
        R(2,1)=zero;
        R(0,2)=zero;
        R(1,2)=zero;
        R(2,2)=uno;

        S(0,0)=Ssmall(0,0);
        S(1,0)=zero;
        S(2,0)=Ssmall(0,1);
        S(0,1)=zero;
        S(1,1)=uno;
        S(2,1)=Ssmall(1,0);
        S(0,2)=zero;
        S(1,2)=zero;
        S(2,2)=Ssmall(1,1);

        T(0,0)=uno;
        T(1,0)=zero;
        T(2,0)=zero;
        T(0,1)=zero;
        T(1,1)=Tsmall(0,0);
        T(2,1)=Tsmall(1,0);
        T(0,2)=zero;
        T(1,2)=Tsmall(0,1);
        T(2,2)=Tsmall(1,1);

        //matrix multiplication S*T = ST
        ST=matrix_multiplication(S,T);
        //matrix multiplication R*ST = X
        X=matrix_multiplication(R,ST);

        ////inverting X
        //std::complex<double> detX = X(0,0)*(X(1,1)*X(2,2)-X(1,2)*X(2,1))-X(0,1)*(X(1,0)*X(2,2)-X(1,2)*X(2,0))+X(0,2)*(X(1,0)*X(2,1)-X(1,1)*X(2,0));
        std::complex<double> detX = det_A(X);

        //std::cout << detX << std::endl;

        // should insure that if det(X) \approx zero we just repeat the process and generate a new X

        double tol = 1e-11;
            
        if (detX.real() > -tol && detX.real() < tol &&
            detX.imag() > -tol && detX.imag() < tol)
        {
            std::cout << "Stuck in loop" << std::endl;
            p = p - 1;
        }
        else{
            normalizeSU3Matrix(X);

            invX = inverse_A(X);

            XSet[2*p]= X;
            XSet[2*p+1]=invX;


        }


    }


}

//parallel
void X_updateSU3(size_t input){


    //for(int p=0; p<NSetXMatrices; p++){
    //    //final matrices
    //    Matrix<rSU,cSU> X;
    //    Matrix<rSU,cSU> invX;
//
//
    //    //X=RST
    //    Matrix<rSU,cSU> R;
    //    Matrix<rSU,cSU> S;
    //    Matrix<rSU,cSU> T;
//
    //    Matrix<2,2> Rsmall;
    //    Matrix<2,2> Ssmall;
    //    Matrix<2,2> Tsmall;
//
    //    Matrix<rSU,cSU> ST;
//
    //    //std::cout << "Matrix declared" << std::endl;
//
//
//
    //    //number required to generate 3 SU(2) matrices, from these we form a SU(3)
//
    //    double r[3];
    //    double r0 = dist(rng1); 
    //    double s[3];
    //    double s0 = dist(rng2); 
    //    double t[3];
    //    double t0 = dist(rng3);
//
    //    r[0]= dist(hotNumbR1);
    //    r[1]= dist(hotNumbR2);
    //    r[2]= dist(hotNumbR3);  
//
    //    double er = dist(hotNumb1Extra);
//
    //    double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    //    r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);
//
//
    //    s[0]= dist(hotNumbS1);
    //    s[1]= dist(hotNumbS2);
    //    s[2]= dist(hotNumbS3);                     
//
    //    double es = dist(hotNumb2Extra);    
//
    //    double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
    //    s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);
//
//
//
    //    t[0]= dist(hotNumbT1);
    //    t[1]= dist(hotNumbT2);
    //    t[2]= dist(hotNumbT3);                    
//
    //    double et = dist(hotNumb3Extra);
    //    double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    //    t0=t0/std::sqrt(t0*t0)*std::sqrt(1-et*et);
  //
    //    for(int i=0; i<3; i++){
    //        s[i]=es*s[i]/sLength;
    //        t[i]=et*t[i]/tLength;
    //        r[i]=er*r[i]/rLength;
//
    //    }
    //    //std::cout << "vectors done" << std::endl;
//
    //    for(int i=0; i<4; i++){
//
    //        if(i==0){
    //            for(int j=0; j<2; j++){
    //                for(int k=0; k<2; k++){
    //                    
//
    //                    Ssmall(j,k)= pauliMatrices[i](j,k)*s0;
    //                    Rsmall(j,k)= pauliMatrices[i](j,k)*r0;
    //                    Tsmall(j,k)= pauliMatrices[i](j,k)*t0;
//
    //                }
    //            } 
    //        }
    //        else{
    //        for(int j=0; j<2; j++){
    //            for(int k=0; k<2; k++){
    //                Ssmall(j,k)+= pauliMatrices[i](j,k)*s[i-1];
    //                Rsmall(j,k)+= pauliMatrices[i](j,k)*r[i-1];
    //                Tsmall(j,k)+= pauliMatrices[i](j,k)*t[i-1];
    //                }
    //            }
    //        }
    //    }       
//
    //    //std::cout << "2D - matrix done" << std::endl;
    //    //fill R,S,T
    //    std::complex<double> uno= {1.0,0.0};
    //    std::complex<double> zero = {0.0,0.0};
//
    //    R(0,0)=Rsmall(0,0);
    //    R(1,0)=Rsmall(1,0);
    //    R(2,0)=zero;
    //    R(0,1)=Rsmall(0,1);
    //    R(1,1)=Rsmall(1,1);
    //    R(2,1)=zero;
    //    R(0,2)=zero;
    //    R(1,2)=zero;
    //    R(2,2)=uno;
//
    //    S(0,0)=Ssmall(0,0);
    //    S(1,0)=zero;
    //    S(2,0)=Ssmall(0,1);
    //    S(0,1)=zero;
    //    S(1,1)=uno;
    //    S(2,1)=Ssmall(1,0);
    //    S(0,2)=zero;
    //    S(1,2)=zero;
    //    S(2,2)=Ssmall(1,1);
//
    //    T(0,0)=uno;
    //    T(1,0)=zero;
    //    T(2,0)=zero;
    //    T(0,1)=zero;
    //    T(1,1)=Tsmall(0,0);
    //    T(2,1)=Tsmall(1,0);
    //    T(0,2)=zero;
    //    T(1,2)=Tsmall(0,1);
    //    T(2,2)=Tsmall(1,1);
//
    //    //matrix multiplication S*T = ST
    //    for(int i=0; i<rSU; i++){
    //        for(int j=0; j<cSU; j++){
    //            std::complex<double> sum;
    //            for(int k = 0; k<rSU; k++){
    //                sum += S(i,k)*T(k,j);
//
    //            }
    //            ST(i,j)= sum;
    //        }
//
    //    }
    //    //matrix multiplication R*ST = X
    //    for(int i=0; i<rSU; i++){
    //        for(int j=0; j<cSU; j++){
    //            std::complex<double> sum;
    //            for(int k = 0; k<rSU; k++){
    //                sum += R(i,k)*ST(k,j);
//
    //            }
    //            X(i,j)= sum;
    //        }
//
    //    }
//
    //    ////inverting X
    //    //std::complex<double> detX = X(0,0)*(X(1,1)*X(2,2)-X(1,2)*X(2,1))-X(0,1)*(X(1,0)*X(2,2)-X(1,2)*X(2,0))+X(0,2)*(X(1,0)*X(2,1)-X(1,1)*X(2,0));
    //    std::complex<double> detX = det_A(X);
//
    //    //std::cout << detX << std::endl;
//
    //    // should insure that if det(X) \approx zero we just repeat the process and generate a new X
//
    //    double tol = 1e-11;
    //        
    //    if (detX.real() > -tol && detX.real() < tol &&
    //        detX.imag() > -tol && detX.imag() < tol)
    //    {
    //        std::cout << "Stuck in loop" << std::endl;
    //        p = p - 1;
    //    }
    //    else{
    //        //normalizes det to 1
    //        // not sure if this needs to be done by dividing over third root //Vincent: ?
    //        for(int i=0; i<rSU; i++){
    //            for(int j=0; j<rSU; j++){
    //                X(i,j)= X(i,j)/pow(detX, 1.0/double(rSU));
    //            }
    //        }
    //        //detX = det_A(X);
    //        //std::cout << detX << std::endl;
//
//
//  //          //the X are already normalized to be detX=1, therefore it becomes redundant to divide by detX
//  //          invX(0,0)= (X(1,1)*X(2,2)-X(1,2)*X(2,1));
//  //          invX(0,1)=-(X(0,1)*X(2,2)-X(0,2)*X(2,1));
//  //          invX(0,2)= (X(0,1)*X(1,2)-X(0,2)*X(1,1));
////
//  //          invX(1,0)=-(X(1,0)*X(2,2)-X(1,2)*X(2,0));
//  //          invX(1,1)=(X(0,0)*X(2,2)-X(0,2)*X(2,0));
//  //          invX(1,2)=-(X(0,0)*X(1,2)-X(0,2)*X(1,0));
////
//  //          invX(2,0)=(X(1,0)*X(2,1)-X(1,1)*X(2,0));
//  //          invX(2,1)=-(X(0,0)*X(2,1)-X(0,1)*X(2,0));
//  //          invX(2,2)=(X(0,0)*X(1,1)-X(0,1)*X(1,0));
////
//
    //        invX = inverse_A(X);
    //        //std::cout << "invertiert" << std::endl;
    //        //save X and invX in our set of matrices
    //        //std::cout << X(0,0) << std::endl;
    //        XSet[2*p]= X;
    //        XSet[2*p+1]=invX;
//
    //        //std::cout << "X generated" << std::endl;
//
    //        //std::cout << "X abgespeichert" << std::endl;
//
    //    }
//
//
    //}

    std::vector<Matrix<rSU,rSU>> bufferSet(XSet.size());
    std::vector<size_t> indices(NSetXMatrices);
    std::iota(indices.begin(), indices.end(),0);
    std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t index){
        Matrix<rSU,cSU> X;
        Matrix<rSU,cSU> invX;

        //X=RST
        Matrix<rSU,cSU> R;
        Matrix<rSU,cSU> S;
        Matrix<rSU,cSU> T;

        Matrix<2,2> Rsmall;
        Matrix<2,2> Ssmall;
        Matrix<2,2> Tsmall;

        Matrix<rSU,cSU> ST;

        std::mt19937_64 hottestNumbR1(index*dhotNumbR1+input);
        std::mt19937_64 hottestNumbR2(index*dhotNumbR2+input);
        std::mt19937_64 hottestNumbR3(index*dhotNumbR3+input);
        std::mt19937_64 hottestNumbR0(index*drng1+input);
        std::mt19937_64 hottestNumbRE(index*dhotNumb1Extra+input);

        std::mt19937_64 hottestNumbS1(index*dhotNumbS1+input);
        std::mt19937_64 hottestNumbS2(index*dhotNumbS2+input);
        std::mt19937_64 hottestNumbS3(index*dhotNumbS3+input);
        std::mt19937_64 hottestNumbS0(index*drng2+input);
        std::mt19937_64 hottestNumbSE(index*dhotNumb2Extra+input);

        std::mt19937_64 hottestNumbT1(index*dhotNumbT1+input);
        std::mt19937_64 hottestNumbT2(index*dhotNumbT2+input);
        std::mt19937_64 hottestNumbT3(index*dhotNumbT3+input);
        std::mt19937_64 hottestNumbT0(index*drng3+input);
        std::mt19937_64 hottestNumbTE(index*dhotNumb3Extra+input);



        double r[3];
        double r0 = dist(hottestNumbR0); 
        double s[3];
        double s0 = dist(hottestNumbS0); 
        double t[3];
        double t0 = dist(hottestNumbT0);

        r[0]= dist(hottestNumbR1);
        r[1]= dist(hottestNumbR2);
        r[2]= dist(hottestNumbR3);  

        double er = dist(hottestNumbSE);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);


        s[0]= dist(hottestNumbS1);
        s[1]= dist(hottestNumbS2);
        s[2]= dist(hottestNumbS3);                     

        double es = dist(hottestNumbSE);    

        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);


        t[0]= dist(hottestNumbT1);
        t[1]= dist(hottestNumbT2);
        t[2]= dist(hottestNumbT3);                    

        double et = dist(hottestNumbRE);
        double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        t0=t0/std::sqrt(t0*t0)*std::sqrt(1-et*et);
        for(int i=0; i<3; i++){
            s[i]=es*s[i]/sLength;
            t[i]=et*t[i]/tLength;
            r[i]=er*r[i]/rLength;

        }

        for(int i=0; i<4; i++){

            if(i==0){
                for(int j=0; j<2; j++){
                    for(int k=0; k<2; k++){
                        

                        Ssmall(j,k)= pauliMatrices[i](j,k)*s0;
                        Rsmall(j,k)= pauliMatrices[i](j,k)*r0;
                        Tsmall(j,k)= pauliMatrices[i](j,k)*t0;

                    }
                } 
            }
            else{
            for(int j=0; j<2; j++){
                for(int k=0; k<2; k++){
                    Ssmall(j,k)+= pauliMatrices[i](j,k)*s[i-1];
                    Rsmall(j,k)+= pauliMatrices[i](j,k)*r[i-1];
                    Tsmall(j,k)+= pauliMatrices[i](j,k)*t[i-1];
                    }
                }
            }
        }       

        //std::cout << "2D - matrix done" << std::endl;
        //fill R,S,T
        std::complex<double> uno= {1.0,0.0};
        std::complex<double> zero = {0.0,0.0};

        R(0,0)=Rsmall(0,0);
        R(1,0)=Rsmall(1,0);
        R(2,0)=zero;
        R(0,1)=Rsmall(0,1);
        R(1,1)=Rsmall(1,1);
        R(2,1)=zero;
        R(0,2)=zero;
        R(1,2)=zero;
        R(2,2)=uno;

        S(0,0)=Ssmall(0,0);
        S(1,0)=zero;
        S(2,0)=Ssmall(0,1);
        S(0,1)=zero;
        S(1,1)=uno;
        S(2,1)=Ssmall(1,0);
        S(0,2)=zero;
        S(1,2)=zero;
        S(2,2)=Ssmall(1,1);

        T(0,0)=uno;
        T(1,0)=zero;
        T(2,0)=zero;
        T(0,1)=zero;
        T(1,1)=Tsmall(0,0);
        T(2,1)=Tsmall(1,0);
        T(0,2)=zero;
        T(1,2)=Tsmall(0,1);
        T(2,2)=Tsmall(1,1);

        //matrix multiplication S*T = ST
        ST =matrix_multiplication(S,T);
        //matrix multiplication R*ST = X
        X=matrix_multiplication(R,ST);
        normalizeSU3Matrix(X);
            

        invX = inverse_A(X);

        bufferSet[2*index]= X;
        bufferSet[2*index+1]=invX;



        });
    std::swap(bufferSet,XSet);


}

//rounding errors causes matrices to potentially digress from det = 1, we correct that from time to time.
void normalizeSU3(std::vector<Matrix<rSU,rSU>>& lattice){

    Matrix<rSU,cSU> U;

    for(int i = 0; i<tAxis; i++){
        for(int j= 0; j<zAxis; j++){
            for(int k=0; k<yAxis; k++){
                for(int l=0; l<xAxis; l++){
                    for(int mu = 0; mu<linksPerSite; mu++){
                        U = lattice[idx(l,k,j,i, mu)];
                        normalizeSU3Matrix(U);
                        lattice[idx(l,k,j,i, mu)]= U;
                }
                }


            }
        }
    }



}


void normalizeSU3Matrix(Matrix<rSU,rSU>& U){


    std::complex<double> sum, vuProduct, vsum;

    std::vector<std::complex<double>> vPrime;
    for(int i = 0; i<rSU; i++){
        sum += std::norm(U(0,i));
}
    sum= std::sqrt(sum);
    for(int i = 0; i<rSU; i++){
    U(0,i)= U(0,i)/sum;
}
    for(int i = 0; i<rSU; i++){
    vuProduct += U(1,i)*std::conj(U(0,i));
}
    for(int i = 0; i<rSU; i++){
    vPrime.push_back(U(1,i)-U(0,i)*(vuProduct));
}
    for(int i = 0; i<rSU; i++){
        vsum += std::norm(vPrime[i]);
}
    vsum= std::sqrt(vsum);

    for(int i = 0; i<rSU; i++){
    U(1,i)= vPrime[i]/vsum;
}
    

    U(2,0)= std::conj(U(0,1))*std::conj(U(1,2))-std::conj(U(0,2))*std::conj(U(1,1));
    U(2,1)= std::conj(U(0,2))*std::conj(U(1,0))-std::conj(U(0,0))*std::conj(U(1,2));
    U(2,2)= std::conj(U(0,0))*std::conj(U(1,1))-std::conj(U(0,1))*std::conj(U(1,0));


    std::complex<double> detU = det_A(U);
    for(int m = 0; m< rSU; m++){
        for(int n=0; n<cSU; n++){
            U(m,n)= U(m,n)/pow(detU, 1.0/double(rSU));
                        }
                    }
                }

    



//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice){

//    for(int i = 0; i<tAxis; i++){
//        for(int j= 0; j<zAxis; j++){
//            for(int k=0; k<yAxis; k++){
//                for(int l=0; l<xAxis; l++){
//                    for(int mu = 0; mu<linksPerSite; mu++){
//                        lattice[idx(l,k,j,i,mu)]=identityMatrix;
//
//                }
//                }
//
//
//
//            }
//        }
//    }
    std::vector<Matrix<rSU,rSU>> bufferLattice(lattice.size());
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);

    std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
        bufferLattice[i]=identityMatrix;
    });
    std::swap(bufferLattice,lattice);

    return;
}

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice){
    //Matrix<rSU,cSU> U;
    //Matrix<rSU,cSU> R;
    //Matrix<rSU,cSU> S;
    //Matrix<rSU,cSU> T;
    //Matrix<2,2> Rsmall;
    //Matrix<2,2> Ssmall;
    //Matrix<2,2> Tsmall;
    //Matrix<rSU,cSU> ST;
//
//
    ////fixed int c and d
    //for(int a = 0; a<tAxis; a++){
    //    for(int b= 0; b<zAxis; b++){
    //        for(int c=0; c<yAxis; c++){
    //            for(int d=0; d<xAxis; d++){
    //                for(int mu = 0; mu<linksPerSite; mu++){
    //                //number required to generate 3 SU(2) matrices, from these we form a SU(3)
//
    //                double r[3];
    //                r[0]= dist(hotNumbR1);
    //                r[1]= dist(hotNumbR2);
    //                r[2]= dist(hotNumbR3); 
    //                double r0 = dist(hotNumb1); 
//
//
    //                double er = hotDist(hotNumb1Extra);
//
    //                double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    //                r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);
//
    //                double s[3];
    //                s[0]= dist(hotNumbS1);
    //                s[1]= dist(hotNumbS2);
    //                s[2]= dist(hotNumbS3);                     
    //                double s0 = dist(hotNumb2); 
//
    //                double es = hotDist(hotNumb2Extra);     
//
    //                double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
    //                s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);
//
//
    //                double t[3];
    //                t[0]= dist(hotNumbT1);
    //                t[1]= dist(hotNumbT2);
    //                t[2]= dist(hotNumbT3);                    
//
    //                double t0 = dist(hotNumb3);
//
    //                double et = hotDist(hotNumb3Extra);
//
    //                double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    //                t0=t0/std::sqrt(t0*t0)*std::sqrt(1-et*et);
//
    //                for(int i=0; i<3; i++){
    //                    s[i]=es*s[i]/sLength;
    //                    t[i]=et*t[i]/tLength;
    //                    r[i]=er*r[i]/rLength;
    //                
    //                }
    //            
    //                for(int i=0; i<4; i++){
    //                
    //                    if(i==0){
    //                        for(int j=0; j<2; j++){
    //                            for(int k=0; k<2; k++){
    //                            
    //                                Ssmall(j,k)= pauliMatrices[i](j,k)*s0;
    //                                Rsmall(j,k)= pauliMatrices[i](j,k)*r0;
    //                                Tsmall(j,k)= pauliMatrices[i](j,k)*t0;
    //                            
    //                            }
    //                        } 
    //                    }
    //                    else{
    //                    for(int j=0; j<2; j++){
    //                        for(int k=0; k<2; k++){
    //                            Ssmall(j,k)+= pauliMatrices[i](j,k)*s[i-1];
    //                            Rsmall(j,k)+= pauliMatrices[i](j,k)*r[i-1];
    //                            Tsmall(j,k)+= pauliMatrices[i](j,k)*t[i-1];
    //                    }
    //                    }  
    //                }     
    //                }               
//
    //                    //fill R,S,T
    //                    std::complex<double> uno= {1.0,0.0};
    //                    std::complex<double> zero = {0.0,0.0};
//
    //                    R(0,0)=Rsmall(0,0);
    //                    R(1,0)=Rsmall(1,0);
    //                    R(2,0)=zero;
    //                    R(0,1)=Rsmall(0,1);
    //                    R(1,1)=Rsmall(1,1);
    //                    R(2,1)=zero;
    //                    R(0,2)=zero;
    //                    R(1,2)=zero;
    //                    R(2,2)=uno;
//
    //                    S(0,0)=Ssmall(0,0);
    //                    S(1,0)=zero;
    //                    S(2,0)=Ssmall(0,1);
    //                    S(0,1)=zero;
    //                    S(1,1)=uno;
    //                    S(2,1)=Ssmall(1,0);
    //                    S(0,2)=zero;
    //                    S(1,2)=zero;
    //                    S(2,2)=Ssmall(1,1);
//
    //                    T(0,0)=uno;
    //                    T(1,0)=zero;
    //                    T(2,0)=zero;
    //                    T(0,1)=zero;
    //                    T(1,1)=Tsmall(0,0);
    //                    T(2,1)=Tsmall(1,0);
    //                    T(0,2)=zero;
    //                    T(1,2)=Tsmall(0,1);
    //                    T(2,2)=Tsmall(1,1);
//
    //                    //matrix multiplication S*T = ST
    //                    for(int i=0; i<rSU; i++){
    //                        for(int j=0; j<cSU; j++){
    //                            std::complex<double> sum;
    //                            for(int k = 0; k<rSU; k++){
    //                                sum += S(i,k)*T(k,j);
    //                            
    //                            }
    //                            ST(i,j)= sum;
    //                        }
    //                    
    //                    }
    //                    //matrix multiplication R*ST = X
    //                    for(int i=0; i<rSU; i++){
    //                        for(int j=0; j<cSU; j++){
    //                            std::complex<double> sum;
    //                            for(int k = 0; k<rSU; k++){
    //                                sum += R(i,k)*ST(k,j);
    //                            
    //                            }
    //                            U(i,j)= sum;
    //                        }
    //                    
    //                    }
    //                    //normalize to det(U)=1
    //                    std::complex<double> detU = det_A(U);
    //                    
    //                    for(int i = 0; i<rSU; i++){
    //                        for(int j = 0; j<cSU; j++){
    //                            U(i,j)=U(i,j)/pow(detU, 1.0/double(rSU));
    //                        }
    //                    }
    //                
    //                    lattice[idx(a,b,c,d, mu)]=U;
    //            }
    //                
    //            }
//
//
    //        }
    //    }
    //}
    std::vector<Matrix<rSU,rSU>> bufferLattice(lattice.size());
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);


    std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){

            
  
        Matrix<rSU,cSU> U;
        Matrix<rSU,cSU> R;
        Matrix<rSU,cSU> S;
        Matrix<rSU,cSU> T;
        Matrix<2,2> Rsmall;
        Matrix<2,2> Ssmall;
        Matrix<2,2> Tsmall;
        Matrix<rSU,cSU> ST;
        
        std::mt19937_64 hottestNumbR1(i*dhotNumbR1);
        std::mt19937_64 hottestNumbR2(i*dhotNumbR2);
        std::mt19937_64 hottestNumbR3(i*dhotNumbR3);
        std::mt19937_64 hottestNumbR0(i*dhotNumb1);
        std::mt19937_64 hottestNumbRE(i*dhotNumb1Extra);

        std::mt19937_64 hottestNumbS1(i*dhotNumbS1);
        std::mt19937_64 hottestNumbS2(i*dhotNumbS2);
        std::mt19937_64 hottestNumbS3(i*dhotNumbS3);
        std::mt19937_64 hottestNumbS0(i*dhotNumb2);
        std::mt19937_64 hottestNumbSE(i*dhotNumb2Extra);

        std::mt19937_64 hottestNumbT1(i*dhotNumbT1);
        std::mt19937_64 hottestNumbT2(i*dhotNumbT2);
        std::mt19937_64 hottestNumbT3(i*dhotNumbT3);
        std::mt19937_64 hottestNumbT0(i*dhotNumb3);
        std::mt19937_64 hottestNumbTE(i*dhotNumb3Extra);




        std::uniform_real_distribution<double> distribution(-0.5,0.5);
        std::uniform_real_distribution<double> hotDistribution(-hotEpsilon,hotEpsilon);



        //number required to generate 3 SU(2) matrices, from these we form a SU(3)
        double r[3];
        r[0]= distribution(hottestNumbR1);
        r[1]= distribution(hottestNumbR2);
        r[2]= distribution(hottestNumbR3); 

        double r0 = distribution(hottestNumbR0); 
        double er = hotDistribution(hottestNumbRE);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);

        double s[3];
        s[0]= distribution(hottestNumbS1);
        s[1]= distribution(hottestNumbS2);
        s[2]= distribution(hottestNumbS3);        

        double s0 = distribution(hottestNumbS0); 
        double es = hotDistribution(hottestNumbSE);  

        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);

        double t[3];
        t[0]= distribution(hottestNumbT1);
        t[1]= distribution(hottestNumbT2);
        t[2]= distribution(hottestNumbT3);   

        double t0 = distribution(hottestNumbT0);

        double et = hotDistribution(hottestNumbTE);

        double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        t0=t0/std::sqrt(t0*t0)*std::sqrt(1-et*et);

        for(int i=0; i<3; i++){
            s[i]=es*s[i]/sLength;
            t[i]=et*t[i]/tLength;
            r[i]=er*r[i]/rLength;
                    
                }
                
        for(int i=0; i<4; i++){
                    
            if(i==0){
                for(int j=0; j<2; j++){
                    for(int k=0; k<2; k++){
                                
                        Ssmall(j,k)= pauliMatrices[i](j,k)*s0;
                        Rsmall(j,k)= pauliMatrices[i](j,k)*r0;
                        Tsmall(j,k)= pauliMatrices[i](j,k)*t0;
                                
                                }
                            } 
                        }
            else{
                for(int j=0; j<2; j++){
                    for(int k=0; k<2; k++){
                        Ssmall(j,k)+= pauliMatrices[i](j,k)*s[i-1];
                        Rsmall(j,k)+= pauliMatrices[i](j,k)*r[i-1];
                        Tsmall(j,k)+= pauliMatrices[i](j,k)*t[i-1];
                        }
                        }  
                    }     
                    }               

        //fill R,S,T
        std::complex<double> uno= {1.0,0.0};
        std::complex<double> zero = {0.0,0.0};

        R(0,0)=Rsmall(0,0);
        R(1,0)=Rsmall(1,0);
        R(2,0)=zero;
        R(0,1)=Rsmall(0,1);
        R(1,1)=Rsmall(1,1);
        R(2,1)=zero;
        R(0,2)=zero;
        R(1,2)=zero;
        R(2,2)=uno;

        S(0,0)=Ssmall(0,0);
        S(1,0)=zero;
        S(2,0)=Ssmall(0,1);
        S(0,1)=zero;
        S(1,1)=uno;
        S(2,1)=Ssmall(1,0);
        S(0,2)=zero;
        S(1,2)=zero;
        S(2,2)=Ssmall(1,1);

        T(0,0)=uno;
        T(1,0)=zero;
        T(2,0)=zero;
        T(0,1)=zero;
        T(1,1)=Tsmall(0,0);
        T(2,1)=Tsmall(1,0);
        T(0,2)=zero;
        T(1,2)=Tsmall(0,1);
        T(2,2)=Tsmall(1,1);

        //matrix multiplication S*T = ST
        for(int i=0; i<rSU; i++){
            for(int j=0; j<cSU; j++){
                std::complex<double> sum;
                for(int k = 0; k<rSU; k++){
                    sum += S(i,k)*T(k,j);
                                
                        }
                    ST(i,j)= sum;
                    }
                        
                    }
        //matrix multiplication R*ST = X
        for(int i=0; i<rSU; i++){
            for(int j=0; j<cSU; j++){
                std::complex<double> sum;
                for(int k = 0; k<rSU; k++){
                    sum += R(i,k)*ST(k,j);
                                
                        }
                        U(i,j)= sum;
                            }
                        
                        }
        //normalize to det(U)=1
        std::complex<double> detU = det_A(U);
                        
        for(int i = 0; i<rSU; i++){
            for(int j = 0; j<cSU; j++){
                U(i,j)=U(i,j)/pow(detU, 1.0/double(rSU));
                            }
                        }
                    
            bufferLattice[i]=U;
                
        
        }
                    
                


            
        
    
    );
    //switch pointers of buffer and original
    std::swap(bufferLattice,lattice);
    return;
}






bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime ,size_t x, size_t y, size_t z, size_t t, size_t mu, const Matrix<rSU,cSU>& A){
    
    bool accept;
    double r;
    double probability;
    
    double SActionDif;
    SActionDif = -beta/(xAxis*yAxis*zAxis*tAxis)*(matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real();
    probability = std::min(1.0, exp(- SActionDif)); // according to my notes -> check in doubt

    r= uniformAcceptReject(acceptReject);
    if(r<=probability){
        accept = true;
    }
    else{
        accept = false;
    }


    return accept;
}

//parallel
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime, const Matrix<rSU,cSU>& A, std::mt19937_64& ActionAccept, std::uniform_real_distribution<double>& Distribution){
    
    bool accept;
    double r;
    double probability;
    
    double SActionDif;
    SActionDif = -beta/(xAxis*yAxis*zAxis*tAxis)*(matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real();
    probability = std::min(1.0, exp(- SActionDif)); // according to my notes -> check in doubt

    r= Distribution(ActionAccept);
    if(r<=probability){
        accept = true;
    }
    else{
        accept = false;
    }


    return accept;
}

Matrix<rSU,cSU> determineA(const std::vector<Matrix<rSU,rSU>>& lattice ,size_t x, size_t y, size_t z, size_t t, size_t mu){
    Matrix<rSU,cSU> ATemp1;
    Matrix<rSU,cSU> ATemp2;
    Matrix<rSU,cSU> ATemp3;
    Matrix<rSU,cSU> A;

    double Sum;
    //enforce periodic boundary condition
    auto bCX = [](size_t i) -> size_t{
        return (i+xAxis) % xAxis;
        };
    auto bCY = [](size_t i) -> size_t{
        return (i+yAxis) % yAxis;
        };
    auto bCZ = [](size_t i) -> size_t{
        return (i+zAxis) % zAxis;
        };
    auto bCT = [](size_t i) -> size_t{
        return (i+tAxis) % tAxis;
        };


    double tP;
    double tM;
    double xP;
    double xM;
    double yP;
    double yM;
    double zP;
    double zM;


    for(int nu=0; nu<linksPerSite; nu++){
        if(nu!=mu){

            if(nu==0){
                //Ich glaube in dieser Schleife ist viel Index-Salat -> Nevermind??
                if(mu==1){

                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xM,yP,z,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                    

                }
                if(mu==2){
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication
                    (matrix_hermitean_conjugate(lattice[idx(xM,y,zP,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication
                    (matrix_hermitean_conjugate(lattice[idx(xM,y,z,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,z,tP,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }
            if(nu==1){

                if(mu==0){

                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,yM,z,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==2){
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yM,zP,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yM,z,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,z,tP,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }
            if(nu==2){


                if(mu==1){

                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yP,zM,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==0){
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,y,zM,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    zP= bCZ(x+1);
                    zM= bCZ(x-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zM,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,z,tP,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }
            if(nu==3){



                if(mu==1){

                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,yP,z,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==2){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,zP,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==0){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(xP,y,z,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }



        }
    }





    return A;
}

Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U){

    int reflect = reflectDist(reflection);

    Matrix<rSU,cSU> UPrime;

    Eigen::Matrix3cd M = translateMatrices(A);
    Eigen::Matrix3cd UEigen = translateMatrices(U);    

    Eigen::Matrix3cd K = M.adjoint()*M;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> es(K);

    Eigen::Matrix3cd D_Sqrt=es.eigenvalues().cwiseSqrt().asDiagonal();
    Eigen::Matrix3cd H = es.eigenvectors()*D_Sqrt*es.eigenvectors().adjoint();

    Eigen::Matrix3cd O=M*H.inverse();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> esH(H);
    Eigen::Matrix3cd V= esH.eigenvectors();
    Eigen::Matrix3cd Vadjoint = V.adjoint();


    Eigen::Matrix3cd Ur = V*UEigen*O*Vadjoint;

    if(reflect = 1){
        Ur(0,1)= -Ur(0,1);
        Ur(1,0)=-Ur(1,0);
        Ur(0,2)=-Ur(0,2);
        Ur(2,0)=-Ur(2,0);

    }
    if(reflect = 2){
        Ur(0,1)= -Ur(0,1);
        Ur(1,0)=-Ur(1,0);
        Ur(1,2)=-Ur(1,2);
        Ur(2,1)=-Ur(2,1);

    }
    if(reflect = 3){
        Ur(0,2)= -Ur(0,2);
        Ur(2,0)=-Ur(2,0);
        Ur(1,2)=-Ur(1,2);
        Ur(2,1)=-Ur(2,1);

    }
    else{
        std::cout << "reflextion failed due to unexpected value" << std::endl;
    }

    Eigen::Matrix3cd UR = Vadjoint*Ur*V*O.adjoint();

    UPrime = retranslateMatrices(UR);

    return UPrime;



}
//parallel
Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U,std::uniform_int_distribution<int>& distribution,std::mt19937_64& random ){
    int reflect = distribution(random);

    Matrix<rSU,cSU> UPrime;

    Eigen::Matrix3cd M = translateMatrices(A);
    Eigen::Matrix3cd UEigen = translateMatrices(U);    

    Eigen::Matrix3cd K = M.adjoint()*M;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> es(K);

    Eigen::Matrix3cd D_Sqrt=es.eigenvalues().cwiseSqrt().asDiagonal();
    Eigen::Matrix3cd H = es.eigenvectors()*D_Sqrt*es.eigenvectors().adjoint();

    Eigen::Matrix3cd O=M*H.inverse();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> esH(H);
    Eigen::Matrix3cd V= esH.eigenvectors();
    Eigen::Matrix3cd Vadjoint = V.adjoint();


    Eigen::Matrix3cd Ur = V*UEigen*O*Vadjoint;

    if(reflect = 1){
        Ur(0,1)= -Ur(0,1);
        Ur(1,0)=-Ur(1,0);
        Ur(0,2)=-Ur(0,2);
        Ur(2,0)=-Ur(2,0);

    }
    if(reflect = 2){
        Ur(0,1)= -Ur(0,1);
        Ur(1,0)=-Ur(1,0);
        Ur(1,2)=-Ur(1,2);
        Ur(2,1)=-Ur(2,1);

    }
    if(reflect = 3){
        Ur(0,2)= -Ur(0,2);
        Ur(2,0)=-Ur(2,0);
        Ur(1,2)=-Ur(1,2);
        Ur(2,1)=-Ur(2,1);

    }
    else{
        std::cout << "reflextion failed due to unexpected value" << std::endl;
    }

    Eigen::Matrix3cd UR = Vadjoint*Ur*V*O.adjoint();

    UPrime = retranslateMatrices(UR);

    return UPrime;
}
void wilsonLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint ){
    Matrix<rSU, cSU> Path1;
    Matrix<rSU, cSU> Path2;

    int dx=(endPoint[0]>startingPoint[0]) ? +1:-1;
    int dy=(endPoint[1]>startingPoint[1]) ? +1:-1;                              
    int dz=(endPoint[2]>startingPoint[2]) ? +1:-1;

    int x= startingPoint[0];
    int y= startingPoint[1];
    int z= startingPoint[2];

    int xDistance=0;
    int yDistance=0;
    int zDistance=0;



    while (x!=endPoint[0])
    {

        Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,0)]);
        Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,0)]);
        loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
        x+=dx;
        xDistance+=1;
        r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));

    }
    while (y!=endPoint[1])
    {
        Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,1)]);
        Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,1)]);
        loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
        y+=dy;
        yDistance+=1;
        r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));                                

    }
    while (z!=endPoint[2])
    {
        Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,0,2)]);
        Path2=matrix_multiplication(Path2,lattice[idx(x,y,z,tAxis-1,2)]);
        loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));
        z+=dz;
        zDistance+=1;
        r.push_back(std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance));  
    }
                            



    }


void polyakovLoop(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<std::complex<double>>& loops,
                            std::vector<double>& r,   const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint ){

    Matrix<rSU, cSU> Path1;
    Matrix<rSU, cSU> Path2;


    int dt=(endPoint[3]>startingPoint[3]) ? +1:-1;
    int x= startingPoint[0];
    int y= startingPoint[1];
    int z= startingPoint[2];
    int t= startingPoint[3];

    int tDistance= 0;

    while(t!=endPoint[3]){
        Path1=matrix_multiplication(Path1, lattice[idx(x,y,z,t,3)]);
        Path2=matrix_multiplication(Path2, lattice[idx(x,y,z,endPoint[3],3)]);
        loops.push_back(matrix_trace(matrix_multiplication(Path2,matrix_hermitean_conjugate(Path1))));

        //save t-distance in r array (for simplicity)
        r.push_back(std::sqrt(tDistance*tDistance));
        tDistance +=1;
        t+=dt;

                                }

                            }




void temporalGauge(std::vector<Matrix<rSU,rSU>>& lattice){

    // only time links are set to identity
    size_t muT= 3; //4; // spatial ones are in {0,1,2}, right? 

    for(int i = 0; i<tAxis; i++){
        for(int j= 0; j<zAxis; j++){
            for(int k=0; k<yAxis; k++){
                for(int l=0; l<xAxis; l++){
                    
                        lattice[idx(l,k,j,i,muT)]=identityMatrix;

                }
            }
        }
    }    

}

void spacialGauge(std::vector<Matrix<rSU,rSU>>& lattice){

    //set all spacial links to identity
    size_t spacialDimensions = 3;
    for(int i = 0; i<tAxis; i++){
        for(int j= 0; j<zAxis; j++){
            for(int k=0; k<yAxis; k++){
                for(int l=0; l<xAxis; l++){
                    //for(int mu= 0; 0<spacialDimensions; mu++)
                    //a loop till the end of time :D
                    for(int mu= 0; mu<spacialDimensions; mu++)
                    
                        lattice[idx(l,k,j,i,mu)]=identityMatrix;

                }
            }
        }
    }    

}