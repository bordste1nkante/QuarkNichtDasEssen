#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"

// this mimics the behaviour of a 4D lattice from our 1D array
double idx(size_t x, size_t y, size_t z, size_t t, size_t mu){
    return 4*(x+ xAxis*(y+ yAxis*(z+zAxis*t)))+mu;
}


//updates the whole set of X matrices
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



        //number required to generate 3 SU(2) matrices, from these we form a SU(3)

        double r[3];
        double r0 = dist(rng); 


        double er = distEpsilon(randNumb);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);

        double s[3];
        double s0 = dist(rng); 

        double es = distEpsilon(randNumb);     
        
        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);


        double t[3];
        double t0 = dist(rng);

        double et = distEpsilon(randNumb);

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
            for(int j=0; j<2; j++){
                for(int k=0; k<2; k++){
                    Ssmall(j,k)= pauliMatrices[i](j,k)*s[i-1];
                    Rsmall(j,k)= pauliMatrices[i](j,k)*r[i-1];
                    Tsmall(j,k)= pauliMatrices[i](j,k)*t[i-1];
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
                X(i,j)= sum;
            }

        }

        //inverting X
        std::complex<double> detX = X(0,0)*(X(1,1)*X(2,2)-X(1,2)*X(2,1))-X(0,1)*(X(1,0)*X(2,2)-X(1,2)*X(2,0))+X(0,2)*(X(1,0)*X(2,1)-X(1,1)*X(2,0));


        // should insure that if det(X) \approx zero we just repeat the process and generate a new X
        double tol = 1e-9;
        if(-tol <detX.real() <tol && -tol<detX.imag()<tol){
            p=p-1;
        }
        else{
            //normalizes det to 1
            for(int i=0; i<rSU; i++){
                for(int j=0; j<rSU; j++){
                    X(i,j)= X(i,j)/detX;
                }
            }
            invX(0,0)= (X(1,1)*X(2,2)-X(1,2)*X(2,1))/detX;
            invX(0,1)=-(X(0,1)*X(2,2)-X(0,2)*X(2,1))/detX;
            invX(0,2)= (X(0,1)*X(1,2)-X(0,2)*X(1,1))/detX;

            invX(1,0)=-(X(1,0)*X(2,2)-X(1,2)*X(2,0))/detX;
            invX(1,1)=(X(0,0)*X(2,2)-X(0,2)*X(2,0))/detX;
            invX(1,2)=-(X(0,0)*X(1,2)-X(0,2)*X(1,0))/detX;

            invX(2,0)=(X(1,0)*X(2,1)-X(1,1)*X(2,0))/detX;
            invX(2,1)=-(X(0,0)*X(2,1)-X(0,1)*X(2,0))/detX;
            invX(2,2)=(X(0,0)*X(1,1)-X(0,1)*X(1,0))/detX;


            //save X and invX in our set of matrices
            XSet[2*p]= X;
            XSet[2*p+1]=invX;

        }


    }
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
                        std::complex<double> detU = U(0,0)*(U(1,1)*U(2,2)-U(1,2)*U(2,1))-U(0,1)*(U(1,0)*U(2,2)-U(1,2)*U(2,0))+U(0,2)*(U(1,0)*U(2,1)-U(1,1)*U(2,0));
                        for(int m = 0; m< rSU; m++){
                            for(int n=0; n<cSU; n++){
                                U(m,n)= U(m,n)/detU;
                        }
                    }
                        lattice[idx(l,k,j,i, mu)]= U;
                }
                }


            }
        }
    }



}

    



//sets all matrices to identity, a potential starting config
void cold_start(std::vector<Matrix<rSU,rSU>>& lattice){

    Matrix<rSU,cSU> identityMatrix;

    for(int i = 0; i<rSU; i++){
        for(int j=0; j<cSU; j++){
            if(i==j){
                std::complex<double> z = {1.0, 0.0};
                identityMatrix(i,j)= z;
            }
            else{

                std::complex<double> z1 = {0.0, 0.0};
                identityMatrix(i,j)=z1;
            }
        }
    }


    for(int i = 0; i<tAxis; i++){
        for(int j= 0; j<zAxis; j++){
            for(int k=0; k<yAxis; k++){
                for(int l=0; l<xAxis; l++){
                    for(int mu = 0; mu<linksPerSite; mu++){
                        lattice[idx(l,k,j,i,mu)]=identityMatrix;

                }
                }



            }
        }
    }

    return;
}

//sets all matrices to random ones, a possible starting config
void hot_start(std::vector<Matrix<rSU,rSU>>& lattice){
    Matrix<rSU,cSU> U;
    Matrix<rSU,cSU> R;
    Matrix<rSU,cSU> S;
    Matrix<rSU,cSU> T;
    Matrix<2,2> Rsmall;
    Matrix<2,2> Ssmall;
    Matrix<2,2> Tsmall;
    Matrix<rSU,cSU> ST;



    for(int a = 0; a<tAxis; a++){
        for(int b= 0; b<zAxis; b++){
            for(int c=0; b<yAxis; c++){
                for(int d=0; b<xAxis; d++){
                    for(int mu = 0; mu<linksPerSite; mu++){
                    //number required to generate 3 SU(2) matrices, from these we form a SU(3)

                    double r[3];
                    double r0 = dist(hotNumb); 


                    double er = hotDist(hotNumb);

                    double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
                    r0=r0/std::sqrt(r0*r0)*std::sqrt(1-er*er);

                    double s[3];
                    double s0 = dist(hotNumb); 

                    double es = hotDist(hotNumb);     

                    double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
                    s0=s0/std::sqrt(s0*s0)*std::sqrt(1-es*es);


                    double t[3];
                    double t0 = dist(hotNumb);

                    double et = hotDist(hotNumb);

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
                        for(int j=0; j<2; j++){
                            for(int k=0; k<2; k++){
                                Ssmall(j,k)= pauliMatrices[i](j,k)*s[i-1];
                                Rsmall(j,k)= pauliMatrices[i](j,k)*r[i-1];
                                Tsmall(j,k)= pauliMatrices[i](j,k)*t[i-1];
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
                    
                        lattice[idx(a,b,c,d, mu)]=U;
                }
                    
                }


            }
        }
    }

    return;

}



bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime ,size_t x, size_t y, size_t z, size_t t, size_t mu){
    
    bool accept;
    double r;
    double probability;
    
    //double SActtionU;
    //double SActionUPrime;
    double SActionDif;
    Matrix<rSU,rSU> ATemp1;
    Matrix<rSU,rSU> ATemp2;
    Matrix<rSU,rSU> ATemp3;
    Matrix<rSU,rSU> A;

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

                if(mu==1){

                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    yP= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,yP,z,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==2){
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,y,zP,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,y,z,t,mu)]), lattice[idx(xM,y,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xM,y,z,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,y,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
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
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,yM,z,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==2){
                    yP= bCX(y+1);
                    yM= bCX(y-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yM,zP,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    yP= bCX(y+1);
                    yM= bCX(y-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yM,z,t,mu)]), lattice[idx(x,yM,z,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yM,z,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yP,z,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
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
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yP,zM,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==0){
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,y,zM,t,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==3){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    zP= bCZ(x+1);
                    zM= bCZ(x-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zM,t,mu)]), lattice[idx(x,y,zM,t,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zM,tP,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zP,t,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,z,tP,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }
            if(nu==3){



                if(mu==1){

                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    y= bCY(y+1);
                    yM= bCY(y-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,yP,z,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,yP,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);


                }
                if(mu==2){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    zP= bCZ(z+1);
                    zM= bCZ(z-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,zP,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(x,y,zP,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }
                if(mu==0){
                    tP= bCT(t+1);
                    tM= bCT(t-1);
                    xP= bCX(x+1);
                    xM= bCX(x-1);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tM,mu)]), lattice[idx(x,y,z,tM,nu)]);
                    ATemp1 = matrix_multiplication(matrix_conjugate(lattice[idx(xP,y,z,tM,nu)]),ATemp1);
                    ATemp2 = matrix_multiplication(matrix_conjugate(lattice[idx(x,y,z,tP,mu)]),matrix_conjugate(lattice[idx(x,y,z,t,nu)]));
                    ATemp2 = matrix_multiplication(lattice[idx(xP,y,z,t,nu)],ATemp2);
                    ATemp3 = matrix_addition(ATemp1,ATemp2);
                    A = matrix_addition(A, ATemp3);
                

                }

            }



        }
    }
    SActionDif = -beta/(xAxis*yAxis*zAxis*tAxis)*(matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real();
    if(SActionDif > 1.0){
        probability=1;
    }
    else{
        probability=SActionDif;
    }
    r= uniformAcceptReject(acceptReject);
    if(r<=probability){
        accept = true;
    }
    else{
        accept = false;
    }





    



    return accept;
}