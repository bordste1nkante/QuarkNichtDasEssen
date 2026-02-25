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



inline size_t shift(size_t coord, int dir, size_t dim) {
    return (coord + dir + dim) % dim; // periodic boundary condition
}

// this mimics the behaviour of a 4D lattice from our 1D array
size_t idx(size_t x, size_t y, size_t z, size_t t, size_t mu){

    return linksPerSite*(x+ xAxis*(y+ yAxis*(z+zAxis*t)))+mu;
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
    size_t t = idx;

    return {mu,x,y,z,t};
}

void normalizeSU2Matrix(Matrix<2,2>& U){
    std::complex<double> sum;

    for(int i = 0; i<2; i++){
        sum += std::norm(U(0,i));
}
    sum= std::sqrt(sum);

    for(int i = 0; i<2; i++){
    U(0,i)= U(0,i)/sum;
}
    U(0,1)= -std::conj(U(1,0));
    U(1,1)=std::conj(U(1,1));
    


    std::complex<double> detU = U(0,0)*U(1,1)-U(1,0)*U(0,1);
    for(int m = 0; m< 2; m++){
        for(int n=0; n<2; n++){
            U(m,n)= U(m,n)/pow(detU, 1.0/double(2));
                        }
                    }
                }





//parallel
void X_updateSU3(size_t input){



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

        std::uniform_real_distribution<double> distEpsilon(-epsilon,epsilon);


        std::mt19937_64 zufall((index+3)*drng1*9731+(input+2)*13063);



        double r[3];

        double r0 = dist(zufall); 
        double s[3];

        double s0 = dist(zufall);
        double t[3];

        double t0 = dist(zufall);

  

        r[0]= dist(zufall);
        r[1]= dist(zufall);
        r[2]= dist(zufall); 


        double er = distEpsilon(zufall);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

        r0 = std::sqrt(std::norm(r0))*std::sqrt(1-er*er);

  

        s[0]= dist(zufall);
        s[1]= dist(zufall);
        s[2]= dist(zufall);   

  
        double es = distEpsilon(zufall);  

        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);

        s0 = std::sqrt(std::norm(s0))*std::sqrt(1-es*es);



        
        t[0]= dist(zufall);
        t[1]= dist(zufall);
        t[2]= dist(zufall); 


        double et = distEpsilon(zufall);
        double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);

        t0 = std::sqrt(std::norm(t0))*std::sqrt(1-et*et);


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
        normalizeSU2Matrix(Ssmall);
        normalizeSU2Matrix(Rsmall);
        normalizeSU2Matrix(Tsmall);


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
        S(2,1)=zero;
        S(0,2)=Ssmall(1,0);
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


        std::mt19937_64 hottestNumbS1(i*dhotNumbS1);


        std::mt19937_64 hottestNumbT1(i*dhotNumbT1);





        std::uniform_real_distribution<double> distribution(-0.5,0.5);
        std::uniform_real_distribution<double> hotDistribution(-hotEpsilon,hotEpsilon);



        //number required to generate 3 SU(2) matrices, from these we form a SU(3)
        double r[3];
        r[0]= distribution(hottestNumbR1);
        r[1]= distribution(hottestNumbR1);
        r[2]= distribution(hottestNumbR1); 

        double r0 = distribution(hottestNumbR1); 
        double er = hotDistribution(hottestNumbR1);

        double rLength = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        r0=std::sqrt(r0*r0)*std::sqrt(1-er*er);

        double s[3];
        s[0]= distribution(hottestNumbS1);
        s[1]= distribution(hottestNumbS1);
        s[2]= distribution(hottestNumbS1);        

        double s0 = distribution(hottestNumbS1); 
        double es = hotDistribution(hottestNumbS1);  

        double sLength = std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        s0=std::sqrt(s0*s0)*std::sqrt(1-es*es);

        double t[3];
        t[0]= distribution(hottestNumbT1);
        t[1]= distribution(hottestNumbT1);
        t[2]= distribution(hottestNumbT1);   

        double t0 = distribution(hottestNumbT1);

        double et = hotDistribution(hottestNumbT1);

        double tLength = std::sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        t0=std::sqrt(t0*t0)*std::sqrt(1-et*et);

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
        ST=matrix_multiplication(S,T);
        //matrix multiplication R*ST = U
        U=matrix_multiplication(R,ST);
        normalizeSU3Matrix(U);
            
        bufferLattice[i]=U;
                    
        }
                    

    );
    //switch pointers of buffer and original
    std::swap(bufferLattice,lattice);
    return;
}




//parallel
bool latticeAction(const std::vector<Matrix<rSU,rSU>>& lattice, const Matrix<rSU,rSU>& U,  const Matrix<rSU,rSU>& UPrime, const Matrix<rSU,cSU>& A, std::mt19937_64& ActionAccept, std::uniform_real_distribution<double>& Distribution){
    
    bool accept;
    double r;
    double probability;
    
    double SActionDif;
    Eigen::Matrix3cd M = translateMatrices(A);
    Eigen::Matrix3cd UEigen = translateMatrices(U);
    Eigen::Matrix3cd UPrimeEigen = translateMatrices(UPrime);


    //Matrix<rSU,cSU> sub = matrix_subtraction(UPrime, U);
    ////std::cout << "sub: "<< SActionDif << std::endl;
//
    //Matrix<rSU,cSU> mult = matrix_multiplication(sub,A);
    //double trace = matrix_trace(mult).real();
    SActionDif = -beta/(rSU)*((UPrimeEigen-UEigen)*M).trace().real(); //*(matrix_trace(matrix_multiplication(matrix_subtraction(UPrime, U),A))).real();
    probability = std::min(1.0, exp(- SActionDif)); 

    //std::cout << "SAction: "<< SActionDif << std::endl;

    r= Distribution(ActionAccept);
    //std::cout << r << std::endl;
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
    Matrix<rSU,cSU> A=zeroMatrix;

    double Sum;


    size_t Pmu_x;
    size_t Mmu_x;
    size_t Pnu_x;
    size_t Mnu_x;
    size_t Pmu_y;
    size_t Mmu_y;
    size_t Pnu_y;
    size_t Mnu_y;
    size_t Pmu_z;
    size_t Mmu_z;
    size_t Pnu_z;
    size_t Mnu_z;
    size_t Pmu_t;
    size_t Mmu_t;
    size_t Pnu_t;
    size_t Mnu_t;

    size_t complicatedx;
    size_t complicatedy;
    size_t complicatedz;
    size_t complicatedt;

    for(int nu=0; nu<linksPerSite; nu++){
        if(nu!=mu){


            Pmu_x=x;
            Mmu_x=x;
            Pnu_x=x;
            Mnu_x=x;
            Pmu_y=y;
            Mmu_y=y;
            Pnu_y=y;
            Mnu_y=y;
            Pmu_z=z;
            Mmu_z=z;
            Pnu_z=z;
            Mnu_z=z;
            Pmu_t=t;
            Mmu_t=t;
            Pnu_t=t;
            Mnu_t=t;

            complicatedx= x;
            complicatedy= y;
            complicatedz= z;
            complicatedt= t;

            switch (mu)
            {
            case 0: Pmu_x = shift(x,1,xAxis); Mmu_x = shift(x,-1,xAxis);complicatedx = shift(x,1,xAxis);break;
            case 1: Pmu_y = shift(y,1,yAxis); Mmu_y = shift(y,-1,yAxis);complicatedy = shift(y,1,yAxis);break;
            case 2: Pmu_z = shift(z,1,zAxis); Mmu_z = shift(z,-1,zAxis);complicatedz = shift(z,1,zAxis);break;
            case 3: Pmu_t = shift(t,1,tAxis); Mmu_t = shift(t,-1,tAxis);complicatedt = shift(t,1,tAxis);break;


            }
            switch (nu)
            {
            case 0: Pnu_x = shift(x,1,xAxis); Mnu_x = shift(x,-1,xAxis);complicatedx = shift(x,-1,xAxis);break;
            case 1: Pnu_y = shift(y,1,yAxis); Mnu_y = shift(y,-1,yAxis);complicatedy = shift(y,-1,yAxis);break;
            case 2: Pnu_z = shift(z,1,zAxis); Mnu_z = shift(z,-1,zAxis);complicatedz = shift(z,-1,zAxis);break;
            case 3: Pnu_t = shift(t,1,tAxis); Mnu_t = shift(t,-1,tAxis);complicatedt = shift(t,-1,tAxis);break;

            }

            ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(Mnu_x,Mnu_y,Mnu_z,Mnu_t,mu)]),lattice[idx(Mnu_x,Mnu_y,Mnu_z,Mnu_t,nu)]);
            ATemp1 = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(complicatedx,complicatedy,complicatedz,complicatedt,nu)]),ATemp1);
            ATemp2=matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(Pnu_x,Pnu_y,Pnu_z,Pnu_t,mu)]), matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
            ATemp2 =matrix_multiplication(lattice[idx(Pmu_x,Pmu_y,Pmu_z,Pmu_t,nu)],ATemp2);
            ATemp3=matrix_addition(ATemp1,ATemp2);
            A=matrix_addition(A,ATemp3);

        }
    }

    return A;
}


//parallel
Matrix<rSU,cSU> overrelaxation(const Matrix<rSU,cSU>& A, const Matrix<rSU,cSU>& U,std::uniform_int_distribution<int>& distribution,std::mt19937_64& random ){
    int reflect = distribution(random);
    Eigen::Matrix3cd R1 = Eigen::Matrix3cd::Identity();
    Eigen::Matrix3cd R2 = Eigen::Matrix3cd::Identity();
    Eigen::Matrix3cd R3 = Eigen::Matrix3cd::Identity();
    
    //R1.diagonal() << -1., -1.,  1.;
    //R2.diagonal() <<  1., -1., -1.;
    //R3.diagonal() << -1.,  1., -1.;

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

    //seems to be the wrong way around in paper
    Eigen::Matrix3cd Ur = Vadjoint*UEigen*O*V;//O*UEigen.adjoint()*O;//V*UEigen*O*Vadjoint;

    Ur = Ur.adjoint().eval();
    if(reflect == 1){
        Ur = R1*Ur*R1;
        //Ur(0,1)=-Ur(0,1);
        //Ur(1,0)=-Ur(1,0);
        //Ur(0,2)=-Ur(0,2);
        //Ur(2,0)=-Ur(2,0);

    }
    else if(reflect == 2){
        Ur=R2*Ur*R2;
        //Ur(0,1)=-Ur(0,1);
        //Ur(1,0)=-Ur(1,0);
        //Ur(1,2)=-Ur(1,2);
        //Ur(2,1)=-Ur(2,1);

    }
    else if(reflect == 3){
        Ur=R3*Ur*R3;
        //Ur(0,2)=-Ur(0,2);
        //Ur(2,0)=-Ur(2,0);
        //Ur(1,2)=-Ur(1,2);
        //Ur(2,1)=-Ur(2,1);

    }
    else{
        std::cout << "reflection failed due to unexpected value" << std::endl;
    }

    Eigen::Matrix3cd UR = V*Ur*Vadjoint*O.adjoint();
    std::cout << ((UEigen)*M).trace() << std::endl;
    std::cout << (UR*M).trace() << std::endl;

    UPrime = retranslateMatrices(UR);//UR
    //std::cout << matrix_trace (matrix_multiplication(U,A)) << std::endl;
    //std::cout << matrix_trace (matrix_multiplication(UPrime,A)) << std::endl;
    
    //std::cout <<"Output:" << ((UEigen-UR)*M).trace() << std::endl;
   // std::cout <<"subEigen:" << ((UEigen-UR)) << std::endl;
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



void plaquette(std::vector<Matrix<rSU,rSU>>& lattice, std::vector<double>& plaquette){
    std::vector<size_t> indices(xAxis*yAxis*zAxis*tAxis);
    std::iota(indices.begin(), indices.end(),0);
    std::vector<double> bufferPlaq (indices.size(), 0.0);


    std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){

    Matrix<rSU,cSU> Plaq = zeroMatrix;

    double tP;
    double tM;
    double xP;
    double xM;
    double yP;
    double yM;
    double zP;
    double zM;

    
    std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(linksPerSite*i);
    size_t mu,x,y,z,t;
    mu= std::get<0>(temp);
    x= std::get<1>(temp);
    y= std::get<2>(temp);
    z= std::get<3>(temp);
    t= std::get<4>(temp);


    for(int mu = 0; mu<linksPerSite; mu++){
        for (int nu = mu+1; nu < linksPerSite; ++nu) {
            // Forward shift coordinates
            size_t x_mu = x, y_mu = y, z_mu = z, t_mu = t;
            size_t x_nu = x, y_nu = y, z_nu = z, t_nu = t;
            switch(mu){ case 0: x_mu = shift(x,1,xAxis); break;
                        case 1: y_mu = shift(y,1,yAxis); break;
                        case 2: z_mu = shift(z,1,zAxis); break;
                        case 3: t_mu = shift(t,1,tAxis); break; }
            switch(nu){ case 0: x_nu = shift(x,1,xAxis); break;
                        case 1: y_nu = shift(y,1,yAxis); break;
                        case 2: z_nu = shift(z,1,zAxis); break;
                        case 3: t_nu = shift(t,1,tAxis); break; }

            Matrix <rSU,cSU> Temp  = matrix_multiplication(matrix_hermitean_conjugate(lattice[idx(x_nu,y_nu,z_nu,t_nu,mu)]),matrix_hermitean_conjugate(lattice[idx(x,y,z,t,nu)]));
            Temp = matrix_multiplication(lattice[idx(x_mu,y_mu,z_mu,t_mu,nu)], Temp);
            Temp= matrix_multiplication(lattice[idx(x,y,z,t,mu)],Temp);
            Plaq = matrix_addition(Plaq,Temp);


            }
        }
    bufferPlaq[i]= 1.0/(double(rSU)*6)* matrix_trace(Plaq).real();

    });
    std::swap(bufferPlaq,plaquette);
}




