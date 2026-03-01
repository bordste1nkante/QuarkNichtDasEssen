#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <execution>
#include <thread>
#include <numeric>
#include <complex>
#include <atomic>
#include <random>
#include <Eigen/Dense>
#include "../header/global.h"
#include "../header/latticeOP.h"
#include "../header/matrixOP.h"
#include "../header/utils.h"

void ensureGroup(H5::H5File& file, const std::string& datasetPath)
{
    auto lastSlash = datasetPath.find_last_of('/');
    if (lastSlash == std::string::npos || lastSlash == 0)
        return;

    std::string parent = datasetPath.substr(0, lastSlash);

    size_t pos = 1;
    while (true) {
        size_t next = parent.find('/', pos);
        std::string group =
            (next == std::string::npos)
            ? parent
            : parent.substr(0, next);

        if (!H5Lexists(file.getId(), group.c_str(), H5P_DEFAULT)) {
            file.createGroup(group);
        }

        if (next == std::string::npos)
            break;

        pos = next + 1;
    }
}


bool Setup_H5(  bool ColdOrHot,  const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps,  
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t observable){

        bool success = false;

        
    try
    {

        std::string fileString = "../h5/"+filenameh5;
        hsize_t xAxish5 = static_cast<hsize_t>(xAxis);
        hsize_t yAxish5 = static_cast<hsize_t>(yAxis);
        hsize_t zAxish5 = static_cast<hsize_t>(zAxis);
        hsize_t tAxish5 = static_cast<hsize_t>(tAxis);



        hsize_t Thermal = static_cast<hsize_t>(numberOfThermalSweeps);
        hsize_t Update = static_cast<hsize_t>(XUpdate);
        hsize_t Configs = static_cast<hsize_t>(NConfigs);
        hsize_t observable = static_cast<hsize_t>(observable);

        std::vector<hsize_t> startingh5;
        std::vector<hsize_t> endh5;
        
        for(int i = 0; i<startingPoint.size(); i++){
            startingh5.push_back(static_cast<hsize_t>(startingPoint[i]));

        }
        for(int i = 0; i<endPoint.size(); i++){
            endh5.push_back(static_cast<hsize_t>(endPoint[i]));

        }



        H5::H5File file(fileString,H5F_ACC_TRUNC);


        H5::Group group =file.createGroup("/metaData");
        H5::Group groupDim =file.createGroup("/metaData/Dimensions");
        H5::Group groupParams =file.createGroup("/metaData/runParams");


        H5::Attribute attrXAxis = groupDim.createAttribute(
            "xAxis",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrXAxis.write(H5::PredType::NATIVE_HSIZE, &xAxish5);

        H5::Attribute attrYAxis = groupDim.createAttribute(
            "yAxis",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrYAxis.write(H5::PredType::NATIVE_HSIZE, &yAxish5);


        H5::Attribute attrZAxis = groupDim.createAttribute(
            "zAxis",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrZAxis.write(H5::PredType::NATIVE_HSIZE, &zAxish5);


        H5::Attribute attrTAxis = groupDim.createAttribute(
            "tAxis",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrTAxis.write(H5::PredType::NATIVE_HSIZE, &tAxish5);


        H5::Attribute attrThermal = groupParams.createAttribute(
            "ThermalSteps",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrThermal.write(H5::PredType::NATIVE_HSIZE, &Thermal);
     


        H5::Attribute attrUpdate = groupParams.createAttribute(
            "XUpdate",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrUpdate.write(H5::PredType::NATIVE_HSIZE, &Update);


        H5::Attribute attrConfigs = groupParams.createAttribute(
            "NConfigurations",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrConfigs.write(H5::PredType::NATIVE_HSIZE, &Configs);


        H5::Attribute attrObservable = groupParams.createAttribute(
            "Collected observable",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrObservable.write(H5::PredType::NATIVE_HSIZE, &observable);

        hsize_t dim[1]={startingPoint.size()};
        H5::DataSpace dataspaceStart(1,dim);


        H5::DataSet start = file.createDataSet(
            "/metaData/StartPosition",
            H5::PredType::NATIVE_HSIZE,
            dataspaceStart
        );

        start.write(startingh5.data(), H5::PredType::NATIVE_HSIZE);




        hsize_t dimEnd[1]={endPoint.size()};
        H5::DataSpace dataspaceEnd(1,dim);


        H5::DataSet end = file.createDataSet(
            "/metaData/EndPosition",
            H5::PredType::NATIVE_HSIZE,
            dataspaceEnd
        );

        end.write(endh5.data(), H5::PredType::NATIVE_HSIZE);


        file.close();

        

        success = true;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    

    return success;


}



bool SaveTune_H5(const double epsilon, const double avgPlaq, const size_t SweepFactor){

        bool success = false;

        
    try
    {

        std::string fileString = "../h5/"+filenameh5;
        hsize_t xAxish5 = static_cast<hsize_t>(xAxis);
        hsize_t yAxish5 = static_cast<hsize_t>(yAxis);
        hsize_t zAxish5 = static_cast<hsize_t>(zAxis);
        hsize_t tAxish5 = static_cast<hsize_t>(tAxis);



        hsize_t Sweep = static_cast<hsize_t>(SweepFactor);
        double epsilonh5 = (epsilon);
        double avgPlaqh5 = (avgPlaq);




        H5::H5File file(fileString,H5F_ACC_RDWR);

        H5::Group groupTune =file.createGroup("/metaData/Tune");


        H5::Attribute Sweeping = groupTune.createAttribute(
            "SweepFactor",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        Sweeping.write(H5::PredType::NATIVE_HSIZE, &Sweep);


        H5::Attribute Epsiloning = groupTune.createAttribute(
            "Epsilon",
            H5::PredType::NATIVE_DOUBLE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        Epsiloning.write(H5::PredType::NATIVE_DOUBLE, &epsilonh5);


        H5::Attribute AvgPlaqing = groupTune.createAttribute(
            "avgPlaq",
            H5::PredType::NATIVE_DOUBLE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        AvgPlaqing.write(H5::PredType::NATIVE_DOUBLE, &avgPlaqh5);

        file.close();

        

        success = true;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    

    return success;


}



bool saveArrayH5(const std::vector<double>& array, std::string dataSetPath){

    bool success= false;
    

    try
    {
        std::string fileString = "../h5/"+filenameh5;
        H5::H5File file(fileString,H5F_ACC_RDWR);

        ensureGroup(file, dataSetPath);

        hsize_t dim[1] = {array.size()};
        H5::DataSpace dataspace(1,dim);

        H5::DataSet Store = file.createDataSet(
            dataSetPath,
            H5::PredType::NATIVE_DOUBLE,
            dataspace
        );
        Store.write(array.data(),H5::PredType::NATIVE_DOUBLE);

        file.close();
        



        success = true;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    

    return success;
}




bool saveArrayH5complex(const std::vector<std::complex<double>>& array, std::string dataSetPath){

    bool success= false;
    

    try
    {
        std::string fileString = "../h5/"+filenameh5;
        H5::H5File file(fileString,H5F_ACC_RDWR);

        ensureGroup(file, dataSetPath);

        std::vector<double> arrayCombined(2*array.size());

        for(size_t i= 0; i<array.size(); i++){
            arrayCombined[2*i]=array[i].real();
            arrayCombined[2*i+1]=array[i].imag();
        }
        hsize_t dim[2] = {array.size(),2};
        H5::DataSpace dataspace(2,dim);

        H5::DataSet Store = file.createDataSet(
            dataSetPath,
            H5::PredType::NATIVE_DOUBLE,
            dataspace
        );
        Store.write(arrayCombined.data(),H5::PredType::NATIVE_DOUBLE);

        file.close();
        



        success = true;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    

    return success;

}


double average(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;  // avoid division by zero
    
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return sum / vec.size();
}

double correlationFunc(const std::vector<double>& Plaqs,const std::vector<double>&OGPlaqs){

    double CX;
    double mean, meanOG, correlationMean;
    std::vector<double> product(Plaqs.size());
    std::vector<double> productOG(OGPlaqs.size());


    mean = average(Plaqs);
    //std::cout << "mean:" << mean << std::endl;
    meanOG = average(OGPlaqs);
    //std::cout << "meanOG:" << meanOG << std::endl;

    //all k or all k and i?
    for(size_t k = 0; k< Plaqs.size(); k++){
        product[k]= Plaqs[k]*OGPlaqs[k];
    }

    correlationMean=average(product);
    //std::cout << "correlation:" << correlationMean << std::endl;
    CX= correlationMean - mean*meanOG;
    return CX;


}

double ThermalTune(std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep){


    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);


    double avg1,avg2, P;
    //account for the fact that multiHits exist
    double PTestSize= PTestSizeGlobal/numberOfMultiHit;

    //ensure that enough Sets are picked
    if(PTestSize < 20){
        PTestSize = 20;
    }
    //std::cout << PTestSize << std::endl;
    std::vector<double> PTest(PTestSize);

    size_t stabilization=0;
    int allcounter = 0;
    while(stabilization < 3){ //stabilization < 5
        int counter = 0;

        while(counter < PTestSize){
        for(int color = 0; color < 2; ++color){
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*12872+(allcounter+2)* 7311);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);

            std::uniform_real_distribution<double> threadAcceptReject(0,1);
            std::uniform_int_distribution<int> threadrefelctDist(1,3);


            Matrix<rSU, cSU> UPrime, X, A, U; 

            U= lattice[i];

            std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(i);
            size_t mu,x,y,z,t;
            mu= std::get<0>(temp);
            x= std::get<1>(temp);
            y= std::get<2>(temp);
            z= std::get<3>(temp);
            t= std::get<4>(temp);
            if( ((x+y+z+t)%2) != color ){
                bufferLattice[i] = U;
                return;
            }
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    //normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,U);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, U, UPrime, A, Threadindexing, threadAcceptReject);

                    if(acceptance==true){
                        U=UPrime;
                        
                        }
                    
                    }


                }
            bufferLattice[i]=U;


        });
        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);
    }
        //Update X matrices
        if(counter%XUpdate ==0 && counter!=0){
            X_updateSU3(counter*drng1+27);
                        }

        std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
        plaquette(lattice, plaquettes);
        //std::cout << average(plaquettes) << std::endl;
        PTest[counter]=average(plaquettes);


        //we saved numberofMultiHit many variables in
        counter++;
    }

    size_t mid = PTest.size() / 2;
    std::vector<double> slice1(PTest.begin(), PTest.begin()+mid );
    std::vector<double> slice2(PTest.begin()+mid,PTest.end());

    avg1 = average(slice1);

    avg2 = average(slice2);


    P = (avg2-avg1)/avg1;
    //std::cout << "avg1: " << avg1 << std::endl;
    //std::cout << "avg2: " << avg2 << std::endl;
    //std::cout << "P: " << P << std::endl;
    //std::cout << "Low: " << changeRateLow << std::endl;
    //std::cout << "High: " << changeRateHigh << std::endl;

    

    if(P < changeRateHigh&& P>changeRateLow){
        stabilization++;
    }
    else{
        stabilization =0;
    }

    allcounter++;
    }
    return avg2;


    }

void epsilonTune( std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep){

    std::atomic<int> acceptanceRate=0;
    std::atomic<int> updates=0;

    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);
    size_t epsilonCounter=0;
    int counter =0;
    while(epsilonCounter<5){
        counter++;
        for(int color = 0; color < 2; ++color){
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*14002+(counter+2)* 12311);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);


            std::uniform_real_distribution<double> threadAcceptReject(0,1);
        

            std::uniform_int_distribution<int> threadrefelctDist(1,3);


            Matrix<rSU, cSU> UPrime, X, A, U; 

            U= lattice[i];

            std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(i);
            size_t mu,x,y,z,t;
            mu= std::get<0>(temp);
            x= std::get<1>(temp);
            y= std::get<2>(temp);
            z= std::get<3>(temp);
            t= std::get<4>(temp);
            if( ((x+y+z+t)%2) != color ){
                bufferLattice[i] = U;
                return;
            }
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);


                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,U);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, U, UPrime, A, Threadindexing, threadAcceptReject);

                    if(acceptance==true){
                        U=UPrime;
                        acceptanceRate.fetch_add(1);
                        }
                    updates.fetch_add(1);
                    }


                }
            bufferLattice[i]=U;

        });
        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);

        }
        //Update X matrices
        if(counter%XUpdate ==0 && counter!=0){
            double rate = double(acceptanceRate.load())/double(updates.load());
            //std::cout << rate << std::endl;
            epsilonCounter++;
            epsilon *= 1 + alpha * (rate - target_rate);
            if(rate > target_rate + rateInterval){
                epsilonCounter = 0;

                }

            if(rate<target_rate - rateInterval){
                epsilonCounter=0;

                }


            if(epsilon>0.7){
                epsilon=0.7;
            }
            //std::cout << epsilonCounter << std::endl;
            //std::cout << epsilon << std::endl;

            acceptanceRate.store(0);
            updates.store(0);
            
            X_updateSU3(counter*counter*drng1+207);
                        }



    }

}


size_t AutoCorrelationTune( std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep, const size_t observable){

    
    std::vector<Matrix<rSU,rSU>> bufferLattice(linksPerSite*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);
    double CX0=0;

    //for plaquette/Wilson
    std::vector<double> copiedplaquettes(xAxis*yAxis*zAxis+tAxis,0.0);
    plaquette(lattice, copiedplaquettes);


    //for polyakov
    std::vector<std::complex<double>> loops;
    std::vector<double> r;
    polyakovLoopSimple(lattice,loops,r);
    std::vector<double> P;
    //get absolute value
    for(int i=0; i < loops.size(); i++){
        P.push_back(std::sqrt(std::norm(loops[i])));
        }

    if(observable == 0 || observable ==1 ||observable ==3){

        CX0 = correlationFunc(copiedplaquettes, copiedplaquettes);
    }
    if(observable==2){

        CX0 = correlationFunc(P,P);
    }
    //std::cout <<"CX0:"<< CX0 << std::endl;

    //that acccounts for CX0/CX0 convention to take the halfs of all
    double integratedCorrelationTime = 0.5;



    bool condition = true;
    int counter = 0;

    while(condition){

        for(int color = 0; color < 2; ++color){

        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*12872+(counter+2)* 7911);
            std::uniform_int_distribution<int> threadIndexDist(0, XSet.size()-1);

            std::uniform_real_distribution<double> threadAcceptReject(0,1);
            std::uniform_int_distribution<int> threadrefelctDist(1,3);


            Matrix<rSU, cSU> UPrime, X, A, U; 

            U= lattice[i];

            std::tuple<size_t, size_t, size_t, size_t, size_t> temp =  ReIdx(i);
            size_t mu,x,y,z,t;
            mu= std::get<0>(temp);
            x= std::get<1>(temp);
            y= std::get<2>(temp);
            z= std::get<3>(temp);
            t= std::get<4>(temp);
            if( ((x+y+z+t)%2) != color ){
                bufferLattice[i] = U;
                return;
            }
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    //normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,U);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, U, UPrime, A, Threadindexing, threadAcceptReject);

                    if(acceptance==true){
                        U=UPrime;
                        
                        }
                    
                    }


                }
            bufferLattice[i]=U;


        });
        //they exchange pointers, so lattice now points to values of buffer and vice versa
        std::swap(lattice, bufferLattice);

        }
        if(observable == 0 || observable ==1 ||observable ==3){
            std::vector<double> plaquettes(xAxis*yAxis*zAxis*tAxis,0.0);
            plaquette(lattice, plaquettes);
            double CX = correlationFunc(plaquettes, copiedplaquettes);
            CX/=CX0;
            //std::cout << CX << std::endl;
            if(CX<=0){
                condition = false;
                break;
            }
            integratedCorrelationTime += CX;

        }
        if(observable==2){
            std::vector<std::complex<double>> loopsy;
            std::vector<double> ry;
            polyakovLoopSimple(lattice,loopsy,ry);
            std::vector<double> Poly;
            for(int i=0; i < loopsy.size(); i++){
                Poly.push_back(std::sqrt(std::norm(loopsy[i])));
                
            }

            double CX = correlationFunc(Poly,P);
            CX/=CX0;
            //std::cout << CX << std::endl;
            if(CX<=0){
                condition = false;
                break;
            }
            std::cout << CX << std::endl;
            integratedCorrelationTime += CX;
        }


        //Update X matrices
        if(counter%XUpdate ==0 && counter!=0){
            X_updateSU3(counter*drng1+27);
                        }  

    }
    //necessary factor from definition
    integratedCorrelationTime *=2;

    //why not just ceil on integrated and then time multihit, isn't the current wrong? it ensure size_t
    size_t temp = static_cast<size_t>(std::ceil(integratedCorrelationTime));
    size_t sweepFactor = temp;

    //size_t sweepFactor = integratedCorrelationTime*numberOfMultiHit;
    return sweepFactor;
    }



Eigen::Matrix3cd translateMatrices(const Matrix<rSU,cSU>& A){


    Eigen::Matrix3cd M;
    for(int i = 0; i<rSU; i++){
        for(int j = 0; j<cSU; j++){
            M(i,j)= A(i,j);

            
        }

    }
    return M;
}


Matrix<rSU,cSU> retranslateMatrices(const Eigen::Matrix3cd& A){

    Matrix<rSU,cSU> M;
    for(int i = 0; i<rSU; i++){
        for(int j = 0; j<cSU; j++){
            M(i,j)= A(i,j);
            
        }

    }
    return M;
}