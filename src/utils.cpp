#include <H5Cpp.h>
#include <stdlib.h>
#include <yaml-cpp/yaml.h> //requires installation of this specific library on Linux: sudo apt install libyaml-cpp-dev
#include <iostream>
#include <cmath>
#include <vector>
#include <execution>
#include <thread>
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


//store important meta data in h5file
bool Setup_H5(  bool ColdOrHot,  const std::vector<size_t>& startingPoint, 
    const std::vector<size_t>& endPoint, 
    const size_t numberOfThermalSweeps,  
    const size_t XUpdate, 
    const size_t NConfigs, 
    const size_t SweepFactor, 
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
        hsize_t Sweep = static_cast<hsize_t>(SweepFactor);
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



        H5::Attribute attrSweep = groupParams.createAttribute(
            "Sweep",
            H5::PredType::NATIVE_HSIZE,
            H5::DataSpace(H5S_SCALAR)
        );
        

        attrSweep.write(H5::PredType::NATIVE_HSIZE, &Sweep);


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




//store an array in h5 file
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



        std::vector<double> arrayCombined(array.size());

        for(size_t i= 0; i<array.size(); i++){
            arrayCombined[2*i]=array[i].real();
            arrayCombined[2*i+1]=array[i].imag();
        }
        hsize_t dim[2] = {arrayCombined.size(),2};
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



void ThermalAndCorrelationTune(std::vector<Matrix<rSU,cSU>>& lattice, 
    const size_t numberOfThermalSweeps, 
    const size_t XUpdate, 
    const size_t numberOfMultiHit,
    const size_t overrelaxationStep){


    std::vector<Matrix<rSU,rSU>> bufferLattice(4*xAxis*yAxis*zAxis*tAxis);
    std::vector<size_t> indices(lattice.size());
    std::iota(indices.begin(), indices.end(),0);

    bool condition = true;
    int counter = 0;
    while(condition){
        counter++;
        std::for_each(std::execution::par, indices.begin(), indices.end(),[&](size_t i){
            thread_local std::mt19937_64 Threadindexing(dindexing*i*12872+(counter+2)* 7311);
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
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,lattice[i]);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, lattice[i], UPrime, A, Threadindexing, threadAcceptReject);

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
    while(epsilonCounter<2){
  //parallelization
        counter++;
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
            A= determineA(lattice, x,y,z,t,mu);


            for(size_t j= 0; j<numberOfMultiHit; j++){


                //in theory automatically accepted
                if(j%overrelaxationStep==0 && j!=0){
                    U = overrelaxation(A, U, threadrefelctDist, Threadindexing);
                    normalizeSU3Matrix(U);

                }
                else{


                    size_t index = threadIndexDist(Threadindexing);
                    X=XSet[index];

                    UPrime = matrix_multiplication(X,lattice[i]);
                    normalizeSU3Matrix(UPrime);
                    bool acceptance = latticeAction(lattice, lattice[i], UPrime, A, Threadindexing, threadAcceptReject);

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
        double rate = double(acceptanceRate.load())/double(updates.load());
        epsilonCounter++;
        if(rate < 0.45){
            epsilonCounter = 0;
            if(rate < 0.35){
                epsilon *=0.8;
            }
            else{
                epsilon *=0.95;
            }
            
        }
        if(rate>0.6){
            epsilonCounter=0;
            if(rate < 0.7){
                epsilon*=1.2;
            }
            else{
                epsilon*=1.05;
            }

        }
        std::cout << epsilonCounter << std::endl;
        std::cout << epsilon << std::endl;

        acceptanceRate.store(0);
        updates.store(0);


    }

}

//translate our matrices to eigen, they have better support
Eigen::Matrix3cd translateMatrices(const Matrix<rSU,cSU>& A){


    Eigen::Matrix3cd M;
    for(int i = 0; i<rSU; i++){
        for(int j = 0; j<cSU; j++){
            M(i,j)= A(i,j);

            
        }

    }
    return M;
}

//translate eigen to our matrices, we use them :(
Matrix<rSU,cSU> retranslateMatrices(const Eigen::Matrix3cd& A){

    Matrix<rSU,cSU> M;
    for(int i = 0; i<rSU; i++){
        for(int j = 0; j<cSU; j++){
            M(i,j)= A(i,j);
            
        }

    }
    return M;
}