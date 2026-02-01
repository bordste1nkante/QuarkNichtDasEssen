cmake_minimum_required(VERSION 3.16)
project(staticquarksVincent LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON )
set(CMAKE_CXX_EXTENSIONS OFF)


add_executable(staticquarks
TestV.cpp
../src/global.cpp
../src/latticeOP.cpp
../src/matrixOP.cpp
../src/observables.cpp
../src/utils.cpp
../src/Finish.cpp

)

target_include_directories(staticquarks PRIVATE ../header)


find_package(yaml-cpp REQUIRED)
target_link_libraries(staticquarks PRIVATE yaml-cpp::yaml-cpp)
find_package(HDF5 REQUIRED COMPONENTS CXX)
target_include_directories(staticquarks PRIVATE HDF5::HDF5)
target_link_libraries(staticquarks PRIVATE HDF5::HDF5)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(staticquarks PRIVATE Eigen3::Eigen)



