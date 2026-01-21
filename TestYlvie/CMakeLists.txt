cmake_minimum_required(VERSION 3.16)
project(staticquarksYlvie LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON )
set(CMAKE_CXX_EXTENSIONS OFF)


add_executable(staticquarks
TestY.cpp
../src/global.cpp
../src/latticeOP.cpp
../src/matrixOP.cpp
../src/observables.cpp

)

target_include_directories(staticquarks PRIVATE ../header)


find_package(yaml-cpp REQUIRED)
target_link_libraries(staticquarks PRIVATE yaml-cpp::yaml-cpp)



