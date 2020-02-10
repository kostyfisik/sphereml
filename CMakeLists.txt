cmake_minimum_required(VERSION 2.8.12)
project(sphereml)

add_subdirectory(pybind11)
pybind11_add_module(sphereml
        src/matrix.cpp
        src/pybind_sphereml.cpp
        src/sphereml.cpp
        src/spfunc.cpp
        src/directivity.cpp)
