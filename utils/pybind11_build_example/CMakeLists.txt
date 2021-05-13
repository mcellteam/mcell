cmake_minimum_required(VERSION 3.14)
project(mcell)

set(PYBIND11_PYTHON_VERSION 3.5)
find_package(Python3 ${PYBIND11_PYTHON_VERSION} EXACT REQUIRED COMPONENTS Interpreter Development)
message("-- Python3 include dirs are ${Python3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

add_subdirectory(${CMAKE_SOURCE_DIR}/../../libs/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/libs/pybind11)

pybind11_add_module(mcell MODULE
    test.cpp
)    
