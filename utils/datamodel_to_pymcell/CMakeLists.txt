cmake_minimum_required(VERSION 2.8.4)

project(datamodel_to_pymcell)

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/deps/)
add_definitions(-std=c++17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -isystem ${CMAKE_SOURCE_DIR}/libs/ -isystem ${CMAKE_SOURCE_DIR}/libs/sparsehash/src")

set(SOURCE_FILES
	pymcell_converter.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
)

add_library(${PROJECT_NAME} ${ALL_SOURCES})

add_executable(test_${PROJECT_NAME}
    ${ALL_SOURCES}
)
