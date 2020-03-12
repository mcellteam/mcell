cmake_minimum_required(VERSION 2.8.4)

project(bng)

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")

set(SOURCE_FILES
	bng_engine.cpp
	complex_species.cpp
	component.cpp
	molecule_instance.cpp
	molecule_species.cpp
	rule.cpp
	
	test.cpp # only for development
)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${INCLUDE_DIRS})

add_definitions(-std=c++17)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
