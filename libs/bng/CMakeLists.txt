cmake_minimum_required(VERSION 2.8.4)

project(bng)

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/deps/)
add_definitions(-std=c++17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -isystem ${CMAKE_SOURCE_DIR}/libs/ -isystem ${CMAKE_SOURCE_DIR}/libs/sparsehash/src")


set(SOURCE_FILES
	ast.cpp
	bng_data.cpp
	bng_engine.cpp
	cplx_instance.cpp
	molecule_instance.cpp
	mol_instance.cpp
	mol_type.cpp
	parser_utils.cpp
	rxn_class.cpp
	rxn_container.cpp
	rxn_rule.cpp
	species.cpp
	species_container.cpp
	semantic_analyzer.cpp
)


# check for bison/flex and set up code gen
find_package(BISON)
find_package(FLEX)
BISON_TARGET(bngl_parser ${CMAKE_CURRENT_SOURCE_DIR}/bngl_parser.y
  ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_parser.cpp)

# the cmake flex plugin requires the flex file to be in the 'deps' directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bngl_scanner.l ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_scanner.l COPYONLY)

FLEX_TARGET(bngl_scanner ${CMAKE_CURRENT_SOURCE_DIR}/bngl_scanner.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_scanner.cpp)


# TODO: add bngl_scanner.l|y dependencies
  
ADD_FLEX_BISON_DEPENDENCY(bngl_scanner bngl_parser bngl_l)

set(ALL_SOURCES
    ${SOURCE_FILES}
    ${BISON_bngl_parser_OUTPUTS}
    ${FLEX_bngl_scanner_OUTPUTS}    
)

add_library(${PROJECT_NAME} ${ALL_SOURCES})

add_executable(test_${PROJECT_NAME}
    ${ALL_SOURCES}
	test.cpp # only for development
)
