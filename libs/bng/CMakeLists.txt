cmake_minimum_required(VERSION 2.8.4)

project(bng)

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/deps/)
add_definitions(-std=c++17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -isystem ${CMAKE_SOURCE_DIR}/libs/")


set(SOURCE_FILES
	bng_engine.cpp
	complex_species.cpp
	component.cpp
	molecule_instance.cpp
	molecule_species.cpp
	rule.cpp
	parser_utils.cpp
	test.cpp # only for development
)


# check for bison/flex and set up code gen
find_package(BISON)
find_package(FLEX)
BISON_TARGET(bngl_parser ${CMAKE_SOURCE_DIR}/bngl_parser.y
  ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_parser.cpp)

# the cmake flex plugin requires the flex file to be in the 'deps' directory
configure_file(${CMAKE_SOURCE_DIR}/bngl_scanner.l ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_scanner.l COPYONLY)

FLEX_TARGET(bngl_scanner ${CMAKE_SOURCE_DIR}/bngl_scanner.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/bngl_scanner.cpp)
  
ADD_FLEX_BISON_DEPENDENCY(bngl_scanner bngl_parser bngl_l)

add_executable(${PROJECT_NAME} 
    ${SOURCE_FILES}
    ${BISON_bngl_parser_OUTPUTS}
    ${FLEX_bngl_scanner_OUTPUTS}    
)
