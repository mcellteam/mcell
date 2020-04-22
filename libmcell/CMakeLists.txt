cmake_minimum_required(VERSION 2.8.4)

project(mcell)


execute_process(
    COMMAND python3 -m pybind11 --includes
    OUTPUT_VARIABLE PYTHON_INCLUDES
)
string(STRIP ${PYTHON_INCLUDES} PYTHON_INCLUDES)

include_directories(${CMAKE_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src4 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include 
    
    ${CMAKE_CURRENT_SOURCE_DIR}/../../mcell_tools/work/build_mcell/deps/
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs  
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_SOURCE_DIR}/../libs/sparsehash/src -Wall -std=c++17 ${PYTHON_INCLUDES}")

set(SOURCE_FILES
    api/bindings.cpp
    api/instantiation_data.cpp
    api/model.cpp
    api/subsystem.cpp

    generated/gen_component_type.cpp
    generated/gen_component_instance.cpp
    generated/gen_molecule_type.cpp
    generated/gen_molecule_instance.cpp
    generated/gen_complex_instance.cpp
    generated/gen_species.cpp

    generated/gen_geometry_object.cpp
    generated/gen_release_site.cpp

    generated/gen_subsystem.cpp
    generated/gen_instantiation_data.cpp
    generated/gen_model.cpp
    
    generated/gen_constants.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
    ${BISON_bngl_parser_OUTPUTS}
    ${FLEX_bngl_scanner_OUTPUTS}    
)

add_library(${PROJECT_NAME} SHARED
    ${ALL_SOURCES}
)
