cmake_minimum_required(VERSION 2.8.4)

project(mcell)

include_directories(${CMAKE_SOURCE_DIR} 
    /usr/include/python3.5m 
    /usr/local/include/python3.5 
    /home/ahusar/.local/include/python3.5m
    
    ${CMAKE_CURRENT_SOURCE_DIR}/../src4 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include 
    
    ${CMAKE_CURRENT_SOURCE_DIR}/../../mcell_tools/work/build_mcell/deps/ 
)

add_definitions(-std=c++17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -isystem ${CMAKE_SOURCE_DIR}/libs/")


set(SOURCE_FILES
    api/bindings.cpp

    generated/gen_geometry_object.cpp
    generated/gen_instantiation_data.cpp
    generated/gen_model.cpp
    generated/gen_release_site.cpp
    generated/gen_species.cpp
    generated/gen_subsystem.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
    ${BISON_bngl_parser_OUTPUTS}
    ${FLEX_bngl_scanner_OUTPUTS}    
)

add_library(${PROJECT_NAME} SHARED
    ${ALL_SOURCES}
)
