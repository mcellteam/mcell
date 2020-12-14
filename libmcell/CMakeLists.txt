cmake_minimum_required(VERSION 2.8.4)

project(libmcell)

option(ENABLE_DOC "Enable build of documentation" OFF)

find_package(Python3 ${PYBIND11_PYTHON_VERSION} EXACT REQUIRED COMPONENTS Interpreter Development)
message("-- Python3 include dirs are ${Python3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src4 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    
    ${CMAKE_CURRENT_SOURCE_DIR}/../../mcell_tools/work/build_mcell/deps/
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs  
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs/jsoncpp/include
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g3 -O0 -isystem ${CMAKE_CURRENT_SOURCE_DIR}/../libs/sparsehash/src -std=c++14 ${PYTHON_INCLUDES}")

set(SOURCE_FILES
    api/bindings.cpp
    
    api/component_type.cpp
    api/component.cpp
    api/elementary_molecule_type.cpp
    api/elementary_molecule.cpp
    api/complex.cpp
    api/species.cpp
    api/reaction_rule.cpp
    api/surface_class.cpp
    
    api/count.cpp
    
    api/subsystem.cpp
    api/instantiation.cpp
    api/observables.cpp
    api/model.cpp
    
    api/geometry_utils.cpp
    api/bngl_utils.cpp
    api/compartment_utils.cpp
    
    api/globals.cpp
    api/mcell4_converter.cpp
    api/callbacks.cpp
    
    api/molecule.cpp
    api/wall.cpp

    generated/gen_component_type.cpp
    generated/gen_component.cpp
    generated/gen_elementary_molecule_type.cpp
    generated/gen_elementary_molecule.cpp
    generated/gen_complex.cpp
    generated/gen_species.cpp
    generated/gen_surface_property.cpp
    generated/gen_surface_class.cpp
    generated/gen_reaction_rule.cpp
    generated/gen_subsystem.cpp

    generated/gen_region.cpp
    generated/gen_initial_surface_release.cpp
    generated/gen_surface_region.cpp
    generated/gen_geometry_object.cpp
    generated/gen_release_pattern.cpp
    generated/gen_molecule_release_info.cpp
    generated/gen_release_site.cpp
    generated/gen_instantiation.cpp
    
    generated/gen_count_term.cpp
    generated/gen_count.cpp
    generated/gen_viz_output.cpp
    generated/gen_observables.cpp
    
    generated/gen_config.cpp
    generated/gen_notifications.cpp
    generated/gen_warnings.cpp
    generated/gen_model.cpp
    
    generated/gen_molecule.cpp
    generated/gen_wall.cpp
    generated/gen_wall_wall_hit_info.cpp
    
    generated/gen_mol_wall_hit_info.cpp
    generated/gen_reaction_info.cpp

    generated/gen_constants.cpp

    generated/gen_geometry_utils.cpp
    generated/gen_bngl_utils.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
    ${BISON_bngl_parser_OUTPUTS}
    ${FLEX_bngl_scanner_OUTPUTS}    
)

add_library(${PROJECT_NAME} STATIC
    ${ALL_SOURCES}
)

add_dependencies(${PROJECT_NAME} version_h)

add_library(libmcell_callbacks_empty STATIC
    api/callbacks_empty.cpp
)

