cmake_minimum_required(VERSION 3.14)

project(mcell4)

include_directories(${CMAKE_SOURCE_DIR} 
    "${CMAKE_CURRENT_SOURCE_DIR}/../src4" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../src" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../include" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../libs"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libmcell"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libs/sparsehash/src/"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libs/jsoncpp/include/"
)

set(SRC4_SOURCE_FILES

    base_event.cpp
    collision_structs.cpp
    defines.cpp
    diffuse_react_event.cpp
    mcell3_world_converter.cpp
    mcell4_iface_for_mcell3.cpp
    molecule.cpp
    partition.cpp
    release_event.cpp
    scheduler.cpp
    count_buffer.cpp
    mol_or_rxn_count_event.cpp
    viz_output_event.cpp
    defragmentation_event.cpp
    sort_mols_by_subpart_event.cpp
    geometry.cpp
    world.cpp
    counted_volumes_util.cpp
    bng_data_to_datamodel_converter.cpp
    
    legacy_callback_info.cpp    
)

add_library(${PROJECT_NAME} STATIC
    ${SRC4_SOURCE_FILES}
)
