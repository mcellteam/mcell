cmake_minimum_required(VERSION 3.14)

project(mcell4)

include_directories(${CMAKE_SOURCE_DIR} 
    "${CMAKE_CURRENT_SOURCE_DIR}/../src4" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../src" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../include" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../libs"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libmcell"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../libbng/"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../libbng/libs/"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../libbng/libs/sparsehash/src/"
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
    clamp_release_event.cpp
    scheduler.cpp
    count_buffer.cpp
    mol_or_rxn_count_event.cpp
    viz_output_event.cpp
    defragmentation_event.cpp
    sort_mols_by_subpart_event.cpp
    rxn_class_cleanup_event.cpp
    species_cleanup_event.cpp
    wall_overlap.cpp
    grid_position.cpp
    geometry.cpp
    memory_limit_checker.cpp
    world.cpp
    simulation_stats.cpp
    simulation_config.cpp
    vtk_utils.cpp
    region_utils.cpp
    bng_data_to_datamodel_converter.cpp
)

add_library(${PROJECT_NAME} STATIC
    ${SRC4_SOURCE_FILES}
)
