cmake_minimum_required(VERSION 3.25...4.3)

project(mcell4 CXX)

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
    mol_order_shuffle_event.cpp
    wall_overlap.cpp
    grid_position.cpp
    region_expr.cpp
    geometry.cpp
    wall.cpp
    memory_limit_checker.cpp
    world.cpp
    simulation_stats.cpp
    simulation_config.cpp
    vtk_utils.cpp
    region_utils.cpp
    bng_data_to_datamodel_converter.cpp
    bngl_exporter.cpp
)

add_library(${PROJECT_NAME} STATIC
    ${SRC4_SOURCE_FILES}
)
add_library(MCell::mcell4 ALIAS ${PROJECT_NAME})

# Usage requirements (PUBLIC dirs propagate to mcell / mcell4_so which link mcell4).
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../libmcell
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libbng
    ${CMAKE_CURRENT_SOURCE_DIR}/../../nfsimCInterface/src)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_BINARY_DIR}/deps)            # generated config.h / version.h
# Third-party bundled headers: SYSTEM, so they never shadow our own headers.
# (sparsehash ships its own config.h, which must NOT win over deps/config.h.)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libbng/libs
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libbng/libs/sparsehash/src)
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs/jsoncpp/include)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

# vtk_utils / world / bngl_exporter include VTK directly. Linking the VTK
# targets pulls in their include directories (no manual paths needed).
target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
