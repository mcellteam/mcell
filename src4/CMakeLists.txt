cmake_minimum_required(VERSION 3.14)

project(mcell4)

# add VTK libs VTK and headers
set(VTK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../mcell_tools/work/build_vtk)
find_package(VTK REQUIRED)

include_directories(${CMAKE_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src4 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include 
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs  
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs/sparsehash/src/
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs/jsoncpp/include/
    
    # VTK - should be set by find_package but it does not work somehow...
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Common/DataModel/
    ${VTK_DIR}/Common/DataModel/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Common/Core/
    ${VTK_DIR}/Common/Core/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Common/Transforms/
    ${VTK_DIR}/Common/Transforms/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Common/Math/
    ${VTK_DIR}/Common/Math/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Filters/Core/
    ${VTK_DIR}/Filters/Core/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Filters/Modeling/
    ${VTK_DIR}/Filters/Modeling/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Common/ExecutionModel/
    ${VTK_DIR}/Common/ExecutionModel/
    ${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/Utilities/KWIML/
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
    geometry.cpp
    world.cpp
    counted_volumes_util.cpp
    bng_data_to_datamodel_converter.cpp
    
    callback_info.cpp    
)

add_library(${PROJECT_NAME} STATIC
    ${SRC4_SOURCE_FILES}
)
