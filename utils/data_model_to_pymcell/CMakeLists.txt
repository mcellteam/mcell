cmake_minimum_required(VERSION 2.8.4)

project(data_model_to_pymcell)

include_directories(
    ${CMAKE_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../.. 
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libs/jsoncpp/include
)

set(CMAKE_CXX_STANDARD 17)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/ ${CMAKE_SOURCE_DIR}/libs/sparsehash/src)

set(SOURCE_FILES
  data_model_to_pymcell.cpp
	mcell4_generator.cpp
	python_generator.cpp
	bngl_generator.cpp
	generator_utils.cpp
	data_model_geometry.cpp
	../../libmcell/api/api_utils.cpp
	../../libmcell/api/python_export_utils.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
)

add_executable(${PROJECT_NAME}
    ${ALL_SOURCES}
)

add_dependencies(data_model_to_pymcell version_h)  

target_link_libraries(${PROJECT_NAME}
    jsoncpp_lib libmcell ${VTK_LIBRARIES}
)