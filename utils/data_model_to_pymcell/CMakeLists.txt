cmake_minimum_required(VERSION 2.8.4)

project(data_model_to_pymcell)

include_directories(
    ${CMAKE_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../.. 
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libs/jsoncpp/include
)

add_definitions(-std=c++14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_SOURCE_DIR}/libs/ -isystem ${CMAKE_SOURCE_DIR}/libs/sparsehash/src")

set(SOURCE_FILES
  data_model_to_pymcell.cpp
	mcell4_generator.cpp
	python_generator.cpp
	bngl_generator.cpp
	generator_utils.cpp
	../../libmcell/api/api_utils.cpp
	../../libmcell/api/python_export_utils.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
)

add_executable(${PROJECT_NAME}
    ${ALL_SOURCES}
)

target_link_libraries(${PROJECT_NAME}
    jsoncpp_lib
)