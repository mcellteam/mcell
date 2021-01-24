cmake_minimum_required(VERSION 2.8.4)

project(bng_analyzer)

include_directories(
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../libs/"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../libs/sparsehash/src/"
)

include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/ ${CMAKE_SOURCE_DIR}/libs/sparsehash/src)

if (NOT WIN32)
  add_definitions(-std=c++14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall ")
endif()

if (NOT TARGET bng)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../libs/bng ${CMAKE_CURRENT_BINARY_DIR}/libs/bng)
endif()

set(SOURCE_FILES
    bng_analyzer.cpp
	nfsim_species_unifier.cpp
)

set(ALL_SOURCES
    ${SOURCE_FILES}
)

add_executable(${PROJECT_NAME}
    ${ALL_SOURCES}
)

target_link_libraries(${PROJECT_NAME}
    bng
)