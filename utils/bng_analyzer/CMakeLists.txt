cmake_minimum_required(VERSION 3.25...4.3)

project(bng_analyzer CXX)

# libbng (and its public include dirs: bng/, sparsehash, bundled libs) is
# normally provided by the parent build; pull it in only for a standalone build.
if (NOT TARGET libbng)
	add_subdirectory(
	  ${CMAKE_CURRENT_SOURCE_DIR}/../../../libbng/bng
	  ${CMAKE_CURRENT_BINARY_DIR}/libs/bng
	)
endif()

add_executable(${PROJECT_NAME}
  bng_analyzer.cpp
	nfsim_species_unifier.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<NOT:$<PLATFORM_ID:Windows>>:-Wall>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(${PROJECT_NAME}
    PRIVATE BNG::bng
)