cmake_minimum_required(VERSION 3.25...4.3)

project(data_model_to_pymcell CXX)

add_executable(${PROJECT_NAME}
  data_model_to_pymcell.cpp
	mcell4_generator.cpp
	python_generator.cpp
	bngl_generator.cpp
	generator_utils.cpp
	data_model_geometry.cpp
	../../libmcell/api/api_utils.cpp
	../../libmcell/api/python_export_utils.cpp
)

add_dependencies(${PROJECT_NAME} version_h)

# src/src4/include, sparsehash and VTK headers arrive transitively from the
# linked libmcell / libbng / VTK targets; only this tool's extras are listed.
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
    ${CMAKE_BINARY_DIR}/deps
    ${CMAKE_SOURCE_DIR}/libs/jsoncpp/include)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

target_link_libraries(${PROJECT_NAME}
    PRIVATE JsonCpp::JsonCpp MCell::libmcell ${VTK_LIBRARIES} BNG::bng
    # libmcell references nanobind's static support-lib symbols; MSVC requires the
    # executable to resolve them, so link nanobind-static on Windows (last, for
    # link-order). Not needed by GNU ld / macOS dynamic_lookup. See mcell/CMakeLists.txt.
    $<$<PLATFORM_ID:Windows>:nanobind-static>
)
