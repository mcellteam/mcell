cmake_minimum_required(VERSION 3.25...4.3)

project(mesh_utils_test CXX)

# Checks the mesh table and the mesh checks (src4/mesh_table.hpp) against the
# expected results stored in mcell_tests/tests/utils_tests/mesh_utils; the test
# set mesh_utils of mcell_tests/test_configs/external.toml runs it. It needs only
# the geometry library.
add_executable(${PROJECT_NAME}
  mesh_utils_test.cpp
)

# New code is C++23 while MCell is built as C++17 (MCELL_CXX_STANDARD); as in
# src4/CMakeLists.txt, the flag is set on the source file so that it comes after
# the project-wide -std option.
if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
  set(MESH_UTILS_TEST_CXX23 /std:c++latest)
else()
  set(MESH_UTILS_TEST_CXX23 -std=gnu++2b)
endif()
set_source_files_properties(mesh_utils_test.cpp PROPERTIES COMPILE_OPTIONS "${MESH_UTILS_TEST_CXX23}")

target_link_libraries(${PROJECT_NAME}
    PRIVATE MCell::geometry
)
