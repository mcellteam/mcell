cmake_minimum_required(VERSION 3.4)
project(mcell)

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -fsanitize=address ")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fsanitize=address ")

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lasan ")

SET(OPTIMIZATION_FLAGS "-O3 -march=core2 -finline-limit=1000 ")

#SET(OPTIMIZATION_FLAGS "-O3")
#SET(OPTIMIZATION_FLAGS "-O2")

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZATION_FLAGS} ${OPTIMIZATION_FLAGS} ${INCLUDE_DIRS} ")
SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11 ")

option(USE_GCOV "Use gcc coverage" OFF) 
option(USE_PROFILING "Use gprof profiling" OFF) 
option(USE_LTO "Use link-time optimization" OFF)
option(PYMCELL "Build also pyMCell" ON) 


if (USE_GCOV)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage ")
endif()

if (USE_PROFILING)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg ")
endif()

if (USE_LTO)
  SET(CMAKE_AR "gcc-ar")
  SET(CMAKE_LD "gcc-ld")
  SET(CMAKE_RANLIB "gcc-ranlib")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto ")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories(
  "${CMAKE_SOURCE_DIR}/src" 
  "${CMAKE_CURRENT_BINARY_DIR}/deps" 
  "${CMAKE_SOURCE_DIR}/nfsimCInterface/src/"
  "${CMAKE_SOURCE_DIR}/include"
)


# move config.h into place. Unix and Windows have different header files
if (UNIX)
  file(COPY ${CMAKE_SOURCE_DIR}/src/config-nix.h DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}/deps)
  file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/deps/config-nix.h
    ${CMAKE_CURRENT_BINARY_DIR}/deps/config.h)
elseif (WIN32)
  file(COPY ${CMAKE_SOURCE_DIR}/src/config-win.h DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}/deps)
  file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/deps/config-win.h
    ${CMAKE_CURRENT_BINARY_DIR}/deps/config.h)
endif()

# check for bison/flex and set up code gen
find_package(BISON)
find_package(FLEX)
BISON_TARGET(mdlParser ${CMAKE_SOURCE_DIR}/src/mdlparse.y
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlparse.c)

file(COPY ${CMAKE_SOURCE_DIR}/src/mdllex.l DESTINATION
  ${CMAKE_CURRENT_BINARY_DIR}/deps)

FLEX_TARGET(mdlScanner ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlex.c COMPILE_FLAGS -Crema)
ADD_FLEX_BISON_DEPENDENCY(mdlScanner mdlParser mdllex_l)

# create version.h target
if (NOT WIN32)
  add_custom_target(
    version_h
    COMMAND ${CMAKE_SOURCE_DIR}/src/version.sh > ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/version.sh
    VERBATIM)
elseif (WIN32)
  add_custom_target(version_h)
  add_custom_command(TARGET version_h COMMAND ${CMAKE_COMMAND} -E
                     copy ${CMAKE_SOURCE_DIR}/appveyor_windows/version.h ${CMAKE_CURRENT_BINARY_DIR}/deps )
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra -Wshadow -Wno-unused-parameter -D_GNU_SOURCE=1 -O2 -std=c11 ${CMAKE_C_FLAGS}" )
set(CMAKE_EXE_LINKER_FLAGS "-lm")

add_subdirectory(nfsim)
add_subdirectory(nfsimCInterface)

set(SRC_FILES
    src/argparse.c
    src/argparse.h
    src/chkpt.c
    src/chkpt.h
    src/config-nix.h
    src/config-win.h
    src/count_util.c
    src/count_util.h
    src/diffuse.c
    src/diffuse.h
    src/diffuse_trimol.c
    src/diffuse_util.c
    src/diffuse_util.h
    src/dyngeom.c
    src/dyngeom.h
    src/dyngeom_lex.c
    src/dyngeom_parse_extras.c
    src/dyngeom_parse_extras.h
    src/dyngeom_yacc.c
    src/grid_util.c
    src/grid_util.h
    src/map_c.cpp
    src/map_c.h
    src/init.c
    src/init.h
    src/isaac64.c
    src/isaac64.h
    src/logging.c
    src/logging.h
    src/mcell.c
    src/mcell_dyngeom.c
    src/mcell_dyngeom.h
    src/mcell_init.c
    src/mcell_init.h
    src/mcell_misc.c
    src/mcell_misc.h
    src/mcell_objects.c
    src/mcell_objects.h
    src/mcell_react_out.c
    src/mcell_react_out.h
    src/mcell_reactions.c
    src/mcell_reactions.h
    src/mcell_release.c
    src/mcell_release.h
    src/mcell_run.c
    src/mcell_run.h
    src/mcell_species.c
    src/mcell_species.h
    src/mcell_structs.h
    src/mcell_surfclass.c
    src/mcell_surfclass.h
    src/mcell_viz.c
    src/mcell_viz.h
    src/mdlparse.y
    src/mdlparse_aux.h
    src/mdlparse_util.c
    src/mdlparse_util.h
    src/mem_util.c
    src/mem_util.h
    src/minrng.c
    src/minrng.h
    src/nfsim_func.c
    src/nfsim_func.h
    src/react.h
    src/react_cond.c
    src/react_nfsim.h
    src/react_outc.c
    src/react_outc_nfsim.c
    src/react_outc_trimol.c
    src/react_output.c
    src/react_output.h
    src/react_trig.c
    src/react_trig_nfsim.c
    src/react_util.c
    src/react_util.h
    src/react_util_nfsim.c
    src/rng.c
    src/rng.h
    src/sched_util.c
    src/sched_util.h
    src/strfunc.c
    src/strfunc.h
    src/sym_table.c
    src/sym_table.h
    src/test_api.c
    src/test_api.h
    src/triangle_overlap.c
    src/triangle_overlap.h
    src/util.c
    src/util.h
    src/vector.c
    src/vector.h
    src/version_info.c
    src/version_info.h
    src/viz_output.c
    src/viz_output.h
    src/vol_util.c
    src/vol_util.h
    src/volume_output.c
    src/volume_output.h
    src/wall_util.c
    src/wall_util.h)

link_directories(${CMAKE_BINARY_DIR}/lib)
if (UNIX)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN/./lib")
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

if (NOT WIN32)
  FIND_PACKAGE(PythonLibs 3 REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  #SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
  SET(CMAKE_SWIG_FLAGS "-ltypemaps;-py3")

  SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES C ON)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
  SWIG_ADD_MODULE(pymcell python
    src/argparse.c
    src/chkpt.c
    src/count_util.c
    src/diffuse.c
    src/diffuse_trimol.c
    src/diffuse_util.c
    src/dyngeom.c
    src/dyngeom_lex.c
    src/dyngeom_parse_extras.c
    src/dyngeom_yacc.c
    src/grid_util.c
    src/map_c.cpp
    src/init.c
    src/isaac64.c
    src/logging.c
    src/mcell_dyngeom.c
    src/mcell_init.c
    src/mcell_misc.c
    src/mcell_objects.c
    src/mcell_react_out.c
    src/mcell_reactions.c
    src/mcell_release.c
    src/mcell_run.c
    src/mcell_species.c
    src/mcell_surfclass.c
    src/mcell_viz.c
    src/mem_util.c
    src/nfsim_func.c
    src/pymcell.i
    src/react_cond.c
    src/react_outc.c
    src/react_outc_nfsim.c
    src/react_outc_trimol.c
    src/react_output.c
    src/react_trig.c
    src/react_trig_nfsim.c
    src/react_util.c
    src/react_util_nfsim.c
    src/rng.c
    src/sched_util.c
    src/strfunc.c
    src/sym_table.c
    src/triangle_overlap.c
    src/util.c
    src/vector.c
    src/version_info.c
    src/viz_output.c
    src/vol_util.c
    src/volume_output.c
    src/wall_util.c)
    SWIG_LINK_LIBRARIES(pymcell nfsim_c NFsim ${PYTHON_LIBRARIES})

  # copy the pyMCell test scripts into place
  file(GLOB PYMCELL_EXAMPLES "${CMAKE_SOURCE_DIR}/src/pymcell_examples/*")
  file(COPY ${PYMCELL_EXAMPLES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

  file(COPY ${CMAKE_SOURCE_DIR}/src/pymcell_unittests.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
endif()

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/src/rules_py/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/bionetgen/bng2/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bng2/)

# build executable
add_executable(mcell
  ${SRC_FILES}
  src/mcell.c
  ${BISON_mdlParser_OUTPUTS}
  ${FLEX_mdlScanner_OUTPUTS})

add_dependencies(mcell version_h)
  
target_link_libraries(mcell nfsim_c_static NFsim_static)
TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
