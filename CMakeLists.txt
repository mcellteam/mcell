cmake_minimum_required(VERSION 3.4)
project(mcell)

# silence UseSWIG warnings
if (${CMAKE_VERSION} VERSION_GREATER "3.13")
  cmake_policy(SET CMP0078 NEW)
endif()

if (${CMAKE_VERSION} VERSION_GREATER "3.13")
  cmake_policy(SET CMP0086 NEW)
endif()

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  SET(USE_SANITIZER_DEFAULT ON)
else()
  SET(USE_SANITIZER_DEFAULT OFF)
endif()

option(USE_GCOV "Use gcc coverage" OFF) 
option(USE_PROFILING "Use gprof profiling" OFF) 
option(USE_LTO "Use link-time optimization" OFF)
option(USE_SANITIZER "Use address sanitizer (ON for debug build, OFF for release build" ${USE_SANITIZER_DEFAULT})
option(PYMCELL "Build also pyMCell" ON) 


if (USE_SANITIZER)
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -fsanitize=address ")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fsanitize=address ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lasan ")
endif()

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

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  SET(OPTIMIZATION_FLAGS "-O3 -g -march=core2 -finline-limit=1000 ")
else()
  SET(OPTIMIZATION_FLAGS "-O0 -g3 ")
endif()

SET(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3 ")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 ")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories(
  "${CMAKE_SOURCE_DIR}/src" 
  "${CMAKE_SOURCE_DIR}/src4" 
  "${CMAKE_CURRENT_BINARY_DIR}/deps" 
  "${CMAKE_SOURCE_DIR}/../nfsimCInterface/src/"
  "${CMAKE_SOURCE_DIR}/include"
  "/nadata/cnl/home/ahusar/tools/anaconda/pkgs/python-3.7.3-h0371630_0/include/python3.7m/" # FIXME: detemine automatically
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
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlparse.cpp)

BISON_TARGET(dynGeomParser ${CMAKE_SOURCE_DIR}/src/dyngeom_parse.y
  ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_yacc.cpp)

# not sure why is the lex file copied
configure_file(${CMAKE_SOURCE_DIR}/src/mdllex.l ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l COPYONLY)

FLEX_TARGET(mdlScanner ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlex.cpp COMPILE_FLAGS -Crema)
ADD_FLEX_BISON_DEPENDENCY(mdlScanner mdlParser mdllex_l)

configure_file(${CMAKE_SOURCE_DIR}/src/dyngeom_lex.l ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.l COPYONLY)

FLEX_TARGET(dynGeomScanner ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.cpp COMPILE_FLAGS -Crema)
ADD_FLEX_BISON_DEPENDENCY(dynGeomScanner dynGeomParser dyngeomlex_l)

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

set(OPTIMIZATION_FLAGS " -march=native -O3")
# FIXME: the -fpermissive arg should be removed
set(WARNING_FLAGS "-Wall -Wextra -Wno-shadow -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable -fpermissive")

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")

# -faligned-new - only for gcc >= 5.x
set(CMAKE_C_FLAGS "-g3  -D_GNU_SOURCE=1 ${WARNING_FLAGS} ${OPTIMIZATION_FLAGS} -std=c11 ${CMAKE_C_FLAGS}" )
set(CMAKE_CXX_FLAGS "-g3 -D_GNU_SOURCE=1 -std=c++14 ${INCLUDE_DIRS} ${WARNING_FLAGS} ${OPTIMIZATION_FLAGS} ${CMAKE_CXX_FLAGS}" )

set(CMAKE_EXE_LINKER_FLAGS "-lm")

add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsim ${CMAKE_CURRENT_BINARY_DIR}/nfsim)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsimCInterface ${CMAKE_CURRENT_BINARY_DIR}/nfsimCInterface)

set(SOURCE_FILES
    src/argparse.c
    src/chkpt.c
    src/count_util.c
    src/diffuse.c
    src/diffuse_trimol.c
    src/diffuse_util.c
    src/dyngeom.c
    #src/dyngeom_lex.c
    src/dyngeom_parse_extras.c
    #src/dyngeom_yacc.c
    
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
    src/dump_state.cpp
    src/mem_util.c
    src/minrng.c
    src/nfsim_func.c
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
    src/test_api.c
    src/triangle_overlap.c
    src/util.c
    src/vector.c
    src/version_info.c
    src/viz_output.c
    src/vol_util.c
    src/volume_output.c
    src/wall_util.c

    src/c_vector.cpp
    
    src4/base_event.cpp
    src4/defines.cpp
    src4/diffuse_react_event.cpp
    src4/mcell3_world_converter.cpp
    src4/molecule.cpp
    src4/partition.cpp
    src4/reaction.cpp
    src4/release_event.cpp
    src4/scheduler.cpp
    src4/species.cpp
    src4/viz_output_event.cpp
    src4/defragmentation_event.cpp
    src4/geometry.cpp
    src4/world.cpp
    src4/world_constants.cpp
    
    src4/callback_info.cpp
)

set(SOURCE_FILES_ONLY_MCELL
    src/mdlparse_util.c
    src/mcell.c
)    


set(SOURCE_FILES_ONLY_PYMCELL
)    


SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILES} PROPERTIES LANGUAGE CXX )
set( SOURCE_FILES ${SOURCE_FILES} src/mdlparse.y src/dyngeom_parse.y )

SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILES_ONLY_MCELL} PROPERTIES LANGUAGE CXX )


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
  SET(CMAKE_SWIG_FLAGS "-c++;-ltypemaps;-py3;-Wall")
  #SET(CMAKE_SWIG_FLAGS "-ltypemaps;-py3")

  SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  #SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES C ON)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES CPLUSPLUS ON)  
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
  
  
  if (${CMAKE_VERSION} VERSION_GREATER "3.8")
    SWIG_ADD_LIBRARY(pymcell 
      TYPE SHARED 
      LANGUAGE python
      SOURCES src/pymcell.i ${SOURCE_FILES} ${SOURCE_FILES_ONLY_PYMCELL} ${BISON_dynGeomParser_OUTPUTS} ${FLEX_dynGeomScanner_OUTPUTS}  
    )
  else()  
    SWIG_ADD_MODULE(pymcell python
      src/pymcell.i
      ${SOURCE_FILES})
  endif()
    
  if (APPLE)
    SWIG_LINK_LIBRARIES(pymcell ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.dylib ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.dylib ${PYTHON_LIBRARIES})
  else()
    SWIG_LINK_LIBRARIES(pymcell ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.so ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.so ${PYTHON_LIBRARIES})
  endif()

  # copy the pyMCell test scripts into place
  file(GLOB PYMCELL_EXAMPLES "${CMAKE_SOURCE_DIR}/src/pymcell_examples/*")
  file(COPY ${PYMCELL_EXAMPLES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
  
  # configure_file copies file even for following runs of make
  configure_file(${CMAKE_SOURCE_DIR}/src/data_model_import.py ${CMAKE_CURRENT_BINARY_DIR}/python/data_model_import.py COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/src/pymcell_unittests.py ${CMAKE_CURRENT_BINARY_DIR}/python/pymcell_unittests.py COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/src/pymcell_helpers.py ${CMAKE_CURRENT_BINARY_DIR}/python/pymcell_helpers.py COPYONLY)  
endif()

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/src/rules_py/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/../bionetgen/bng2/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bng2/)

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  ${SOURCE_FILES_ONLY_MCELL}
  ${BISON_mdlParser_OUTPUTS}
  ${BISON_dynGeomParser_OUTPUTS}
  ${FLEX_mdlScanner_OUTPUTS}
  ${FLEX_dynGeomScanner_OUTPUTS}  
)

add_dependencies(mcell version_h)
add_dependencies(pymcell version_h nfsim_c NFsim)
  
target_link_libraries(mcell nfsim_c_static NFsim_static)
TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
