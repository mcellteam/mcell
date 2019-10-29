cmake_minimum_required(VERSION 3.14)
project(mcell)

#TODO: unify case for cmake commands, e.g. set vs SET


# Using default gcc 8.3.0 provided in Debian 10
# causes a segfault in nfSIM, however different gcc is ok
if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 8.3.0)
  message(FATAL_ERROR 
    "There is an error in nfSim, possibly due to a compiler issue, with gcc-8.3.0. "
    "Use a different compiler with CC=... and CXX=... "
  )
endif()

# silence UseSWIG warnings
if (${CMAKE_VERSION} VERSION_GREATER "3.13")
  cmake_policy(SET CMP0078 NEW)
endif()

if (${CMAKE_VERSION} VERSION_GREATER "3.13")
  cmake_policy(SET CMP0086 NEW)
endif()

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 ")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 ")

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  SET(USE_SANITIZER_DEFAULT ON)
else()
  SET(USE_SANITIZER_DEFAULT OFF)
  SET(OPTIMIZATION_FLAGS " -O3 -march=core2 -finline-limit=1000 ")
endif()

set(WARNING_FLAGS "-Wall -Wextra  -Wno-return-type -Wno-shadow -Wno-sign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable")

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")
SET(COMMON_C_CXX_FLAGS "${WARNING_FLAGS} -g3 ${OPTIMIZATION_FLAGS} ${INCLUDE_DIRS} -D_GNU_SOURCE=1")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_CXX_FLAGS}  -std=c11 ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_C_CXX_FLAGS} -std=c++14")

set(CMAKE_EXE_LINKER_FLAGS "-lm")

if (CYGWIN)
    # cygwin app run in cygwin environment is somehow unable to load cygstdc++-6.dll
    # although it works fine when ran from 'cmd' 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")
endif()


option(USE_GCOV "Use gcc coverage" OFF) 
option(USE_PROFILING "Use gprof profiling" OFF) 
option(USE_LTO "Use link-time optimization" OFF)
option(USE_SANITIZER "Use address sanitizer (ON for debug build, OFF for release build" ${USE_SANITIZER_DEFAULT})
option(PYMCELL "Build also pyMCell" ON) 


if (USE_SANITIZER)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address ")
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

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories(
  "${CMAKE_SOURCE_DIR}/src" 
  "${CMAKE_CURRENT_BINARY_DIR}/deps" 
  "${CMAKE_SOURCE_DIR}/../nfsimCInterface/src/"
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

add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsim ${CMAKE_CURRENT_BINARY_DIR}/nfsim)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsimCInterface ${CMAKE_CURRENT_BINARY_DIR}/nfsimCInterface)

set(SOURCE_FILES
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
    src/wall_util.h
)

set(SOURCE_FILES_ONLY_MCELL
    src/mdlparse_util.c
    src/mdlparse_util.h
)    

link_directories(${CMAKE_BINARY_DIR}/lib)
if (UNIX)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN/./lib")
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

if (NOT WIN32 AND NOT CYGWIN)
  FIND_PACKAGE(PythonLibs 3 REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  #SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
  SET(CMAKE_SWIG_FLAGS "-ltypemaps;-py3")

  SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES C ON)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
 
  
  if (${CMAKE_VERSION} VERSION_GREATER "3.8")
    SWIG_ADD_LIBRARY(pymcell 
      TYPE SHARED 
      LANGUAGE python
      SOURCES src/pymcell.i ${SOURCE_FILES}
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
  
  # configure_file copies file even for subsequent runs of make
  configure_file(${CMAKE_SOURCE_DIR}/src/data_model_import.py ${CMAKE_CURRENT_BINARY_DIR}/python/data_model_import.py COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/src/pymcell_unittests.py ${CMAKE_CURRENT_BINARY_DIR}/python/pymcell_unittests.py COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/src/pymcell_helpers.py ${CMAKE_CURRENT_BINARY_DIR}/python/pymcell_helpers.py COPYONLY)  
endif()


configure_file(${CMAKE_SOURCE_DIR}/src/pymcell_helpers.py ${CMAKE_CURRENT_BINARY_DIR}/python/pymcell_helpers.py COPYONLY)

# copy tools for MDLr to MDL
# TODO: list all files using configure_file
configure_file(${CMAKE_SOURCE_DIR}/src/rules_py/nfsim_python.py ${CMAKE_CURRENT_BINARY_DIR}/nfsim_python.py COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/src/rules_py/mdlr2mdl.py ${CMAKE_CURRENT_BINARY_DIR}/mdlr2mdl.py COPYONLY)

file(COPY ${CMAKE_SOURCE_DIR}/src/rules_py/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/../bionetgen/bng2/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bng2/)

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  ${SOURCE_FILES_ONLY_MCELL}
  src/mcell.c
  ${BISON_mdlParser_OUTPUTS}
  ${FLEX_mdlScanner_OUTPUTS})

add_dependencies(mcell version_h)

if (NOT WIN32 AND NOT CYGWIN)
    add_dependencies(pymcell version_h nfsim_c NFsim)
endif()
  
target_link_libraries(mcell nfsim_c_static NFsim_static)
TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
