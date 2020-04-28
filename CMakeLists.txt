cmake_minimum_required(VERSION 3.14)
project(mcell)

find_package(Python3 COMPONENTS Interpreter Development)
message("-- Python3 include dirs are ${Python3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})


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

# add VTK libs VTK and headers
set(VTK_DIR ${CMAKE_SOURCE_DIR}/../mcell_tools/work/build_vtk)
find_package(VTK REQUIRED)

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 ")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 ")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  SET(OPTIMIZATION_FLAGS " -O3 -march=core2 -finline-limit=10000 ")
  
  # must not be used for pymcell (at least for now)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie ")
endif()

set(WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable ")

set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")
SET(COMMON_C_CXX_FLAGS "${WARNING_FLAGS} -fPIC -g3 ${OPTIMIZATION_FLAGS} ${INCLUDE_DIRS} -D_GNU_SOURCE=1 ")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_CXX_FLAGS}  -std=c11 ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_C_CXX_FLAGS} -std=c++14 -fpermissive")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -L ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp/src/lib_json")

if (CYGWIN)
    # cygwin app run in cygwin environment is somehow unable to load cygstdc++-6.dll
    # although it works fine when ran from 'cmd' 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")
endif()


option(ENABLE_LTO "Use link-time optimization" OFF)
option(ENABLE_GCOV "Use gcc coverage" OFF) 
option(ENABLE_GPROF "Use gprof profiling" OFF)
option(ENABLE_GPERFTOOLS "Build with gperftools profiling" OFF)  
option(ENABLE_SANITIZER "Use address sanitizer" OFF)
option(PYMCELL "Build also PyMCell" OFF) 


if (ENABLE_LTO)
  SET(CMAKE_AR "gcc-ar")
  SET(CMAKE_LD "gcc-ld")
  SET(CMAKE_RANLIB "gcc-ranlib")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto ")
endif()

if (ENABLE_GCOV)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage ")
endif()

if (ENABLE_GPROF)
  # gcc 6.3 generates PIC code by default, this does not work with gprof, -no-pie disables it
  SET(COMMON_PROFILING_FLAGS "-pg -no-pie -fno-omit-frame-pointer -finline-limit=20")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_PROFILING_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_PROFILING_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COMMON_PROFILING_FLAGS}")
endif()

if (ENABLE_GPERFTOOLS)
  # TODO: build and install gperftools, for now it is done manually
  # then to run profiling 
  # export CPUPROFILE=./prof.out 
  # mcell
  # pprof --text <app> ./prof.out # text output
  # pprof --gv <app> ./prof.out # graphical output, requires gv installed
  # (http://euccas.github.io/blog/20170827/cpu-profiling-tools-on-linux.html)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITHGPERFTOOLS -I${CMAKE_CURRENT_BINARY_DIR}/libs/")
  SET(CMAKE_EXE_LINKER_FLAGS 
  	"${CMAKE_EXE_LINKER_FLAGS} -L${CMAKE_CURRENT_BINARY_DIR}/libs/install_gperftools/lib/ -L${CMAKE_CURRENT_SOURCE_DIR}/../../src4_bng/other/install_gperftools/lib/ -lprofiler")

endif()
	
if (ENABLE_SANITIZER)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lasan ")
endif()



set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories(
  "${CMAKE_SOURCE_DIR}/src" 
  "${CMAKE_SOURCE_DIR}/src4"   
  "${CMAKE_CURRENT_BINARY_DIR}/deps" 
  "${CMAKE_SOURCE_DIR}/../nfsimCInterface/src/"
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/libs/sparsehash/src/"
  "${CMAKE_SOURCE_DIR}/libs/jsoncpp/include/"
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

# the cmake flex plugin requires the flex file to be in the 'deps' directory
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

add_subdirectory(${CMAKE_SOURCE_DIR}/libs/bng ${CMAKE_CURRENT_BINARY_DIR}/libs/bng)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsim ${CMAKE_CURRENT_BINARY_DIR}/nfsim)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsimCInterface ${CMAKE_CURRENT_BINARY_DIR}/nfsimCInterface)
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/jsoncpp ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp)

set(SOURCE_FILES
    src/argparse.c
    src/chkpt.c
    src/count_util.c
    src/diffuse.c
    src/diffuse_trimol.c
    src/diffuse_util.c
    src/dyngeom.c
    src/dyngeom_parse_extras.c
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
    
    src/dump_state.cpp
    src/mdl2datamodel.cpp
    src/c_vector.cpp
    
    src4/base_event.cpp
    src4/collision_structs.cpp
    src4/defines.cpp
    src4/diffuse_react_event.cpp
    src4/mcell3_world_converter.cpp
    src4/molecule.cpp
    src4/partition.cpp
    src4/release_event.cpp
    src4/scheduler.cpp
    src4/count_buffer.cpp
    src4/mol_count_event.cpp
    src4/viz_output_event.cpp
    src4/defragmentation_event.cpp
    src4/geometry.cpp
    src4/world.cpp
    src4/counted_volumes_util.cpp
    
    src4/callback_info.cpp    
)

set(SOURCE_FILES_ONLY_MCELL
    src/mcell.c
    src/mdlparse_util.c
)    

SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILES} PROPERTIES LANGUAGE CXX )
set( SOURCE_FILES ${SOURCE_FILES} src/mdlparse.y src/dyngeom_parse.y )

SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILES_ONLY_MCELL} PROPERTIES LANGUAGE CXX )

link_directories(${CMAKE_BINARY_DIR}/lib)

# fallback to a different swig version:
if (NOT /usr/bin/swig2.0)
  set(SWIG_EXECUTABLE /usr/bin/swig)
endif()

if (UNIX)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN/./lib")
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

if (PYMCELL AND NOT WIN32 AND NOT CYGWIN)

  #SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
  SET(CMAKE_SWIG_FLAGS "-c++;-ltypemaps;-py3")

  SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
 
  
  SWIG_ADD_LIBRARY(pymcell 
    TYPE SHARED 
    LANGUAGE python
    SOURCES src/pymcell.i ${SOURCE_FILES}  ${BISON_dynGeomParser_OUTPUTS} ${FLEX_dynGeomScanner_OUTPUTS}  
  )
  
  set(PYMCELL_ARCHIVES ${CMAKE_CURRENT_BINARY_DIR}/libs/bng/libbng.a)
  
  if (APPLE)
    SWIG_LINK_LIBRARIES(pymcell ${PYMCELL_ARCHIVES} ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.dylib ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.dylib ${PYTHON_LIBRARIES})
  else()
    SWIG_LINK_LIBRARIES(pymcell  ${PYMCELL_ARCHIVES}  ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp/src/lib_json/libjsoncpp.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.so ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.so ${PYTHON_LIBRARIES})
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
SET(MCELL_RULES_DIR ${CMAKE_SOURCE_DIR}/src/rules_py)

configure_file(${MCELL_RULES_DIR}/grammar_definition.py ${CMAKE_CURRENT_BINARY_DIR}/grammar_definition.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/mcell3r.py ${CMAKE_CURRENT_BINARY_DIR}/mcell3r.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/mdlr2mdl.py ${CMAKE_CURRENT_BINARY_DIR}/mdlr2mdl.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/nfsim_python.py ${CMAKE_CURRENT_BINARY_DIR}/nfsim_python.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/postprocess_mcell3r.py ${CMAKE_CURRENT_BINARY_DIR}/postprocess_mcell3r.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/read_bngxml.py ${CMAKE_CURRENT_BINARY_DIR}/read_bngxml.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/read_mdl.py ${CMAKE_CURRENT_BINARY_DIR}/read_mdl.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/small_structures.py ${CMAKE_CURRENT_BINARY_DIR}/small_structures.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/split_bngxml.py ${CMAKE_CURRENT_BINARY_DIR}/split_bngxml.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/write_bngxmle.py ${CMAKE_CURRENT_BINARY_DIR}/write_bngxmle.py COPYONLY)
configure_file(${MCELL_RULES_DIR}/write_mdl.py ${CMAKE_CURRENT_BINARY_DIR}/write_mdl.py COPYONLY)

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/../bionetgen/bng2/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bng2/)

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  ${SOURCE_FILES_ONLY_MCELL}
  ${BISON_mdlParser_OUTPUTS}
  ${BISON_dynGeomParser_OUTPUTS}  
  ${FLEX_mdlScanner_OUTPUTS}
  ${FLEX_dynGeomScanner_OUTPUTS})

add_dependencies(mcell version_h)

if (PYMCELL AND NOT WIN32 AND NOT CYGWIN)
    add_dependencies(pymcell version_h bng nfsim_c NFsim jsoncpp_lib)
    target_link_libraries(pymcell ${VTK_LIBRARIES})
endif()
  
target_link_libraries(mcell bng nfsim_c_static NFsim_static jsoncpp_lib ${VTK_LIBRARIES})


TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
