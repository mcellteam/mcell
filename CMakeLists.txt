cmake_minimum_required(VERSION 3.14)
project(mcell)

set(PYBIND11_PYTHON_VERSION 3.5)
find_package(Python3 ${PYBIND11_PYTHON_VERSION} EXACT REQUIRED COMPONENTS Interpreter Development)
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
  # when -march=native or -march=haswell is enabled, there is some error in MCell3 wall 
  # overlap detection, can be checked with mcell_tests_private/benchmarks/mdl/B4110_neuropil_no_custom_timestep 
  # seems to be related to AVX2, maybe gcc issue
  # changing to sandybridge gives 2-6% in performance (in executed instructions), however 
  # MCell does not work afterwards on machines such as Mac Mini 2018 with Intel i5 
  SET(OPTIMIZATION_FLAGS " -O3 -g -march=core2 -finline-limit=10000 ")
  
  # must not be used for pymcell (at least for now)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie ")
endif()

set(WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable ")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-ignored-optimization-argument")
endif()

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
option(PYMCELL "Build also PyMCell" ON)


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

# jsoncpp_test cannot be built with gperftools linking options
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/jsoncpp ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp)

if (ENABLE_GPERFTOOLS)
  # build and install gperftools, for now it is done manually
  # usign script build_gperftools.sh
  # then to run profiling 
  # CPUPROFILE=./prof.out mcell
  # pprof --text <app> ./prof.out # text output
  # pprof --gv <app> ./prof.out # graphical output, requires gv installed
  # (http://euccas.github.io/blog/20170827/cpu-profiling-tools-on-linux.html)
  # pybind11 ignores CMAKE_SHARED_LINKER_FLAGS so we must put the -L args to CXX flags as well 
  SET(CMAKE_INCLUDE_GPERF_FLAGS "-L${CMAKE_CURRENT_BINARY_DIR}/libs/install_gperftools/lib/ -L${CMAKE_CURRENT_SOURCE_DIR}/build/install_gperftools/lib/")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITHGPERFTOOLS -I${CMAKE_CURRENT_BINARY_DIR}/libs/ ${CMAKE_INCLUDE_GPERF_FLAGS}")
  SET(CMAKE_LD_GPERF_FLAGS "${CMAKE_INCLUDE_GPERF_FLAGS} -lprofiler -ltcmalloc")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LD_GPERF_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_LD_GPERF_FLAGS}")
  SET(GPERFTOOLS_LIB "profiler")
  SET(GPERFTOOLS_MCELLSO_DEP "tcmalloc")
else()
  SET(GPERFTOOLS_LIB "")
  SET(GPERFTOOLS_MCELLSO_DEP "")  
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
  "${CMAKE_SOURCE_DIR}/libs"  
  "${CMAKE_SOURCE_DIR}/libmcell"  
      
   # VTK - should be set by find_package but it does not work somehow...
   "${CMAKE_SOURCE_DIR}/../VTK/Common/DataModel/"
   "${VTK_DIR}/Common/DataModel/"
   "${CMAKE_SOURCE_DIR}/../VTK/Common/Core/"
   "${VTK_DIR}/Common/Core/"
   "${CMAKE_SOURCE_DIR}/../VTK/Common/Transforms/"
   "${VTK_DIR}/Common/Transforms/"
   "${CMAKE_SOURCE_DIR}/../VTK/Common/Math/"
   "${VTK_DIR}/Common/Math/"
   "${CMAKE_SOURCE_DIR}/../VTK/Filters/Core/"
   "${VTK_DIR}/Filters/Core/"
   "${CMAKE_SOURCE_DIR}/../VTK/Filters/Modeling/"
   "${VTK_DIR}/Filters/Modeling/"
   "${CMAKE_SOURCE_DIR}/../VTK/Common/ExecutionModel/"
   "${VTK_DIR}/Common/ExecutionModel/"
   "${CMAKE_SOURCE_DIR}/../VTK/Utilities/KWIML/"
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

# linking on centos7 with 7.3.1 fails with undefined reference to `casc::detail::asc_NodeBase::~asc_NodeBase()'
set(PYBIND11_LTO_CXX_FLAGS "")
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/libs/pybind11)
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/bng ${CMAKE_CURRENT_BINARY_DIR}/libs/bng)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsim ${CMAKE_CURRENT_BINARY_DIR}/nfsim)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsimCInterface ${CMAKE_CURRENT_BINARY_DIR}/nfsimCInterface)

add_subdirectory(${CMAKE_SOURCE_DIR}/libmcell ${CMAKE_CURRENT_BINARY_DIR}/libmcell)
add_subdirectory(${CMAKE_SOURCE_DIR}/src4 ${CMAKE_CURRENT_BINARY_DIR}/src4)

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
    #src/isaac64.c
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
    #src/rng.c
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
    src/edge_util.c

    src/bng_util.cpp
    src/dump_state.cpp
    src/c_vector.cpp
)

set(SOURCE_FILES_ONLY_MCELL
    src/mcell.c
    src/mdlparse_util.c
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

if (PYMCELL AND NOT WIN32 AND NOT CYGWIN)

  #SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
  SET(CMAKE_SWIG_FLAGS "-c++;-ltypemaps;-py3")

  SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
  #SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
 
  
  SWIG_ADD_LIBRARY(pymcell 
    TYPE SHARED 
    LANGUAGE python
    SOURCES src/pymcell.i ${SOURCE_FILES}  ${BISON_dynGeomParser_OUTPUTS} ${FLEX_dynGeomScanner_OUTPUTS}  
  )
  
  set(PYMCELL_ARCHIVES ${CMAKE_CURRENT_BINARY_DIR}/libs/bng/libbng.a ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp/src/lib_json/libjsoncpp.a)
  
  if (APPLE)
    SWIG_LINK_LIBRARIES(pymcell 
        ${PYMCELL_ARCHIVES}         
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.dylib 
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.dylib 
        ${PYTHON_LIBRARIES}
    )
  else()
    SWIG_LINK_LIBRARIES(pymcell 
        -Wl,--whole-archive
        ${PYMCELL_ARCHIVES} 
        -Wl,--no-whole-archive 
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.so 
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.so 
        ${PYTHON_LIBRARIES}
    )
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

SET(MCELL_UTILS_DIR ${CMAKE_SOURCE_DIR}/utils)
SET(BUILD_UTILS_DIR ${CMAKE_CURRENT_BINARY_DIR}/utils)
configure_file(${MCELL_UTILS_DIR}/mcell4_runner/mcell4_runner.py ${BUILD_UTILS_DIR}/mcell4_runner.py COPYONLY)
configure_file(${MCELL_UTILS_DIR}/plotter/plotter.py ${BUILD_UTILS_DIR}/plotter.py COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/libmcell/generated/mcell.pyi ${BUILD_UTILS_DIR}/mcell.pyi COPYONLY)

# copy tools for MDLr to MDL
# TODO: list all files using configure_file
file(COPY ${CMAKE_SOURCE_DIR}/../bionetgen/bng2/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bng2/)

pybind11_add_module(mcell4_so MODULE  
    # few files from mcell3 that are used
    # TODO: get rid of all deps except for the random generator
    src/logging.c
    src/mem_util.c
    src/strfunc.c
    src/util.c
    src/edge_util.c
    src/sym_table.c
    src/vector.c
    src/dump_state.cpp
    src/bng_util.cpp
            
    libmcell/api/exports.cpp
)

add_dependencies(mcell4_so libmcell mcell4 bng jsoncpp_lib ${GPERFTOOLS_MCELLSO_DEP})
target_link_libraries(mcell4_so 
    PUBLIC libmcell
    PRIVATE mcell4 bng jsoncpp_lib nauty ${VTK_LIBRARIES} ${GPERFTOOLS_LIB} ${GPERFTOOLS_MCELLSO_DEP}
)
set_property(TARGET mcell4_so PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(mcell4_so PROPERTIES PREFIX "")
set_target_properties(mcell4_so PROPERTIES OUTPUT_NAME "mcell")
set_target_properties(mcell4_so PROPERTIES SUFFIX ".so") # must be .so on MacOS as well, TODO fix for Windows

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  ${SOURCE_FILES_ONLY_MCELL}
  ${BISON_mdlParser_OUTPUTS}
  ${BISON_dynGeomParser_OUTPUTS}  
  ${FLEX_mdlScanner_OUTPUTS}
  ${FLEX_dynGeomScanner_OUTPUTS})

TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
add_dependencies(mcell version_h mcell4)  
target_link_libraries(mcell mcell4 libmcell_callbacks_empty bng nfsim_c_static NFsim_static jsoncpp_lib nauty ${VTK_LIBRARIES})

if (PYMCELL AND NOT WIN32 AND NOT CYGWIN)
    add_dependencies(pymcell mcell4 version_h bng nfsim_c NFsim jsoncpp_lib)
    target_link_libraries(pymcell mcell4 nauty ${VTK_LIBRARIES})
    set_target_properties(pymcell PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)	
endif()

# utils

add_subdirectory(${CMAKE_SOURCE_DIR}/utils/data_model_to_pymcell)
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/bng_analyzer)
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/bng_analyzer_py)
