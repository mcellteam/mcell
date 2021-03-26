cmake_minimum_required(VERSION 3.14)
project(mcell)

set(PYBIND11_PYTHON_VERSION 3.5)
find_package(Python3 ${PYBIND11_PYTHON_VERSION} EXACT REQUIRED COMPONENTS Interpreter Development)
message("-- Python3 include dirs are ${Python3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")


#TODO: unify case for cmake commands, e.g. set vs SET

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")
  add_definitions(-DYY_NO_UNISTD_H)
endif(MSVC)

# set default locations of libraries and executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin )
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib )
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/archives )
# Second, for multi-config builds (e.g. msvc)
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_CURRENT_BINARY_DIR}/bin )
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_CURRENT_BINARY_DIR}/lib )
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_CURRENT_BINARY_DIR}/archives  )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

#set(CMAKE_CXX_EXTENSIONS OFF)

# Using default gcc 8.3.0 provided in Debian 10
# causes a segfault in nfSIM, however different gcc is ok
if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 8.3.0)
  message(FATAL_ERROR 
    "There is an error in nfSim, possibly due to a compiler issue, with gcc-8.3.0. "
    "Use a different compiler with CC=... and CXX=... "
  )
endif()


# add VTK libs VTK and headers
set(VTK_DIR ${CMAKE_SOURCE_DIR}/../mcell_tools/work/build_vtk)
find_package(VTK REQUIRED)

# jsoncpp_test cannot be built with gperftools linking options
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/jsoncpp ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp)


if (NOT WIN32) 
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 ")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 ")

  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # when -march=native or -march=haswell is enabled, there is some error in MCell3 wall 
    # overlap detection, can be checked with mcell_tests_private/benchmarks/mdl/B4110_neuropil_no_custom_timestep 
    # seems to be related to AVX2, maybe gcc issue
    # changing to sandybridge gives 2-6% in performance (in executed instructions), however 
    # MCell does not work afterwards on machines such as Mac Mini 2018 with Intel i5 
    set(OPTIMIZATION_FLAGS " -O3 -g -march=core2 -finline-limit=10000 ")
      
      # must not be used for pymcell (at least for now)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie ")
  endif()


  set(WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable -Wno-class-memaccess")

  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-ignored-optimization-argument")
  endif()

  set(INCLUDE_DIRS " -isystem ${CMAKE_SOURCE_DIR}/libs/")
  set(COMMON_C_CXX_FLAGS "${WARNING_FLAGS} -fPIC -g3 ${OPTIMIZATION_FLAGS} ${INCLUDE_DIRS} -D_GNU_SOURCE=1 ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_CXX_FLAGS}  -std=c11 ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_C_CXX_FLAGS} -std=c++17 -fpermissive ")

  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -L ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp/src/lib_json")

  # special performance and debug options are not supported on Windows 
  option(ENABLE_LTO "Use link-time optimization" OFF)
  option(ENABLE_GCOV "Use gcc coverage" OFF) 
  option(ENABLE_GPROF "Use gprof profiling" OFF)
  option(ENABLE_GPERFTOOLS "Build with gperftools profiling" OFF)  
  option(ENABLE_SANITIZER "Use address sanitizer" OFF)

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
    # build and install gperftools, for now it is done manually
    # using script build_gperftools.sh
    # then to run profiling 
    # CPUPROFILE=./prof.out mcell
    # pprof --text <app> ./prof.out # text output
    # pprof --gv <app> ./prof.out # graphical output, requires gv installed
    # (http://euccas.github.io/blog/20170827/cpu-profiling-tools-on-linux.html)
    # pybind11 ignores CMAKE_SHARED_LINKER_FLAGS so we must put the -L args to CXX flags as well 
    SET(CMAKE_INCLUDE_GPERF_FLAGS "-L${CMAKE_CURRENT_BINARY_DIR}/libs/install_gperftools/lib/ -L${CMAKE_CURRENT_SOURCE_DIR}/build/install_gperftools/lib/")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITHGPERFTOOLS -I${CMAKE_CURRENT_BINARY_DIR}/libs/ ${CMAKE_INCLUDE_GPERF_FLAGS}")
    SET(CMAKE_LD_GPERF_FLAGS "${CMAKE_INCLUDE_GPERF_FLAGS} -lprofiler ")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LD_GPERF_FLAGS}")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_LD_GPERF_FLAGS}")
    SET(GPERFTOOLS_LIB "profiler")
  else()
    SET(GPERFTOOLS_LIB "")
  endif()
    
  if (ENABLE_SANITIZER)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address ")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lasan ")
  endif()

endif() # if (NOT WIN32) 


# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories(
  "${CMAKE_SOURCE_DIR}/src" 
  "${CMAKE_SOURCE_DIR}/src4"   all
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
  configure_file(${CMAKE_SOURCE_DIR}/src/config-nix.h ${CMAKE_CURRENT_BINARY_DIR}/deps/config.h COPYONLY)
elseif (WIN32)
  configure_file(${CMAKE_SOURCE_DIR}/src/config-win.h ${CMAKE_CURRENT_BINARY_DIR}/deps/config.h COPYONLY)
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
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlex.cpp COMPILE_FLAGS -Cema)
ADD_FLEX_BISON_DEPENDENCY(mdlScanner mdlParser mdllex_l)

configure_file(${CMAKE_SOURCE_DIR}/src/dyngeom_lex.l ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.l COPYONLY)

FLEX_TARGET(dynGeomScanner ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/dyngeom_lex.cpp COMPILE_FLAGS -Cema)
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

if (WIN32)
  # not compiling for now
  # add_subdirectory(${CMAKE_SOURCE_DIR}/libs/libunistd ${CMAKE_CURRENT_BINARY_DIR}/libs/libunistd)
  # include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/libunistd/unistd)
endif()

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

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

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

# for some reason MSVC searches for python library when linking utils although they should not need it
get_filename_component(PYTHON_LIB_DIR ${Python3_LIBRARY_RELEASE} DIRECTORY)
link_directories(${PYTHON_LIB_DIR})

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  ${SOURCE_FILES_ONLY_MCELL}
  ${BISON_mdlParser_OUTPUTS}
  ${BISON_dynGeomParser_OUTPUTS}  
  ${FLEX_mdlScanner_OUTPUTS}
  ${FLEX_dynGeomScanner_OUTPUTS})

# mcell must be for backwards compatibility located in the root directory
set_target_properties(mcell PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${CMAKE_BINARY_DIR}"
)

TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
add_dependencies(mcell version_h mcell4)  
target_link_libraries(mcell mcell4 libmcell_dummy bng nfsim_c_static NFsim_static jsoncpp_lib nauty ${VTK_LIBRARIES})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils/data_model_to_pymcell)
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/bng_analyzer) # used in MCell4 vs MCell4 validation
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/bng_analyzer_py)

if (NOT WIN32 AND NOT APPLE)
    set(STDC_FS stdc++fs)
endif()

add_dependencies(mcell4_so libmcell mcell4 bng jsoncpp_lib)
target_link_libraries(mcell4_so 
    PUBLIC libmcell
    PRIVATE mcell4 bng jsoncpp_lib nauty ${VTK_LIBRARIES} ${GPERFTOOLS_LIB} ${GPERFTOOLS_MCELLSO_DEP} ${STDC_FS} ${WIN_PYTHON}
    PUBLIC libmcell  # linking dependency issue, we must include limbcell once again to the list of libraries 
)

set_property(TARGET mcell4_so PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(mcell4_so PROPERTIES PREFIX "")
set_target_properties(mcell4_so PROPERTIES OUTPUT_NAME "mcell")
if (WIN32)
    set_target_properties(mcell4_so PROPERTIES SUFFIX ".pyd") 
else()
    set_target_properties(mcell4_so PROPERTIES SUFFIX ".so") # must be .so on MacOS as well
endif()
