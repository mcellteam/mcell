cmake_minimum_required(VERSION 2.8)
project(mcell C)

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories("${CMAKE_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}/deps")

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

# we need to hack around the current mdllex.l using sed which is targeted
# toward autotools' ylwrap and won't work with flex proper
#set_source_files_properties(mdllex_l PROPERTIES GENERATED true)
find_program (SED_TOOL NAMES sed)
if (NOT SED_TOOL)
  message( FATAL_ERROR "sed not found" )
endif (NOT SED_TOOL)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l
  COMMAND ${SED_TOOL} -e "/^%option outfile/ d" ${CMAKE_SOURCE_DIR}/src/mdllex.l > ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l
  DEPENDS ${CMAKE_SOURCE_DIR}/src/mdllex.l
  VERBATIM)

FLEX_TARGET(mdlScanner ${CMAKE_CURRENT_BINARY_DIR}/deps/mdllex.l
  ${CMAKE_CURRENT_BINARY_DIR}/deps/mdlex.c COMPILE_FLAGS -Crema)
ADD_FLEX_BISON_DEPENDENCY(mdlScanner mdlParser mdllex_l)

# create version.h target
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
  COMMAND ${CMAKE_SOURCE_DIR}/src/version.sh > ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
  DEPENDS ${CMAKE_SOURCE_DIR}/src/version.sh
  VERBATIM)

# check for needed libraries
find_library(M_LIB m)

set(CMAKE_C_FLAGS "-Wall -Wextra -Wshadow -Wno-unused-parameter -D_GNU_SOURCE=1 -O2 -std=c99 ${CMAKE_C_FLAGS}" )
set(CMAKE_EXE_LINKER_FLAGS ${M_LIB})

set(SOURCE_FILES
    src/api_test.c
    src/api_test.h
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
    src/grid_util.c
    src/grid_util.h
    src/init.c
    src/init.h
    src/isaac64.c
    src/isaac64.h
    src/logging.c
    src/logging.h
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
    src/react.h
    src/react_cond.c
    src/react_outc.c
    src/react_outc_trimol.c
    src/react_output.c
    src/react_output.h
    src/react_trig.c
    src/react_util.c
    src/react_util.h
    src/rng.c
    src/rng.h
    src/sched_util.c
    src/sched_util.h
    src/strfunc.c
    src/strfunc.h
    src/sym_table.c
    src/sym_table.h
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

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs 3.5 REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

#SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
SET(CMAKE_SWIG_FLAGS "-ltypemaps;-py3")

SET_SOURCE_FILES_PROPERTIES(src/mcellSwig.i PROPERTIES C ON)
SET_SOURCE_FILES_PROPERTIES(src/mcellSwig.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(mcellSwig python
  src/mcellSwig.i src/mcell_react_out.c src/mcell_misc.c src/util.c
  src/argparse.c src/mcell_objects.c src/mcell_reactions.c src/mcell_release.c
  src/mcell_species.c src/mcell_viz.c src/mcell_surfclass.c src/strfunc.c
  src/mcell_init.c src/chkpt.c src/init.c src/react_output.c src/version_info.c
  src/sched_util.c src/wall_util.c src/mem_util.c src/react_trig.c
  src/react_cond.c src/count_util.c src/vol_util.c src/diffuse.c src/vector.c
  src/diffuse_util.c src/sym_table.c src/triangle_overlap.c src/logging.c
  src/grid_util.c src/react_outc.c src/react_util.c src/react_outc_trimol.c
  src/diffuse_trimol.c src/isaac64.c src/rng.c src/viz_output.c src/mcell_run.c
  src/volume_output.c)
SWIG_LINK_LIBRARIES(mcellSwig ${PYTHON_LIBRARIES})

# copy the pyMCell test script into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test_api.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# build executable
add_executable(mcell
  ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
  ${SOURCE_FILES}
  src/mcell.c
  ${BISON_mdlParser_OUTPUTS}
  ${FLEX_mdlScanner_OUTPUTS})
target_link_libraries(mcell ${M_LIB})

# build little libmcell test
#add_executable(libmcell_test
#  ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
#  ${SOURCE_FILES}
#  src/libmcell_test.c
#  ${BISON_mdlParser_OUTPUTS}
#  ${FLEX_mdlScanner_OUTPUTS})
#target_link_libraries(libmcell_test ${M_LIB})
