cmake_minimum_required(VERSION 3.4)
project(mcell)

add_definitions(-DPACKAGE_BUGREPORT="mcell@salk.edu")

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# directories holding flex/bison files and out of source includes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps)
include_directories("${CMAKE_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}/deps" "${CMAKE_SOURCE_DIR}/include")

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
add_custom_target(
  version_h
  COMMAND ${CMAKE_SOURCE_DIR}/src/version.sh > ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
  DEPENDS ${CMAKE_SOURCE_DIR}/src/version.sh
  VERBATIM)

# check for needed libraries
find_library(M_LIB m)

set(CMAKE_C_FLAGS "-Wall -Wextra -Wshadow -Wno-unused-parameter -D_GNU_SOURCE=1 -O2 -std=c11 ${CMAKE_C_FLAGS}" )
set(CMAKE_EXE_LINKER_FLAGS ${M_LIB})

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
    src/hashmap.c
    src/hashmap.h
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
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs 3 REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/src")

#SET(CMAKE_SWIG_FLAGS "-debug-typemap;-debug-tmsearch;-debug-tmused;-ltypemaps")
SET(CMAKE_SWIG_FLAGS "-ltypemaps;-py3")

SET(CMAKE_SWIG_OUTDIR "${CMAKE_BINARY_DIR}/python")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES C ON)
SET_SOURCE_FILES_PROPERTIES(src/pymcell.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_LIBRARY(pymcell 
  LANGUAGE python
  SOURCES
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
  src/hashmap.c
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
SWIG_LINK_LIBRARIES(pymcell ${CMAKE_CURRENT_BINARY_DIR}/lib/libnfsim_c.so ${CMAKE_CURRENT_BINARY_DIR}/lib/libNFsim.so ${PYTHON_LIBRARIES})

# copy the pyMCell test script into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/lipidraft.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/pymcell_unittests.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/test_dm_read.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/test_new_api.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/test_api.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
file(COPY ${CMAKE_SOURCE_DIR}/src/torus.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python)
file(COPY ${CMAKE_SOURCE_DIR}/src/box.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python)

# copy the pyMCell Boolean logic test script into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test_bool.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

# copy the pyMCell list release test script into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test_release_list.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

# copy the pyMCell Release pattern test script into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test_release_pattern.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

# copy the mcell neuron interface test into place
file(COPY ${CMAKE_SOURCE_DIR}/src/test_mcell_neuron_python.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

# copy the Jupyter notebook into place
file(COPY ${CMAKE_SOURCE_DIR}/src/pymcell_test.ipynb DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)

# copy tools for MDLr to MDL
file(COPY ${CMAKE_SOURCE_DIR}/src/rules_py/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# build executable
add_executable(mcell
  ${SOURCE_FILES}
  src/mcell.c
  ${BISON_mdlParser_OUTPUTS}
  ${FLEX_mdlScanner_OUTPUTS})
target_link_libraries(mcell ${M_LIB} nfsim_c NFsim)

#set(CMAKE_EXE_LINKER_FLAGS "${M_LIB} -static-libgcc -static-libstdc++")

#add_executable(mcell_static
#  ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
#  ${SOURCE_FILES}
#  ${BISON_mdlParser_OUTPUTS}
#  ${FLEX_mdlScanner_OUTPUTS})
#target_link_libraries(mcell_static ${M_LIB} ${CMAKE_SOURCE_DIR}/lib/libnfsim_c_static.a ${CMAKE_SOURCE_DIR}/lib/libNFsim_static.a)

target_link_libraries(mcell ${M_LIB})
TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)

# build little libmcell test
#add_executable(libmcell_test
#  ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
#  ${SOURCE_FILES}
#  src/libmcell_test.c
#  ${BISON_mdlParser_OUTPUTS}
#  ${FLEX_mdlScanner_OUTPUTS})
#target_link_libraries(libmcell_test ${M_LIB})

# build nfsim and nfsimCInterface before trying to build MCell
add_custom_target(
  build_nfsim
  COMMAND python ${CMAKE_SOURCE_DIR}/requirements.py) 
add_dependencies(mcell build_nfsim)
add_dependencies(_pymcell build_nfsim)
add_dependencies(mcell version_h)
