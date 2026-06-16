cmake_minimum_required(VERSION 3.25...4.3)
project(mcell)

if ("${PYTHON_VERSION}" STREQUAL "")
  set(PYTHON_VERSION 3.13)
endif()
string(REPLACE "\." "" PYTHON_VERSION_NO_DOT ${PYTHON_VERSION})
set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})

find_package(Python3 ${PYBIND11_PYTHON_VERSION} EXACT REQUIRED COMPONENTS Interpreter Development)
message("-- Python3 include dirs are ${Python3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

# All project-wide flags, definitions, the C++ standard and the platform/compiler
# branches that used to be appended to CMAKE_CXX_FLAGS are centralised in the
# mcell_project_options INTERFACE target (defined below) and consumed by every
# first-party and bundled target.

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
# Request only the modules MCell actually uses (the rest come transitively):
#   CommonCore        - vtkSmartPointer
#   CommonDataModel   - vtkPolyData / vtkPointData / vtkTriangle
#   CommonTransforms  - vtkTransform
#   FiltersCore       - vtkCleanPolyData / vtkTriangleFilter / vtkMassProperties / ...
#   FiltersModeling   - vtkCollisionDetectionFilter / vtkSelectEnclosedPoints
#   FiltersSources    - vtkPlatonicSolidSource
#   IOExport          - vtkOBJExporter
#   RenderingCore     - vtkRenderer / vtkRenderWindow / vtkPolyDataMapper
# After this, ${VTK_LIBRARIES} expands to exactly these VTK:: targets, which
# carry their own include dirs (no manual include paths needed).
find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  CommonTransforms
  FiltersCore
  FiltersModeling
  FiltersSources
  IOExport
  RenderingCore)

# jsoncpp 1.9.8 defaults the BUILD_SHARED_LIBS option to ON. Because jsoncpp is
# pulled in via add_subdirectory(), that cache variable would leak into the rest
# of the project and turn libraries that do not specify STATIC/SHARED explicitly
# (e.g. nauty in libbng) into shared libraries, which then fail to link
# (undefined reference to updatecan). jsoncpp 1.9.3 defaulted this to OFF, so
# force it OFF before adding jsoncpp to preserve the previous static-link behavior.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# jsoncpp is added before mcell's global -fPIC flags are configured below, so its
# targets do not inherit them. jsoncpp 1.9.3 set POSITION_INDEPENDENT_CODE ON on
# its static target explicitly; jsoncpp 1.9.8 does not do this for the static
# target. Because the static jsoncpp is linked into the shared module mcell.so, it
# must be position independent, otherwise linking fails with "relocation ... can
# not be used when making a shared object; recompile with -fPIC". Force PIC here.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# jsoncpp_test cannot be built with gperftools linking options
add_subdirectory(${CMAKE_SOURCE_DIR}/libs/jsoncpp ${CMAKE_CURRENT_BINARY_DIR}/libs/jsoncpp)

# jsoncpp 1.9.3 exposes a single static target named "jsoncpp_lib", while
# jsoncpp 1.9.8 names the static target "jsoncpp_static" (jsoncpp_lib only exists
# when BUILD_SHARED_LIBS is ON). Pick whichever exists so both versions work.
if(TARGET jsoncpp_lib)
  set(JSONCPP_LIB jsoncpp_lib)
else()
  set(JSONCPP_LIB jsoncpp_static)
endif()
# Give jsoncpp a stable namespaced name so the rest of the tree links
# JsonCpp::JsonCpp regardless of which underlying target the vendored version
# exposes (jsoncpp_lib vs jsoncpp_static).
if(NOT TARGET JsonCpp::JsonCpp)
  add_library(JsonCpp::JsonCpp ALIAS ${JSONCPP_LIB})
endif()

# ===========================================================================
# Project-wide build options as a single INTERFACE target.
#
# Everything that used to be appended to CMAKE_C_FLAGS / CMAKE_CXX_FLAGS lives
# here as generator expressions (per language / config / platform / compiler).
# Every first-party target AND the bundled siblings link this PRIVATE (see the
# loop after the subdirectories), so the whole build shares one flag set.
# ===========================================================================
add_library(mcell_project_options INTERFACE)

# Single project-wide C++ standard knob. Bump to 20 or 23 to experiment; CMake
# emits -std=gnu++<N> automatically once <N> exceeds the compiler default. The
# "max" of this and any lower per-target floor wins, so the bundled nfsim/libbng
# get raised to the same standard. (Note: a value the compiler does not yet
# advertise -- e.g. 26 on GCC 12 -- is silently ignored, leaving the default.)
set(MCELL_CXX_STANDARD 17 CACHE STRING "C++ standard for the whole MCell build (17/20/23)")
# The feature gives a configure-time guarantee the compiler supports the standard
# (and is the floor that overrides the siblings' lower cxx_std_* requirements).
target_compile_features(mcell_project_options INTERFACE cxx_std_${MCELL_CXX_STANDARD})
# ...and we *also* pass -std=gnu++<N> explicitly so the standard is always visible
# on the compile command line, even when it matches the compiler default (GCC's
# default is already gnu++17, so the feature alone would emit no flag). We use the
# gnu++ dialect (extensions on) to preserve the historical behaviour the legacy
# C-compiled-as-C++ code relies on.
target_compile_options(mcell_project_options INTERFACE
    $<$<AND:$<NOT:$<PLATFORM_ID:Windows>>,$<COMPILE_LANGUAGE:CXX>>:-std=gnu++${MCELL_CXX_STANDARD}>)

target_compile_definitions(mcell_project_options INTERFACE
    PACKAGE_BUGREPORT="mcell@salk.edu"
    PYTHON_VERSION=${PYTHON_VERSION_NO_DOT}
    $<$<NOT:$<PLATFORM_ID:Windows>>:_GNU_SOURCE=1>
    $<$<PLATFORM_ID:Windows>:YY_NO_UNISTD_H>)

# Common GCC/Clang warning set (kept identical to the historical global flags).
set(MCELL_NIX_WARNINGS
    -Wall -Wextra
    -Wno-unused-parameter -Wno-unused-function
    -Wno-unused-variable -Wno-deprecated-declarations)

target_compile_options(mcell_project_options INTERFACE
    # warnings (non-Windows)
    $<$<NOT:$<PLATFORM_ID:Windows>>:${MCELL_NIX_WARNINGS}>
    $<$<AND:$<NOT:$<PLATFORM_ID:Windows>>,$<CXX_COMPILER_ID:Clang>>:-Wno-ignored-optimization-argument>
    # -fpermissive is C++-only and required by the legacy C-compiled-as-C++ code
    $<$<AND:$<NOT:$<PLATFORM_ID:Windows>>,$<COMPILE_LANGUAGE:CXX>>:-fpermissive>
    # the few genuine C sources are C11 (matches the historical -std=c11)
    $<$<AND:$<NOT:$<PLATFORM_ID:Windows>>,$<COMPILE_LANGUAGE:C>>:-std=c11>
    # Release optimisation (historically Release-only)
    $<$<AND:$<NOT:$<PLATFORM_ID:Windows>>,$<CONFIG:Release>>:-O3;-g>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Release>>:-finline-limit=10000>
    # MSVC equivalents of the old if(MSVC) block.
    $<$<CXX_COMPILER_ID:MSVC>:/EHsc>
    # The C++ standard on MSVC comes from the cxx_std_<N> feature above (=>
    # /std:c++17). We deliberately do NOT force /std:c++latest: on modern MSVC it
    # implies /permissive- -> /Zc:strictStrings, which rejects the legacy
    # string-literal-to-char* code in the vendored nauty (and nfsim). Turning
    # strict strings off is the MSVC counterpart of the -fpermissive we pass to
    # GCC/Clang for the same legacy code.
    $<$<CXX_COMPILER_ID:MSVC>:/Zc:strictStrings->
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/O2>)

# The bundled mcell/libs headers were exposed via a global -isystem; keep that
# as a SYSTEM usage requirement so consumers' warnings stay quiet.
target_include_directories(mcell_project_options SYSTEM INTERFACE
    ${CMAKE_SOURCE_DIR}/libs)

# Debuggable Debug builds, cross-platform:
#   * GCC (Linux) and Clang/AppleClang (macOS): -O0 disables optimisation so
#     every local variable stays live and stepping is linear, and -g emits the
#     symbols gdb/lldb need. We deliberately use -O0 rather than -Og: under
#     Clang/lldb (macOS) -Og still optimises enough to drop variables and make
#     lldb warn "compiled with optimization - variables may not be available".
#     The cost is a slower Debug binary, which is the expected Debug tradeoff;
#     override with -DCMAKE_CXX_FLAGS_DEBUG="-Og -g" if you want a faster one.
#   * MSVC (Windows): CMake's default CMAKE_*_FLAGS_DEBUG already carries
#     /Zi /Od /RTC1 -- a full debug build -- so we leave it untouched. (NOT MSVC
#     also excludes clang-cl, which uses the MSVC-style driver, not -O0/-g.)
if (NOT MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
endif()

if (NOT WIN32)
  # SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 ")
  # SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 ")

  # Release-only LINK option: the mcell executable is linked -no-pie (avoids an
  # AVX2 / MCell3 wall-overlap miscompile; see git history). It is a link option,
  # so it stays on CMAKE_EXE_LINKER_FLAGS rather than mcell_project_options.
  # (Compile-side optimisation -O3/-g/-finline-limit now lives in that target.)
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie ")
  endif()

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
    # or
    # LD_PRELOAD='/mnt/code/src4/mcell/build/install_gperftools/lib/libprofiler.so' CPUPROFILE=./prof.out python model.py
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

# ---------------------------------------------------------------------------
# Build-flag strings for the `mcell --fullversion` banner.
#
# These are assembled here, where the real project flags are defined, and
# emitted to deps/version_flags.h via file(GENERATE) -- the same mechanism on
# Linux, macOS and Windows. This replaces the old fragile shell flag-splitting
# (and the separate -CFlags/-LdFlags plumbing into version.sh / version.ps1),
# and also fixes MCELL_LDFLAGS, which used to be empty on Linux/macOS.
# ---------------------------------------------------------------------------
if (WIN32)
  set(_banner_cflags /EHsc /std:c++${MCELL_CXX_STANDARD} /Zc:strictStrings-
      -DPYTHON_VERSION=${PYTHON_VERSION_NO_DOT} /O2)
else()
  set(_banner_cflags ${MCELL_NIX_WARNINGS}
      -fpermissive -std=gnu++${MCELL_CXX_STANDARD} -D_GNU_SOURCE=1
      -DPYTHON_VERSION=${PYTHON_VERSION_NO_DOT})
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    list(APPEND _banner_cflags -O3 -g)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      list(APPEND _banner_cflags -finline-limit=10000)
    endif()
  elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Mirror the per-config optimisation flags CMake actually injects for Debug
    # (CMAKE_CXX_FLAGS_DEBUG is set to "-Og -g" for GCC/Clang above). Without
    # this the banner would under-report a Debug build.
    separate_arguments(_banner_debug_flags UNIX_COMMAND "${CMAKE_CXX_FLAGS_DEBUG}")
    list(APPEND _banner_cflags ${_banner_debug_flags})
  endif()
endif()
list(JOIN _banner_cflags " " MCELL_BANNER_CFLAGS)
string(STRIP "${CMAKE_EXE_LINKER_FLAGS}" MCELL_BANNER_LDFLAGS)

file(GENERATE
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deps/version_flags.h
  CONTENT "/* Auto-generated by CMake -- do not edit. Build flags for --fullversion. */
#pragma once
#define MCELL_CFLAGS  \"${MCELL_BANNER_CFLAGS}\"
#define MCELL_LDFLAGS \"${MCELL_BANNER_LDFLAGS}\"
#define MCELL_LFLAGS  \"-Cema\"
#define MCELL_YFLAGS  \"\"
")
# Include directories are now expressed as per-target usage requirements:
# - each library (libbng, nfsim, mcell4, libmcell, ...) exports its own headers,
# - VTK headers come from linking the VTK:: targets,
# so the mcell executable and the mcell4_so module only need to spell out the
# handful of dirs unique to their own translation units (see APP_INCLUDE_DIRS
# near the target definitions below).

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



# create version.h target (build environment / tool / git info only; the build
# flags now live in version_flags.h generated above via file(GENERATE)).
if (NOT WIN32)
  add_custom_target(
    version_h
    COMMAND ${CMAKE_SOURCE_DIR}/src/version.sh > ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/version.sh
    VERBATIM)
elseif (WIN32)
  add_custom_target(
    version_h
    COMMAND powershell -ExecutionPolicy Bypass -File ${CMAKE_SOURCE_DIR}/src/version.ps1 -OutputFile ${CMAKE_CURRENT_BINARY_DIR}/deps/version.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/version.ps1
    VERBATIM)
endif()


# linking on centos7 with 7.3.1 fails with undefined reference to `casc::detail::asc_NodeBase::~asc_NodeBase()'
set(PYBIND11_LTO_CXX_FLAGS "")

if (WIN32)
  # not compiling for now
  # add_subdirectory(${CMAKE_SOURCE_DIR}/libs/libunistd ${CMAKE_CURRENT_BINARY_DIR}/libs/libunistd)
  # include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/libunistd/unistd)
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/libs/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/libs/pybind11)
add_subdirectory(${CMAKE_SOURCE_DIR}/../libbng/ ${CMAKE_CURRENT_BINARY_DIR}/libs/libbng)

add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsim ${CMAKE_CURRENT_BINARY_DIR}/nfsim)
add_subdirectory(${CMAKE_SOURCE_DIR}/../nfsimCInterface ${CMAKE_CURRENT_BINARY_DIR}/nfsimCInterface)

add_subdirectory(${CMAKE_SOURCE_DIR}/libmcell ${CMAKE_CURRENT_BINARY_DIR}/libmcell)
add_subdirectory(${CMAKE_SOURCE_DIR}/src4 ${CMAKE_CURRENT_BINARY_DIR}/src4)

# set default locations of libraries and executables, must be done after add_subdirectory setups
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

# may be empty
link_directories(${SYSTEM_PYTHON_LIB_DIR})

if (UNIX)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN/./lib")
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# Include dirs needed by the mcell executable's and the mcell4_so module's OWN
# translation units (most others arrive transitively from the linked libraries).
# The generated config.h / version.h and the flex/bison output live under deps.
set(APP_INCLUDE_DIRS
    ${CMAKE_CURRENT_BINARY_DIR}/deps
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src4
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/libmcell
    ${CMAKE_SOURCE_DIR}/libs
    ${CMAKE_SOURCE_DIR}/libs/jsoncpp/include
    ${CMAKE_SOURCE_DIR}/../nfsimCInterface/src
    ${CMAKE_SOURCE_DIR}/../libbng
    ${CMAKE_SOURCE_DIR}/../libbng/libs
    ${CMAKE_SOURCE_DIR}/../libbng/libs/sparsehash/src)

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

if (NOT WIN32 AND NOT APPLE)
    set(STDC_FS stdc++fs)
endif()

target_include_directories(mcell PRIVATE ${APP_INCLUDE_DIRS})
TARGET_COMPILE_DEFINITIONS(mcell PRIVATE NOSWIG=1)
add_dependencies(mcell version_h mcell4)
target_link_libraries(mcell
    PRIVATE
    MCell::mcell4 MCell::libmcell_dummy BNG::bng NFsimC::NFsimC NFsim::NFsim JsonCpp::JsonCpp Nauty::nauty ${VTK_LIBRARIES}
    $<$<NOT:$<PLATFORM_ID:Darwin>>:Python3::Python>
    $<$<PLATFORM_ID:Darwin>:pybind11::python_link_helper>
)

add_subdirectory(${CMAKE_SOURCE_DIR}/utils/data_model_to_pymcell)
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/)

add_dependencies(mcell4_so libmcell mcell4 libbng ${JSONCPP_LIB})
target_link_libraries(mcell4_so
    PUBLIC MCell::libmcell
    PRIVATE MCell::mcell4 BNG::bng JsonCpp::JsonCpp Nauty::nauty ${VTK_LIBRARIES} ${GPERFTOOLS_LIB} ${GPERFTOOLS_MCELLSO_DEP} ${STDC_FS}
    PUBLIC MCell::libmcell  # linking dependency issue, we must include libmcell once again to the list of libraries
)

target_include_directories(mcell4_so PRIVATE ${APP_INCLUDE_DIRS})
set_property(TARGET mcell4_so PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(mcell4_so PROPERTIES PREFIX "")
set_target_properties(mcell4_so PROPERTIES OUTPUT_NAME "mcell")
if (WIN32)
    set_target_properties(mcell4_so PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(mcell4_so PROPERTIES SUFFIX ".so") # must be .so on MacOS as well
endif()

# Apply the shared project options to every target in this build -- the
# first-party targets and the bundled siblings (nfsim/libbng/nauty), which used
# to receive these flags through the global CMAKE_*_FLAGS. Linking is PRIVATE so
# the flags do not leak through public interfaces. Guarded by if(TARGET) so the
# Windows-only / optional variants are skipped cleanly.
foreach(_t IN ITEMS
    mcell4 libmcell libmcell_dummy mcell mcell4_so
    data_model_to_pymcell bng_analyzer
    libbng parser_tester_libbng nauty
    NFsim_static NFsim_static_fpic NFsim NFsim_EXE
    nfsim_c nfsim_c_static)
  if (TARGET ${_t})
    target_link_libraries(${_t} PRIVATE mcell_project_options)
  endif()
endforeach()

# nfsim (NFinput.cpp) and its C interface use Boost (algorithm/string, range).
# The only copy in this source tree is the one bundled in libbng/libs; the
# first-party targets (mcell4/libmcell/...) already pick it up via their own
# include dirs, but the nfsim family used to receive it solely from mcell's old
# global include_directories. Wire the bundled Boost to them explicitly here.
# (Linux has a system Boost that masked this regression; macOS does not.)
foreach(_t IN ITEMS
    NFsim_static NFsim_static_fpic NFsim NFsim_EXE
    nfsim_c nfsim_c_static)
  if (TARGET ${_t})
    target_include_directories(${_t} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/../libbng/libs)
  endif()
endforeach()

# VTK 9 object factories register via per-target autoinit headers. Apply it to
# every target that compiles VTK-using translation units so factory-backed
# objects (e.g. renderers) resolve at runtime. MODULES is the set we requested.
vtk_module_autoinit(
  TARGETS mcell4 libmcell data_model_to_pymcell mcell mcell4_so
  MODULES ${VTK_LIBRARIES})
