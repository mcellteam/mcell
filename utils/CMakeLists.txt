cmake_minimum_required(VERSION 2.8.4)

SET(BUILD_UTILS_DIR ${CMAKE_BINARY_DIR}/utils)

add_subdirectory(bng_analyzer)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bng_analyzer_py/bng_analysis_example.py ${BUILD_UTILS_DIR}/bng_analysis_example.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bng_analyzer_py/viz_reader.py ${BUILD_UTILS_DIR}/viz_reader.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bng_species_reader_py/bng_species_reader_example.py ${BUILD_UTILS_DIR}/bng_species_reader_example.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bng_species_reader_py/species_reader.py ${BUILD_UTILS_DIR}/species_reader.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bngl/generate_bngl_function.py ${BUILD_UTILS_DIR}/generate_bngl_function.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mcell4_runner/bng_runner.py ${BUILD_UTILS_DIR}/bng_runner.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mcell4_runner/mcell4_runner.py ${BUILD_UTILS_DIR}/mcell4_runner.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mcell4_runner/mcell3_runner.py ${BUILD_UTILS_DIR}/mcell3_runner.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mcell4_runner/run_remotely.sh ${BUILD_UTILS_DIR}/run_remotely.sh COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_averages_into_images.py ${BUILD_UTILS_DIR}/plot_averages_into_images.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_single_run.py ${BUILD_UTILS_DIR}/plot_single_run.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_multicolumn_dat_file.py ${BUILD_UTILS_DIR}/plot_multicolumn_dat_file.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/visualize/visualize.sh ${BUILD_UTILS_DIR}/visualize.sh COPYONLY)

