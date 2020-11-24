cmake_minimum_required(VERSION 2.8.4)

SET(BUILD_UTILS_DIR ${CMAKE_BINARY_DIR}/utils)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bng_analysis_example.py ${BUILD_UTILS_DIR}/bng_analysis_example.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/viz_reader.py ${BUILD_UTILS_DIR}/viz_reader.py COPYONLY)
