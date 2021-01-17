cmake_minimum_required(VERSION 2.8.4)

project(nauty)

if (NOT WIN32)
  set(CMAKE_CXX_FLAGS 
	"${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-sign-compare -Wno-write-strings")
endif()	

# compiling only sources needed for function Trace 
set(SOURCE_FILES
	nauty.c
	traces.c
	schreier.c
	nausparse.c
	naututil.c
	naurng.c
	gtools.c
	nautil.c
)

# compile everything as C++

SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILES} PROPERTIES LANGUAGE CXX )

add_library(${PROJECT_NAME} ${SOURCE_FILES})