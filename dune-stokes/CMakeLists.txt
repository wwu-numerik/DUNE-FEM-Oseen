PROJECT(dune-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)
#----------------------------------------------------------------------------------------------------
# General Settings
#----------------------------------------------------------------------------------------------------
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
SET( CMAKE_CXX_COMPILER "g++-4.2" )
SET( CMAKE_COLOR_MAKEFILE ON)
SET( CMAKE_VERBOSE_MAKEFILE ON)

# Add Definitions, Compiler-Switches, etc.: -Wall -O2 -g3 ...
ADD_DEFINITIONS(-Wall -O0 -DPOLORDER=1 )

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune-stokes
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"DESCRIPTION" )

SET( MODULES common istl grid fem stuff )

FOREACH( MODULE ${MODULES} )
    INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR}/../dune-${MODULE} )
    LINK_DIRECTORIES(${dune-stokes_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/.libs )
ENDFOREACH(MODULE)
INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR} )
LINK_DIRECTORIES( "/usr/lib" )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
SET(dune-stokesSrc 
	src/dune_stokes.cc
)

ADD_EXECUTABLE(dune_stokes ${dune-stokesSrc} )

TARGET_LINK_LIBRARIES(dune_stokes "common" "fem" "grid" ${BLAS_LIB} )
