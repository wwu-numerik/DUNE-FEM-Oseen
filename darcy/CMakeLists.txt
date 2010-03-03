PROJECT(darcy)

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

EXECUTE_PROCESS(
	COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
	OUTPUT_VARIABLE COMMIT
	ERROR_VARIABLE shell_error
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

#----------------------------------------------------------------------------------------------------
# Options, that can be changed by the User
#----------------------------------------------------------------------------------------------------
# Add Definitions, Compiler-Switches, etc.: -Wall -O2 -g3 ...
SET( CXX_FLAGS
	"-Wall -O3 -fomit-frame-pointer -funroll-loops" CACHE STRING
 	"CXX_FLAGS" )

SET( BLAS_LIB
	"blas" CACHE STRING
	"BLAS_LIB" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( GRIDTYPE
	"ALUGRID_COMFORM" CACHE STRING
	"GRIDTYPE" )

SET( MACRO_POLORDER
	"1" CACHE STRING
	"MACRO_POLORDER" )

SET( MACRO_PROBLEM
	"MACRO_PROBLEM" CACHE STRING
	"MACRO_PROBLEM" )

SET( MACRO_CG_SOLVERTYPE
	"CGInverseOp" CACHE STRING
	"MACRO_CG_SOLVERTYPE" )

SET( MICRO_GRIDTYPE
	"ALUGRID_COMFORM" CACHE STRING
	"MICRO_GRIDTYPE" )

SET( MICRO_POLORDER
	"1" CACHE STRING
	"MICRO_POLORDER" )

SET( MICRO_PROBLEM
	"MICRO_PROBLEM" CACHE STRING
	"MICRO_PROBLEM" )

SET( MICRO_INNER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"MICRO_INNER_CG_SOLVERTYPE" )

SET( MICRO_OUTER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"MICRO_OUTER_CG_SOLVERTYPE" )

SET( ALUGRID_BASE_PATH
	"/data/dune_work/private/f_albr01/modules/alugrid-1.14" CACHE STRING
	"ALUGRID_BASE_PATH" )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ADD_DEFINITIONS(-DGRIDDIM=${GRIDDIM} -D${GRIDTYPE})
ADD_DEFINITIONS(-DMACRO_POLORDER=${MACRO_POLORDER} -DMACRO_PROBLEM=${MACRO_PROBLEM})
ADD_DEFINITIONS(-DMACRO_CG_SOLVERTYPE=${MACRO_CG_SOLVERTYPE})
ADD_DEFINITIONS(-DMICRO_GRIDTYPE=${MICRO_GRIDTYPE})
ADD_DEFINITIONS(-DMICRO_POLORDER=${MICRO_POLORDER} -DMICRO_PROBLEM=${MICRO_PROBLEM})
ADD_DEFINITIONS(-DMICRO_INNER_CG_SOLVERTYPE=${MICRO_INNER_CG_SOLVERTYPE} -DMICRO_OUTER_CG_SOLVERTYPE=${MICRO_OUTER_CG_SOLVERTYPE} -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME -DUSE_ALTERNATIVE_SOLVER)
ADD_DEFINITIONS(
	-I${ALUGRID_BASE_PATH}/include
	-I${ALUGRID_BASE_PATH}/include/serial
	-I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

#needs to come first to include correct config.h
INCLUDE_DIRECTORIES( ${darcy_SOURCE_DIR/../dune-stokes} )

SET( MODULES common istl grid fem stuff stokes )

FOREACH( MODULE ${MODULES} )
    INCLUDE_DIRECTORIES( ${darcy_SOURCE_DIR}/../../dune-${MODULE} )
    LINK_DIRECTORIES(${darcy_SOURCE_DIR}/../../dune-${MODULE}/${MODULE}/.libs )
ENDFOREACH(MODULE)
LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" )
ADD_DEFINITIONS( ${CXX_FLAGS} )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
SET(darcySrc
	darcy.cc
)

ADD_EXECUTABLE(darcy ${darcySrc} )

TARGET_LINK_LIBRARIES(darcy "common" "fem" "grid" ${BLAS_LIB} "alugrid" )
