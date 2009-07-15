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
SET( CXX_FLAGS "-Wall -O3 -fomit-frame-pointer -funroll-loops" CACHE STRING
 	"" )

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune-stokes
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )

SET( GRIDTYPE
	"UGGRID" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( POLORDER
	"2" CACHE STRING
	"POLORDER" )

SET( INNER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"INNER_CG_SOLVERTYPE" )
	
SET( OUTER_CG_SOLVERTYP
	"OEMCGOp" CACHE STRING
	"OUTER_CG_SOLVERTYP" )

SET( PROBLEM
	"COCKBURN_PROBLEM" CACHE STRING
	"PROBLEM" )
	
SET( ALUGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ALUGrid-1.1_Parallel" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_DEFINITIONS( -I${UGGRID_BASE_PATH}/include -DENABLE_UG )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

ADD_DEFINITIONS(-DHAVE_CONFIG_H -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME -DUSE_ALTERNATIVE_SOLVER )
ADD_DEFINITIONS(-D${GRIDTYPE} -DGRIDDIM=${GRIDDIM} -DPOLORDER=${POLORDER} -D${PROBLEM})
ADD_DEFINITIONS(-DINNER_CG_SOLVERTYPE=${INNER_CG_SOLVERTYPE} -DOUTER_CG_SOLVERTYP=${OUTER_CG_SOLVERTYP} )
ADD_DEFINITIONS( 
	-I${ALUGRID_BASE_PATH}/include 
	-I${ALUGRID_BASE_PATH}/include/serial 
	-I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )


SET( MODULES common istl grid fem stuff )

FOREACH( MODULE ${MODULES} )
    INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR}/../dune-${MODULE} )
    LINK_DIRECTORIES(${dune-stokes_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/.libs )
ENDFOREACH(MODULE)
INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR} )
LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
SET(dune-stokesSrc 
	src/dune_stokes.cc
)

ADD_EXECUTABLE(dune_stokes ${dune-stokesSrc} )

TARGET_LINK_LIBRARIES(dune_stokes "common" "fem" "grid" ${BLAS_LIB} "alugrid" ${UGLIBS} )
