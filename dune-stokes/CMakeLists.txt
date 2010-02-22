PROJECT(dune-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)
#----------------------------------------------------------------------------------------------------
# General Settings
#----------------------------------------------------------------------------------------------------
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
SET( CMAKE_COLOR_MAKEFILE ON)

EXECUTE_PROCESS(
	COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
	OUTPUT_VARIABLE COMMIT
	ERROR_VARIABLE shell_error
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune-stokes
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )

SET( GRIDTYPE
	"ALUGRID_SIMPLEX" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( POLORDER
	"2" CACHE STRING
	"Sigma POLORDER" )

SET( PRESSUREPOLORDER
	${POLORDER} CACHE STRING
	"Pressure POLORDER" )

SET( VELOCITYPOLORDER
	${POLORDER} CACHE STRING
	"Velocity POLORDER" )

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
	"/share/dune/modules/modules_x86_64/alugrid/alugrid-1.14" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET( GRAPE_PATH
	"/share/dune/Modules/modules_x86_64/grape" CACHE STRING
	"GRAPE_PATH" )

SET( CXX_FLAGS
	"-Wall -O3 -fomit-frame-pointer -funroll-loops -g -ggdb -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME" CACHE STRING
	"")

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_DEFINITIONS( -I${UGGRID_BASE_PATH}/include -DENABLE_UG )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

ADD_DEFINITIONS(-DHAVE_CONFIG_H )
ADD_DEFINITIONS(-D${GRIDTYPE} -DGRIDDIM=${GRIDDIM} -D${PROBLEM} -DCOMMIT="${COMMIT}" )
ADD_DEFINITIONS(-DPOLORDER=${POLORDER} -DPRESSURE_POLORDER=${PRESSUREPOLORDER} -DVELOCITY_POLORDER=${VELOCITYPOLORDER} )
ADD_DEFINITIONS(-DINNER_CG_SOLVERTYPE=${INNER_CG_SOLVERTYPE} -DOUTER_CG_SOLVERTYP=${OUTER_CG_SOLVERTYP} )
ADD_DEFINITIONS( -I${GRAPE_PATH}
	-I${ALUGRID_BASE_PATH}/include
	-I${ALUGRID_BASE_PATH}/include/serial
	-I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

#needs to come first to include correct config.h
INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR} )

SET( MODULES common istl grid fem stuff )

FOREACH( MODULE ${MODULES} )
    INCLUDE_DIRECTORIES( ${dune-stokes_SOURCE_DIR}/../dune-${MODULE} )
    LINK_DIRECTORIES(${dune-stokes_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/.libs )
ENDFOREACH(MODULE)
LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" "${GRAPE_PATH}" )
ADD_DEFINITIONS( ${CXX_FLAGS} )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
SET(dune-stokesSrc
	src/dune_stokes.cc
)

FILE( GLOB header "src/*.hh" )
FILE( GLOB stokes "stokes/*.hh" )
FILE( GLOB stuff "../dune-stuff/stuff/*.hh" )

set_source_files_properties( ${stuff} ${stokes} ${header} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_EXECUTABLE(dune_stokes ${dune-stokesSrc} ${header} ${stokes} ${stuff} )

TARGET_LINK_LIBRARIES(dune_stokes "common" "fem" "grid" ${BLAS_LIB} "alugrid" ${UGLIBS} "gr" "GL")

ADD_EXECUTABLE(minimal src/minimal.cc ${header} ${stokes} )
TARGET_LINK_LIBRARIES(minimal "common" "fem" "grid" ${BLAS_LIB} "alugrid" ${UGLIBS} "gr" "GL")

ADD_EXECUTABLE(visualisation src/visualisation.cc ${header} ${stokes} )
TARGET_LINK_LIBRARIES(visualisation "common" "fem" "grid" ${BLAS_LIB} "alugrid" ${UGLIBS} "gr" "GL")
