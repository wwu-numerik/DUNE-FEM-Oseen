PROJECT(dune-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(UseDoxygen)

FIND_PACKAGE(Boost 1.35.0 REQUIRED)
INCLUDE_DIR(${Boost_INCLUDE_DIR})
FIND_PACKAGE( PkgConfig )
pkg_check_modules( CCGNU REQUIRED libccgnu2 )
ADD_CXX_FLAGS( "${CCGNU_CFLAGS}" )

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

SET( OUTER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"OUTER_CG_SOLVERTYP" )

SET( PROBLEM
	"COCKBURN_PROBLEM" CACHE STRING
	"PROBLEM" )

SET( DIRICHLET_DATA 
	"DirichletData" CACHE STRING
	"DIRICHLET_DATA implementation for default dirichlet traits" )

SET( ALUGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ALUGrid-1.22" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET( GRAPE_PATH
	"/share/dune/Modules/modules_x86_64/grape" CACHE STRING
	"GRAPE_PATH" )

SET( CXX_FLAGS
	-Wall -O3 -fomit-frame-pointer -DDNDEBUG- funroll-loops -g -ggdb -DENABLE_ADAPTIVE=1 -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME CACHE STRING
	"")

SET( ENABLE_PARALELL
	0 CACHE BOOL
	"Enable parallel features" )

SET ( METIS_DIR
	"/share/dune/Modules/modules_x86_64/ParMetis-3.1.1" CACHE STRING
	"metis toplevel directory" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

#needs to come first to include correct config.h
INCLUDE_DIR( ${dune-stokes_BINARY_DIR} ${dune-stokes_SOURCE_DIR} )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -I${UGGRID_BASE_PATH}/include -DENABLE_UG )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

INCLUDE_DIR( ${GRAPE_PATH}
	${ALUGRID_BASE_PATH}/include
	${ALUGRID_BASE_PATH}/include/serial
	${ALUGRID_BASE_PATH}/include/duneinterface )

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	ADD_CXX_FLAGS( -I${ALUGRID_BASE_PATH}/include/parallel
					-DMPIPP_H -DENABLE_MPI=1
					)
ENDIF ( ENABLE_PARALELL )

ADD_DUNE_MODULES( common istl grid fem stuff )

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" "${GRAPE_PATH}" ${CCGNU_LIBRARY_DIRS} )
ADD_CXX_FLAGS( ${CXX_FLAGS} )

FILE( GLOB header "src/*.hh" )
FILE( GLOB_RECURSE stokes "stokes/*.hh" )
FILE( GLOB stuff "../dune-stuff/stuff/*.hh" )

set_source_files_properties( ${DUNE_HEADERS} ${stuff} ${stokes} ${header} PROPERTIES HEADER_FILE_ONLY 1 )

set( COMMON_LIBS "common" "fem" "grid" ${BLAS_LIB} "alugrid" ${UGLIBS}  "GL" ${ParaLIBS} )
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${DUNE_HEADERS} )

ADD_EXECUTABLE(dune_stokes src/dune_stokes.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(dune_stokes ${COMMON_LIBS} )

ADD_EXECUTABLE(scaling_test src/scaling_test.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(scaling_test ${COMMON_LIBS} )

ADD_DUNE_EXECUTABLE(minimal src/minimal.cc ${COMMON_HEADER} )
# TARGET_LINK_LIBRARIES(minimal ${COMMON_LIBS} )

ADD_EXECUTABLE(visualisation src/visualisation.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(visualisation ${COMMON_LIBS} )

HEADERCHECK( ${header} ${stokes} )