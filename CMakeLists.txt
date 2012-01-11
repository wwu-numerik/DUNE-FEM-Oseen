PROJECT(dune-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(UseDoxygen)

SET (CMAKE_VERBOSE_MAKEFILE "on" )
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

SET( USE_ISTL
	0 CACHE STRING
	"Use fem-wrapped istl matrix" )
		
SET( GRIDTYPE
	"ALUGRID_SIMPLEX" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( POLORDER
	"1" CACHE STRING
	"Sigma POLORDER" )

SET( PRESSUREPOLORDER
	${POLORDER} CACHE STRING
	"Pressure POLORDER" )

SET( VELOCITYPOLORDER
	${POLORDER} CACHE STRING
	"Velocity POLORDER" )

SET( INNER_SOLVER
	"CG" CACHE STRING
	"INNER_CG_SOLVERTYPE" )

SET( OUTER_SOLVER
	"CG" CACHE STRING
	"OUTER_CG_SOLVERTYP" )

SET( PROBLEM_NAMESPACE
	"StokesProblems::Cockburn" CACHE STRING
	"PROBLEM_NAMESPACE" )

#sadly this produces a drop down list only in the qt gui
SET_PROPERTY(CACHE PROBLEM_NAMESPACE PROPERTY STRINGS 
			"StokesProblems::Simple"
			"StokesProblems::Cockburn"
			"StokesProblems::Generalized"
			"StokesProblems::Constant"
			"StokesProblems::Rotate"
			"StokesProblems::Aorta"
			"StokesProblems::TimeDisc"
	)

SET( ALUGRID_BASE_PATH
	"/opt/dune/modules/alugrid/alugrid-1.50" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET( ALBERTA_BASE_PATH
	"/opt/dune/modules/alberta/alberta-2.0" CACHE STRING
	"ALBERTA_BASE_PATH" )

SET( GRAPE_PATH
	"/share/dune/Modules/modules_x86_64/grape" CACHE STRING
	"GRAPE_PATH" )

SET( ENABLE_PARALELL
	0 CACHE BOOL
	"Enable parallel features" )
	
SET( ENABLE_OMP
	0 CACHE BOOL
	"Enable openmp features" )
	
SET ( METIS_DIR
	"/share/dune/Modules/modules_x86_64/ParMetis-3.1.1" CACHE STRING
	"metis toplevel directory" )

#needs to come first to include correct config.h
INCLUDE_DIR( ${dune-stokes_BINARY_DIR}  ${dune-stokes_SOURCE_DIR} )
INCLUDE_SYS_DIR( ${dune-stokes_SOURCE_DIR}/stokes/oemsolver/ )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -DENABLE_UG )
	INCLUDE_SYS_DIR( ${UGGRID_BASE_PATH}/include )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

SET( ENABLE_ALUGRID OFF )
IF( ${GRIDTYPE} MATCHES "ALUGRID" )
	INCLUDE_SYS_DIR( 	${GRAPE_PATH} 
						${ALUGRID_BASE_PATH}/include 
						${ALUGRID_BASE_PATH}/include/serial 
						${ALUGRID_BASE_PATH}/include/duneinterface )
	LINK_DIRECTORIES( "${ALUGRID_BASE_PATH}/lib"  )
	SET( ALULIBS "alugrid" )
	SET( ENABLE_ALUGRID ON )
	ADD_CXX_FLAGS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

IF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )
	INCLUDE_SYS_DIR( 	${ALBERTA_BASE_PATH}/include/alberta )
	LINK_DIRECTORIES( "${ALBERTA_BASE_PATH}/lib"  )
	SET( ALULIBS "alberta_${GRIDDIM}d" "alberta_util" )
	ADD_CXX_FLAGS( "-DENABLE_ALBERTA" )
ENDIF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ADD_CXX_FLAGS( -DMPIPP_H -DENABLE_MPI=1 )
ENDIF ( ENABLE_PARALELL )

SET( USE_OMP "0" )
IF( ENABLE_OMP )
	ADD_CXX_FLAGS( -fopenmp -DUSE_OMP=1)
	LIST( APPEND ParaLIBS "-fopenmp" )
ELSE( ENABLE_OMP )
	ADD_CXX_FLAGS( -DUSE_OMP=0)
ENDIF( ENABLE_OMP )

IF( USE_ISTL )
	SET( OUTER_CG_SOLVERTYPE "ISTL${OUTER_SOLVER}Op" )
	SET( INNER_CG_SOLVERTYPE "ISTL${INNER_SOLVER}Op" )
ELSE( USE_ISTL )
	SET( OUTER_CG_SOLVERTYPE "OEM${OUTER_SOLVER}Op" )
	SET( INNER_CG_SOLVERTYPE "OEM${INNER_SOLVER}Op" )
ENDIF( USE_ISTL )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

ADD_DUNE_MODULES( common istl grid fem geometry )
ADD_MY_MODULES( stuff )

LINK_DIRECTORIES( "/usr/lib" "${GRAPE_PATH}" ${CCGNU_LIBRARY_DIRS} )

FILE( GLOB header "src/*.hh" )
FILE( GLOB_RECURSE stokes "stokes/*.hh" )

set( COMMON_LIBS "dunecommon" "dunegeometry"  "dunefem" "dunegrid" ${BLAS_LIB} ${ALULIBS} ${UGLIBS}  ${ParaLIBS} ${CCGNU_LIBRARIES} "boost_filesystem" "boost_system")
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_EXECUTABLE(dune_stokes src/dune_stokes.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(dune_stokes ${COMMON_LIBS} )

#legacy dbg binary
#ADD_EXECUTABLE(scaling_test src/scaling_test.cc ${COMMON_HEADER} )
#TARGET_LINK_LIBRARIES(scaling_test ${COMMON_LIBS} )

# ADD_DUNE_EXECUTABLE(minimal src/minimal.cc ${COMMON_HEADER} ${COMMON_LIBS} )
ADD_EXECUTABLE(minimal src/minimal.cc ${COMMON_HEADER}  )
TARGET_LINK_LIBRARIES(minimal ${COMMON_LIBS} )

ADD_EXECUTABLE(visualisation src/visualisation.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(visualisation ${COMMON_LIBS} )

HEADERCHECK( ${header} ${stokes} )
ADD_CPPCHECK( src/dune_stokes.cc ${COMMON_HEADER} )
ADD_DEFINITIONS( ${MY_CXX_FLAGS} )

SET( PROBLEMS Simple Cockburn Generalized  TimeDisc )
SET( INCOMPLETEPROBLEMS Rotate Constant )
SET( PROBLEMS3D Aorta )


FOREACH( problem ${PROBLEMS} )
	SET( targetName so_${problem} )
	ADD_EXECUTABLE(${targetName} src/dune_stokes.cc ${COMMON_HEADER} )
	TARGET_LINK_LIBRARIES( ${targetName} ${COMMON_LIBS} )
	set_target_properties(${targetName} PROPERTIES COMPILE_FLAGS -DPROBLEM_NAMESPACE=StokesProblems::${problem} )
ENDFOREACH( problem )