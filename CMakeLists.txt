PROJECT(dune-fem-oseen)

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
# Options, that can be changed be the User in order to customise dune-fem-oseen
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )
		
SET( GRIDTYPE
	"YASPGRID" CACHE STRING
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
INCLUDE_DIR( ${dune-fem-oseen_BINARY_DIR}  ${dune-fem-oseen_SOURCE_DIR} )
INCLUDE_SYS_DIR( ${dune-fem-oseen_SOURCE_DIR}/oseen/oemsolver/ )
ADD_DEFINITIONS( ${CUSTOM_FLAGS} )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_DEFINITIONS( -DENABLE_UG )
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
	ADD_DEFINITIONS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

IF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )
	INCLUDE_SYS_DIR( 	${ALBERTA_BASE_PATH}/include/alberta )
	LINK_DIRECTORIES( "${ALBERTA_BASE_PATH}/lib"  )
	SET( ALULIBS "alberta_${GRIDDIM}d" "alberta_util" )
	ADD_DEFINITIONS( "-DENABLE_ALBERTA" )
ENDIF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ADD_DEFINITIONS( -DMPIPP_H -DENABLE_MPI=1 )
ENDIF ( ENABLE_PARALELL )

SET( USE_OMP "0" )
IF( ENABLE_OMP )
	ADD_DEFINITIONS( -fopenmp -DUSE_OMP=1)
	LIST( APPEND ParaLIBS "-fopenmp" )
ELSE( ENABLE_OMP )
	ADD_DEFINITIONS( -DUSE_OMP=0)
ENDIF( ENABLE_OMP )

SET( OUTER_CG_SOLVERTYPE "OEM${OUTER_SOLVER}Op" )
SET( INNER_CG_SOLVERTYPE "OEM${INNER_SOLVER}Op" )


CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

ADD_DUNE_MODULES(common istl grid fem geometry localfunctions stuff spgrid)

LINK_DIRECTORIES( "/usr/lib" "${GRAPE_PATH}" ${CCGNU_LIBRARY_DIRS} )

FILE( GLOB header "src/*.hh" )
FILE( GLOB_RECURSE oseen "dune/fem/oseen/*.h*" )

set( COMMON_LIBS "dunestuff" "dunecommon" "dunegeometry"  "dunefem" "dunegrid"  
    ${BLAS_LIB} ${ALULIBS} ${UGLIBS}  ${ParaLIBS} ${CCGNU_LIBRARIES}  ${BOOST_LIBS})
set( COMMON_HEADER ${header} ${oseen} ${stuff} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_EXECUTABLE(oseen src/dune_stokes.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(oseen ${COMMON_LIBS} )


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
