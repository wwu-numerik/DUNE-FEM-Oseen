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
SET( CMAKE_EXE_LINKER_FLAGS "-static")

# Add Definitions, Compiler-Switches, etc.: -Wall -O2 -g3 ...
ADD_DEFINITIONS(-Wall -Wno-strict-aliasing -O0 -DPOLORDER=1 )
#----------------------------------------------------------------------------------------------------
# Load needed Modules
#----------------------------------------------------------------------------------------------------
# PKG-Config
FIND_PACKAGE( PkgConfig )

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune-stokes
#----------------------------------------------------------------------------------------------------
SET( OPTION_TORRENT_SYSTEM
	1 CACHE BOOL
	"Enables the builtin map/mod p2p downloader (requires libtorrent-rasterbar library)" )
SET( OPTION_SOUND
	1 CACHE BOOL
	"Enables sound notification functionality (requires SDL and SDL_sound libraries)" )
SET( OPTION_TRANSLATION_SUPPORT
	1 CACHE BOOL
	"Enables translation support to the programs and adds facilities for helping translators (requires GNU Gettext)" )

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


#ADD_CUSTOM_TARGET(test_susynclib ALL COMMAND tools/test-susynclib.awk src/springunitsynclib.cpp WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

#----------------------------------------------------------------------------------------------------
# Build target defintions
#----------------------------------------------------------------------------------------------------

# Here we define the executable dune-stokes ( or on Windows dune-stokes.exe )

#auto-remove whitespaces before/after lib paths
#cmake_policy(SET CMP0004 OLD)


ADD_EXECUTABLE(dune_stokes ${dune-stokesSrc} )

TARGET_LINK_LIBRARIES(dune_stokes "common" "fem" "grid" "cblas" "f77blas" "atlas")
