# IGNORE_HEADERS

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) 
noinst_PROGRAMS = dgtest datadisp 
#euler_SOURCES = euler.cc
dgtest_SOURCES = dgtest.cc
datadisp_SOURCES = datadisp.cc

if BUILD_DOCS
SUBDIRS = doc
endif

GRIDTYPE = ALUGRID_SIMPLEX
GRIDDIM = 2
POLORDER = 1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE)
DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

#CXXFLAGS = -O3 -DNDEBUG -funroll-loops -finline-functions --param max-inline-insns-single=4000 --param large-function-growth=4500 --param inline-unit-growth=4000 -ffast-math -fomit-frame-pointer -march=opteron -msse3 -mfpmath=sse # -Wall
#CXXFLAGS = -g -Wall
