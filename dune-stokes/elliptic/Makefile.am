# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  

#Default values:
GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = UGGRID
GRIDDIM = 3
POLORDER = 1
#PROBLEM = POISSON
PROBLEM = AORTA

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

noinst_PROGRAMS = elliptic
elliptic_SOURCES = elliptic.cc

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -D$(PROBLEM) -Wall -O0 \
	--param max-inline-insns-single=1000 --param large-function-growth=1500 \
	-I../../dune-stuff/  -DDEBUG -g2
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -Wall -g --param max-inline-insns-single=1000 --param large-function-growth=1500

include $(top_srcdir)/am/global-rules
