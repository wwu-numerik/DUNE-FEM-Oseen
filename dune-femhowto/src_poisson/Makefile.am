noinst_PROGRAMS = lagrange
lagrange_SOURCES = lagrange.cc

# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
#Default values:
#GRIDTYPE=YASPGRID
#GRIDTYPE=ALBERTAGRID
GRIDTYPE=ALUGRID_SIMPLEX
GRIDDIM=2
POLORDER=1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

AUTOMAKE_OPTIONS = foreign 1.5

# additional file to clean 
CLEANFILES = dgfparser.log unitcube*.dgf.* solution-*.vtu grid*.vtu manager.*.log
