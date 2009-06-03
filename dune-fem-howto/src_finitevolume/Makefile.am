GRIDTYPE=ALUGRID_SIMPLEX
GRIDDIM=3

CXX=g++-4.2

# programs to compile 
noinst_PROGRAMS = finitevolume

finitevolume_SOURCES = finitevolume.cc 
finitevolume_CXXFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)
finitevolume_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) 

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used

CLEANFILES = dgfparser.log unitcube*.dgf.* *.pvtu solution*.vtu grid*.vtu solution-*.vtp ParaView*.pvs manager.*.log

CXXFLAGS = -g -Wall -DNDEBUG
