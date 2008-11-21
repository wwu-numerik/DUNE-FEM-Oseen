# programs to compile 
noinst_PROGRAMS = dune_hello_world gettingstarted dgfparser

# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
GRIDTYPE=ALUGRID_CUBE
GRIDDIM=3

dune_hello_world_SOURCES = dune_hello_world.cc 
dune_hello_world_CXXFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)
dune_hello_world_LDADD = $(ALL_PKG_LDFLAGS) $(MPI_LDFLAGS) 

gettingstarted_SOURCES = gettingstarted.cc
gettingstarted_CXXFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)
gettingstarted_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) 

dgfparser_SOURCES = dgfparser.cc 
dgfparser_CXXFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)
dgfparser_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) 

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used

CLEANFILES = dgfparser.log unitcube*.dgf.* *.pvtu solution-*.vtu grid*.vtu manager.*.log
