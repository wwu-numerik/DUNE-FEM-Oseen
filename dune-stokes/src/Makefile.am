
# own defines
GRIDDIM=2
POLORDER=1
GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=SGRID

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE)

noinst_PROGRAMS = dune_stokes

dune_stokes_SOURCES = dune_stokes.cc

dune_stokes_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(AM_CPPFLAGS)
dune_stokes_LDADD = $(MPI_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) \
$(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GRAPE_LDFLAGS) $(GRAPE_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
