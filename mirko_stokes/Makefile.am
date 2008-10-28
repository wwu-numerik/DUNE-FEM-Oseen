# $Id: duneproject 5104 2008-03-06 16:24:04Z robertk $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)

noinst_PROGRAMS = ldg_stokes

ldg_stokes_SOURCES = stokes.cc

#LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
#ldg_stokes_CXXFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS)


ldg_stokes_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS) \ 
	-DALBERTAGRID -DGRIDDIM=2
ldg_stokes_LDADD = $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(GRAPE_LDFLAGS) $(ALL_PKG_LDFLAGS) \
	$(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) -lgr -lglut -lGL -lGLU -lblas

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
