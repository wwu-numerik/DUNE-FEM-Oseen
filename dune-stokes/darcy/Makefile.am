# #### global ######################################################################################

# grid dimension
# don't take ALUGRID_CONFORM if griddim != 2
GRIDDIM=2

# grid type
GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_CONFORM
#GRIDTYPE=YASPGRID
#GRIDTYPE=SGRID
#GRIDTYPE=UGGRID
#GRIDTYPE=ALBERTAGRID

# #### end global ##################################################################################


# #### macro problem ###############################################################################

# polorder
MACRO_POLORDER=1

# problem
MACRO_PROBLEM=MACRO_PROBLEM

# solver
MACRO_CG_SOLVERTYPE=CGInverseOp

# #### end macro problem ###########################################################################


# #### micro problem ###############################################################################

# grid
MICRO_GRIDTYPE=UGGrid

# polorder
MICRO_POLORDER=1

# problem
MICRO_PROBLEM=MICRO_PROBLEM

# solver
MICRO_INNER_CG_SOLVERTYPE=OEMCGOp
MICRO_OUTER_CG_SOLVERTYPE=OEMCGOp

# #### end micro problem ###########################################################################


AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) $(ALBERTA_CPPFLAGS)\
    -DMACRO_POLORDER=$(MACRO_POLORDER) -D$(MACRO_PROBLEM) -DMACRO_CG_SOLVERTYPE=$(MACRO_CG_SOLVERTYPE) \
    -DMICRO_GRIDTYPE=$(MICRO_GRIDTYPE) -DMICRO_POLORDER=$(MICRO_POLORDER) -D$(MICRO_PROBLEM) -DINNER_CG_SOLVERTYPE=$(MICRO_INNER_CG_SOLVERTYPE) -DOUTER_CG_SOLVERTYP=$(MICRO_OUTER_CG_SOLVERTYPE) \
    -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME -DUSE_ALTERNATIVE_SOLVER -DALBERTAGRID \
    -I../../dune-stuff -g -O0 -DDEBUG

noinst_PROGRAMS = darcy

darcy_SOURCES = darcy.cc

darcy_CXXFLAGS = $(MPI_CPPFLAGS)  $(AM_CPPFLAGS)
darcy_LDADD = $(MPI_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
    $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GRAPE_LDFLAGS) \
    $(GRAPE_LIBS) $(BLAS_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
