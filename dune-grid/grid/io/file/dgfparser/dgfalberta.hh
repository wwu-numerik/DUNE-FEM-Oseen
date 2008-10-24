#ifndef DUNE_DGFPARSERALBERTA_HH
#define DUNE_DGFPARSERALBERTA_HH

// only include if ALBERTA is used 
#if defined ENABLE_ALBERTA 
#include <dune/grid/albertagrid.hh>
#include "dgfparser.hh"
namespace Dune {
template <int dim,int dimworld>
class MacroGrid::Impl<AlbertaGrid<dim,dimworld> > {
  typedef MPIHelper::MPICommunicator MPICommunicatorType;
public:
  static AlbertaGrid<dim,dimworld>* generate(MacroGrid& mg,
					     const char* filename, 
					     MPICommunicatorType MPICOMM = MPIHelper::getCommunicator());
};

template <int dimworld>
struct DGFGridInfo< AlbertaGrid<dimworld,dimworld> > {
    static int refineStepsForHalf() {return dimworld;}
    static double refineWeight() {return 0.5;}
};
}
#include "dgfalberta.cc"
#endif

#endif
