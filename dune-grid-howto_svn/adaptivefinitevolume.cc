#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class

// checks for defined gridtype and inlcudes appropriate dgfparser implementation 
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> 

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/common/mpihelper.hh> // include mpi helper class

#include "vtkout.hh"
#include "unitcube.hh" 
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"
#include "finitevolumeadapt.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (Dune::GeometryType gt)
  {
	if (gt.dim()==dim) return true;
	return false;
  }
}; 

template<class G>
void gnuplot (G& grid, std::vector<double>& c)
{
  // first we extract the dimensions of the grid  
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  // such a type is exported by every grid implementation
  typedef typename G::ctype ct;

  // the grid has an iterator providing the access to
  // all elements (better codim 0 entities) which are leafs
  // of the refinement tree.
  // Note the use of the typename keyword and the traits class
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;

  // make a mapper for codim 0 entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout> 
	mapper(grid);

  // iterate through all entities of codim 0 at the leafs
  int count = 0;
  for (ElementLeafIterator it = grid.template leafbegin<0>(); 
	   it!=grid.template leafend<0>(); ++it)
	{
	  Dune::GeometryType gt = it->type();
	  const Dune::FieldVector<ct,dim>& 
		local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
	  Dune::FieldVector<ct,dimworld> 
		global = it->geometry().global(local);
	  std::cout << global[0] << " " << c[mapper.map(*it)] << std::endl;
	  count++;
	}
}


template<class G>
void timeloop (G& grid, double tend, int lmin, int lmax)
{
  // make a mapper for codim 0 entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout> 
	mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);
  for (int i=grid.maxLevel(); i<lmax; i++)
	{
	  if (grid.maxLevel()>=lmax) break;
	  finitevolumeadapt(grid,mapper,c,lmin,lmax,0);
	  initialize(grid,mapper,c);
	}

  // write initial data 
  vtkout(grid,c,"concentration",0);

  // variables for time, timestep etc.
  double dt, t=0;
  double saveStep = 0.1;
  const double saveInterval = 0.1;
  int counter = 0;
  int k = 0;

  std::cout << "s=" << grid.size(0) << " k=" << k << " t=" << t << std::endl;
  while (t<tend)
	{
    // augment time step counter 
    ++k;

    // apply finite volume scheme 
	  evolve(grid,mapper,c,t,dt);
    
    // augment time 
	  t += dt;

    // check if data should be written 
    if (t >= saveStep) 
    {
      // write data 
      vtkout(grid,c,"concentration",counter);

      // increase counter and saveStep for next interval 
      saveStep += saveInterval;
      ++counter;
    }

    // print info about time, timestep size and counter  
    std::cout << "s=" << grid.size(0) << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
    
    // for unstructured grids call adaptation algorithm 
    if( Dune :: Capabilities :: IsUnstructured<G> :: v )
    {
 	    finitevolumeadapt(grid,mapper,c,lmin,lmax,k);
    }
	}

  // write last time step 
  vtkout(grid,c,"concentration",counter);
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv)
{ 
  // initialize MPI, finalize is done automatically on exit 
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    using namespace Dune;

    // use unitcube from grids 
    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    GridPtr<GridType> gridPtr( dgfFileName.str() );

    // grid reference 
    GridType& grid = *gridPtr;

    // minimal allowed level during refinement 
    int minLevel = 2 * DGFGridInfo<GridType>::refineStepsForHalf();

    // refine grid until upper limit of level 
    grid.globalRefine(minLevel); 

    // maximal allowed level during refinement 
    int maxLevel = minLevel + 3 * DGFGridInfo<GridType>::refineStepsForHalf();

    // do time loop until end time 0.5 
    timeloop(grid, 0.5, minLevel, maxLevel);
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}