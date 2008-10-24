// Id: gettingstarted.cc 2006-11-03 gersbach

// Dune includes
#include<config.h>               // file constructed by ./configure script 
#include<dune/grid/sgrid.hh>     // load sgrid definition 
#include <dune/grid/common/gridpart.hh>

//- local inlcudes
#include "../../space/lagrangespace.hh"

using namespace Dune;

const int dim = 2;
typedef Dune::SGrid<dim,dim> GridType;

//! the index set we are using
//typedef LeafGridPart<GridType> GridPartType;
typedef LevelGridPart < GridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
typedef FunctionSpace < double , double, dim , 1 > FuncSpace;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace < FuncSpace , GridPartType , 1 > FuncSpaceType ;

int main (int argc, char **argv)
{
  
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }
  int maxLevel = atoi( argv[1] );

  //! construct a grid
  Dune::FieldVector<int,dim> N(3); 
  Dune::FieldVector<GridType::ctype,dim> L(-1.0); 
  Dune::FieldVector<GridType::ctype,dim> H(1.0); 
  GridType grid(N,L,H); 
  grid.globalRefine(maxLevel);
  
  for (int level=0; level < maxLevel; ++level ) 
  {
    //! create GridPart of level and FunctionSpace
    GridPartType part(grid, level);
    FuncSpaceType linFuncSpace ( part );

    //! print some information
    std::cout << "\nSpace is continuous? " << linFuncSpace.continuous();
    std::cout << "\nSpace has " << linFuncSpace.size() << " number of DOFs. \n";
  } 
}


