#include"config.h"               // know what grids are present
#include<iostream>               // for input/output to shell
#include<fstream>                // for input/output to files
#include<vector>                 // STL vector class
#include <dune/common/mpihelper.hh> // include mpi helper class 

#include <dune/grid/common/gridinfo.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh> // VTK output routines

// checks for defined gridtype and includes appropriate dgfparser implementation
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

// for manual grid definition
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>


void dgfTest ()
{
  using namespace Dune;

  // use unitcube from grids 
  std::stringstream dgfFileName;
  dgfFileName << "./unitcube2.dgf";

  std::cout << "Try to open " << dgfFileName.str() << std::endl;

  // create grid
  typedef YaspGrid<2,2> GridType;
  GridPtr<GridType> gridPtr ( dgfFileName.str() );

  // print info 
  gridinfo( *gridPtr );
  std::cout << std::endl;

  // write grid in VTK format
  VTKWriter<GridType> vtkwriter(*gridPtr);
  vtkwriter.write("grid1",Dune::VTKOptions::ascii);
}

void dgfGridType () 
{
  using namespace Dune;

  // use unitcube from grids 
  std::stringstream dgfFileName;
  // GridType is defined by dgfgridtype.hh
  dgfFileName << "./unitcube" << GridType :: dimension << ".dgf";

  std::cout << "Try to open " << dgfFileName.str() << std::endl;

  // create grid pointer
  GridPtr<GridType> gridPtr( dgfFileName.str() );

  // grid reference 
  GridType& grid = *gridPtr;

  // half grid width 4 times 
  int level = 4 * DGFGridInfo<GridType>::refineStepsForHalf();

  // refine grid until upper limit of level 
  grid.globalRefine(level);

  // print info 
  gridinfo( *gridPtr );
  std::cout << std::endl;

  // write grid in VTK format
  VTKWriter<GridType> vtkwriter(*gridPtr);
  vtkwriter.write("grid2",Dune::VTKOptions::ascii);
}

//===============================================================
// main rountine 
//===============================================================
int main (int argc , char ** argv)
{ 
  // initialize MPI, finalize is done automatically on exit 
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    using namespace Dune;

    // use manual typedef
    dgfTest();

    // use grid type definition 
    dgfGridType();
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
