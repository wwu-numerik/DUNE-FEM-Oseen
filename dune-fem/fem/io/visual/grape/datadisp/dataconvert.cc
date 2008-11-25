//************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//************************************************************
#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#if HAVE_MPI == 1 
<<<<<<< HEAD:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc
#warning "Visualization does not work in parallel" 
=======
#error "Visualization only works without MPI" 
>>>>>>> reverts all post 1.1.1 commits except adding a new build target for suse10.3:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc
#endif 

#include <dune/common/misc.hh>
#include <dune/common/exceptions.hh>
using namespace Dune;

// include definition of grid type 
<<<<<<< HEAD:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc
// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
=======
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
>>>>>>> reverts all post 1.1.1 commits except adding a new build target for suse10.3:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc

// include data reading 
#include <dune/fem/io/visual/grape/datadisp/printhelp.cc>

// uses readtuple data instead of readiotupledata.
#include <dune/fem/io/visual/grape/datadisp/readtupledata.cc>
#include <dune/fem/io/visual/grape/datadisp/readioparams.cc> 
<<<<<<< HEAD:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc
#include <dune/fem/misc/mpimanager.hh>

int main(int argc, char **argv)
{
  MPIManager::initialize(argc,argv);
  try {
    Parameter::append(argc,argv);
=======

int main(int argc, char **argv)
{
  try {
>>>>>>> reverts all post 1.1.1 commits except adding a new build target for suse10.3:dune-fem/fem/io/visual/grape/datadisp/dataconvert.cc
    if (argc < 2)
    {
      print_help(argv[0]);
      return(0);
    }   

    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help"))
    {
      print_help(argv[0]);
      return(0);
    }
    return readParameterList(argc,argv,false);
  }
  catch (Dune::Exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  return 0;
}
