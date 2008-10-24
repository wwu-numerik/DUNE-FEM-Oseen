// Dune includes
#include <config.h>           // file constructed by ./configure script
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/common/gridinfo.hh> // definition of gridinfo
#include <dune/common/mpihelper.hh> // include mpi helper class

template <int dim> 
struct Info 
{
  static void print() 
  {
    std::cout << "**************************************************" << std::endl;
    std::cout << "Print grid infos for dim = " << dim << std::endl;
    std::cout << "**************************************************" << std::endl;

    // make a grid
    typedef typename Dune::SGrid<dim,dim> GridType; 
    Dune::FieldVector<int,dim> N(1); 
    Dune::FieldVector<typename GridType::ctype,dim> L(0.0);
    Dune::FieldVector<typename GridType::ctype,dim> H(1.0);
    GridType grid(N,L,H);

    // print some information about the grid
    std::cout << "gridinfo:" << std::endl; 
    Dune::gridinfo(grid);
    std::cout << std::endl; 

    // print level list 
    std::cout << "gridlevellist:" << std::endl; 
    Dune::gridlevellist(grid,0,"LevelList ");
    std::cout << std::endl; 

    // print leaf list  
    std::cout << "gridleaflist:" << std::endl; 
    Dune::gridleaflist(grid,"LeafList ");

    std::cout << std::endl << std::endl; 
  }
};

// main routine 
int main(int argc, char **argv)
{
  // start try/catch block to get error messages from dune
  try{
    // initialize MPI, finalize is done automatically on exit 
    Dune::MPIHelper::instance(argc,argv);

    // print grid info for dim = 1 
    Info<1>::print();

    // print grid info for dim = 2 
    Info<2>::print();

    // print grid info for dim = 3 
    Info<3>::print();
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
