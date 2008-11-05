/** \file dune_stokes.cc
    \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
//#define GRIDTYPE ALUGRID_SIMPLEX
#include <iostream>
#include <memory> //not including this gives error of undefined autopointer in dgfparser.hh
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid
//#include <dune/grid/io/file/dgfparser/dgfparser.hh> // for the grid
//#include <dune/grid/utility/gridtype.hh>

//#include "traits.hh"
#include "parametercontainer.hh"
#include "logging.hh"
#include "problem.hh"

/**
 *  \brief main function
 *
 *  \c more main function
 *
 *  \attention attention
 *
 *  \param argc number of arguments from command line
 *  \param argv array of arguments from command line
 **/
int main( int argc, char** argv )
{
  try{
    //Maybe initialize Mpi
    //Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    Logger().Create(
        Logging::LOG_CONSOLE |
        Logging::LOG_FILE |
        Logging::LOG_ERR |
        Logging::LOG_DEBUG |
        Logging::LOG_INFO );

    ParameterContainer& parameters = Parameters();
    if ( !( parameters.ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
    if ( !( parameters.SetUp() ) ) {
        return 1;
    }
    else {
        //parameters.SetPolOrder( POLORDER );
        parameters.Print( std::cout );
    }





    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
//    Dune::GridPtr<GridType> gridptr( "grid.dgf" );
//    parameters.SetGridDimension( GridType::dimensionworld );
//
//    /* ********************************************************************** *
//     * initialize the analytical problem                                      *
//     * ********************************************************************** */
//    Velocity< GridType::dimensionworld > velocity;
//    Pressure< GridType::dimensionworld > pressure;
//    Force< GridType::dimensionworld > force;
//    DirichletData< GridType::dimensionworld > dirichletData;
//
//    Dune::FieldVector< double, 2 > x, y;
//    x[0] = 1.0;
//    x[1] = 1.0;
//    y = velocity( x );
//    std::cout << "\ny1: " << y[0]
//        << "\ny2: " << y[1] << std::endl;


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
