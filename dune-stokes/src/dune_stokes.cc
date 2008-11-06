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
#include "velocity.hh"
#include "pressure.hh"
#include "force.hh"


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
    ParameterContainer& parameters = Parameters();
    if ( !( parameters.ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
    if ( !( parameters.SetUp() ) ) {
        return 1;
    }
    else {
        parameters.SetGridDimension( GridType::dimensionworld );
        parameters.SetPolOrder( POLORDER );
        parameters.Print( std::cout );
    }

    Logger().Create(
        Logging::LOG_CONSOLE |
        Logging::LOG_FILE |
        Logging::LOG_ERR |
        Logging::LOG_DEBUG |
        Logging::LOG_INFO );

    Logging::LogStream& infoStream = Logger().Info();
    //Logging::LogStream& debugStream = Logger().Dbg();
    //Logging::LogStream& errorStream = Logger().Err();

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\ninitialising the grid..." << std::endl;
    Dune::GridPtr<GridType> gridptr( parameters.DgfFilename() );
    infoStream << "...done." << std::endl;

    /* ********************************************************************** *
     * initialize the analytical problem                                      *
     * ********************************************************************** */
    infoStream << "\ninitializing the analytical problem..." << std::endl;
    Velocity< GridType::dimensionworld > velocity;
    velocity.TestMe();
    Pressure< GridType::dimensionworld > pressure;
    pressure.TestMe();
    Force< GridType::dimensionworld > force;
    force.TestMe();
    DirichletData< GridType::dimensionworld > dirichletData;
    dirichletData.TestMe();
    infoStream << "\n...done." << std::endl;


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
