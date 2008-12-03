/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
//#define GRIDTYPE ALUGRID_SIMPLEX
#include <iostream>
//#include <memory> //not including this gives error of undefined autopointer in dgfparser.hh
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid
//#include <dune/grid/io/file/dgfparser/dgfparser.hh> // for the grid
//#include <dune/grid/utility/gridtype.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/saddlepoint_inverse_operator.hh>
#include <dune/stokes/stokespass.hh>

#include "parametercontainer.hh"
#include "logging.hh"
#include "problem.hh"
#include "postprocessing.hh"
#include "profiler.hh"

/**
 *  \brief  main function
 *
 *  \attention  attention
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main( int argc, char** argv )
{
  try{
    //Maybe initialize Mpi
    #if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
    #else
        typedef Dune::CollectiveCommunication<double > CollectiveCommunication;
    #endif
    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc, argv);
    CollectiveCommunication mpicomm ( mpihelper.getCommunicator() );

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

    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;

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
    infoStream << "\ninitialising grid..." << std::endl;
    typedef Dune::LeafGridPart< GridType >
        GridPartType;
    Dune::GridPtr< GridType > gridPtr( parameters.DgfFilename() );
    GridPartType gridPart( *gridPtr );
    infoStream << "...done." << std::endl;

    /* ********************************************************************** *
     * initialize model                                                       *
     * ********************************************************************** */
    infoStream << "\ninitialising model..." << std::endl;

    typedef Dune::DiscreteStokesModelDefault<
                Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    gridDim,
                    polOrder > >
        StokesModelType;

    StokesModelType stokesModel;

    StokesModelType::DiscreteVelocityFunctionSpaceType velocitySpace( gridPart );
    StokesModelType::DiscreteSigmaFunctionSpaceType sigmaSpace( gridPart );
    StokesModelType::DiscretePressureFunctionSpaceType pressureSpace( gridPart );

    infoStream << "...done." << std::endl;


    /* ********************************************************************** *
     * initialize model (and profiler example)                                *
     * ********************************************************************** */
    infoStream << "\ninitialising passes..." << std::endl;
    profiler().Reset( 1 ); //prepare for one single run of code
    profiler().StartTiming( "Problem/Postprocessing" );
        StartPassType;
    StartPassType startPass;
    postProcessor.save( *gridPtr, exactPressure, exactVelocity ); //dummy params, should be computed solutions );

    profiler().StopTiming( "Problem/Postprocessing" );
    profiler().Output( mpicomm, 0, exactPressure.size() );

    typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
        StokesPassType;

    StokesPassType stokesPass(  startPass,
                                velocitySpace,
                                sigmaSpace,
                                pressureSpace,
                                stokesModel,
                                gridPart );

    //StokesPassType::
    StokesPassType::DomainType uDummy( "uDummy", velocitySpace );
    stokesPass( uDummy, uDummy );


    infoStream << "...done." << std::endl;

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
