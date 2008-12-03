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
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction

#include <dune/stokes/discretestokesmodelinterface.hh>

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
    Logging::LogStream& debugStream = Logger().Dbg();
    Logging::LogStream& errorStream = Logger().Err();

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
     * initialize function spaces and functions                               *
     * ********************************************************************** */
    infoStream << "\ninitialising functions..." << std::endl;
    // velocity
    typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
        VelocityFunctionSpaceType;
    typedef Dune::DiscontinuousGalerkinSpace<   VelocityFunctionSpaceType,
                                                GridPartType,
                                                polOrder >
        DiscreteVelocityFunctionSpaceType;
    DiscreteVelocityFunctionSpaceType velocitySpace( gridPart );
    typedef Dune::AdaptiveDiscreteFunction< DiscreteVelocityFunctionSpaceType >
        DiscreteVelocityFunctionType;
    DiscreteVelocityFunctionType exactVelocity( "exact_velocity",
                                                velocitySpace );
    exactVelocity.clear();
    //sigma
    typedef Dune::MatrixFunctionSpace< double, double, gridDim, gridDim, gridDim >
        SigmaFunctionSpaceType;
    typedef Dune::DiscontinuousGalerkinSpace<   SigmaFunctionSpaceType,
                                                GridPartType,
                                                polOrder >
        DiscreteSigmaFunctionSpaceType;
    DiscreteSigmaFunctionSpaceType sigmaSpace( gridPart );
    typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
        DiscreteSigmaFunctionType;
    DiscreteSigmaFunctionType exactSigma(   "exact_sigma",
                                            sigmaSpace );
    exactSigma.clear();
    // pressure
    typedef Dune::FunctionSpace< double, double, gridDim, 1 >
        PressureFunctionSpaceType;
    typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
                                                GridPartType,
                                                polOrder >
        DiscretePressureFunctionSpaceType;
    DiscretePressureFunctionSpaceType pressureSpace( gridPart );
    typedef Dune::AdaptiveDiscreteFunction< DiscretePressureFunctionSpaceType >
        DiscretePressureFunctionType;
    DiscretePressureFunctionType exactPressure( "exact_pressure",
                                                pressureSpace );
    exactPressure.clear();
    // right hand side
    DiscreteVelocityFunctionType righthandSide( "rhs",
                                                velocitySpace );
    righthandSide.clear();
    // dirichlet data
    DiscreteVelocityFunctionType dirichletData( "g_D",
                                                velocitySpace );
    dirichletData.clear();

    infoStream << "...done." << std::endl;

    /* ********************************************************************** *
     * initialize model (and profiler example)                                *
     * ********************************************************************** */
    infoStream << "\ninitialising model..." << std::endl;
    profiler().Reset( 1 ); //prepare for one single run of code
    profiler().StartTiming( "Problem/Postprocessing" );
    typedef Problem< ProblemTraits< gridDim, VelocityFunctionSpaceType, PressureFunctionSpaceType > >
        Problemtype;
    Problemtype problem( parameters.viscosity(), velocitySpace, pressureSpace );
    //problem.testMe();
    typedef PostProcessor< Problemtype, GridPartType, DiscreteVelocityFunctionType, DiscretePressureFunctionType >
        PostProcessorType;
    PostProcessorType postProcessor( problem, gridPart, velocitySpace, pressureSpace );
    postProcessor.save( *gridPtr, exactPressure, exactVelocity ); //dummy params, should be computed solutions );

    infoStream << "...done." << std::endl;

    profiler().StopTiming( "Problem/Postprocessing" );
    profiler().Output( mpicomm, 0, exactPressure.size() );

    return 0;
  } //end try
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    return -1;
  }
  catch (std::exception &f){
    std::cerr << "stdlib reported error: " << f.what() << std::endl;
    return -2;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return -3;
  } //end main
}
