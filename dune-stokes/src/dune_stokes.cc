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
#include <memory> //not including this gives error of undefined autopointer in dgfparser.hh
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
     * initialize model                                                       *
     * ********************************************************************** */
    infoStream << "\ninitialising model..." << std::endl;
    typedef Problem< ProblemTraits< gridDim, VelocityFunctionSpaceType, PressureFunctionSpaceType > >
        Problemtype;
    Problemtype problem( parameters.viscosity(), velocitySpace, pressureSpace );
    //problem.testMe();
    typedef PostProcessor< Problemtype, GridPartType, DiscreteVelocityFunctionType, DiscretePressureFunctionType >
        PostProcessorType;
    PostProcessorType postProcessor( problem, gridPart, velocitySpace, pressureSpace );
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
