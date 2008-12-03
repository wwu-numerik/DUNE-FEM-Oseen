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
#include <dune/fem/misc/eoc.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/saddlepoint_inverse_operator.hh>
#include <dune/stokes/stokespass.hh>

#include "parametercontainer.hh"
#include "logging.hh"
#include "problem.hh"
#include "postprocessing.hh"
#include "profiler.hh"
#include "stuff.hh"

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication<double > CollectiveCommunication;
#endif

typedef std::vector< std::vector<double> > L2ErrorVector;
typedef std::vector<std::string> ErrorColumnHeaders;
const std::string errheaders[] = { "Velocity L2 Error", "Pressure L2 Error" };
const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );

int singleRun( CollectiveCommunication mpicomm, Dune::GridPtr< GridType > gridPtr,
                L2ErrorVector& l2_errors );

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

    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc, argv);
    CollectiveCommunication mpicomm ( mpihelper.getCommunicator() );

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
    if ( !(  Parameters().SetUp() ) ) {
        return 2;
    }
    else {
        Parameters().SetGridDimension( GridType::dimensionworld );
        Parameters().SetPolOrder( POLORDER );
        Parameters().Print( std::cout );
    }

    Logger().Create(
        Logging::LOG_CONSOLE |
        Logging::LOG_FILE |
        Logging::LOG_ERR |
        Logging::LOG_DEBUG |
        Logging::LOG_INFO );

    L2ErrorVector l2_errors;
    ErrorColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;

    Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
    int err = singleRun( mpicomm, gridPtr, l2_errors );

    err += chdir( "data" );
    Dune::EocOutput eoc_output ( "eoc", "info" );
    Stuff::TexOutput texOP;
    eoc_output.printInput( *gridPtr, texOP );

    eoc_output.printTexAddError( gridPtr->size(0), l2_errors[0], errorColumnHeaders, 150, 0 );
    eoc_output.printTexEnd( 650 );

    return err;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

int singleRun( CollectiveCommunication mpicomm, Dune::GridPtr< GridType > gridPtr,
                L2ErrorVector& l2_errors )
{
    Logging::LogStream& infoStream = Logger().Info();
    ParameterContainer& parameters = Parameters();

    typedef Dune::LeafGridPart< GridType >
        GridPartType;
    GridPartType gridPart( *gridPtr );

    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;

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

    typedef Problem< ProblemTraits< gridDim, StokesModelType::VelocityFunctionSpaceType, StokesModelType::PressureFunctionSpaceType > >
        Problemtype;
    Problemtype problem( parameters.viscosity(), velocitySpace, pressureSpace );
    //problem.testMe();

    typedef Dune::StartPass< StokesModelType::DiscreteVelocityFunctionType, -1 >
        StartPassType;
    StartPassType startPass;

    typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
        StokesPassType;

    StokesPassType stokesPass(  startPass,
                                velocitySpace,
                                sigmaSpace,
                                pressureSpace,
                                stokesModel,
                                gridPart );

    StokesPassType::DomainType uDummy( "uDummy", velocitySpace );
    stokesPass( uDummy, uDummy );

    typedef PostProcessor<  Problemtype, GridPartType,
                            StokesModelType::DiscreteVelocityFunctionType,
                            StokesModelType::DiscretePressureFunctionType >
        PostProcessorType;

    PostProcessorType postProcessor( problem, gridPart, velocitySpace, pressureSpace );

    StokesModelType::DiscretePressureFunctionType pDummy ( "pDummy", pressureSpace );
    postProcessor.save( *gridPtr, pDummy, uDummy ); //dummy params, should be computed solutions );
    l2_errors.push_back( postProcessor.getError() );

    profiler().StopTiming( "Problem/Postprocessing" );
    profiler().Output( mpicomm, 0, uDummy.size() );


    infoStream << "...done." << std::endl;

    return 0;
}
