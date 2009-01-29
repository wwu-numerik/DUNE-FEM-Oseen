/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#define POLORDER 0

//#define SIMPLE_PROBLEM
#define CONSTANT_PROBLEM
#define NLOG

#include <iostream>
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
//#include <dune/fem/grid/gridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/eoc.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/stokes/discretestokesfunctionspacewrapper.hh>
#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>

#include "analyticaldata.hh"
#include "parametercontainer.hh"
#include "logging.hh"
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

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 46
    Logger().Create( Parameters().getParam( "loglevel", 62 ),
                     Parameters().getParam( "logfile", std::string("dune_stokes") ) );

    L2ErrorVector l2_errors;
    ErrorColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;

    Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
    long el = gridPtr->size(0);
    profiler().Reset( 1 ); //prepare for one single run
    int err = singleRun( mpicomm, gridPtr, l2_errors );
    profiler().NextRun( l2_errors[0][0] ); //finish this run

    err += chdir( "data" );
    Dune::EocOutput eoc_output ( "eoc", "info" );
    Stuff::TexOutput texOP;
    eoc_output.printInput( *gridPtr, texOP );

    long prof_time = profiler().Output( mpicomm, 0, el );
    eoc_output.printTexAddError( el, l2_errors[0], errorColumnHeaders, 150, 0 );
    eoc_output.printTexEnd( prof_time );

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

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\ninitialising grid..." << std::endl;

    typedef Dune::LeafGridPart< GridType >
        GridPartType;
    GridPartType gridPart( *gridPtr );
    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;
    const double viscosity = Parameters().getParam( "viscosity", 1.0 );;


    infoStream << "...done." << std::endl;
    /* ********************************************************************** *
     * initialize problem                                                     *
     * ********************************************************************** */
    infoStream << "\ninitialising problem..." << std::endl;

    typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
        VelocityFunctionSpaceType;
    VelocityFunctionSpaceType velocitySpace;

    typedef Force< VelocityFunctionSpaceType >
        AnalyticalForceType;
    AnalyticalForceType analyticalForce( viscosity , velocitySpace );

    typedef DirichletData< VelocityFunctionSpaceType >
        AnalyticalDirichletDataType;
    AnalyticalDirichletDataType analyticalDirichletData( velocitySpace );

    infoStream << "...done." << std::endl;
    /* ********************************************************************** *
     * initialize model (and profiler example)                                *
     * ********************************************************************** */
    infoStream << "\ninitialising model..." << std::endl;
    typedef Dune::DiscreteStokesModelDefault<
                Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    AnalyticalForceType,
                    AnalyticalDirichletDataType,
                    gridDim,
                    polOrder > >
        StokesModelImpType;

    Dune::GridWidthProvider< GridType > gw ( *gridPtr );

    double grid_width = gw.gridWidth();
    infoStream << " \n max grid width: " << grid_width << std::endl;

    Dune::FieldVector< double, gridDim > ones( 1.0 );
    Dune::FieldVector< double, gridDim > vec_1_h( 1 - grid_width );
    //ones /= grid_width;
    Dune::FieldVector< double, gridDim > zeros( 0.0 );
    double h_fac = Parameters().getParam( "h-factor", 1.0 );

    StokesModelImpType stokesModel( h_fac / ( grid_width) ,
                                    vec_1_h,
                                    h_fac * grid_width,
                                    vec_1_h,
                                    analyticalForce,
                                    analyticalDirichletData,
                                    viscosity  );

    typedef Dune::DiscreteStokesModelInterface<
                Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    AnalyticalForceType,
                    AnalyticalDirichletDataType,
                    gridDim,
                    polOrder > >
        StokesModelType;

    typedef StokesModelType::DiscreteVelocityFunctionSpaceType
        DiscreteVelocityFunctionSpaceType;

    typedef StokesModelType::DiscretePressureFunctionSpaceType
        DiscretePressureFunctionSpaceType;

    typedef Dune::DiscreteStokesFunctionSpaceWrapper< Dune::DiscreteStokesFunctionSpaceWrapperTraits<
                DiscreteVelocityFunctionSpaceType,
                DiscretePressureFunctionSpaceType > >
        DiscreteStokesFunctionSpaceWrapperType;

    DiscreteStokesFunctionSpaceWrapperType discreteStokesFunctionSpaceWrapper( gridPart );

    typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
                DiscreteStokesFunctionSpaceWrapperType,
                StokesModelType::DiscreteVelocityFunctionType,
                StokesModelType::DiscretePressureFunctionType > >
        DiscreteStokesFunctionWrapperType;

    DiscreteStokesFunctionWrapperType discreteStokesFunctionWrapper(    "wrapped",
                                                                        discreteStokesFunctionSpaceWrapper );
    DiscreteStokesFunctionWrapperType discreteStokesFunctionWrapper2(    "wrapped2",
                                                                        discreteStokesFunctionSpaceWrapper );

    infoStream << "...done." << std::endl;
    /* ********************************************************************** *
     * initialize passes                                                      *
     * ********************************************************************** */
    infoStream << "\ninitialising passes..." << std::endl;


    typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
        StartPassType;
    StartPassType startPass;

    typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
        StokesPassType;
    StokesPassType stokesPass(  startPass,
                                stokesModel,
                                gridPart,
                                discreteStokesFunctionSpaceWrapper );

    infoStream << "...done." << std::endl;

    infoStream << "\nstarting pass..." << std::endl;
    discreteStokesFunctionWrapper.discretePressure().clear();
    discreteStokesFunctionWrapper.discreteVelocity().clear();
    stokesPass.apply( discreteStokesFunctionWrapper, discreteStokesFunctionWrapper );

    infoStream << "...done." << std::endl;

    /* ********************************************************************** *
     * Problem postprocessing (with profiler example)                                 *
     * ********************************************************************** */
    infoStream << "postprocesing" << std::endl;

    profiler().StartTiming( "Problem/Postprocessing" );

    typedef Problem< gridDim, DiscreteStokesFunctionWrapperType >
        ProblemType;
    ProblemType problem( viscosity , discreteStokesFunctionWrapper );

    typedef PostProcessor< StokesPassType, ProblemType >
        PostProcessorType;

    PostProcessorType postProcessor( discreteStokesFunctionSpaceWrapper, problem );

    postProcessor.save( *gridPtr, discreteStokesFunctionWrapper ); //dummy params, should be computed solutions );
    l2_errors.push_back( postProcessor.getError() );

    profiler().StopTiming( "Problem/Postprocessing" );

    infoStream << "...done." << std::endl;

    return 0;
}
