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
//#define NLOG

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

template < class GridPartType >
int singleRun(  CollectiveCommunication mpicomm,
                Dune::GridPtr< GridType > gridPtr,
                GridPartType& gridPart,
                L2ErrorVector& l2_errors, int pow1, int pow2  );
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

//
//
    long el = 0;//gridPtr->size(0);
    profiler().Reset( 9 ); //prepare for one single run
    int err = 0;
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
            Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
            el = gridPtr->size(0);
            typedef Dune::AdaptiveLeafGridPart< GridType >
                GridPartType;
            GridPartType gridPart( *gridPtr );
            std::string ff = "matlab__pow1_" + Stuff::toString( i ) + "_pow2_" + Stuff::toString( j );
            Logger().SetPrefix( ff );
            err += singleRun( mpicomm, gridPtr, gridPart, l2_errors, i, j );
            profiler().NextRun( l2_errors[0][0] ); //finish this run
        }
    }


    err += chdir( "data" );
//    Dune::EocOutput eoc_output ( "eoc", "info" );
//    Stuff::TexOutput texOP;
//    eoc_output.printInput( *gridPtr, texOP );

    long prof_time = profiler().Output( mpicomm, 0, el );
//    eoc_output.printTexAddError( el, l2_errors[0], errorColumnHeaders, 150, 0 );
//    eoc_output.printTexEnd( prof_time );

    return err;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

template < class GridPartType >
int singleRun(  CollectiveCommunication mpicomm,
                Dune::GridPtr< GridType > gridPtr,
                GridPartType& gridPart,
                L2ErrorVector& l2_errors, int pow1, int pow2  )
{
    Logging::LogStream& infoStream = Logger().Info();
    ParameterContainer& parameters = Parameters();

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\ninitialising grid..." << std::endl;

    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;
    const double viscosity = Parameters().getParam( "viscosity", 1.0 );


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

    StokesModelImpType stokesModel( std::pow( grid_width, pow1 ),
                                    zeros,
                                    std::pow( grid_width, pow2 ),
                                    zeros,
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

    typedef typename StokesModelType::DiscreteVelocityFunctionSpaceType
        DiscreteVelocityFunctionSpaceType;

    typedef typename StokesModelType::DiscretePressureFunctionSpaceType
        DiscretePressureFunctionSpaceType;

    typedef Dune::DiscreteStokesFunctionSpaceWrapper< Dune::DiscreteStokesFunctionSpaceWrapperTraits<
                DiscreteVelocityFunctionSpaceType,
                DiscretePressureFunctionSpaceType > >
        DiscreteStokesFunctionSpaceWrapperType;

    DiscreteStokesFunctionSpaceWrapperType discreteStokesFunctionSpaceWrapper( gridPart );

    typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
                DiscreteStokesFunctionSpaceWrapperType,
                typename StokesModelType::DiscreteVelocityFunctionType,
                typename StokesModelType::DiscretePressureFunctionType > >
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
