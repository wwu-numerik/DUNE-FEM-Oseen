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
//#define CONSTANT_PROBLEM
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
#include <dune/fem/misc/femeoc.hh>
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
typedef std::vector<std::string> ColumnHeaders;
const std::string errheaders[] = { "h", "tris","C11","C12","D11","D12","Velocity L2 Error", "Pressure L2 Error" };
const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );

struct RunInfo
{
    std::vector< double > L2Errors;
    double grid_width;
    int refine_level;
    double run_time;
    long codim0;
    double c11,d11,c12,d12;
};

template < class GridPartType >
RunInfo singleRun(  CollectiveCommunication mpicomm,
                Dune::GridPtr< GridType > gridPtr,
                GridPartType& gridPart,
                int pow1, int pow2, int pow3, int pow4  );
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
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;

    int err = 0;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::TexOutput< RunInfo > texwriter( errorColumnHeaders );

    profiler().Reset( 9 ); //prepare for 9 single runs

    for ( int ref = 0; ref < 6 ; ++ref ) {
        for ( int i = 0, num = 0; i < 3; ++i ) {
            for ( int j = 0; j < 3; ++j,++num ) {
                Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
                gridPtr->globalRefine( ref );
                typedef Dune::AdaptiveLeafGridPart< GridType >
                    GridPartType;
                GridPartType gridPart( *gridPtr );
                std::string ff = "matlab__pow1_" + Stuff::toString( i ) + "_pow2_" + Stuff::toString( j );
                Logger().SetPrefix( ff );
                RunInfo info = singleRun( mpicomm, gridPtr, gridPart, i, j, -2, -2 );
                l2_errors.push_back( info.L2Errors );
    //            profiler().NextRun( info.L2Errors ); //finish this run

                eoc_output.setErrors( idx,info.L2Errors );
                eoc_output.write( info.grid_width, info.codim0 ,num,num);
                texwriter.setInfo( info );
                eoc_output.write( texwriter );

            }
        }
    }

//    Stuff::TexOutput texOP;
//    eoc_output.printInput( *gridPtr, texOP );
//
//    long prof_time = profiler().Output( mpicomm, 0, el );


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
RunInfo singleRun(  CollectiveCommunication mpicomm,
                Dune::GridPtr< GridType > gridPtr,
                GridPartType& gridPart,
                int pow1, int pow2, int pow3, int pow4  )
{
    Logging::LogStream& infoStream = Logger().Info();
    ParameterContainer& parameters = Parameters();
    RunInfo info;

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\ninitialising grid..." << std::endl;

    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;
    const double viscosity = Parameters().getParam( "viscosity", 1.0 );
    info.codim0 = gridPtr->size( 0 );
    info.codim0 = gridPart.grid().size( 0 );



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
    info.grid_width = grid_width;

    typedef Dune::FieldVector< double, gridDim >
        ConstVec;

    const double minpow = -2;
    double c11 = pow1 > minpow ? std::pow( grid_width, pow1 ) :0;
    double d11 = pow2 > minpow ? std::pow( grid_width, pow2 ) :0;
    double c12 = pow3 > minpow ? std::pow( grid_width, pow3 ) :0;
    double d12 = pow4 > minpow ? std::pow( grid_width, pow4 ) :0;
    StokesModelImpType stokesModel( c11,
                                    ConstVec ( c12 ),
                                    d11,
                                    ConstVec ( d12 ),
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

    postProcessor.save( *gridPtr, discreteStokesFunctionWrapper );
    info.L2Errors = postProcessor.getError();
    info.c11 = c11;
    info.c12 = c12;
    info.d11 = d11;
    info.d12 = d12;

    profiler().StopTiming( "Problem/Postprocessing" );

    infoStream << "...done." << std::endl;

    return info;
}

