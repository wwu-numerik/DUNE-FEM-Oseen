/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//#define POLORDER 1

//#define SIMPLE_PROBLEM
//#define CONSTANT_PROBLEM
//#define ROTATE_PROBLEM
//#define NLOG

#include <iostream>
#include <cmath>
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

#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/postprocessing.hh>
#include <dune/stuff/profiler.hh>

#include "analyticaldata.hh"
#include "velocity.hh"
#include "pressure.hh"
#include "problem.hh"

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

typedef std::vector< std::vector<double> > L2ErrorVector;
typedef std::vector<std::string> ColumnHeaders;
const std::string errheaders[] = { "h", "el't","runtime","C11","C12","D11","D12","Velocity", "Pressure" };
const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );

struct RunInfo
{
    std::vector< double > L2Errors;
    double grid_width;
    int refine_level;
    double run_time;
    long codim0;
    int polorder_velocity;
    int polorder_pressure;
    int polorder_sigma;
    double c11,d11,c12,d12;
    bool bfg;
    std::string gridname;
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
        std::cerr << "\nUsage: " << argv[0] << " parameterfile \n" << "\t --- OR --- \n";
        std::cerr << "\nUsage: " << argv[0] << " paramfile:"<<"file" << " more-opts:val ..." << std::endl;
        Parameters().PrintParameterSpecs( std::cerr );
        std::cerr << std::endl;
        return 2;
    }
    else {
        Parameters().SetGridDimension( GridType::dimensionworld );
        Parameters().SetPolOrder( POLORDER );
//        Parameters().Print( std::cout );
    }

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    Logger().Create( Parameters().getParam( "loglevel", 62 ),
                     Parameters().getParam( "logfile", std::string("dune_stokes") ) );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    // coloumn headers for eoc table output
    L2ErrorVector l2_errors;
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;

    int err = 0;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::TexOutput< RunInfo > texwriter( errorColumnHeaders );

    profiler().Reset( 9 ); //prepare for 9 single runs

    //the outermost loop counts from 0 to this
    //also sets the refinelevel in non-variation context
    int maxref = Parameters().getParam( "maxref", 0 );
    int minref = Parameters().getParam( "minref", 0 );

    // the min andmax exponent that will be used in stabilizing terms
    // ie c11 = ( h_^minpow ) .. ( h^maxpow )
    // pow = -2 is interpreted as 0
    // in non-variation context minpow is used for c12 exp and maxpow for d12,
    //      c11 and D11 are set to zero
    int maxpow = Parameters().getParam( "maxpow", 2 );
    //set minref == maxref to get only a single run in non variation part
    int minpow = Parameters().getParam( "minpow", -2 );



    if ( Parameters().getParam( "multirun", true ) ) {
        /** all four stab parameters are permutated in [ minpow ; maxpow ]
            inside an outer loop that increments the grid's refine level
        **/
        for ( int ref = minref; ref <= maxref; ++ref ) {
            int i,j,k,l;
            i = j = k = l = maxpow - 1;
//            for ( int i = minpow; i < maxpow; ++i ) {
//                for ( int j = minpow; j < maxpow; ++j ) {
//                    for ( int k = minpow; k < maxpow; ++k ) {
//                        for ( int l = minpow; l < maxpow; ++l ) {
                            if ( per_run_log_target ) { //sets unique per run filename if requested
                                std::string ff = "matlab__pow1_" + Stuff::toString( i ) + "_pow2_" + Stuff::toString( j );
                                Logger().SetPrefix( ff );
                            }
                            Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
                            gridPtr->globalRefine( Dune::DGFGridInfo< GridType >::refineStepsForHalf()* ref );
                            typedef Dune::AdaptiveLeafGridPart< GridType >
                                GridPartType;
                            GridPartType gridPart( *gridPtr );
                            Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
                            //do some matlab magic to suppress errors and display a little info baout current params
                            matlabLogStream <<std::endl<< "\nclear;\ntry\ntic;warning off all;" << std::endl;
                            matlabLogStream << "disp(' c/d-11/12 h powers: "<< i << " "
                                    << j << " " << k << " "
                                    << l << "') \n" << std::endl;
                            //actual work
                            profiler().StartTiming( "SingleRun" );
                            RunInfo info = singleRun( mpicomm, gridPtr, gridPart, i, j, k, l );
                            profiler().StopTiming( "SingleRun" );
                            info.run_time = profiler().GetTiming( "SingleRun" );
                            // new line and closing try/catch in m-file
                            matlabLogStream << "\ncatch\ndisp('errors here');\nend\ntoc\ndisp(' ');\n" << std::endl;
                            l2_errors.push_back( info.L2Errors );
//    			            profiler().NextRun( info.L2Errors ); //finish this run

                            //push errors to eoc-outputter class
                            eoc_output.setErrors( idx,info.L2Errors );
                            texwriter.setInfo( info );
                            bool lastrun = ( //this test is somewhat stupid, make it smart!!
                                ( ref >= ( maxref  ) ) &&
                                ( j + k + l + i >= 4 * ( maxpow - 1 ) ) );
                            //the writer needs to know if it should close the table etc.
                            eoc_output.write( texwriter, lastrun );
//                        }
//                    }
//                }
//            }
        }
    }
    else { //we don't do any variation here, just one simple run, no eoc, nothing

        for ( int ref = minref; ref <= maxref; ++ref ) {
            Logger().SetPrefix( "dune_stokes_ref_"+Stuff::toString(ref) );
            Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
            gridPtr->globalRefine( Dune::DGFGridInfo< GridType >::refineStepsForHalf()*ref );
            typedef Dune::AdaptiveLeafGridPart< GridType >
                GridPartType;
            GridPartType gridPart( *gridPtr );
            profiler().StartTiming( "SingleRun" );
            RunInfo info = singleRun( mpicomm, gridPtr, gridPart, minpow, maxpow, -9, -9 );
            profiler().StopTiming( "SingleRun" );
            info.run_time = profiler().GetTiming( "SingleRun" );
            l2_errors.push_back( info.L2Errors );
            eoc_output.setErrors( idx,info.L2Errors );
            texwriter.setInfo( info );
            eoc_output.write( texwriter, ( ref == maxref ) );
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
                    const int pow1,
                    const int pow2,
                    const int pow3,
                    const int pow4  )
{
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    ParameterContainer& parameters = Parameters();
    RunInfo info;

    debugStream << "\nsingleRun( pow1: " << pow1 << ","
                << "\n           pow2: " << pow2 << ","
                << "\n           pow3: " << pow3 << ","
                << "\n           pow4: " << pow4 << " )" << std::endl;

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\n- initialising grid" << std::endl;

    const int gridDim = GridType::dimensionworld;
    info.codim0 = gridPtr->size( 0 );
    info.codim0 = gridPart.grid().size( 0 );
    Dune::GridWidthProvider< GridType > gw ( *gridPtr );
    double grid_width = gw.gridWidth();
//    infoStream << "  - max grid width: " << grid_width << std::endl;
    info.grid_width = grid_width;

    /* ********************************************************************** *
     * initialize problem                                                     *
     * ********************************************************************** */
    infoStream << "\n- initialising problem" << std::endl;

    const int polOrder = POLORDER;
    debugStream << "  - polOrder: " << polOrder << std::endl;
    const double viscosity = Parameters().getParam( "viscosity", 1.0 );
    debugStream << "  - viscosity: " << viscosity << std::endl;

    // analytical data
    typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
        VelocityFunctionSpaceType;
    VelocityFunctionSpaceType velocitySpace;

    typedef Force< VelocityFunctionSpaceType >
        AnalyticalForceType;
    AnalyticalForceType analyticalForce( viscosity , velocitySpace );

    typedef DirichletData< VelocityFunctionSpaceType >
        AnalyticalDirichletDataType;
    AnalyticalDirichletDataType analyticalDirichletData( velocitySpace );

    // model traits
    typedef Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    AnalyticalForceType,
                    AnalyticalDirichletDataType,
                    gridDim,
                    polOrder >
        StokesModelTraitsImp;
    typedef Dune::DiscreteStokesModelDefault< StokesModelTraitsImp >
        StokesModelImpType;

    // determine pows
//    typedef Dune::FieldVector< double, gridDim >
//        ConstVec;
//    const double minpow = -9; // all less or equal treatet as zero
    const int c11 = pow1;// > minpow ? std::pow( grid_width, pow1 ) : 0;
    const int d11 = pow2;// > minpow ? std::pow( grid_width, pow2 ) : 0;
    const int c12 = pow3;// > minpow ? std::pow( grid_width, pow3 ) : 0;
    const int d12 = pow4;// > minpow ? std::pow( grid_width, pow4 ) : 0;

//    debugStream << "  - flux constants" << std::endl
//                << "    C_11: " << c11 << std::endl
//                << "    C_12: " << c12 << std::endl
//                << "    D_11: " << d11 << std::endl
//                << "    D_12: " << d12 << std::endl;

    // model
    StokesModelImpType stokesModel( c11,
                                    c12,
                                    d11,
                                    d12,
                                    analyticalForce,
                                    analyticalDirichletData,
                                    viscosity  );

    // treat as interface
    typedef Dune::DiscreteStokesModelInterface< StokesModelTraitsImp >
        StokesModelType;

    // function wrapper for the solutions
    typedef typename StokesModelTraitsImp::DiscreteStokesFunctionSpaceWrapperType
        DiscreteStokesFunctionSpaceWrapperType;

    DiscreteStokesFunctionSpaceWrapperType
        discreteStokesFunctionSpaceWrapper( gridPart );

    typedef typename StokesModelTraitsImp::DiscreteStokesFunctionWrapperType
        DiscreteStokesFunctionWrapperType;

    DiscreteStokesFunctionWrapperType
        discreteStokesFunctionWrapper(  "wrapped",
                                        discreteStokesFunctionSpaceWrapper );

    // some info about the analytical data
    typedef typename StokesModelType::DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType
        DiscreteAnalyticalDataFunctionSpaceType;
    DiscreteAnalyticalDataFunctionSpaceType
        discreteAnalyticalDataFunctionSpace( gridPart );
    typedef typename StokesModelType::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
        DiscreteAnalyticalDataFunctionType;
    DiscreteAnalyticalDataFunctionType
        discreteAnalyticalForce(    "discrete_analytical_force",
                                    discreteAnalyticalDataFunctionSpace );
    DiscreteAnalyticalDataFunctionType
        discreteAnalyticalDirichletData(    "discrete analytical dirichlet data",
                                            discreteAnalyticalDataFunctionSpace );
    typedef Dune::L2Projection< double,
                                double,
                                AnalyticalForceType,
                                DiscreteAnalyticalDataFunctionType >
        AnalyticalForceProjectionType;
    AnalyticalForceProjectionType analyticalForceProjection( 0 );
    analyticalForceProjection( analyticalForce, discreteAnalyticalForce );
    typedef Dune::L2Projection< double,
                                double,
                                AnalyticalDirichletDataType,
                                DiscreteAnalyticalDataFunctionType >
        AnalyticalDirichletDataProjectionType;
    AnalyticalDirichletDataProjectionType analyticalDirichletDataProjection( 0 );
    analyticalDirichletDataProjection(  analyticalDirichletData,
                                        discreteAnalyticalDirichletData );
    double analyticalForceMin = 0.0;
    double analyticalForceMax = 0.0;
    Stuff::getMinMaxOfDiscreteFunction( discreteAnalyticalForce,
                                        analyticalForceMin,
                                        analyticalForceMax );
    debugStream << "  - force" << std::endl
                << "    min: " << std::sqrt( 2.0 ) * analyticalForceMin << std::endl
                << "    max: " << std::sqrt( 2.0 ) * analyticalForceMax << std::endl;
    double analyticalDirichletDataMin = 0.0;
    double analyticalDirichletDataMax = 0.0;
    Stuff::getMinMaxOfDiscreteFunction( discreteAnalyticalDirichletData,
                                        analyticalDirichletDataMin,
                                        analyticalDirichletDataMax );
    debugStream << "  - dirichlet data" << std::endl
                << "    min: " << std::sqrt( 2.0 ) * analyticalDirichletDataMin << std::endl
                << "    max: " << std::sqrt( 2.0 ) * analyticalDirichletDataMax << std::endl;

    // exact solution, for testing
    typedef Velocity< VelocityTraits< gridDim, VelocityFunctionSpaceType > >
        AnalyticalExactVelocityType;
    AnalyticalExactVelocityType analyticalExactVelocity( velocitySpace );
    typedef typename StokesModelType::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
        DiscreteVelocityFunctionType;
    typedef Dune::L2Projection< double,
                                double,
                                AnalyticalExactVelocityType,
                                DiscreteVelocityFunctionType >
        AnalyticalVelocityProjectionType;
    AnalyticalVelocityProjectionType analyticalVelocityProjection( 0 );
    DiscreteVelocityFunctionType discreteExactVelocity( "discrete exact velocity",
                                                        discreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );
    analyticalVelocityProjection(   analyticalExactVelocity,
                                    discreteExactVelocity );

    typedef Dune::FunctionSpace< double, double, gridDim, 1 >
        PressureFunctionSpaceType;
    PressureFunctionSpaceType pressureSpace;
    typedef Pressure< PressureTraits< gridDim, PressureFunctionSpaceType > >
        AnalyticalExactPressureType;
    AnalyticalExactPressureType analyticalExactPressure( pressureSpace );
    typedef typename StokesModelType::DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
        DiscretePressureFunctionType;
    typedef Dune::L2Projection< double,
                                double,
                                AnalyticalExactPressureType,
                                DiscretePressureFunctionType >
        AnalyticalPressureProjectionType;
    AnalyticalPressureProjectionType analyticalPressureProjection( 0 );
    DiscretePressureFunctionType discreteExactPressure( "discrete exact pressure",
                                                        discreteStokesFunctionSpaceWrapper.discretePressureSpace() );
    analyticalPressureProjection(   analyticalExactPressure,
                                    discreteExactPressure );

    DiscreteStokesFunctionWrapperType discreteExactSolutions(   discreteStokesFunctionSpaceWrapper,
                                                                discreteExactVelocity,
                                                                discreteExactPressure );

    /* ********************************************************************** *
     * initialize passes                                                      *
     * ********************************************************************** */
    infoStream << "\n- starting pass" << std::endl;

    typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
        StartPassType;
    StartPassType startPass;

    typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
        StokesPassType;
    StokesPassType stokesPass(  startPass,
                                stokesModel,
                                gridPart,
                                discreteStokesFunctionSpaceWrapper );

    discreteStokesFunctionWrapper.discretePressure().clear();
    discreteStokesFunctionWrapper.discreteVelocity().clear();

    stokesPass.apply( discreteExactSolutions, discreteStokesFunctionWrapper );

    /* ********************************************************************** *
     * Problem postprocessing (with profiler example)                         *
     * ********************************************************************** */
    infoStream << "\n- postprocesing" << std::endl;

//    typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
//        DiscreteVelocityFunctionType;
    DiscreteVelocityFunctionType& computedVelocity = discreteStokesFunctionWrapper.discreteVelocity();
//    typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
//        DiscretePressureFunctionType;
    DiscretePressureFunctionType& computedPressure = discreteStokesFunctionWrapper.discretePressure();
//Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
//Stuff::printDiscreteFunctionMatlabStyle( computedPressure, "p", matlabLogStream );
//Stuff::printDiscreteFunctionMatlabStyle( computedVelocity, "u", matlabLogStream );

    double computedVelocityMin = 0.0;
    double computedVelocityMax = 0.0;
    Stuff::getMinMaxOfDiscreteFunction( computedVelocity,
                                        computedVelocityMin,
                                        computedVelocityMax );
    debugStream << "  - computed velocity" << std::endl
                << "    min: " << std::sqrt( 2.0 ) * computedVelocityMin << std::endl
                << "    max: " << std::sqrt( 2.0 ) * computedVelocityMax << std::endl;

    double computedPressureMin = 0.0;
    double computedPressureMax = 0.0;

    Stuff::getMinMaxOfDiscreteFunction( computedPressure,
                                        computedPressureMin,
                                        computedPressureMax );
    debugStream << "  - computed pressure" << std::endl
                << "    min: " << std::sqrt( 2.0 ) * computedPressureMin << std::endl
                << "    max: " << std::sqrt( 2.0 ) * computedPressureMax << std::endl
                << std::endl;


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
    info.bfg = Parameters().getParam( "do-bfg", true );
    info.gridname = gridPart.grid().name();

    info.polorder_pressure = StokesModelTraitsImp::pressureSpaceOrder;
    info.polorder_sigma = StokesModelTraitsImp::sigmaSpaceOrder;
    info.polorder_velocity = StokesModelTraitsImp::velocitySpaceOrder;

    profiler().StopTiming( "Problem/Postprocessing" );

    return info;
}

