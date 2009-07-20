/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(UGGRID) && defined(DEBUG)
    #warning ("UGGRID in debug mode is likely to produce a segfault")
#endif

#include <vector>
#include <string>

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


RunInfo singleRun(  CollectiveCommunication& mpicomm,
                    int refine_level,
                    Dune::StabilizationCoefficients stabil_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients()  );


void BfgRun( CollectiveCommunication& mpicomm );
void RefineRun( CollectiveCommunication& mpicomm );
void StabRun( CollectiveCommunication& mpicomm );

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

    int err = 0;

    profiler().Reset( 9 ); //prepare for 9 single runs


    //a little trickery so felix doesn't scream at me for breaking any of his scripts/parameterfiles
    const int runtype = Parameters().getParam( "runtype", -1 ) != -1  ? Parameters().getParam( "runtype", -1 ) : Parameters().getParam( "multirun", true );
    switch ( runtype ) {
        case 1: {
            /** all four stab parameters are permutated in [ minpow ; maxpow ]
                inside an outer loop that increments the grid's refine level
            **/
            StabRun( mpicomm );
            break;
        }
        default:
        case 0: { //we don't do any variation here
            RefineRun( mpicomm );
            break;
        }
        case 2: { //bfg-tau variance run on set reflevel
            BfgRun( mpicomm );
            break;
        }
    } // end case
    //profiler output in current format is somewhat meaningless
    profiler().Output( mpicomm, 0, 0 );

    return err;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

void RefineRun( CollectiveCommunication& mpicomm )
{
    // coloumn headers for eoc table output
    const std::string errheaders[] = { "h", "el't","runtime","Velocity", "Pressure" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
    L2ErrorVector l2_errors;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::RefineOutput eoc_texwriter( errorColumnHeaders );

    int maxref = Parameters().getParam( "maxref", 0 );
    int minref = Parameters().getParam( "minref", 0 );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    for ( int ref = minref; ref <= maxref; ++ref ) {
        if ( per_run_log_target )
            Logger().SetPrefix( "dune_stokes_ref_"+Stuff::toString(ref) );
        const int refine_level = Dune::DGFGridInfo< GridType >::refineStepsForHalf()* ref;
        RunInfo info = singleRun( mpicomm, refine_level );
        l2_errors.push_back( info.L2Errors );
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );
        eoc_output.write( eoc_texwriter, ( ref >= maxref ) );
        profiler().NextRun( info.L2Errors[0] ); //finish this run
    }
}

void StabRun( CollectiveCommunication& mpicomm )
{
    // coloumn headers for eoc table output
    const std::string errheaders[] = { "h", "el't","runtime","C11","C12","D11","D12","Velocity", "Pressure" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
    L2ErrorVector l2_errors;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::EocOutput eoc_texwriter( errorColumnHeaders );

    // the min andmax exponent that will be used in stabilizing terms
    // ie c11 = ( h_^minpow ) .. ( h^maxpow )
    // pow = -2 is interpreted as 0
    // in non-variation context minpow is used for c12 exp and maxpow for d12,
    //      c11 and D11 are set to zero
    int maxpow = Parameters().getParam( "maxpow", -1 );
    int minpow = Parameters().getParam( "minpow", 1 );
    double minfactor = Parameters().getParam( "minfactor", 0.0 );
    double maxfactor = Parameters().getParam( "maxfactor", 2.0 );
    double incfactor = Parameters().getParam( "incfactor", 0.5 );

    int ref = Parameters().getParam( "minref", 0 );

    const std::vector<std::string> coefficient_names = Dune::StabilizationCoefficients::getCoefficientNames();
    std::vector<std::string>::const_iterator name_it = coefficient_names.begin();

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    for ( int pow = minpow; pow <= maxpow; ++pow ) {
        for ( double factor = minfactor; factor <= maxfactor; factor+=incfactor ) {
            for ( ; name_it != coefficient_names.end(); ++name_it ) {
                if ( per_run_log_target ) { //sets unique per run filename if requested
                    std::string ff = "matlab__pow1_" + Stuff::toString( pow ) + "_pow2_" + Stuff::toString( factor ); //! not correct
                    Logger().SetPrefix( ff );
                }

                Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
                //do some matlab magic to suppress errors and display a little info baout current params
                matlabLogStream <<std::endl<< "\nclear;\ntry\ntic;warning off all;" << std::endl;
                matlabLogStream << "disp(' c/d-11/12 h powers: ";
                //coeff.print( matlabLogStream );
                matlabLogStream << "') \n" << std::endl;

                //actual work
                const int refine_level = Dune::DGFGridInfo< GridType >::refineStepsForHalf()* ref;
//                Dune::StabilizationCoefficients coeff =
                RunInfo info = singleRun( mpicomm, refine_level );

                // new line and closing try/catch in m-file
                matlabLogStream << "\ncatch\ndisp('errors here');\nend\ntoc\ndisp(' ');\n" << std::endl;
                l2_errors.push_back( info.L2Errors );

                //push errors to eoc-outputter class
                eoc_output.setErrors( idx,info.L2Errors );
                eoc_texwriter.setInfo( info );
                bool lastrun = ( pow >= maxpow && !(factor < maxfactor) && name_it == coefficient_names.end() );
                //the writer needs to know if it should close the table etc.
                eoc_output.write( eoc_texwriter, lastrun );

                profiler().NextRun( info.L2Errors[0] ); //finish this run
            }
        }
    }
}

void BfgRun( CollectiveCommunication& mpicomm )
{
    //first up a non-bfg run for reference
    const int refine_level = Dune::DGFGridInfo< GridType >::refineStepsForHalf()* Parameters().getParam( "minref", 0 );
    Parameters().setParam( "do-bfg", false );
    RunInfo nobfg_info = singleRun( mpicomm, refine_level );

    L2ErrorVector l2_errors;
    const std::string bfgheaders[] = { "h", "el't","runtime","$\\tau$","avg inner","min inner","max inner","total outer","max inner acc.","Velocity", "Pressure" };
    const unsigned int num_bfgheaders = sizeof ( bfgheaders ) / sizeof ( bfgheaders[0] );
    ColumnHeaders bfgColumnHeaders ( bfgheaders, bfgheaders + num_bfgheaders ) ;
    l2_errors.push_back( nobfg_info.L2Errors );

    Dune::FemEoc& bfg_output = Dune::FemEoc::instance( );
    bfg_output.initialize( "data","eoc-file", "bfg-desc", "eoc-template.tex" );
    size_t idx = bfg_output.addEntry( bfgColumnHeaders );
    Stuff::BfgOutput bfg_texwriter( bfgColumnHeaders, nobfg_info );
    bfg_output.setErrors( idx,nobfg_info.L2Errors );
    bfg_texwriter.setInfo( nobfg_info );
    bfg_output.write( bfg_texwriter, false );
    Parameters().setParam( "do-bfg", true );

    const double start_tau = Parameters().getParam( "bfg-tau-start", 0.01 ) ;
    const double stop_tau = Parameters().getParam( "bfg-tau-stop", 0.4 ) ;
    const double tau_inc = Parameters().getParam( "bfg-tau-inc", 0.025 ) ;

    for ( double tau = start_tau; tau < stop_tau; tau += tau_inc ) {
        Parameters().setParam( "bfg-tau", tau );

        RunInfo info = singleRun( mpicomm, refine_level );
        l2_errors.push_back( info.L2Errors );
        bfg_output.setErrors( idx,info.L2Errors );
        bfg_texwriter.setInfo( info );
        bfg_output.write( bfg_texwriter, !( tau + tau_inc < stop_tau ) );
        profiler().NextRun( info.L2Errors[0] ); //finish this run
    }
}

RunInfo singleRun(  CollectiveCommunication& mpicomm,
                    int refine_level,
                    Dune::StabilizationCoefficients stabil_coeff )
{
    profiler().StartTiming( "SingleRun" );
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    RunInfo info;

    debugStream << "\nsingleRun( ";
    stabil_coeff.print( debugStream );
    debugStream << " )" << std::endl;

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\n- initialising grid" << std::endl;
    Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename() );
    gridPtr->globalRefine( refine_level );
    typedef Dune::AdaptiveLeafGridPart< GridType >
        GridPartType;
    GridPartType gridPart( *gridPtr );
    const int gridDim = GridType::dimensionworld;
    info.codim0 = gridPtr->size( 0 );
    info.codim0 = gridPart.grid().size( 0 );
    Dune::GridWidthProvider< GridType > gw ( *gridPtr );
    double grid_width = gw.gridWidth();
    infoStream << "  - max grid width: " << grid_width << std::endl;
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


    // model traits
    typedef Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    Force,
                    DirichletData,
                    gridDim,
                    polOrder >
        StokesModelTraitsImp;
    typedef Dune::DiscreteStokesModelDefault< StokesModelTraitsImp >
        StokesModelImpType;

    // treat as interface
    typedef Dune::DiscreteStokesModelInterface< StokesModelTraitsImp >
        StokesModelType;

    // function wrapper for the solutions
    typedef StokesModelTraitsImp::DiscreteStokesFunctionSpaceWrapperType
        DiscreteStokesFunctionSpaceWrapperType;

    DiscreteStokesFunctionSpaceWrapperType
        discreteStokesFunctionSpaceWrapper( gridPart );

    typedef StokesModelTraitsImp::DiscreteStokesFunctionWrapperType
        DiscreteStokesFunctionWrapperType;

    DiscreteStokesFunctionWrapperType
        computedSolutions(  "computed_",
                                        discreteStokesFunctionSpaceWrapper );

    DiscreteStokesFunctionWrapperType initArgToPass( "init_", discreteStokesFunctionSpaceWrapper );

     typedef StokesModelTraitsImp::AnalyticalForceType
         AnalyticalForceType;
     AnalyticalForceType analyticalForce( viscosity , discreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );

     typedef StokesModelTraitsImp::AnalyticalDirichletDataType
         AnalyticalDirichletDataType;
     AnalyticalDirichletDataType analyticalDirichletData( discreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );

    StokesModelImpType stokesModel( stabil_coeff,
                                    analyticalForce,
                                    analyticalDirichletData,
                                    viscosity  );

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

    computedSolutions.discretePressure().clear();
    computedSolutions.discreteVelocity().clear();

    profiler().StartTiming( "Pass:apply" );
    stokesPass.apply( initArgToPass, computedSolutions );
    profiler().StopTiming( "Pass:apply" );
    info.run_time = profiler().GetTiming( "Pass:apply" );
    stokesPass.getRuninfo( info );

    /* ********************************************************************** *
     * Problem postprocessing
     * ********************************************************************** */
    infoStream << "\n- postprocesing" << std::endl;


    profiler().StartTiming( "Problem/Postprocessing" );

#ifndef COCKBURN_PROBLEM //bool tpl-param toggles ana-soltion output in post-proc
    typedef Problem< gridDim, DiscreteStokesFunctionWrapperType, false >
        ProblemType;
#else
    typedef Problem< gridDim, DiscreteStokesFunctionWrapperType, true >
        ProblemType;
#endif
    ProblemType problem( viscosity , computedSolutions );

    typedef PostProcessor< StokesPassType, ProblemType >
        PostProcessorType;

    PostProcessorType postProcessor( discreteStokesFunctionSpaceWrapper, problem );

    postProcessor.save( *gridPtr, computedSolutions, refine_level );
    info.L2Errors = postProcessor.getError();
    info.c11 = stabil_coeff.Power( "C11" );
    info.c12 = stabil_coeff.Power( "C12" );
    info.d11 = stabil_coeff.Power( "D11" );
    info.d12 = stabil_coeff.Power( "D12" );
    info.bfg = Parameters().getParam( "do-bfg", true );
    info.gridname = gridPart.grid().name();

    info.polorder_pressure = StokesModelTraitsImp::pressureSpaceOrder;
    info.polorder_sigma = StokesModelTraitsImp::sigmaSpaceOrder;
    info.polorder_velocity = StokesModelTraitsImp::velocitySpaceOrder;

    info.solver_accuracy = Parameters().getParam( "absLimit", 1e-4 );
    info.bfg_tau = Parameters().getParam( "bfg-tau", 0.1 );

    profiler().StopTiming( "Problem/Postprocessing" );
    profiler().StopTiming( "SingleRun" );

    return info;
}

