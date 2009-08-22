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

#if ! defined(POLORDER)
    #define POLORDER 0
    #warning ("using default polorder 0 for all spaces")
#endif

#if ! defined(PRESSURE_POLORDER)
    #define PRESSURE_POLORDER POLORDER
#endif

#if ! defined(VELOCITY_POLORDER)
    #define VELOCITY_POLORDER POLORDER
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

//! return type of getXXX_Permutations()
typedef std::vector<Dune::StabilizationCoefficients>
    CoeffVector;

//! used in all runs to store L2 errors across runs, but i forgot why...
typedef std::vector< RunInfo >
    RunInfoVector;

//! the strings used for column headers in tex output
typedef std::vector<std::string>
    ColumnHeaders;


/** \brief one single application of the discretisation and solver

    \param  mpicomm
            mostly useless atm, but mandatory
    \param  refine_level_factor
            integer to be multiplied by Dune::DGFGridInfo< GridType >::refineStepsForHalf()
            to get the used refine level for the constructed grid
    \param  stabil_coeff
            the set of coefficients to be used in the run. Default is used in all run types but StabRun().

**/
RunInfo singleRun(  CollectiveCommunication& mpicomm,
                    int refine_level_factor,
                    Dune::StabilizationCoefficients stabil_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients()  );

//! multiple runs with bfg tau set in [bfg-tau-start : bfg-tau-inc : bfg-tau-stop] and everything else on default
void BfgRun     ( CollectiveCommunication& mpicomm );
//! multiple runs with refine level in [minref*refineStepsForHalf : maxref*refineStepsForHalf] and everything else on default
void RefineRun  ( CollectiveCommunication& mpicomm );
//! LOTS of runs where any combination of stabilisation coefficients getXXX_Permutations() yields is used ( currently hardcoded to getC_power_Permutations )
void StabRun    ( CollectiveCommunication& mpicomm );
//! multiple runs with  set in [accuracy_start :  *=accuracy_factor : accuracy_stop] and everything else on default
void AccuracyRun( CollectiveCommunication& mpicomm );

//! brute force all permutations
CoeffVector getAll_Permutations();
//! get only permutations for C_11 and C_12
CoeffVector getC_Permutations();
//! get only the permutations in power for C_11 and C_12
CoeffVector getC_power_Permutations();

/**
 *  \brief  main function
 *
 *  ParameterContainer and Logger setup, select run type from parameter file and go
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
    if ( argc < 2 ) {
        std::cerr << "\nUsage: " << argv[0] << " parameterfile \n" << "\n\t --- OR --- \n";
        std::cerr << "\nUsage: " << argv[0] << " paramfile:"<<"file" << " more-opts:val ..." << std::endl;
        Parameters().PrintParameterSpecs( std::cerr );
        std::cerr << std::endl;
        return 2;
    }
    if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }



    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
                     Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
                     Parameters().getParam( "fem.io.logdir",    std::string(),              useLogger )
                    );

    int err = 0;

    //a little trickery so felix doesn't scream at me for breaking any of his scripts/parameterfiles
    const int runtype = Parameters().getParam( "runtype", -1 ) != -1  ? Parameters().getParam( "runtype", -1 ) : Parameters().getParam( "multirun", true );
    switch ( runtype ) {
        case 1: {
            StabRun( mpicomm );
            break;
        }
        default:
        case 0: {
            RefineRun( mpicomm );
            break;
        }
        case 2: {
            BfgRun( mpicomm );
            break;
        }
        case 3: {
            AccuracyRun( mpicomm );
            break;
        }
    } // end case

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
    Logger().Info() << "starting refine run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","runtime","Velocity", "Pressure" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
    RunInfoVector run_infos;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::RefineOutput eoc_texwriter( errorColumnHeaders );

    int maxref = Parameters().getParam( "maxref", 0 );
    int minref = Parameters().getParam( "minref", 0 );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    profiler().Reset( maxref - minref + 1 );
    for ( int ref = minref; ref <= maxref; ++ref ) {
        if ( per_run_log_target )
            Logger().SetPrefix( "dune_stokes_ref_"+Stuff::toString(ref) );

        RunInfo info = singleRun( mpicomm, ref );
        run_infos.push_back( info );
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );
        eoc_output.write( eoc_texwriter, ( ref >= maxref ) );
        profiler().NextRun(); //finish this run
    }
    profiler().Output( mpicomm, run_infos );
}

void AccuracyRun( CollectiveCommunication& mpicomm )
{
    Logger().Info() << "starting accurracy run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","runtime (s)","avg. inner iter.","inner acc.","total outer iter.","outer acc.","Velocity", "Pressure" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
    RunInfoVector run_infos;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::AccurracyOutput  eoc_texwriter( errorColumnHeaders );

    double  accurracy_start  = Parameters().getParam( "accurracy_start", 10e-3 );
    int     accurracy_steps  = Parameters().getParam( "accurracy_steps", 5 );
    double  accurracy_factor = Parameters().getParam( "accurracy_factor", 10e-3 );

    Logger().Dbg() << " accurracy_start: " <<  accurracy_start
            << " accurracy_steps: " <<  accurracy_steps
            << " accurracy_factor: " <<  accurracy_factor << std::endl;

    int ref = Parameters().getParam( "minref", 0 );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    int numruns = std::pow( accurracy_steps, 2.0 );
    profiler().Reset( numruns );
    for ( int i = 0; i < accurracy_steps; ++i ) {
        for ( int j = 0; j < accurracy_steps; ++j ) {

            if ( per_run_log_target )
                Logger().SetPrefix( "dune_stokes_ref_"+Stuff::toString(ref) );

            double inner_acc = accurracy_start * std::pow( accurracy_factor, double( j ) );
            double outer_acc = accurracy_start * std::pow( accurracy_factor, double( i ) );
            Parameters().setParam( "absLimit", outer_acc );
            Parameters().setParam( "inner_absLimit", inner_acc );
            RunInfo info = singleRun( mpicomm, ref );
            run_infos.push_back( info );
            eoc_output.setErrors( idx,info.L2Errors );
            eoc_texwriter.setInfo( info );
            numruns--;
            eoc_output.write( eoc_texwriter, !numruns );
            profiler().NextRun(); //finish this run

            Logger().Info() << numruns << " runs remaining" << std::endl;
//            if ( numruns<=0 ) //yeah, stop conditions are a bit flaky
//                break;
        }
//        if ( numruns<=0 )
//            break;
    }
    profiler().Output( mpicomm, run_infos );
}

void StabRun( CollectiveCommunication& mpicomm )
{
    Logger().Info() << "starting stab run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","runtime","C11","C12","D11","D12","Velocity", "Pressure" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
    RunInfoVector run_infos;
    Dune::FemEoc& eoc_output = Dune::FemEoc::instance( );
    eoc_output.initialize( "data","eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    Stuff::EocOutput eoc_texwriter( errorColumnHeaders );

    int ref = Parameters().getParam( "minref", 0 );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = Parameters().getParam( "per-run-log-target", true );

    CoeffVector coeff_vector = getC_power_Permutations();

    Logger().Info().Resume();
    Logger().Info() << "beginning " << coeff_vector.size() << " runs" << std::endl ;
    profiler().Reset( coeff_vector.size() );
    CoeffVector::const_iterator it = coeff_vector.begin();
    for ( ; it != coeff_vector.end(); ++it ) {
        RunInfo info = singleRun( mpicomm, ref, *it );
        run_infos.push_back( info );

        //push errors to eoc-outputter class
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );

        //the writer needs to know if it should close the table etc.
        eoc_output.write( eoc_texwriter, ( (it+1) == coeff_vector.end() ) );

        profiler().NextRun(); //finish this run
    }
    Logger().Info().Resume();
    Logger().Info() << "completed " << coeff_vector.size() << " runs" << std::endl ;

    profiler().Output( mpicomm, run_infos );
}

void BfgRun( CollectiveCommunication& mpicomm )
{
    //first up a non-bfg run for reference
    const int refine_level_factor = Parameters().getParam( "minref", 0 );
    Parameters().setParam( "do-bfg", false );
    RunInfo nobfg_info = singleRun( mpicomm, refine_level_factor );

    RunInfoVector run_infos;
    const std::string bfgheaders[] = { "h", "el't","runtime","$\\tau$","avg inner","min inner","max inner","total outer","max inner acc.","Velocity", "Pressure" };
    const unsigned int num_bfgheaders = sizeof ( bfgheaders ) / sizeof ( bfgheaders[0] );
    ColumnHeaders bfgColumnHeaders ( bfgheaders, bfgheaders + num_bfgheaders ) ;
    run_infos.push_back( nobfg_info );

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

    profiler().Reset( ( stop_tau - start_tau ) / tau_inc + 1 );

    for ( double tau = start_tau; tau < stop_tau; tau += tau_inc ) {
        Parameters().setParam( "bfg-tau", tau );

        RunInfo info = singleRun( mpicomm, refine_level_factor );
        run_infos.push_back( info );
        bfg_output.setErrors( idx,info.L2Errors );
        bfg_texwriter.setInfo( info );
        bfg_output.write( bfg_texwriter, !( tau + tau_inc < stop_tau ) );
        profiler().NextRun(); //finish this run
    }
    profiler().Output( mpicomm, run_infos );
}

RunInfo singleRun(  CollectiveCommunication& mpicomm,
                    int refine_level_factor,
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
    const int gridDim = GridType::dimensionworld;
    Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename( gridDim ) );
    const int refine_level = refine_level_factor * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
    gridPtr->globalRefine( refine_level );
    typedef Dune::AdaptiveLeafGridPart< GridType >
        GridPartType;
    GridPartType gridPart( *gridPtr );
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
                    polOrder,
                    VELOCITY_POLORDER,
                    PRESSURE_POLORDER >
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

    profiler().StartTiming( "Pass -- APPLY" );
    stokesPass.apply( initArgToPass, computedSolutions );
    profiler().StopTiming( "Pass -- APPLY" );
    info.run_time = profiler().GetTiming( "Pass -- APPLY" );
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
    typedef Dune::StabilizationCoefficients::ValueType
        Pair;
    info.c11 = Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
    info.c12 = Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
    info.d11 = Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
    info.d12 = Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
    info.bfg = Parameters().getParam( "do-bfg", true );
    info.gridname = gridPart.grid().name();
    info.refine_level = refine_level;

    info.polorder_pressure = StokesModelTraitsImp::pressureSpaceOrder;
    info.polorder_sigma = StokesModelTraitsImp::sigmaSpaceOrder;
    info.polorder_velocity = StokesModelTraitsImp::velocitySpaceOrder;

    info.solver_accuracy = Parameters().getParam( "absLimit", 1e-4 );
    info.inner_solver_accuracy = Parameters().getParam( "inner_absLimit", 1e-4 );
    info.bfg_tau = Parameters().getParam( "bfg-tau", 0.1 );

    profiler().StopTiming( "Problem/Postprocessing" );
    profiler().StopTiming( "SingleRun" );

    return info;
}

CoeffVector getAll_Permutations() {
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

    CoeffVector coeff_vector;
    for ( int c11_pow = minpow; c11_pow <= maxpow; ++c11_pow ) {
        for ( double c11_factor = minfactor; c11_factor <= maxfactor; c11_factor+=incfactor ) {
            for ( int c12_pow = minpow; c12_pow <= maxpow; ++c12_pow ) {
                for ( double c12_factor = minfactor; c12_factor <= maxfactor; c12_factor+=incfactor ) {
                    for ( int d11_pow = minpow; d11_pow <= maxpow; ++d11_pow ) {
                        for ( double d11_factor = minfactor; d11_factor <= maxfactor; d11_factor+=incfactor ) {
                            for ( int d12_pow = minpow; d12_pow <= maxpow; ++d12_pow ) {
                                for ( double d12_factor = minfactor; d12_factor <= maxfactor; d12_factor+=incfactor ) {
                                    Dune::StabilizationCoefficients c ( c11_pow, c12_pow, d11_pow, d12_pow,
                                                                        c11_factor, c12_factor, d11_factor, d12_factor );
                                    coeff_vector.push_back( c );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return coeff_vector;
}

CoeffVector getC_Permutations(){
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

    typedef Dune::StabilizationCoefficients::ValueType::first_type
        PowerType;
    typedef Dune::StabilizationCoefficients::ValueType::second_type
        FactorType;

    const PowerType invalid_power = Dune::StabilizationCoefficients::invalid_power;
    const FactorType invalid_factor = Dune::StabilizationCoefficients::invalid_factor;

    CoeffVector coeff_vector;
    for ( PowerType c11_pow = minpow; c11_pow <= maxpow; ++c11_pow ) {
        for ( FactorType c11_factor = minfactor; c11_factor <= maxfactor; c11_factor+=incfactor ) {
            for ( PowerType c12_pow = minpow; c12_pow <= maxpow; ++c12_pow ) {
                for ( FactorType c12_factor = minfactor; c12_factor <= maxfactor; c12_factor+=incfactor ) {
                    Dune::StabilizationCoefficients c ( c11_pow, c12_pow, invalid_power, invalid_power,
                                                        c11_factor, c12_factor, invalid_factor, invalid_factor );
                    coeff_vector.push_back( c );
                }
            }
        }
    }
    return coeff_vector;
}

CoeffVector getC_power_Permutations(){
        // the min andmax exponent that will be used in stabilizing terms
    // ie c11 = ( h_^minpow ) .. ( h^maxpow )
    // pow = -2 is interpreted as 0
    // in non-variation context minpow is used for c12 exp and maxpow for d12,
    //      c11 and D11 are set to zero
    int maxpow = Parameters().getParam( "maxpow", -1 );
    int minpow = Parameters().getParam( "minpow", 1 );
    double minfactor = Parameters().getParam( "minfactor", 0.0 );

    typedef Dune::StabilizationCoefficients::ValueType::first_type
        PowerType;
    typedef Dune::StabilizationCoefficients::ValueType::second_type
        FactorType;

    const PowerType invalid_power = Dune::StabilizationCoefficients::invalid_power;
    const FactorType invalid_factor = Dune::StabilizationCoefficients::invalid_factor;

    CoeffVector coeff_vector;
    for ( PowerType c11_pow = minpow; c11_pow <= maxpow; ++c11_pow ) {
            for ( PowerType c12_pow = minpow; c12_pow <= maxpow; ++c12_pow ) {
             Dune::StabilizationCoefficients c ( c11_pow, c12_pow, invalid_power, invalid_power,
                                                minfactor, minfactor, invalid_factor, invalid_factor );
            coeff_vector.push_back( c );
        }
    }
    return coeff_vector;
}
