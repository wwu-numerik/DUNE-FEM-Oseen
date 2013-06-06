/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#include "cmake_config.h"

#include <cstdio>
#include <vector>
#include <string>

#include <iostream>
#include <cmath>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/mpimanager.hh>


#include <dune/grid/utility/gridtype.hh>
typedef Dune::GridSelector::GridType
    GridType;

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/oseen/functionspacewrapper.hh>
#include <dune/fem/oseen/modelinterface.hh>
#include <dune/fem/oseen/pass.hh>
#include <dune/fem/oseen/boundarydata.hh>
#include <dune/fem/oseen/runinfo.hh>

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/fem/oseen/postprocessing.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/signals.hh>
#include <dune/fem/oseen/tex.hh>

#include <dune/fem/oseen/problems.hh>
#include <dune/stuff/fem/femeoc.hh>

#include "analyticaldata.hh"
#include "velocity.hh"
#include "pressure.hh"
#include "estimator.hh"

#ifndef COMMIT
    #define COMMIT "undefined"
#endif

static const std::string commit_string (COMMIT);

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

//! return type of getXXX_Permutations()
typedef std::vector<Dune::StabilizationCoefficients>
    CoeffVector;

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
DSC::RunInfo singleRun(CollectiveCommunication& mpicomm,
                    int refine_level_factor,
                    const Dune::StabilizationCoefficients &stabil_coeff );

//! multiple runs with bfg tau set in [bfg-tau-start : bfg-tau-inc : bfg-tau-stop] and everything else on default
void BfgRun     ( CollectiveCommunication& mpicomm );
//! multiple runs with refine level in [minref*refineStepsForHalf : maxref*refineStepsForHalf] and everything else on default
void RefineRun  ( CollectiveCommunication& mpicomm );
//! LOTS of runs where any combination of stabilisation coefficients getXXX_Permutations() yields is used ( currently hardcoded to getC_power_Permutations )
void StabRun    ( CollectiveCommunication& mpicomm );
//! multiple runs with  set in [accuracy_start :  *=accuracy_factor : accuracy_stop] and everything else on default
void AccuracyRun( CollectiveCommunication& mpicomm );
//! multiple runs with  set in [accuracy_start :  *=accuracy_factor : accuracy_stop] and everything else on default (only outer acc is varied)
void AccuracyRunOuter( CollectiveCommunication& mpicomm );

//! display last computed pressure/velo with grape
int display( int argc, char** argv );

//! brute force all permutations
CoeffVector getAll_Permutations();
//! get only permutations for C_11 and C_12
CoeffVector getC_Permutations();
//! get only the permutations in power for C_11 and C_12
CoeffVector getC_power_Permutations();

//! output alert for neg. EOC
void eocCheck( const DSC::RunInfoVector& runInfos );

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
    DSC::installSignalHandler();
  try{

	Dune::MPIManager::initialize(argc, argv);
	//assert( Dune::Capabilities::isParallel< GridType >::v );
    CollectiveCommunication mpicomm;//Dune::MPIManager::helper().getCommunicator() );

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    if ( argc < 2 ) {
        std::cerr << "\nUsage: " << argv[0] << " parameterfile \n" << "\n\t --- OR --- \n";
        std::cerr << "\nUsage: " << argv[0] << " paramfile:"<<"file" << " more-opts:val ..." << std::endl;
        std::cerr << std::endl;
        return 2;
    }

    DSC_CONFIG.readCommandLine(argc, argv);

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    DSC_LOG.create( DSC_CONFIG_GETB( "loglevel",         62,                         useLogger ),
                     DSC_CONFIG_GETB( "logfile",          std::string("dune_stokes"), useLogger ),
                     DSC_CONFIG_GETB( "fem.io.datadir",   std::string(),              useLogger )
                    );

    int err = 0;

	const int runtype = DSC_CONFIG_GET( "runtype", 5 );
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
        case 4: {
            AccuracyRunOuter( mpicomm );
            break;
        }
        case 5: {
            DSC_PROFILER.reset( 1 );
			DSC::RunInfoVector rf;
			Dune::StabilizationCoefficients st = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
			st.FactorFromParams( "C11" );
			st.FactorFromParams( "C12" );
			st.FactorFromParams( "D11" );
			st.FactorFromParams( "D12" );
			rf.push_back(singleRun( mpicomm, DSC_CONFIG_GET( "minref", 0 ), st ) );
            DSC_PROFILER.outputTimings();
			DSC::dumpRunInfoVectorToFile( rf );
            break;
        }
    } // end case

    DSC_LOG_ERROR << "\nRun from: " << commit_string << std::endl;
    return err;
  }



  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch ( std::bad_alloc& b ) {
      std::cerr << "Memory allocation failed: " << b.what() ;
      DSC_LOG_INFO.resume();
      DSC::meminfo( DSC_LOG_INFO );
  }
  catch ( std::runtime_error& a ) {
      std::cerr << "Runtime error:\n" << a.what() << std::endl ;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

void RefineRun( CollectiveCommunication& mpicomm )
{
#if !(ENABLE_ADAPTIVE)
	throw std::runtime_error("refine runs don't work with adaptation disabled");
#endif
    DSC_LOG_INFO << "starting refine run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","Laufzeit (s)","Geschwindigkeit", "Druck" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
	DSC::RunInfoVector run_infos;
    auto& eoc_output = DSFe::FemEoc::instance( );
	eoc_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    DSC::RefineOutput eoc_texwriter( errorColumnHeaders );

    int minref = DSC_CONFIG_GET( "minref", 0 );
	// ensures maxref>=minref
	const int maxref = DSC::clamp( DSC_CONFIG_GET( "maxref", 0 ), minref, DSC_CONFIG_GET( "maxref", 0 ) );

	Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
	if ( DSC_CONFIG_GET( "custom_stabilization_coefficients", false ) ) {
		stab_coeff.FactorFromParams( "D11" );
	}
    // setting this to true will give each run a unique logilfe name
	bool per_run_log_target = DSC_CONFIG_GET( "per-run-log-target", true );

    DSC_PROFILER.reset( maxref - minref + 1 );
    for ( int ref = minref; ref <= maxref; ++ref ) {
        if ( per_run_log_target )
            DSC_LOG.setPrefix( "dune_stokes_ref_"+DSC::toString(ref) );

		DSC::RunInfo info = singleRun( mpicomm, ref, stab_coeff );
        run_infos.push_back( info );
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );
        eoc_output.write( eoc_texwriter, ( ref >= maxref ) );
        DSC_PROFILER.nextRun(); //finish this run
    }
	run_infos[0].cumulative_run_time = run_infos[0].run_time;
	for ( size_t i = 1; i < run_infos.size(); ++i )
			run_infos[i].cumulative_run_time = run_infos[i-1].cumulative_run_time + run_infos[i].run_time;
    DSC_PROFILER.outputTimings();
	DSC::dumpRunInfoVectorToFile( run_infos );
	eocCheck( run_infos );
}

void AccuracyRun( CollectiveCommunication& mpicomm )
{
    DSC_LOG_INFO << "starting accurracy run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","Laufzeit (s)","\\o{} Iter. (innen)","Genauigkeit (innen)","\\# Iter. (aussen)","Genauigkeit (aussen)","Geschwindigkeit", "Druck" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
	DSC::RunInfoVector run_infos;
    auto& eoc_output = DSFe::FemEoc::instance( );
	eoc_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    DSC::AccurracyOutput  eoc_texwriter( errorColumnHeaders );

    double  accurracy_start  = DSC_CONFIG_GET( "accurracy_start", 10e-3 );
    int     accurracy_steps  = DSC_CONFIG_GET( "accurracy_steps", 5 );
    double  accurracy_factor = DSC_CONFIG_GET( "accurracy_factor", 10e-3 );

    DSC_LOG_DEBUG << " accurracy_start: " <<  accurracy_start
            << " accurracy_steps: " <<  accurracy_steps
            << " accurracy_factor: " <<  accurracy_factor << std::endl;

    int ref = DSC_CONFIG_GET( "minref", 0 );

	int numruns = int( std::pow( accurracy_steps, 2.0 ) );
    DSC_PROFILER.reset( numruns );
    for ( int i = 0; i < accurracy_steps; ++i ) {
        for ( int j = 0; j < accurracy_steps; ++j ) {

            double inner_acc = accurracy_start * std::pow( accurracy_factor, double( j ) );
            double outer_acc = accurracy_start * std::pow( accurracy_factor, double( i ) );
            DSC_CONFIG.set( "absLimit", outer_acc );
            DSC_CONFIG.set( "inner_absLimit", inner_acc );
			Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
			DSC::RunInfo info = singleRun( mpicomm, ref, stab_coeff );
            run_infos.push_back( info );
            eoc_output.setErrors( idx,info.L2Errors );
            eoc_texwriter.setInfo( info );
            numruns--;
            eoc_output.write( eoc_texwriter, !numruns );
            DSC_PROFILER.nextRun(); //finish this run

            DSC_LOG_INFO << numruns << " runs remaining" << std::endl;
        }
    }
    DSC_PROFILER.outputTimings();
	DSC::dumpRunInfoVectorToFile( run_infos );
}

void AccuracyRunOuter( CollectiveCommunication& mpicomm )
{
    DSC_LOG_INFO << "starting accurracyOuter run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","Laufzeit (s)","\\o{} Iter. (innen)","\\# Iter. (aussen)","Genauigkeit (aussen)","Geschwindigkeit", "Druck" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
	DSC::RunInfoVector run_infos;
    auto& eoc_output = DSFe::FemEoc::instance( );
	eoc_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    DSC::AccurracyOutputOuter  eoc_texwriter( errorColumnHeaders );

    double  accurracy_start  = DSC_CONFIG_GET( "accurracy_start", 10e-3 );
    int     accurracy_steps  = DSC_CONFIG_GET( "accurracy_steps", 5 );
    double  accurracy_factor = DSC_CONFIG_GET( "accurracy_factor", 10e-3 );

    DSC_LOG_DEBUG << " accurracy_start: " <<  accurracy_start
            << " accurracy_steps: " <<  accurracy_steps
            << " accurracy_factor: " <<  accurracy_factor << std::endl;

    int ref = DSC_CONFIG_GET( "minref", 0 );

    // setting this to true will give each run a unique logilfe name
    bool per_run_log_target = DSC_CONFIG_GET( "per-run-log-target", true );

    int numruns = accurracy_steps;
    DSC_PROFILER.reset( numruns );
    for ( int i = 0; i < accurracy_steps; ++i ) {
        if ( per_run_log_target )
            DSC_LOG.setPrefix( "dune_stokes_ref_"+DSC::toString(ref) );

        double outer_acc = accurracy_start * std::pow( accurracy_factor, double( i ) );
        double inner_acc = outer_acc;
        DSC_CONFIG.set( "absLimit", outer_acc );
        DSC_CONFIG.set( "inner_absLimit", inner_acc );
		Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
		DSC::RunInfo info = singleRun( mpicomm, ref, stab_coeff );
        run_infos.push_back( info );
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );
        numruns--;
        eoc_output.write( eoc_texwriter, !numruns );
        DSC_PROFILER.nextRun(); //finish this run

        DSC_LOG_INFO << numruns << " runs remaining" << std::endl;
    }
    DSC_PROFILER.outputTimings();
	DSC::dumpRunInfoVectorToFile( run_infos );
}

void StabRun( CollectiveCommunication& mpicomm )
{
    DSC_LOG_INFO << "starting stab run " << std::endl;
    // column headers for eoc table output
    const std::string errheaders[] = { "h", "el't","Laufzeit (s)","C11","C12","D11","D12","Geschwindigkeit", "Druck" };
    const unsigned int num_errheaders = sizeof ( errheaders ) / sizeof ( errheaders[0] );
    ColumnHeaders errorColumnHeaders ( errheaders, errheaders + num_errheaders ) ;
	DSC::RunInfoVector run_infos;
    auto& eoc_output = DSFe::FemEoc::instance( );
	eoc_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "eoc-desc", "eoc-template.tex" );
    size_t idx = eoc_output.addEntry( errorColumnHeaders );
    DSC::EocOutput eoc_texwriter( errorColumnHeaders );

    int ref = DSC_CONFIG_GET( "minref", 0 );

	// setting this to true will give each run a unique logfile name
    bool per_run_log_target = DSC_CONFIG_GET( "per-run-log-target", true );

    CoeffVector coeff_vector = getC_power_Permutations();

    DSC_LOG_INFO.resume();
    DSC_LOG_INFO << "beginning " << coeff_vector.size() << " runs" << std::endl ;
    DSC_PROFILER.reset( coeff_vector.size() );
	CoeffVector::iterator it = coeff_vector.begin();
	for ( unsigned i = 0; it != coeff_vector.end(); ++it,++i ) {
		DSC::RunInfo info = singleRun( mpicomm, ref, *it );
        run_infos.push_back( info );

        //push errors to eoc-outputter class
        eoc_output.setErrors( idx,info.L2Errors );
        eoc_texwriter.setInfo( info );

        //the writer needs to know if it should close the table etc.
        eoc_output.write( eoc_texwriter, ( (it+1) == coeff_vector.end() ) );

        DSC_PROFILER.nextRun(); //finish this run
		if ( per_run_log_target )
            DSC_LOG.setPrefix( "dune_stokes_permutation_"+DSC::toString(i) );
    }
    DSC_LOG_INFO.resume();
    DSC_LOG_INFO << "completed " << coeff_vector.size() << " runs" << std::endl ;

    DSC_PROFILER.outputTimings();
	DSC::dumpRunInfoVectorToFile( run_infos );
}

void BfgRun( CollectiveCommunication& mpicomm )
{
    //first up a non-bfg run for reference
    const int refine_level_factor = DSC_CONFIG_GET( "minref", 0 );
    DSC_CONFIG.set( "do-bfg", false );
    const Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
	DSC::RunInfo nobfg_info = singleRun( mpicomm, refine_level_factor, stab_coeff );
	nobfg_info.bfg_tau =  std::numeric_limits<double>::quiet_NaN();

	DSC::RunInfoVector run_infos;
    const std::string bfgheaders[] = { "h", "el't","Laufzeit (s)","$\\tau$","\\o{} Iter. (i)","min \\# Iter. (i)","max \\# Iter. (i)","\\# Iter. (a)","min. Genau. (i)","Geschwindigkeit", "Druck" };
    const unsigned int num_bfgheaders = sizeof ( bfgheaders ) / sizeof ( bfgheaders[0] );
    ColumnHeaders bfgColumnHeaders ( bfgheaders, bfgheaders + num_bfgheaders ) ;
    run_infos.push_back( nobfg_info );

    auto& bfg_output = DSFe::FemEoc::instance( );
	bfg_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "bfg-desc", "eoc-template.tex" );
    size_t idx = bfg_output.addEntry( bfgColumnHeaders );
    DSC::BfgOutput bfg_texwriter( bfgColumnHeaders, nobfg_info );
    bfg_output.setErrors( idx,nobfg_info.L2Errors );
    bfg_texwriter.setInfo( nobfg_info );
    bfg_output.write( bfg_texwriter, false );
    DSC_CONFIG.set( "do-bfg", true );

    const double start_tau = DSC_CONFIG_GET( "bfg-tau-start", 0.01 ) ;
    const double stop_tau = DSC_CONFIG_GET( "bfg-tau-stop", 0.4 ) ;
    const double tau_inc = DSC_CONFIG_GET( "bfg-tau-inc", 0.025 ) ;

    DSC_PROFILER.reset( int( ( stop_tau - start_tau ) / tau_inc + 1 ) );

    for ( double tau = start_tau; tau < stop_tau; tau += tau_inc ) {
        DSC_CONFIG.set( "bfg-tau", tau );
		DSC::RunInfo info = singleRun( mpicomm, refine_level_factor, stab_coeff );
        run_infos.push_back( info );
        bfg_output.setErrors( idx,info.L2Errors );
        bfg_texwriter.setInfo( info );
        bfg_output.write( bfg_texwriter, !( tau + tau_inc < stop_tau ) );
        DSC_PROFILER.nextRun(); //finish this run
    }
    DSC_PROFILER.outputTimings();
	DSC::dumpRunInfoVectorToFile( run_infos );
}

DSC::RunInfo singleRun(  CollectiveCommunication& /*mpicomm*/,
                    int refine_level_factor,
                    const Dune::StabilizationCoefficients& stabil_coeff )
{
    DSC_PROFILER.startTiming( "SingleRun" );
    auto& infoStream = DSC_LOG_INFO;
    auto& debugStream = DSC_LOG_DEBUG;
//	stabil_coeff.Add( "E12", 0.0 );
	DSC::RunInfo info;

    debugStream << "\nsingleRun( ";
    stabil_coeff.print( debugStream );
    debugStream << " )" << std::endl;

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
	infoStream << boost::format("\n- initialising grid (level %d)") % refine_level_factor << std::endl;
    const int gridDim = GridType::dimensionworld;
    static Dune::GridPtr< GridType > gridPtr( DSC_CONFIG.get( "dgf_file", "unitsquare.dgf" ) );
    static bool firstRun = true;
    int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	static int last_refine_level = refine_level;
    if ( firstRun && refine_level_factor > 0 ) { //since we have a couple of local statics, only do this once, further refinement done in estimator
        gridPtr->globalRefine( refine_level );
    }

    /* ********************************************************************** *
     * initialize problem                                                     *
     * ********************************************************************** */
    infoStream << "\n- initialising problem" << std::endl;

    const int polOrder = POLORDER;
    debugStream << "  - polOrder: " << polOrder << std::endl;
    const double viscosity = DSC_CONFIG_GET( "viscosity", 1.0 );
	info.viscosity = viscosity;
    debugStream << "  - viscosity: " << viscosity << std::endl;
	const double alpha = DSC_CONFIG_GET( "alpha", 0.0 );
	info.alpha = alpha;

    // analytical data

    // model traits
	#if 0 //defined( AORTA_PROBLEM )
    typedef Dune::DiscreteOseenModelDefaultTraits<
                    GridType,
					PROBLEM_NAMESPACE::Force,
					Dune::GeometryBasedBoundaryFunctionTraits<VariableDirichletData,FirstOrderBoundaryShapeFunction>,
                    gridDim,
                    polOrder,
                    VELOCITY_POLORDER,
                    PRESSURE_POLORDER >
        StokesModelTraitsImp;
	#else
    typedef Dune::DiscreteOseenModelDefaultTraits<
                    GridType,
					PROBLEM_NAMESPACE::Force,
					DefaultDirichletDataTraits<PROBLEM_NAMESPACE::DirichletData>,
                    gridDim,
                    polOrder,
                    VELOCITY_POLORDER,
                    PRESSURE_POLORDER >
        StokesModelTraitsImp;
	#endif
    static typename StokesModelTraitsImp::GridPartType gridPart( *gridPtr );

    typedef Dune::DiscreteOseenModelDefault< StokesModelTraitsImp >
        StokesModelImpType;

    // treat as interface
    typedef Dune::DiscreteOseenModelInterface< StokesModelTraitsImp >
        StokesModelType;

    // function wrapper for the solutions
    typedef StokesModelTraitsImp::DiscreteOseenFunctionSpaceWrapperType
        DiscreteOseenFunctionSpaceWrapperType;

    static DiscreteOseenFunctionSpaceWrapperType
        discreteStokesFunctionSpaceWrapper( gridPart );

    typedef StokesModelTraitsImp::DiscreteOseenFunctionWrapperType
        DiscreteOseenFunctionWrapperType;

    static DiscreteOseenFunctionWrapperType
        computedSolutions(  "computed_",
                            discreteStokesFunctionSpaceWrapper,
                            gridPart );
	DiscreteOseenFunctionWrapperType
		dummyFunctions(  "dummy_",
							discreteStokesFunctionSpaceWrapper,
							gridPart );
#if ENABLE_ADAPTIVE
    if ( !firstRun ) {
        Dune::Estimator<DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType>
            estimator ( computedSolutions.discretePressure() );
		for ( int i = refine_level - last_refine_level; i > 0; --i ) {
			estimator.mark( 0.0 /*dummy*/ ); //simpler would be to use real weights in mark(), but alas, that doesn't work as advertised
			computedSolutions.adapt();
        }

        if ( DSC_CONFIG_GET( "clear_u" , true ) )
            computedSolutions.discreteVelocity().clear();
        if ( DSC_CONFIG_GET( "clear_p" , true ) )
            computedSolutions.discretePressure().clear();
    }
#endif
	last_refine_level = refine_level;

	info.codim0 = gridPtr->size( 0 );
	const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
    infoStream << "  - max grid width: " << grid_width << std::endl;
    info.grid_width = grid_width;

    typedef StokesModelTraitsImp::AnalyticalForceType
        AnalyticalForceType;
	AnalyticalForceType analyticalForce( viscosity, discreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), alpha );

    typedef StokesModelTraitsImp::AnalyticalDirichletDataType
        AnalyticalDirichletDataType;
	AnalyticalDirichletDataType analyticalDirichletData =
			StokesModelTraitsImp::AnalyticalDirichletDataTraitsImplementation::getInstance( discreteStokesFunctionSpaceWrapper );

    typedef Dune::OseenPass< StokesModelImpType >
	    OseenPassType;

	{
		typedef StokesProblems::Container< gridDim, DiscreteOseenFunctionWrapperType>
			ProblemType;
	//	ProblemType problem( viscosity , computedSolutions, analyticalDirichletData );

		typedef PostProcessor< OseenPassType, ProblemType >
			PostProcessorType;

	//	PostProcessorType ( discreteStokesFunctionSpaceWrapper, problem ).save( *gridPtr, computedSolutions, refine_level );
	}

    StokesModelImpType stokesModel( stabil_coeff,
                                    analyticalForce,
                                    analyticalDirichletData,
									viscosity, /*viscosity*/
									alpha, /*alpha*/
									0.0,/*convection_scale_factor*/
									1.0 /*pressure_gradient_scale_factor*/ );

    /* ********************************************************************** *
     * initialize passes                                                      *
     * ********************************************************************** */
    infoStream << "\n- starting pass" << std::endl;


    OseenPassType stokesPass(  stokesModel,
                                gridPart,
								discreteStokesFunctionSpaceWrapper,
								dummyFunctions.discreteVelocity(),
								false );

    PROBLEM_NAMESPACE::SetupCheck check;
    if ( !check( gridPart, stokesPass, stokesModel, computedSolutions ) )
        DUNE_THROW( Dune::InvalidStateException, check.error() );
    DSC_PROFILER.startTiming( "Pass -- APPLY" );
	stokesPass.printInfo();
    auto last_wrapper ( computedSolutions );
    stokesPass.apply( last_wrapper, computedSolutions);
    DSC_PROFILER.stopTiming( "Pass -- APPLY" );
    info.run_time = DSC_PROFILER.getTiming( "Pass -- APPLY" );
    stokesPass.getRuninfo( info );

    /* ********************************************************************** *
     * Problem postprocessing
     * ********************************************************************** */
    infoStream << "\n- postprocesing" << std::endl;
    DSC_PROFILER.startTiming( "Problem/Postprocessing" );

	if ( !DSC_CONFIG_GET( "disableSolver", false ) )
	{
		typedef StokesProblems::Container< gridDim, DiscreteOseenFunctionWrapperType>
			ProblemType;
		ProblemType problem( viscosity , computedSolutions, analyticalDirichletData );

		typedef PostProcessor< OseenPassType, ProblemType >
			PostProcessorType;

		PostProcessorType postProcessor( discreteStokesFunctionSpaceWrapper, problem );

        if ( DSC_CONFIG_GET( "save_solutions", true ) )
            postProcessor.save( *gridPtr, computedSolutions, refine_level );
        else
            postProcessor.calcError( computedSolutions );
		info.L2Errors = postProcessor.getError();
	}
	typedef Dune::StabilizationCoefficients::ValueType
		Pair;
	info.c11 = Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
	info.c12 = Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
	info.d11 = Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
	info.d12 = Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
    info.bfg = DSC_CONFIG_GET( "do-bfg", true );
    //TODO GRIDNAME
//    info.gridname = gridPart.grid().name();
    info.refine_level = refine_level;

    info.polorder_pressure = StokesModelTraitsImp::pressureSpaceOrder;
    info.polorder_sigma = StokesModelTraitsImp::sigmaSpaceOrder;
    info.polorder_velocity = StokesModelTraitsImp::velocitySpaceOrder;

    info.solver_accuracy = DSC_CONFIG_GET( "absLimit", 1e-4 );
    info.inner_solver_accuracy = DSC_CONFIG_GET( "inner_absLimit", 1e-4 );
    info.bfg_tau = DSC_CONFIG_GET( "bfg-tau", 0.1 );

	info.problemIdentifier = PROBLEM_NAMESPACE::identifier;

    DSC_PROFILER.stopTiming( "Problem/Postprocessing" );
    DSC_PROFILER.stopTiming( "SingleRun" );

    firstRun = false;

    return info;
}

void eocCheck( const DSC::RunInfoVector& runInfos )
{
	bool ups = false;
	DSC::RunInfoVector::const_iterator it = runInfos.begin();
	DSC::RunInfo last = *it;
	++it;
	for ( ; it != runInfos.end(); ++it ) {
		ups = ( last.L2Errors[0] < it->L2Errors[0]
			|| last.L2Errors[1] < it->L2Errors[1] );
		last = *it;
	}
	if ( ups ) {
        DSC_LOG_ERROR 	<< 	"----------------------------------------------------------\n"
						<<	"-                                                        -\n"
						<<	"-                  negative EOC                          -\n"
						<<	"-                                                        -\n"
						<< 	"----------------------------------------------------------\n"
						<< std::endl;
	}
}

CoeffVector getAll_Permutations() {
        // the min andmax exponent that will be used in stabilizing terms
    // ie c11 = ( h_^minpow ) .. ( h^maxpow )
    // pow = -2 is interpreted as 0
    // in non-variation context minpow is used for c12 exp and maxpow for d12,
    //      c11 and D11 are set to zero
    int maxpow = DSC_CONFIG_GET( "maxpow", -1 );
    int minpow = DSC_CONFIG_GET( "minpow", 1 );
    double minfactor = DSC_CONFIG_GET( "minfactor", 0.0 );
    double maxfactor = DSC_CONFIG_GET( "maxfactor", 2.0 );
    double incfactor = DSC_CONFIG_GET( "incfactor", 0.5 );

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
    int maxpow = DSC_CONFIG_GET( "maxpow", -1 );
    int minpow = DSC_CONFIG_GET( "minpow", 1 );
    double minfactor = DSC_CONFIG_GET( "minfactor", 0.0 );
    double maxfactor = DSC_CONFIG_GET( "maxfactor", 2.0 );
    double incfactor = DSC_CONFIG_GET( "incfactor", 0.5 );

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
    int maxpow = DSC_CONFIG_GET( "maxpow", -1 );
    int minpow = DSC_CONFIG_GET( "minpow", 1 );
    double minfactor = DSC_CONFIG_GET( "minfactor", 0.0 );

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


/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

