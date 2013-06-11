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
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/oseen/functionspacewrapper.hh>
#include <dune/fem/oseen/modeldefault.hh>

#include <dune/fem/oseen/ldg_method.hh>
#include <dune/fem/oseen/boundarydata.hh>
#include <dune/fem/oseen/runinfo.hh>
#include <dune/fem/oseen/tex.hh>
#include <dune/fem/oseen/postprocessing.hh>
#include <dune/fem/oseen/problems.hh>

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/signals.hh>
#include <dune/stuff/fem/femeoc.hh>

#include <dune/grid/utility/gridtype.hh>

#ifndef COMMIT
    #define COMMIT "undefined"
#endif

typedef Dune::GridSelector::GridType
    GridType;

static const std::string commit_string (COMMIT);

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

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

//! multiple runs with refine level in [minref*refineStepsForHalf : maxref*refineStepsForHalf] and everything else on default
void RefineRun  ( CollectiveCommunication& mpicomm );

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
    CollectiveCommunication mpicomm;// = Dune::MPIManager::helper().getCommunicator();

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    if ( argc < 2 ) {
        std::cerr << "\nUsage: " << argv[0] << " parameterfile \n" << "\n\t --- OR --- \n";
        std::cerr << "\nUsage: " << argv[0] << " paramfile:"<<"file" << " more-opts:val ..." << std::endl << std::endl;
        return 2;
    }

    DSC_CONFIG.readCommandLine(argc, argv);

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    DSC_LOG.create( DSC_CONFIG_GETB( "loglevel",         62,                          useLogger ),
                     DSC_CONFIG_GETB( "logfile",          std::string("dune_stokes"), useLogger ),
                     DSC_CONFIG_GETB( "fem.io.datadir",   std::string(),              useLogger )
                    );

    if (DSC_CONFIG_GET( "runtype", 5 ) == 5)
    {
        const auto minref = DSC_CONFIG_GET( "minref", 0 );
        DSC_CONFIG.set("maxref", minref);
    }
    RefineRun( mpicomm );
    DSC_LOG_ERROR << "\nRun from: " << commit_string << std::endl;
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
  return 0;
}

void RefineRun( CollectiveCommunication& mpicomm )
{
    DSC_LOG_INFO << "starting refine run " << std::endl;
    // column headers for eoc table output
    const ColumnHeaders errorColumnHeaders = { "h", "el't","Laufzeit (s)","Geschwindigkeit", "Druck" };
	DSC::RunInfoVector run_infos;
    auto& eoc_output = DSFe::FemEoc::instance( );
	eoc_output.initialize( DSC_CONFIG_GET("fem.io.datadir", std::string("data") ),"eoc-file", "eoc-desc", "eoc-template.tex" );
    const size_t idx = eoc_output.addEntry( errorColumnHeaders );
    DSC::RefineOutput eoc_texwriter( errorColumnHeaders );

    const int minref = DSC_CONFIG_GET( "minref", 0 );
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
    if ( firstRun && refine_level_factor > 0 ) {
        //since we have a couple of local statics, only do this once, further refinement done in estimator
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
    if ( !firstRun ) {
        for ( int i = refine_level - last_refine_level; i > 0; --i )
        {
            //simpler would be to use real weights in mark(), but alas, that doesn't work as advertised
            for(const auto& entity :computedSolutions.discretePressure().space())
                gridPart.grid().mark(1, entity);
            computedSolutions.adapt();
        }
        if ( DSC_CONFIG_GET( "clear_u" , true ) )
            computedSolutions.discreteVelocity().clear();
        if ( DSC_CONFIG_GET( "clear_p" , true ) )
            computedSolutions.discretePressure().clear();
    }
    gridPart.grid().loadBalance();

	last_refine_level = refine_level;

	info.codim0 = gridPtr->size( 0 );
	const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
    infoStream << "  - max grid width: " << grid_width << std::endl;
    info.grid_width = grid_width;

    typedef StokesModelTraitsImp::AnalyticalForceType
        AnalyticalForceType;
    AnalyticalForceType analyticalForce( viscosity, alpha );

    typedef StokesModelTraitsImp::AnalyticalDirichletDataType
        AnalyticalDirichletDataType;
	AnalyticalDirichletDataType analyticalDirichletData =
			StokesModelTraitsImp::AnalyticalDirichletDataTraitsImplementation::getInstance( discreteStokesFunctionSpaceWrapper );

    typedef Dune::OseenLDGMethod< StokesModelImpType >
        OseenLDGMethodType;

	{
		typedef StokesProblems::Container< gridDim, DiscreteOseenFunctionWrapperType>
			ProblemType;
        ProblemType problem( viscosity , analyticalDirichletData );
        typedef PostProcessor< OseenLDGMethodType, ProblemType >
			PostProcessorType;
        PostProcessorType ( discreteStokesFunctionSpaceWrapper, problem ).save( *gridPtr, computedSolutions, refine_level );
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


    OseenLDGMethodType oseenLDG(  stokesModel,
                                gridPart,
								discreteStokesFunctionSpaceWrapper,
								dummyFunctions.discreteVelocity(),
								false );

    PROBLEM_NAMESPACE::SetupCheck check;
    if ( !check( gridPart, oseenLDG, stokesModel, computedSolutions ) )
        DUNE_THROW( Dune::InvalidStateException, check.error() );
    DSC_PROFILER.startTiming( "Pass -- APPLY" );
    auto last_wrapper ( computedSolutions );
    oseenLDG.apply( last_wrapper, computedSolutions);
    DSC_PROFILER.stopTiming( "Pass -- APPLY" );
    info.run_time = DSC_PROFILER.getTiming( "Pass -- APPLY" );
    oseenLDG.getRuninfo( info );

    /* ********************************************************************** *
     * Problem postprocessing
     * ********************************************************************** */
    infoStream << "\n- postprocesing" << std::endl;
    DSC_PROFILER.startTiming( "Problem/Postprocessing" );

	if ( !DSC_CONFIG_GET( "disableSolver", false ) )
	{
		typedef StokesProblems::Container< gridDim, DiscreteOseenFunctionWrapperType>
			ProblemType;
        ProblemType problem( viscosity , analyticalDirichletData );

        typedef PostProcessor< OseenLDGMethodType, ProblemType >
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

