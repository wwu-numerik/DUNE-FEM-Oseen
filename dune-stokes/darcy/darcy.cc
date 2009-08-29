/**
 *  \file   darcy.cc
 *
 *  \todo   doc
 **/

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#if defined(UGGRID) && defined(DEBUG)
    #warning ("UGGRID in debug mode is likely to produce a segfault")
#endif

#if ! defined(MACRO_POLORDER)
    #define MACRO_POLORDER 1
    #warning ("MACRO_POLORDER undefined, defaulting to 1")
#endif

#if ! defined(MICRO_POLORDER)
    #define MICRO_POLORDER 1
    #warning ("MICRO_POLORDER undefined, defaulting to 1")
#endif

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
//#if HAVE_GRAPE
//#include <dune/grid/io/visual/grapedatadisplay.hh>
//#endif

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/operator/feop.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/elementintegrators.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/stuff/logging.hh>
//#include <dune/stuff/printing.hh>
//#include <dune/stuff/misc.hh>
//#include <dune/stuff/parametercontainer.hh>
//#include <dune/stuff/postprocessing.hh>
//#include <dune/stuff/profiler.hh>

//#include <dune/stokes/discretestokesfunctionspacewrapper.hh>
//#include <dune/stokes/discretestokesmodelinterface.hh>
//#include <dune/stokes/stokespass.hh>
//#include "analyticaldata.hh"
//#include "velocity.hh"
//#include "pressure.hh"
//#include "problem.hh"

#include "elementintegrators.hh"
#include "ellipticmodels.hh"

#if ENABLE_MPI
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
    typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

// forward
void singleRun( const int refineLevel );

/**
 *  \brief  main function
 *
 *          ParameterContainer and Logger setup
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main( int argc, char** argv )
{
    try {

        Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc, argv);
        CollectiveCommunication mpicomm ( mpihelper.getCommunicator() );

        if ( argc < 2 ) {
            std::cerr << "\nUsage: " << argv[0] << " parameterfile" << std::endl;
            return 1;
        }
        else {
            Dune::Parameter::append( argv[1] );
            Dune::Parameter::append( argc, argv );

            // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
            //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
            Logger().Create( Dune::Parameter::getValue( "loglevel", 62 ),
                             Dune::Parameter::getValue( "logfile", std::string("darcy") ),
                             Dune::Parameter::getValue( "fem.io.logdir", std::string() )
                           );

            const int refineLevel = 1;

            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();

            infoStream << "This is " << argv[0] << "." << std::endl;
            debugStream << "Doing singleRun() with refinelevel " << refineLevel << "." << std::endl;

            singleRun( refineLevel );

            return 0;
        }
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}

/**
 *  \brief  single run
 *
 *  \todo   doc
 **/
void singleRun( const int refineLevel )
{
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();

    debugStream << "This is singleRun( " << refineLevel << " )." << std::endl;

    debugStream << "Initialising macro model..." << std::endl;

    /*
     * macro grid
     */

    const int gridDim = GridType::dimensionworld;

    std::string macroGridFile( "macro_grid_2d.dgf" );
    if ( gridDim == 2 ) {
        macroGridFile = Dune::Parameter::getValue( "macro_grid_2d", std::string("macro_grid_2d.dgf") );
    }
    else if ( gridDim == 3 ) {
        assert( !"Darcy only implemented in 2D!" );
        macroGridFile = Dune::Parameter::getValue( "macro_grid_3d", std::string("macro_grid_3d.dgf") );
    }
    else {
        assert( !"Darcy only implemented in 2D and 3D!" );
    }
    Dune::GridPtr< GridType > macroGridPointer( macroGridFile );

    macroGridPointer->globalRefine( refineLevel * Dune::DGFGridInfo< GridType >::refineStepsForHalf() );

    typedef Dune::AdaptiveLeafGridPart< GridType >
        MacroGridPartType;

    MacroGridPartType macroGridPart( *macroGridPointer );

    infoStream << "Initialised the grid (with " << macroGridPart.grid().size( 0 ) << " elements)." << std::endl;

    /*
     * macro function space etc
     */

    const int macroPolOrder = MACRO_POLORDER;

    typedef Dune::FunctionSpace< double, double, gridDim, 1 >
        MacroFunctionSpaceType;

    typedef Dune::LagrangeDiscreteFunctionSpace< MacroFunctionSpaceType, MacroGridPartType, macroPolOrder >
        MacroDiscreteFunctionSpaceType;

    MacroDiscreteFunctionSpaceType macroDiscreteFunctionSpace( macroGridPart );

    typedef Dune::AdaptiveDiscreteFunction< MacroDiscreteFunctionSpaceType >
        MacroDiscreteFunctionType;

    MacroDiscreteFunctionType macroRightHandSide( "macro_rhs", macroDiscreteFunctionSpace );
    macroRightHandSide.clear();

    MacroDiscreteFunctionType macroPressure( "macro_pressure", macroDiscreteFunctionSpace );
    macroPressure.clear();

    /*
     * macro model
     */

    typedef Darcy::DarcyModel< MacroFunctionSpaceType >
        MacroModelType;

    MacroModelType macroModel;

    infoStream << "Initialised macro model." << std::endl;

    debugStream << "Initialising macro elliptic operator..." << std::endl;

    /*
     * macro integrators and system matrix
     */

    typedef Dune::DefaultMatrixElementIntegratorTraits< MacroDiscreteFunctionType, 100 >
        MacroElementIntegratorTraitsType;

    typedef Dune::SimpleElementMatrixIntegrator< Dune::DefaultMatrixElementIntegratorTraits< MacroDiscreteFunctionType, 100 >, MacroModelType >
        MacroElementMatrixIntegratorType;

    MacroElementMatrixIntegratorType macroElementMatrixIntegrator( macroModel, macroDiscreteFunctionSpace, 0 );

    typedef Darcy::DarcyElementRhsIntegrator< MacroElementIntegratorTraitsType, MacroModelType, MacroDiscreteFunctionSpaceType >
        MacroElementRhsIntegratorType;

    MacroElementRhsIntegratorType macroElementRhsIntegrator( macroModel, macroDiscreteFunctionSpace, 0 );

    typedef Dune::RhsAssembler< MacroElementRhsIntegratorType >
        MacroRhsAssemblerType;

    MacroRhsAssemblerType macroRhsAssembler( macroElementRhsIntegrator, 0 );

    typedef Dune::SparseRowMatrix< double >
        MacroSystemMatrixType;

    /*
     * macro elliptic operator
     */

    typedef Dune::FEOp< MacroSystemMatrixType, MacroElementMatrixIntegratorType >
        MacroEllipticOperatorType;

    MacroEllipticOperatorType macroEllipticOperator(    macroElementMatrixIntegrator,
                                                        MacroEllipticOperatorType::ASSEMBLED,
                                                        200,
                                                        0 );

    assert( macroEllipticOperator.systemMatrix().checkConsistency() );

    infoStream << "Initialised macro elliptic operator." << std::endl;

    /*
     * assembling
     */

    debugStream << "Assembling macro system matrix and right hand side..." << std::endl;

    macroEllipticOperator.assemble();

    macroRhsAssembler.assemble( macroRightHandSide );

    infoStream << "Assembled macro system matrix and right hand side." << std::endl;

    /*
     * solver
     */

    debugStream << "Solving macro system for the pressure (" << macroDiscreteFunctionSpace.size() <<" unknowns)..." << std::endl;

    typedef Dune::MACRO_CG_SOLVERTYPE< MacroDiscreteFunctionType, MacroEllipticOperatorType >
        MacroInverseOperatorType;

    MacroInverseOperatorType macroInverseOperator(  macroEllipticOperator,
                                                    12345.67890,
                                                    Dune::Parameter::getValue( "macro_solver_accuracy", 1e-10 ),
                                                    20000,
                                                    Dune::Parameter::getValue( "macro_solver_verbosity", 0 ) );

    macroInverseOperator( macroRightHandSide, macroPressure );

    infoStream << "Solved macro system for the pressure." << std::endl;

    /*
     * output
     */

    debugStream << "Writing output..." << std::endl;

    typedef Dune::VTKIO< MacroGridPartType >
        MacroVTKWriterType;

    MacroVTKWriterType macroVtkWriter( macroGridPart );

    macroVtkWriter.addVertexData( macroPressure );
    macroVtkWriter.write( "data/macroPressure" );
    macroVtkWriter.clear();

    macroVtkWriter.addVertexData( macroRightHandSide );
    macroVtkWriter.write( "data/macroRightHandSide" );
    macroVtkWriter.clear();

    infoStream << "Output written. Have a nice Day." << std::endl;
}
