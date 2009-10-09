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

//replacing gridtype.hh
#ifndef ENABLE_ALUGRID
    #define ENABLE_ALUGRID
#endif
#ifndef HAVE_ALUGRID
    #define HAVE_ALUGRID
#endif
#ifndef ALUGRID_SIMPLEX
    #define ALUGRID_SIMPLEX
#endif
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
const int gridDim = GRIDDIM;
typedef Dune::ALUSimplexGrid< gridDim, gridDim >
    GridType;

const int polOrder = POLORDER;

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/misc.hh>

#if ENABLE_MPI
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
    typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

// forward
void computePermeabilityTensor( const std::string microSolutionsFilenamePrefix,
                                const int microSolutionsRefineLevel );

/**
 *  \brief  main function
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
            return 2;
        }
        if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
            return 1;
        }
        else {
            // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
            //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
            Logger().Create( Dune::Parameter::getValue( "loglevel", 62 ),
                             Dune::Parameter::getValue( "logfile", std::string("darcy") ),
                             Dune::Parameter::getValue( "fem.io.logdir", std::string("log") )
                           );

            computePermeabilityTensor(  Dune::Parameter::getValue( "microSolutionsFilenamePrefix", std::string("") ),
                                        Dune::Parameter::getValue( "microSolutionsRefineLevel", 0 ) );


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

void computePermeabilityTensor( const std::string microSolutionsFilenamePrefix,
                                const int microSolutionsRefineLevel )
{
    Logging::LogStream& info = Logger().Info();
    info.Flush();
    Logging::LogStream& debug = Logger().Dbg();
    debug.Flush();

    // get filenames
    std::string microVelocityDgfFilenameX( microSolutionsFilenamePrefix + "_X_ref_" + Stuff::toString( microSolutionsRefineLevel ) + ".dgf" );
    std::string microVelocityDdfFilenameX( microSolutionsFilenamePrefix + "_X_ref_" + Stuff::toString( microSolutionsRefineLevel ) + ".ddf" );
    std::string microVelocityDgfFilenameY( microSolutionsFilenamePrefix + "_Y_ref_" + Stuff::toString( microSolutionsRefineLevel ) + ".dgf" );
    std::string microVelocityDdfFilenameY( microSolutionsFilenamePrefix + "_Y_ref_" + Stuff::toString( microSolutionsRefineLevel ) + ".ddf" );

    // make sure gridtype and refinement are the same
    const int refineLevel = Stuff::readRefineLevelFromDGF( microVelocityDgfFilenameX );
//    debug << "refineLevel: " << refineLevel << std::endl;
    assert( refineLevel == Stuff::readRefineLevelFromDGF( microVelocityDgfFilenameY ) );

    const std::string gridType( Stuff::readGridTypeFromDGF( microVelocityDgfFilenameX ) );
//    debug << "gridType: " << gridType << std::endl;
    assert( !gridType.compare( Stuff::readGridTypeFromDGF( microVelocityDgfFilenameY ) ) );
    assert( !gridType.compare( "ALUGRID_SIMPLEX" ) );

    //make grid and stuff
    typedef Dune::AdaptiveLeafGridPart< GridType >
        GridPartType;

    typedef Dune::GridPtr< GridType >
        GridPointerType;

    GridPointerType gridPointer( microVelocityDgfFilenameX );
    gridPointer->globalRefine( refineLevel );
    GridPartType gridPart( *gridPointer );

    typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
        FunctionSpaceType;

    typedef FunctionSpaceType::DomainType
        DomainType;

    typedef Dune::DiscontinuousGalerkinSpace<   FunctionSpaceType,
                                                GridPartType,
                                                polOrder >
        DiscreteFunctionSpaceType;

    typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
        DiscreteFunctionType;

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

    DiscreteFunctionType microVelocityX( "microVelocityX", discreteFunctionSpace );
    microVelocityX.read_ascii( microVelocityDdfFilenameX );

    DiscreteFunctionType microVelocityY( "microVelocityY", discreteFunctionSpace );
    microVelocityY.read_ascii( microVelocityDdfFilenameY );

    typedef GridPartType::Codim< 0 >::IteratorType
        EntityIteratorType;

    typedef GridType::Codim< 0 >::Entity
        EntityType;

    typedef EntityType::Geometry
        EntityGeometryType;

    typedef Dune::CachingQuadrature< GridPartType, 0 >
        VolumeQuadratureType;

    typedef DiscreteFunctionType::LocalFunctionType
        LocalFunctionType;

    typedef LocalFunctionType::JacobianRangeType
        JacobianRangeType;

    typedef Dune::FieldMatrix< EntityGeometryType::ctype,
                                        EntityGeometryType::coorddimension,
                                        EntityGeometryType::mydimension >
        JacobianInverseTransposedType;

    double permeabilityTensor_0_0( 0.0 );
    double permeabilityTensor_0_1( 0.0 );
    double permeabilityTensor_1_0( 0.0 );
    double permeabilityTensor_1_1( 0.0 );

    // do gridwalk
    EntityIteratorType entityIteratorEnd = discreteFunctionSpace.end();
    for (   EntityIteratorType entityIterator = discreteFunctionSpace.begin();
            entityIterator != entityIteratorEnd;
            ++entityIterator ) {

        // get entity and geometry
        const EntityType& entity = *entityIterator;
        const EntityGeometryType& geometry = entity.geometry();

        // get local functions
        LocalFunctionType localFunctionX = microVelocityX.localFunction( entity );
        LocalFunctionType localFunctionY = microVelocityY.localFunction( entity );

        // quadrature
        VolumeQuadratureType quadrature( entity, ( 2 * discreteFunctionSpace.order() ) + 1 );

        // do loop over all quadrature points
        for ( int quad = 0; quad < quadrature.nop(); ++quad ) {

            // quadrature point in reference coordinates
            const DomainType x = quadrature.point( quad );
            // in worlsd coordinates
            const DomainType xWorld = geometry.global( x );

            // get the integration factor
            const double integrationFactor = geometry.integrationElement( x );
            // get the quadrature weight
            const double quadratureWeight = quadrature.weight( quad );

            // do integration over unit square only
            if ( !( ( xWorld[0] < 0.0 ) || ( xWorld[0] > 1.0 ) ) ) {
                if ( !( ( xWorld[1] < 0.0 ) || ( xWorld[1] > 1.0 ) ) ) {
                    JacobianRangeType gradient_x_untransposed( 0.0 );
                    localFunctionX.jacobian( x, gradient_x_untransposed );
                    const JacobianInverseTransposedType jacobianInverseTransposed = geometry.jacobianInverseTransposed( x );
                    JacobianRangeType gradient_x( 0.0 );
                    jacobianInverseTransposed.mv( gradient_x_untransposed[0], gradient_x[0] );
                    jacobianInverseTransposed.mv( gradient_x_untransposed[1], gradient_x[1] );
                    JacobianRangeType gradient_y_untransposed( 0.0 );
                    localFunctionY.jacobian( x, gradient_y_untransposed );
                    JacobianRangeType gradient_y( 0.0 );
                    jacobianInverseTransposed.mv( gradient_y_untransposed[0], gradient_y[0] );
                    jacobianInverseTransposed.mv( gradient_y_untransposed[1], gradient_y[1] );
//                    permeabilityTensor_0_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x, gradient_x );
//                    permeabilityTensor_0_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x, gradient_y );
//                    permeabilityTensor_1_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y, gradient_x );
//                    permeabilityTensor_1_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y, gradient_y );
                    permeabilityTensor_0_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x_untransposed, gradient_x_untransposed );
                    permeabilityTensor_0_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x_untransposed, gradient_y_untransposed );
                    permeabilityTensor_1_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y_untransposed, gradient_x_untransposed );
                    permeabilityTensor_1_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y_untransposed, gradient_y_untransposed );
                }
            } // done integration over unit square only
        } // done loop over all quadrature points
    } // done gridwalk

    std::cout << "permeabilitytensor:" << std::endl;
    std::cout << "\t" << permeabilityTensor_0_0 << " \t" << permeabilityTensor_0_1 << std::endl;
    std::cout << "\t" << permeabilityTensor_1_0 << " \t" << permeabilityTensor_1_1 << std::endl;

}

