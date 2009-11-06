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
//typedef Dune::ALUConformGrid< gridDim, gridDim >
//    GridType;

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
#include <dune/stuff/printing.hh>

#if ENABLE_MPI
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
    typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

// forward
void computePermeabilityTensor( const std::string microSolutionsFilenamePrefix,
                                const int microSolutionsRefineLevel );

int writeComsolFile(    const std::string comsolFilename,
                        const double a_0_0,
                        const double a_0_1,
                        const double a_1_0,
                        const double a_1_1 );

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
    debug << "refineLevel: " << refineLevel << std::endl;
    assert( refineLevel == Stuff::readRefineLevelFromDGF( microVelocityDgfFilenameY ) );

    const std::string gridType( Stuff::readGridTypeFromDGF( microVelocityDgfFilenameX ) );
    debug << "gridType: " << gridType << std::endl;
    assert( !gridType.compare( Stuff::readGridTypeFromDGF( microVelocityDgfFilenameY ) ) );
//    assert( !gridType.compare( "ALUGRID_CONFORM" ) );

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

    DomainType velocity_x_min( 0.0 );
    DomainType velocity_x_max( 0.0 );
    DomainType velocity_y_min( 0.0 );
    DomainType velocity_y_max( 0.0 );

    JacobianRangeType velocity_x_gradient_min( 0.0 );
    JacobianRangeType velocity_x_gradient_max( 0.0 );
    JacobianRangeType velocity_y_gradient_min( 0.0 );
    JacobianRangeType velocity_y_gradient_max( 0.0 );

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
                    // do calculation of min and max of velocities
                    DomainType velocity_x( 0.0 );
                    DomainType velocity_y( 0.0 );
                    localFunctionX.evaluate( x, velocity_x );
                    localFunctionY.evaluate( x, velocity_y );
                    const double velocity_x_norm = velocity_x.two_norm();
                    const double velocity_y_norm = velocity_y.two_norm();
                    velocity_x_min[0] = velocity_x[0] < velocity_x_min[0] ? velocity_x[0] : velocity_x_min[0];
                    velocity_x_max[0] = velocity_x[0] > velocity_x_max[0] ? velocity_x[0] : velocity_x_max[0];
                    velocity_x_min[1] = velocity_x[1] < velocity_x_min[1] ? velocity_x[1] : velocity_x_min[1];
                    velocity_x_max[1] = velocity_x[1] > velocity_x_max[1] ? velocity_x[1] : velocity_x_max[1];
                    velocity_y_min[0] = velocity_y[0] < velocity_y_min[0] ? velocity_y[0] : velocity_y_min[0];
                    velocity_y_max[0] = velocity_y[0] > velocity_y_max[0] ? velocity_y[0] : velocity_y_max[0];
                    velocity_y_min[1] = velocity_y[1] < velocity_y_min[1] ? velocity_y[1] : velocity_y_min[1];
                    velocity_y_max[1] = velocity_y[1] > velocity_y_max[1] ? velocity_y[1] : velocity_y_max[1];

                    // do calculation of permeabilitytensor
                    JacobianRangeType gradient_x( 0.0 );
                    localFunctionX.jacobian( x, gradient_x );
                    JacobianRangeType gradient_y( 0.0 );
                    localFunctionY.jacobian( x, gradient_y );
                    permeabilityTensor_0_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x, gradient_x );
                    permeabilityTensor_0_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_x, gradient_y );
                    permeabilityTensor_1_0 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y, gradient_x );
                    permeabilityTensor_1_1 += integrationFactor * quadratureWeight * Stuff::colonProduct( gradient_y, gradient_y );

                    //do calculation of min and max of velocity gradients
                    velocity_x_gradient_min[0][0] = gradient_x[0][0] < velocity_x_gradient_min[0][0] ? gradient_x[0][0] : velocity_x_gradient_min[0][0];
                    velocity_x_gradient_min[0][1] = gradient_x[0][1] < velocity_x_gradient_min[0][1] ? gradient_x[0][1] : velocity_x_gradient_min[0][1];
                    velocity_x_gradient_min[1][0] = gradient_x[1][0] < velocity_x_gradient_min[1][0] ? gradient_x[1][0] : velocity_x_gradient_min[1][0];
                    velocity_x_gradient_min[1][1] = gradient_x[1][1] < velocity_x_gradient_min[1][1] ? gradient_x[1][1] : velocity_x_gradient_min[1][1];
                    velocity_x_gradient_max[0][0] = gradient_x[0][0] > velocity_x_gradient_max[0][0] ? gradient_x[0][0] : velocity_x_gradient_max[0][0];
                    velocity_x_gradient_max[0][1] = gradient_x[0][1] > velocity_x_gradient_max[0][1] ? gradient_x[0][1] : velocity_x_gradient_max[0][1];
                    velocity_x_gradient_max[1][0] = gradient_x[1][0] > velocity_x_gradient_max[1][0] ? gradient_x[1][0] : velocity_x_gradient_max[1][0];
                    velocity_x_gradient_max[1][1] = gradient_x[1][1] > velocity_x_gradient_max[1][1] ? gradient_x[1][1] : velocity_x_gradient_max[1][1];
                    velocity_y_gradient_min[0][0] = gradient_y[0][0] < velocity_y_gradient_min[0][0] ? gradient_y[0][0] : velocity_y_gradient_min[0][0];
                    velocity_y_gradient_min[0][1] = gradient_y[0][1] < velocity_y_gradient_min[0][1] ? gradient_y[0][1] : velocity_y_gradient_min[0][1];
                    velocity_y_gradient_min[1][0] = gradient_y[1][0] < velocity_y_gradient_min[1][0] ? gradient_y[1][0] : velocity_y_gradient_min[1][0];
                    velocity_y_gradient_min[1][1] = gradient_y[1][1] < velocity_y_gradient_min[1][1] ? gradient_y[1][1] : velocity_y_gradient_min[1][1];
                    velocity_y_gradient_max[0][0] = gradient_y[0][0] > velocity_y_gradient_max[0][0] ? gradient_y[0][0] : velocity_y_gradient_max[0][0];
                    velocity_y_gradient_max[0][1] = gradient_y[0][1] > velocity_y_gradient_max[0][1] ? gradient_y[0][1] : velocity_y_gradient_max[0][1];
                    velocity_y_gradient_max[1][0] = gradient_y[1][0] > velocity_y_gradient_max[1][0] ? gradient_y[1][0] : velocity_y_gradient_max[1][0];
                    velocity_y_gradient_max[1][1] = gradient_y[1][1] > velocity_y_gradient_max[1][1] ? gradient_y[1][1] : velocity_y_gradient_max[1][1];


                }
            } // done integration over unit square only
        } // done loop over all quadrature points
    } // done gridwalk


    Stuff::printFieldVector( velocity_x_min, "velocity_x_min", std::cout );
    Stuff::printFieldVector( velocity_x_max, "velocity_x_max", std::cout );
    Stuff::printFieldVector( velocity_y_min, "velocity_y_min", std::cout );
    Stuff::printFieldVector( velocity_y_max, "velocity_y_max", std::cout );
    Stuff::printFieldMatrix( velocity_x_gradient_min, "velocity_x_gradient_min", std::cout );
    Stuff::printFieldMatrix( velocity_x_gradient_max, "velocity_x_gradient_max", std::cout );
    Stuff::printFieldMatrix( velocity_y_gradient_min, "velocity_y_gradient_min", std::cout );
    Stuff::printFieldMatrix( velocity_y_gradient_max, "velocity_y_gradient_max", std::cout );
    std::cout << std::endl;

    std::cout << "permeabilitytensor:" << std::endl;
    std::cout << "\t" << permeabilityTensor_0_0 << " \t" << permeabilityTensor_0_1 << std::endl;
    std::cout << "\t" << permeabilityTensor_1_0 << " \t" << permeabilityTensor_1_1 << std::endl;

    std::string comsolFilename( Dune::Parameter::getValue( "comsolFilename", std::string("standardComsolFilename.m") ) );
    if ( !writeComsolFile(  comsolFilename,
                            permeabilityTensor_0_0,
                            permeabilityTensor_0_1,
                            permeabilityTensor_1_0,
                            permeabilityTensor_1_1 ) ) {
        std::cout << "comsol file written to " << comsolFilename << std::endl;
    }
    else {
        std::cerr << "Error: could not write to comsol file " << comsolFilename << std::endl;
    }

}

int writeComsolFile(   const std::string comsolFilename,
                        const double a_0_0,
                        const double a_0_1,
                        const double a_1_0,
                        const double a_1_1 )
{
    double rightHandSide = a_0_1 - a_1_0;
    std::string g_n_1( "'-1.0*" + Stuff::toString( a_0_0 ) + "*y + 1.0*" + Stuff::toString( a_0_1 ) + "*x'" );
    std::string g_n_2( "'-1.0*" + Stuff::toString( a_0_1 ) + "*y + 1.0* " + Stuff::toString( a_1_1 ) + "*x'" );
    std::string g_n_3( "'1.0*" + Stuff::toString( a_1_0 ) + "*y - 1.0* " + Stuff::toString( a_1_1 ) + "*x'" );
    std::string g_n_4( "'1.0*" + Stuff::toString( a_0_0 ) + "*y - 1.0* " + Stuff::toString( a_0_1 ) + "*x'" );
    std::string one_over_mu( Stuff::toString( 1.0 / Dune::Parameter::getValue( "micro_viscosity", 1.0 ) ) );

    int errorState( 0 );

    std::ofstream writeFile( comsolFilename.c_str() );
    if ( writeFile.is_open() ) {
        writeFile   << "% COMSOL Multiphysics Model M-file" << std::endl
                    << "% Generated by darcy_permeabilitytensor.cc" << std::endl
                    << "% by Felix Albrecht (felix.albrecht@math.uni-muenster.de)" << std::endl
                    << "%" << std::endl
                    << "%permeabilitytensor: " << a_0_0 << " " << a_0_1 << std::endl
                    << "%                    " << a_1_0 << " " << a_1_1 << std::endl
                    << std::endl
                    << "flclear fem" << std::endl
                    << std::endl
                    << "% COMSOL version" << std::endl
                    << "clear vrsn" << std::endl
                    << "vrsn.name = 'COMSOL 3.4';" << std::endl
                    << "vrsn.ext = '';" << std::endl
                    << "vrsn.major = 0;" << std::endl
                    << "vrsn.build = 248;" << std::endl
                    << "vrsn.rcs = '$Name:  $';" << std::endl
                    << "vrsn.date = '$Date: 2007/10/10 16:07:51 $';" << std::endl
                    << "fem.version = vrsn;" << std::endl
                    << std::endl
                    << "% Geometry" << std::endl
                    << "g1=rect2(2,2,'base','corner','pos',[-1,-1]);" << std::endl
                    << std::endl
                    << "% Analyzed geometry" << std::endl
                    << "clear s" << std::endl
                    << "s.objs={g1};" << std::endl
                    << "s.name={'R1'};" << std::endl
                    << "s.tags={'g1'};" << std::endl
                    << std::endl
                    << "fem.draw=struct('s',s);" << std::endl
                    << "fem.geom=geomcsg(fem);" << std::endl
                    << std::endl
                    << "% Initialize mesh" << std::endl
                    << "fem.mesh=meshinit(fem, ..." << std::endl
                    << "                  'hauto',5);" << std::endl
                    << std::endl
                    << "% Refine mesh" << std::endl
                    << "fem.mesh=meshrefine(fem, ..." << std::endl
                    << "                    'mcase',0, ..." << std::endl
                    << "                    'rmethod','regular');" << std::endl
                    << "fem.mesh=meshrefine(fem, ..." << std::endl
                    << "                    'mcase',0, ..." << std::endl
                    << "                    'rmethod','regular');" << std::endl
                    << std::endl
                    << "% (Default values are not included)" << std::endl
                    << std::endl
                    << "% Application mode 1" << std::endl
                    << "clear appl" << std::endl
                    << "appl.mode.class = 'Poisson';" << std::endl
                    << "appl.dim = {'p'};" << std::endl
                    << "appl.assignsuffix = '_poeq';" << std::endl
                    << "clear prop" << std::endl
                    << "prop.elemdefault='Lag1';" << std::endl
                    << "appl.prop = prop;" << std::endl
                    << "clear bnd" << std::endl
                    << "bnd.type = 'neu';" << std::endl
                    << "bnd.g = {" << g_n_4 << "," << g_n_1 << "," << g_n_2 << "," << g_n_3 << "};" << std::endl
                    << "bnd.ind = [2,3,4,1];" << std::endl
                    << "appl.bnd = bnd;" << std::endl
                    << "clear equ" << std::endl
                    << "equ.c = {{{" << a_0_0 << "," << a_0_1 << ";" << a_1_0 << "," << a_1_1 << "}}};" << std::endl
                    << "equ.f = " << rightHandSide << ";" << std::endl
                    << "equ.ind = [1];" << std::endl
                    << "appl.equ = equ;" << std::endl
                    << "fem.appl{1} = appl;" << std::endl
                    << "fem.frame = {'ref'};" << std::endl
                    << "fem.border = 1;" << std::endl
                    << "clear units;" << std::endl
                    << "units.basesystem = 'SI';" << std::endl
                    << "fem.units = units;" << std::endl
                    << std::endl
                    << "% ODE Settings" << std::endl
                    << "clear ode" << std::endl
                    << "clear units;" << std::endl
                    << "units.basesystem = 'SI';" << std::endl
                    << "ode.units = units;" << std::endl
                    << "fem.ode=ode;" << std::endl
                    << "% Multiphysics" << std::endl
                    << "fem=multiphysics(fem);" << std::endl
                    << std::endl
                    << "% Extend mesh" << std::endl
                    << "fem.xmesh=meshextend(fem);" << std::endl
                    << std::endl
                    << "% Solve problem" << std::endl
                    << "fem.sol=femstatic(fem, ..." << std::endl
                    << "                  'solcomp',{'p'}, ..." << std::endl
                    << "                  'outcomp',{'p'});" << std::endl
                    << std::endl
                    << "% eigener kram" << std::endl
                    << "[x,y] = meshgrid(-1.0:0.01:1.0,-1.0:0.01:1.0);" << std::endl
                    << "points = [x(:)';y(:)'];" << std::endl
                    << "pressure = postinterp(fem,'p',points);" << std::endl
                    << "pressure = reshape(pressure,size(x));" << std::endl
                    << "pressure_gradient_x = postinterp(fem,'px',points);" << std::endl
                    << "pressure_gradient_x = reshape(pressure_gradient_x,size(x));" << std::endl
                    << "pressure_gradient_y = postinterp(fem,'py',points);" << std::endl
                    << "pressure_gradient_y = reshape(pressure_gradient_y,size(x));" << std::endl
                    << std::endl
                    << "velocity_x = " << one_over_mu << ".*" << a_0_0 << ".*(y - 1.0.*pressure_gradient_x) + 1.0.*" << a_0_1 << ".*(-1.0.*x - 1.0.*pressure_gradient_y);" << std::endl
                    << "velocity_y = " << one_over_mu << ".*" << a_1_0 << ".*(y - 1.0.*pressure_gradient_x) + 1.0.*" << a_1_1 << ".*(-1.0.*x - 1.0.*pressure_gradient_y);" << std::endl
                    << std::endl
                    << "save( 'saved_x_" << comsolFilename << "', 'x' );" << std::endl
                    << "save( 'saved_y_" << comsolFilename << "', 'y' );" << std::endl
                    << "save( 'saved_pressure_" << comsolFilename << "', 'pressure' );" << std::endl
                    << std::endl;

        writeFile.flush();
        writeFile.close();

        return 0;
    }
    else {
        ++errorState;
        std::cerr << "Error: could not write to " << comsolFilename << ".m!" << std::endl;
        return errorState;
    }

    return errorState;
}































