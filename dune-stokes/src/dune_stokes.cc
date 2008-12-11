/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <memory> //not including this gives error of undefined autopointer in dgfparser.hh
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction

#include <dune/stokes/discretestokesmodelinterface.hh>
//#include <dune/stokes/saddlepoint_inverse_operator.hh>
#include <dune/stokes/stokespass.hh>

#include "parametercontainer.hh"
#include "logging.hh"
//#include "problem.hh"
//#include "postprocessing.hh"
#include "analyticaldata.hh"

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
    //Maybe initialize Mpi
    //Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */

    ParameterContainer& parameters = Parameters();
    if ( !( parameters.ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
    if ( !( parameters.SetUp() ) ) {
        return 1;
    }
    else {
        parameters.SetGridDimension( GridType::dimensionworld );
        parameters.SetPolOrder( POLORDER );
        parameters.Print( std::cout );
    }

    const int gridDim = GridType::dimensionworld;
    const int polOrder = POLORDER;

    Logger().Create(
        Logging::LOG_CONSOLE |
        Logging::LOG_FILE |
        Logging::LOG_ERR |
        Logging::LOG_DEBUG |
        Logging::LOG_INFO );

    Logging::LogStream& infoStream = Logger().Info();
    //Logging::LogStream& debugStream = Logger().Dbg();
    //Logging::LogStream& errorStream = Logger().Err();

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
    infoStream << "\ninitialising grid..." << std::endl;

    typedef Dune::LeafGridPart< GridType >
        GridPartType;
    Dune::GridPtr< GridType > gridPtr( parameters.DgfFilename() );
    GridPartType gridPart( *gridPtr );

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
    AnalyticalForceType analyticalForce( 0.5, velocitySpace );

    typedef DirichletData< VelocityFunctionSpaceType >
        AnalyticalDirichletDataType;
    AnalyticalDirichletDataType analyticalDirichletData( velocitySpace );

    infoStream << "...done." << std::endl;
    /* ********************************************************************** *
     * initialize model                                                       *
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

    Dune::FieldVector< double, gridDim > ones( 1.0 );
    StokesModelImpType stokesModel( 1.0,
                                    ones,
                                    1.0,
                                    ones,
                                    analyticalForce,
                                    analyticalDirichletData );

    typedef Dune::DiscreteStokesModelInterface<
                Dune::DiscreteStokesModelDefaultTraits<
                    GridPartType,
                    AnalyticalForceType,
                    AnalyticalDirichletDataType,
                    gridDim,
                    polOrder > >
        StokesModelType;

    StokesModelType::DiscreteVelocityFunctionSpaceType discreteVelocitySpace( gridPart );
    StokesModelType::DiscreteSigmaFunctionSpaceType discreteSigmaSpace( gridPart );
    StokesModelType::DiscretePressureFunctionSpaceType discretePressureSpace( gridPart );

    typedef StokesModelType::Traits::DiscreteFunctionSpacePair
        DiscreteFunctionSpacePair;

    typedef Dune::Pair < GridPartType, GridPartType >
        GridPartTypePair;

    GridPartTypePair gridPartTypePair( gridPart, gridPart );

    DiscreteFunctionSpacePair discreteFunctionSpacePair( gridPartTypePair );

    typedef Dune::DiscreteFunctionPair< Dune::DiscreteFunctionPairTraits <
                StokesModelType::DiscreteVelocityFunctionType,
                StokesModelType::DiscretePressureFunctionType,
                DiscreteFunctionSpacePair > >
        DiscreteFunctionPairType;

    Dune::Pair< const std::string, const std::string > namePair( "velocity", "pressure" );

    DiscreteFunctionPairType discreteFunctionPair( namePair, discreteFunctionSpacePair );


    infoStream << "...done." << std::endl;
    /* ********************************************************************** *
     * initialize passes                                                      *
     * ********************************************************************** */
    infoStream << "\ninitialising passes..." << std::endl;

    typedef Dune::StartPass< DiscreteFunctionPairType, -1 >
        StartPassType;
    StartPassType startPass;

    typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
        StokesPassType;
    StokesPassType stokesPass(  startPass,
                                discreteVelocitySpace,
                                discreteSigmaSpace,
                                discretePressureSpace,
                                stokesModel,
                                gridPart );

    StokesPassType::DomainType dummyDomain( namePair, discreteFunctionSpacePair );
    StokesPassType::RangeType dummyRange( namePair, discreteFunctionSpacePair );
    stokesPass.apply( dummyDomain, dummyRange );

    infoStream << "...done." << std::endl;

    return 0;
  }
  catch ( Dune::Exception &e ){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch ( ... ){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
