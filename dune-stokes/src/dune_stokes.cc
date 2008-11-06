/** \file dune_stokes.cc
    \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
//#define GRIDTYPE ALUGRID_SIMPLEX
#include <iostream>
#include <memory> //not including this gives error of undefined autopointer in dgfparser.hh
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid
//#include <dune/grid/io/file/dgfparser/dgfparser.hh> // for the grid
//#include <dune/grid/utility/gridtype.hh>

#include "traits.hh"
#include "parametercontainer.hh"
#include "parameterhandler.hh"
#include "logging.hh"

#include "saddlepoint_inverse_operator.hh"
#include <stokes/stokespass.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/gridpart.hh>

/**
 *  \brief main function
 *
 *  \c more main function
 *
 *  \attention attention
 *
 *  \param argc number of arguments from command line
 *  \param argv array of arguments from command line
 **/
int main( int argc, char** argv )
{
  try{
    //Maybe initialize Mpi
    //Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    ParameterContainer parameters( argc, argv );
    if ( !( parameters.ReadCommandLine() ) ) {
        return 1;
    }
    if ( !( parameters.SetUp() ) ) {
        return 1;
    }
    else {
        parameters.Print( std::cout );
    }
    Logger().Create( Logging::LOG_CONSOLE | Logging::LOG_FILE | Logging::LOG_ERR | Logging::LOG_INFO );

    /* ********************************************************************** *
     * initialize the grid                                                    *
     * ********************************************************************** */
     using Dune::GridPtr;
    GridPtr<GridType> gridptr( "grid.dgf" );

    Logging::LogStream& myStream = Logger().Err();
    myStream << "hgude" << " pudge";
    myStream << std::endl ;
    myStream << std::setw(12) << std::setprecision(8) << 6.786968789659698697 ;
    myStream << std::endl ;
    Logging::LogStream& myStream2 = Logger().Dbg();
    myStream2 << "\ndebugout" << std::endl;
    Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_DEBUG | Logging::LOG_FILE );
    myStream2 << "\ndebugout22" << std::endl;

    int newStreamID = Logger().AddStream( Logging::LOG_CONSOLE );
    Logging::LogStream& blah = Logger().GetStream( newStreamID );
    blah << "blah" << std::endl;

    const int polOrd = 2;
    const int dim = 2;


    typedef Dune::FunctionSpace<double, double, dim, dim>
        FunctionSpaceType;
    typedef FunctionSpaceType::DomainFieldType
        DomainFieldType;
    typedef FunctionSpaceType::RangeFieldType
        RangeFieldType;
    typedef Dune::FunctionSpace<DomainFieldType,RangeFieldType,dim,dim>
        PressureSpaceType;
    typedef Dune::LeafGridPart<GridType>
        GridPartType;
    typedef Dune::DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd >
        DiscreteFunctionSpaceType;
    typedef Dune::AdaptiveDiscreteFunction<DiscreteFunctionSpaceType>
        DiscreteFunctionType;
    typedef Dune::StokesPass< DiscreteFunctionType, DiscreteFunctionType >
        PassType;
    typedef Dune::OEMCGOp <DiscreteFunctionType, PassType >
        InverseOperatorType;
    typedef Dune::SaddlepointInverseOperator< PassType, InverseOperatorType >
        SaddlepointInverseOperatortype;

    PassType pass;
    InverseOperatorType aufSolver(pass,1e-10,1e-10,5000,false);
    GridPartType gridpart ( *gridptr );
    DiscreteFunctionSpaceType pressurespc ( gridpart  );
    DiscreteFunctionSpaceType velospace ( gridpart  );
    SaddlepointInverseOperatortype invOp(pass,1e-8,1e-8,5000,1,aufSolver,pressurespc,velospace);


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
