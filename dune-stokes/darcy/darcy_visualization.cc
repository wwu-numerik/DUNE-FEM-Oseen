/**
 *  \file   darcy_visualization.cc
 *
 *  \todo   doc
 **/

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#ifndef POLORDER
    #define POLORDER 1
    #warning ("POLORDER undefined, defaulting to 1")
#endif
const int polOrder = POLORDER;

#ifndef GRIDTYPE
    #define GRIDTYPE ALUGRID_SIMPLEX
    #warning ("GRIDTYPE undefined, defaulting to ALUGRID_SIMPLEX")
#endif
#ifndef GRIDDIM
    #define GRIDDIM 2
    #warning ("GRIDDIM undefined, defaulting to 2")
#endif
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
const int gridDim = GRIDDIM;

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/functions.hh>

#if ENABLE_MPI
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
    typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif


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

            typedef Dune::AdaptiveLeafGridPart< GridType >
                GridPartType;

            typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
                FunctionSpaceType;

            typedef Dune::DiscontinuousGalerkinSpace<   FunctionSpaceType,
                                                        GridPartType,
                                                        polOrder >
                DiscreteFunctionSpaceType;

            typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
                DiscreteFunctionType;

            Stuff::loadDiscreteFunction< DiscreteFunctionType >( Dune::Parameter::getValue( "loadDiscreteFunctionFilenamePrefix", std::string("") ), Logger().Err() );

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

