/** \file dune_stokes.cc
    \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include "dune/common/mpihelper.hh" // An initializer of MPI
#include "dune/common/exceptions.hh" // We use exceptions

#include "traits.hh"
#include "parameterhandler.hh"
#include "logging.hh"
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
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    //Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    Logger().Create( LogStream::ALL );

    ParameterHandler pm ( "test.param" ) ;
    if ( pm.Ok() ) {
        //pm.Print( std::cout );
        pm.Print( Logger().Min() );
    }
    else {
        return 1;
    }
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
