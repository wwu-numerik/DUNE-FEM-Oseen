#include "cmake_config.h"

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/progressbar.hh>
//do whatever you like to this file to test out simple and small stuff

int main( int argc, char** argv ) {
	Dune::MPIManager::initialize(argc, argv);
	if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
	const bool useLogger = true;
    Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
                     Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
					 Parameters().getParam( "fem.io.datadir",   std::string("data"),        useLogger ),
					 Parameters().getParam( "fem.io.logdir",    std::string("logs"),        useLogger )
                    );

	std::vector<double> dummies = Parameters().getList( "list", 0.9 );
	Stuff::MinMaxAvg<double> nums ( dummies );
	nums.output( std::cout );
	nums( 	4.0 );
	nums.output( std::cout );
	nums( 12.0 );
	nums.output( std::cout );
	nums( 	4.0 );
	nums.output( std::cout );

	unsigned k = 10;
	Stuff::SimpleProgressBar<> pbar(k);
	 for(unsigned i=0; i<k+1; i++,++pbar) {
		 sleep( 1 );
	  }
	 	++pbar;

	return 0;
}
