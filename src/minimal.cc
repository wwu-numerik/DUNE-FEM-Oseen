#include "cmake_config.h"

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/math.hh>

//do whatever you like to this file to test out simple and small stuff
#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>

#include <dune/istl/io.hh>

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

    typedef Dune::FieldMatrix<double,2,2> M;
    Dune::BCRSMatrix<M> B(4,4,Dune::BCRSMatrix<M>::random);

    // initially set row size for each row
    B.setrowsize(0,1);
    B.setrowsize(3,4);
    B.setrowsize(2,1);
    B.setrowsize(1,1);
    // increase row size for row 2
    B.incrementrowsize(2);

    // finalize row setup phase
    B.endrowsizes();

    // add column entries to rows
    B.addindex(0,0);
    B.addindex(3,1);
    B.addindex(2,2);
    B.addindex(1,1);
    B.addindex(2,0);
    B.addindex(3,2);
    B.addindex(3,0);
    B.addindex(3,3);

    // finalize column setup phase
    B.endindices();

    // set entries using the random access operator
    B[0][0] = 1;
    B[1][1] = 2;
    B[2][0] = 3;
    B[2][2] = 4;
    B[3][1] = 5;
    B[3][2] = 6;
    B[3][0] = 7;
    B[3][3] = 8;

     Dune::printSparseMatrix( std:: cout, B, "title", "row");

	return 0;
}
