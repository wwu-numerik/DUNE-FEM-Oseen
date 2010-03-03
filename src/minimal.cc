#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/misc.hh>
//do whatever you like to this file to test out simple and small stuff

int main( int argc, char** argv ) {
	if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
	const bool useLogger = true;
    Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
                     Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
                     Parameters().getParam( "fem.io.logdir",    std::string(),              useLogger )
                    );

	std::vector<double> dummies = Parameters().getList( "list", 0.9 );
	Stuff::MinMaxAvg<double> nums ( dummies );
	nums.output( std::cout );
	nums.push( 	4.0 );
	nums.output( std::cout );
	nums.push( 12.0 );
	nums.output( std::cout );
	nums.push( 	4.0 );
	nums.output( std::cout );

	return 0;
}
