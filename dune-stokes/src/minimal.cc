#include <dune/stuff/parametercontainer.hh>

//do whatever you like to this file to test out simple and small stuff

int main( int argc, char** argv ) {
	if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
        return 1;
    }
	std::vector<double> dummies = Parameters().getList( "list", 0.9 );
	for ( std::vector<double>::const_iterator it = dummies.begin(); it != dummies.end(); ++it ) {
		std::cout << *it << "\n";
	}
	return 0;
}
