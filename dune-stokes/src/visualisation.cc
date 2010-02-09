// $Id$

#include <config.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <dune/common/mpihelper.hh> // include mpi helper class
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>

#include"elementdata.hh"

#ifdef GRIDDIM
const int dimGrid = GRIDDIM;
#else
const int dimGrid = 3;
#endif

class FunctorBase {
	public:
		FunctorBase ( const std::string filename )
			: filename_( filename ) {}
		const std::string filename() const { return filename_; }
	protected:
		const std::string filename_;

};

class VolumeFunctor : public FunctorBase {
	public:
		VolumeFunctor ( const std::string filename )
			: FunctorBase( filename ) {}

		template <class Entity>
		double operator() ( const Entity& ent ) const
		{
			return ent.geometry().volume();
		}
};

class BoundaryFunctor : public FunctorBase {
	public:
		BoundaryFunctor ( const std::string filename )
			: FunctorBase( filename ) {}

		template <class Entity>
		double operator() ( const Entity& entity ) const
		{
			double ret( 0.0 );
			int numberOfBoundarySegments( 0 );
			bool isOnBoundary = false;
			typedef typename Entity::LeafIntersectionIterator
			IntersectionIteratorType;
			IntersectionIteratorType endIntersection = entity.ileafend();
			for (	IntersectionIteratorType intersection = entity.ileafbegin();
					intersection != endIntersection;
					++intersection ) {
				if ( !intersection.neighbor() && intersection.boundary() ) {
					isOnBoundary = true;
					numberOfBoundarySegments += 1;
					ret += double( intersection.boundaryId() );
				}
			}
			if ( isOnBoundary ) {
				ret /= double( numberOfBoundarySegments );
			}
			return ret;
		}
};

class AreaMarker : public FunctorBase {

	public:

		AreaMarker( const std::string filename )
			: FunctorBase( filename ) {}

		template <class Entity>
		double operator() ( const Entity& entity ) const
		{
			typedef typename Entity::Geometry
			EntityGeometryType;

			typedef Dune::FieldVector< typename EntityGeometryType::ctype, EntityGeometryType::coorddimension>
			DomainType;

			const EntityGeometryType& geometry = entity.geometry();

			DomainType baryCenter( 0.0 );

			for ( int corner = 0; corner < geometry.corners(); ++corner ) {
				baryCenter += geometry[ corner ];
			}
			baryCenter /= geometry.corners();

			double ret( 0.0 );

			if ( !( ( baryCenter[0] < 0.0 ) || ( baryCenter[0] > 1.0 ) ) ) { // only in unit square
				if ( !( ( baryCenter[1] < 0.0 ) || ( baryCenter[1] > 1.0 ) ) ) {
					ret = 1.0;
				}
			}
			return ret;
		}
};

class GeometryFunctor : public FunctorBase {
	public:
		GeometryFunctor ( const std::string filename )
			: FunctorBase( filename ) {}

		template <class Entity>
		double operator() ( const Entity& ent ) const
		{
			const typename Entity::Geometry& geo = ent.geometry();
			double vol = geo.volume();
			if ( vol < 0 ) {
				std::cout << std::setiosflags( std::ios::fixed ) << std::setprecision( 6 ) << std::setw( 8 );
				//std::cout.showpoint();
				for ( int i = 0; i < geo.corners(); ++i ) {
					std::cout << geo[i] << "\t\t" ;
				}
				std::cout << std::endl;
			}
			return vol;
		}
};

//! supply functor
template<class Grid>
void dowork ( Grid& grid, int refSteps = 0 )
{
	std::string outputDir = Parameters().getParam( "visualisationOutputDir", std::string("visualisation") );
	Stuff::testCreateDirectory( outputDir );
	// make function objects
	BoundaryFunctor boundaryFunctor( outputDir + std::string("/boundaryFunctor") );
	AreaMarker areaMarker( outputDir + std::string("/areaMarker") );
	GeometryFunctor geometryFunctor( outputDir + std::string("/geometryFunctor") );
	VolumeFunctor volumeFunctor( outputDir + std::string("/volumeFunctor") );
	// refine the grid
	grid.globalRefine( refSteps );

	// call the visualization functions
	elementdata( grid, areaMarker );
	elementdata( grid, boundaryFunctor );
	elementdata( grid, volumeFunctor );
	elementdata( grid, geometryFunctor );
}
using namespace Dune;

int main( int argc, char **argv )
{
	// initialize MPI, finalize is done automatically on exit
	MPIHelper& mpihelper = Dune::MPIHelper::instance( argc, argv );

	// start try/catch block to get error messages from dune
	try
	{
		if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
			return 1;
		}

		const bool useLogger = false;
		Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
						 Parameters().getParam( "logfile",          std::string( "dune_stokes" ), useLogger ),
						 Parameters().getParam( "fem.io.logdir",    std::string(),              useLogger )
					   );

		Parameter::append( argc,argv );                         /*@\label{dg:param0}@*/
		if ( argc == 2 ) {
			Parameter::append( argv[1] );
		} else {
			Parameter::append( "parameter" );                     /*@\label{dg:paramfile}@*/
		}                                                       /*@\label{dg:param1}@*/

		GridPtr<GridType> gridptr ( Parameters().DgfFilename( dimGrid ) );
		int refineLevel( 0 );
		Parameter::get( "minref", refineLevel );
		dowork( *gridptr, refineLevel );

	}
	catch ( std::exception & e ) {
		std::cout << "STL ERROR: " << e.what() << std::endl;
		return 1;
	}
	catch ( Dune::Exception & e ) {
		std::cout << "DUNE ERROR: " << e.what() << std::endl;
		return 1;
	}
	catch ( ... ) {
		std::cout << "Unknown ERROR" << std::endl;
		return 1;
	}

	// done
	return 0;
}
