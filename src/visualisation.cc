#include "cmake_config.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/tex.hh>
#include <dune/stuff/tikzgrid.hh>

#include"elementdata.hh"

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

class ProcessIdFunctor : public FunctorBase {
	public:
		ProcessIdFunctor ( const std::string filename, Dune::MPIHelper& mpiHelper )
			: FunctorBase( filename ),
			mpiHelper_( mpiHelper )
		{}

		template <class Entity>
		double operator() ( const Entity& /*ent*/ ) const
		{
			return mpiHelper_.rank();
		}

	protected:
		Dune::MPIHelper& mpiHelper_;
};
#include <dune/grid/alugrid/common/intersectioniteratorwrapper.hh>
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
				if ( !intersection->neighbor() && intersection->boundary() ) {
					isOnBoundary = true;
					numberOfBoundarySegments += 1;
					ret += double( intersection->boundaryId() );
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
				baryCenter += geometry.corner( corner );
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
					std::cout << geo.corner( i ) << "\t\t" ;
				}
				std::cout << std::endl;
			}
			return vol;
		}
};

//! supply functor
template<class Grid>
void dowork ( Grid& grid, int refSteps, Dune::MPIHelper& mpiHelper )
{
	std::string outputDir = Parameters().getParam( "visualisationOutputDir",
					Parameters().getParam("fem.io.datadir", std::string("data") ) + std::string("/visualisation") );
	// make function objects
	BoundaryFunctor boundaryFunctor( outputDir + std::string("/boundaryFunctor") );
	AreaMarker areaMarker( outputDir + std::string("/areaMarker") );
	GeometryFunctor geometryFunctor( outputDir + std::string("/geometryFunctor") );
	ProcessIdFunctor processIdFunctor( outputDir + std::string("/ProcessIdFunctor"), mpiHelper );
//	 refine the grid

//	 call the visualization functions
//	elementdata( grid, areaMarker );
	for ( int i = 0; i < Parameters().getParam( "minref", 0 ); ++i )
	    grid.globalRefine( 1 );
	elementdata( grid, boundaryFunctor );

//	elementdata( grid, geometryFunctor );
//	elementdata( grid, processIdFunctor );
	std::ofstream file( Stuff::getFileinDatadir( "grids.tex" ).c_str() );
//	file << "\\documentclass{article}\n"
//			"\\usepackage{tikz}\n"
//			"\\usetikzlibrary{calc,intersections}\n"
//			"\\pagestyle{empty}\n"
//			"\\begin{document}\n";

//	for ( int i = 0; i < Parameters().getParam( "maxref", 0 ); ++i )
//	{
//		grid.globalRefine( 1 );
//		Stuff::Tex::PgfGrid<Grid> pgfGrid( grid );
//		pgfGrid.output( (boost::format("grid_level%d.tex") % i).str() );
//		file << boost::format("\\begin{tikzpicture}\n\\input{grid_level%d}\n\\end{tikzpicture}\n") % i ;
//	}
//	file << "\\end{document}\n";
    Stuff::Tex::RefineSeriesPgfGrid<Grid> pgfGrid( grid );
    pgfGrid.output( "series.tex", Parameters().getParam( "maxref", 3 ), !Parameters().getParam( "standalone_tex", true ) );
    Stuff::Tex::StackedPgfGrid<Grid> pgfGrid2( grid );
    pgfGrid2.output( "stacked.tex", Parameters().getParam( "maxref", 3 ), !Parameters().getParam( "standalone_tex", true ) );
//	{
//	VolumeFunctor volumeFunctor( outputDir + std::string("/volumeFunctor") );
//	std::cout << Stuff::GridDimensions<Grid>( grid );
//	elementdata( grid, volumeFunctor );
//	}
//	Stuff::EnforceMaximumEntityVolume( grid, 1.2 );
//	{
//	VolumeFunctor volumeFunctor( outputDir + std::string("/volumeFunctorNew") );
//	std::cout << Stuff::GridDimensions<Grid>( grid );
//	elementdata( grid, volumeFunctor );
//	}

}

typedef Dune::GridSelector::GridType
    GridType;

#if ENABLE_MPI
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
    #define DUNE_GET_MPI_COMM Dune::MPIManager::helper().getCommunicator()
#else
    typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
    #define DUNE_GET_MPI_COMM 0.0
#endif

int main( int argc, char **argv )
{
	// initialize MPI, finalize is done automatically on exit
	Dune::MPIManager::initialize(argc, argv);
    CollectiveCommunication mpicomm ( DUNE_GET_MPI_COMM );
	// start try/catch block to get error messages from dune
	try
	{
		if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
			return 1;
		}

		const bool useLogger = false;
		Logger().Create( Parameters().getParam( "loglevel",         62,								useLogger ),
						 Parameters().getParam( "logfile",          std::string( "dune_stokes" ),	useLogger ),
						 Parameters().getParam( "fem.io.datadir",   std::string("data"),			useLogger ),
						 Parameters().getParam( "fem.io.logdir",    std::string(),					useLogger )
					   );

        Dune::GridPtr<GridType> gridptr ( Parameters().DgfFilename( GridType::dimensionworld ) );
		gridptr->loadBalance();
		int refineLevel = Parameters().getParam( "minref", 0 );
		dowork( *gridptr, refineLevel, Dune::MPIManager::helper() );

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
