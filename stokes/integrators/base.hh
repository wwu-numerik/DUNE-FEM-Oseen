#ifndef DUNE_STOKES_INTEGRATORS_BASE_HH
#define DUNE_STOKES_INTEGRATORS_BASE_HH

#include <dune/stuff/grid.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/profiler.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class SigmaJacobianRangeType, class VelocityRangeType >
	static SigmaJacobianRangeType
			prepareVelocityRangeTypeForSigmaDivergence( const VelocityRangeType& arg )
	{
		SigmaJacobianRangeType ret( 0.0 );
		assert( arg.dim() == ret[0].dim() );
		for ( int i = 0; i < int(arg.dim()) ; ++i ) {
			for ( int j = 0; j < int(arg.dim()); ++j ) {
				VelocityRangeType row( 0.0 );
				row[ j ] = arg[ i ];
				ret[ i * arg.dim() + j ] = row;
			}
		}
		return ret;
	}

	template < class Traits >
	static typename Traits::VelocityJacobianRangeType preparePressureRangeTypeForVelocityDivergence(
												const typename Traits::PressureRangeType& arg )
	{
		typename Traits::VelocityJacobianRangeType ret( 0.0 );
		for ( unsigned int i = 0; i < ret[0].dim(); ++i ) {
			typename Traits::VelocityRangeType row( 0.0 );
			row[ i ] = arg;
			ret[ i ] = row;
		}
		return ret;
	}

	template < class Traits, class IntegratorTuple >
	class Coordinator
	{
	protected:
		typedef Coordinator< Traits, IntegratorTuple >
			CoordinatorType;

		const typename Traits::DiscreteModelType&					discrete_model_;
		const typename Traits::GridPartType&						grid_part_;
		const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space_;
		const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space_;
		const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space_;

	public:
		typedef typename Traits::ElementCoordinateType
			ElementCoordinateType;
		typedef typename Traits::SigmaRangeType
			SigmaRangeType;
		typedef typename Traits::VelocityRangeType
			VelocityRangeType;
		typedef typename Traits::PressureRangeType
			PressureRangeType;
		typedef typename Traits::VelocityJacobianRangeType
			VelocityJacobianRangeType;
		typedef typename Traits::PressureJacobianRangeType
			PressureJacobianRangeType;
		typedef typename Traits::DiscreteSigmaFunctionSpaceType::BaseFunctionSetType
				SigmaBaseFunctionSetType;
		typedef typename Traits::DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
			VelocityBaseFunctionSetType;
		typedef typename Traits::DiscretePressureFunctionSpaceType::BaseFunctionSetType
			PressureBaseFunctionSetType;

		Coordinator(const typename Traits::DiscreteModelType&					discrete_model,
						const typename Traits::GridPartType&						grid_part,
						const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space,
						const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space,
						const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space)
					:discrete_model_(discrete_model),
					grid_part_(grid_part),
					velocity_space_(velocity_space),
					pressure_space_(pressure_space),
					sigma_space_(sigma_space)
		{}

		//! just to avoid overly long argument lists
		struct InfoContainerVolume {
			const typename Traits::EntityType& entity;
			const typename Traits::EntityType::Geometry& geometry;

			const SigmaBaseFunctionSetType
					sigma_basefunction_set_element;
			const VelocityBaseFunctionSetType
					velocity_basefunction_set_element;
			const PressureBaseFunctionSetType
					pressure_basefunction_set_element;
			const int numSigmaBaseFunctionsElement;
			const int numVelocityBaseFunctionsElement;
			const int numPressureBaseFunctionsElement;
			const typename Traits::VolumeQuadratureType volumeQuadratureElement;
			const typename Traits::DiscreteModelType&	discrete_model;
			const double eps;
			const double viscosity;
			const double convection_scaling;
			const double pressure_gradient_scaling;
			//! generalized stokes alpha
			const double alpha;
			InfoContainerVolume(const CoordinatorType& interface,
								const typename Traits::EntityType& ent,
								const typename Traits::DiscreteModelType& discrete_modelIn )
				: entity( ent ),
				  geometry( entity.geometry() ),
				  sigma_basefunction_set_element( interface.sigma_space_.baseFunctionSet( entity ) ),
				  velocity_basefunction_set_element( interface.velocity_space_.baseFunctionSet( entity ) ),
				  pressure_basefunction_set_element( interface.pressure_space_.baseFunctionSet( entity ) ),
				  numSigmaBaseFunctionsElement( sigma_basefunction_set_element.numBaseFunctions() ),
				  numVelocityBaseFunctionsElement( velocity_basefunction_set_element.numBaseFunctions() ),
				  numPressureBaseFunctionsElement( pressure_basefunction_set_element.numBaseFunctions() ),
				  volumeQuadratureElement( entity, ( 4 * Traits::pressureSpaceOrder ) + 1 ),
				  discrete_model( discrete_modelIn ),
				  eps( Parameters().getParam( "eps", 1.0e-14 ) ),
				  viscosity( discrete_modelIn.viscosity() ),
				  convection_scaling( discrete_modelIn.convection_scaling() ),
				  pressure_gradient_scaling( discrete_modelIn.pressure_gradient_scaling() ),
				  alpha( discrete_modelIn.alpha() )
			{}
		};
		struct InfoContainerFace : public InfoContainerVolume {
			const typename Traits::IntersectionIteratorType::Intersection& intersection;
			const typename Traits::IntersectionIteratorType::Geometry& intersectionGeometry;
			const typename Traits::FaceQuadratureType faceQuadratureElement;
			const double lengthOfIntersection;
			const StabilizationCoefficients& stabil_coeff;
			const double C_11;
			const double D_11;
			typename Traits::VelocityRangeType D_12;

			InfoContainerFace (const CoordinatorType& interface,
								const typename Traits::EntityType& ent,
							   const typename Traits::IntersectionIteratorType::Intersection& inter,
								const typename Traits::DiscreteModelType& discrete_modelIn )
				:InfoContainerVolume( interface, ent, discrete_modelIn ),
				  intersection( inter ),
				  intersectionGeometry( intersection.intersectionGlobal() ),
				  faceQuadratureElement( interface.sigma_space_.gridPart(),
																  intersection,
																  ( 4 * Traits::pressureSpaceOrder ) + 1,
																  Traits::FaceQuadratureType::INSIDE ),
				  lengthOfIntersection( Stuff::getLenghtOfIntersection( intersection ) ),
				  stabil_coeff( discrete_modelIn.getStabilizationCoefficients() ),
				  C_11( stabil_coeff.Factor("C11") * std::pow( lengthOfIntersection, stabil_coeff.Power("C11") ) ),
				  D_11( stabil_coeff.Factor("D11") * std::pow( lengthOfIntersection, stabil_coeff.Power("D11") ) ),
				  D_12( 1 )
			{
				D_12 /= D_12.two_norm();
				D_12 *= stabil_coeff.Factor("D12");
			}
		};

		struct InfoContainerInteriorFace : public InfoContainerFace {
			const typename Traits::EntityType& neighbour;
			const SigmaBaseFunctionSetType
					sigma_basefunction_set_neighbour;
			const VelocityBaseFunctionSetType
					velocity_basefunction_set_neighbour;
			const PressureBaseFunctionSetType
					pressure_basefunction_set_neighbour;
			const int numSigmaBaseFunctionsNeighbour;
			const int numVelocityBaseFunctionsNeighbour;
			const int numPressureBaseFunctionsNeighbour;
			const typename Traits::FaceQuadratureType faceQuadratureNeighbour;

			InfoContainerInteriorFace (const CoordinatorType& interface,
								const typename Traits::EntityType& ent,
							   const typename Traits::EntityType& nei,
							   const typename Traits::IntersectionIteratorType::Intersection& inter,
								const typename Traits::DiscreteModelType& discrete_modelIn )
				:InfoContainerFace( interface, ent, inter, discrete_modelIn ),
				  neighbour( nei ),
				  sigma_basefunction_set_neighbour( interface.sigma_space_.baseFunctionSet( neighbour ) ),
				  velocity_basefunction_set_neighbour( interface.velocity_space_.baseFunctionSet( neighbour ) ),
				  pressure_basefunction_set_neighbour( interface.pressure_space_.baseFunctionSet( neighbour ) ),
				  numSigmaBaseFunctionsNeighbour( sigma_basefunction_set_neighbour.numBaseFunctions() ),
				  numVelocityBaseFunctionsNeighbour( velocity_basefunction_set_neighbour.numBaseFunctions() ),
				  numPressureBaseFunctionsNeighbour( pressure_basefunction_set_neighbour.numBaseFunctions() ),
				  faceQuadratureNeighbour( interface.sigma_space_.gridPart(),
																  inter,
																  ( 4 * Traits::pressureSpaceOrder ) + 1,
																  Traits::FaceQuadratureType::OUTSIDE )
			{
				//some integration logic depends on this
				assert( InfoContainerFace::faceQuadratureElement.nop() == faceQuadratureNeighbour.nop() );
			}
		};

		struct ApplyVolume {
			template < class IntegratorType >
			static void apply( IntegratorType& integrator, const InfoContainerVolume& info )
			{
				integrator.applyVolume( info );
			}
		};

		struct ApplyInteriorFace {
			template < class IntegratorType >
			static void apply( IntegratorType& integrator, const InfoContainerInteriorFace& info )
			{
				integrator.applyInteriorFace( info );
			}
		};

		struct ApplyBoundaryFace {
			template < class IntegratorType >
			static void apply( IntegratorType& integrator, const InfoContainerFace& info )
			{
				integrator.applyBoundaryFace( info );
			}
		};

		//! A slightly modified version of fem/misc/utility:ForEachValue
		template < class FunctorType, class InfoType >
		class ForEachIntegrator {
		public:
			//! \brief Constructor
			//! \param tuple The tuple which we want to process.
			ForEachIntegrator(IntegratorTuple& tuple, const InfoType& info)
				: tuple_(tuple),
				  info_(info)
			{
				apply( tuple_ );
			}
		private:
			//! Specialisation for the last element
			template <class Head>
			void apply(Pair<Head, Nil>& last) {
				FunctorType::apply( last.first(), info_ );
			}

			//! Specialisation for a standard tuple element
			template <class Head, class Tail>
			void apply(Pair<Head, Tail>& pair) {
				FunctorType::apply( pair.first(), info_ );
				apply(pair.second());
			}
		private:
			IntegratorTuple& tuple_;
			const InfoType& info_;
		};

		void apply ( IntegratorTuple& integrator_tuple ) const
		{
			Profiler::ScopedTiming assembler_time("assembler");
			typedef typename GridType::LeafGridView
				GridView;
			typedef typename GridView::template Codim<0>::
					template Partition<Dune::All_Partition>::Iterator
				LeafIterator;

			const GridView gridView = grid_part_.grid().leafView();

			const LeafIterator endit = gridView.template end<0,Dune::All_Partition>();
			for ( LeafIterator entityIt = gridView.template begin<0,Dune::All_Partition>();
				 entityIt!=endit;
				 ++entityIt)
			{
				const typename Traits::EntityType& entity = *entityIt;
				InfoContainerVolume info( *this, entity, discrete_model_ );
				ForEachIntegrator<ApplyVolume,InfoContainerVolume>( integrator_tuple, info );

				// walk the intersections
				const typename Traits::IntersectionIteratorType intItEnd = gridView.iend( entity );
				for (   typename Traits::IntersectionIteratorType intIt = gridView.ibegin( entity );
						intIt != intItEnd;
						++intIt )
				{
					const typename Traits::IntersectionIteratorType::Intersection& intersection = *intIt;

					// if we are inside the grid
					if ( intersection.neighbor() && !intersection.boundary() )
					{
						//! DO NOT TRY TO DEREF outside() DIRECTLY
						const typename Traits::IntersectionIteratorType::EntityPointer neighbourPtr = intersection.outside();
						InfoContainerInteriorFace info( *this, entity, *neighbourPtr, intersection, discrete_model_ );
						ForEachIntegrator<ApplyInteriorFace,InfoContainerInteriorFace>( integrator_tuple, info );
					}
					else if ( !intersection.neighbor() && intersection.boundary() )
					{
						InfoContainerFace info( *this, entity, intersection, discrete_model_ );
						ForEachIntegrator<ApplyBoundaryFace,InfoContainerFace>( integrator_tuple, info );
					}
				}
			}
		}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_BASE_HH
