#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include <dune/stuff/grid.hh>
#include <dune/stuff/misc.hh>
#include <dune/fem/misc/femtuples.hh>

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

	template < class Traits, class MatrixIntegratorTuple >
	class MatrixInterface
	{
	protected:
		typedef MatrixInterface < Traits, MatrixIntegratorTuple >
			MatrixInterfaceType;

		// VelocityRangeType is expected to be a FieldVector,
		// SigmaJacobianRangeType to be a Matrixmapping and
		// SigmaJacobianRangeType[i] to be a FieldVector
		//! \todo   doc me

		const typename Traits::DiscreteModelType&					discrete_model_;
		const typename Traits::GridPartType&						grid_part_;
		const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space_;
		const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space_;
		const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space_;

		const double eps_,viscosity_;

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

		MatrixInterface(const typename Traits::DiscreteModelType&					discrete_model,
						const typename Traits::GridPartType&						grid_part,
						const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space,
						const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space,
						const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space)
					:discrete_model_(discrete_model),
					grid_part_(grid_part),
					velocity_space_(velocity_space),
					pressure_space_(pressure_space),
					sigma_space_(sigma_space),
					eps_( Parameters().getParam( "eps", 1.0e-14 ) ),
					viscosity_( discrete_model_.viscosity() )
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
			InfoContainerVolume(const MatrixInterfaceType& interface,
								const typename Traits::EntityType& ent,
								const typename Traits::DiscreteModelType& discrete_modelIn )
				:entity( ent ),
				  geometry( entity.geometry() ),
				  sigma_basefunction_set_element( interface.sigma_space_.baseFunctionSet( entity ) ),
				  velocity_basefunction_set_element( interface.velocity_space_.baseFunctionSet( entity ) ),
				  pressure_basefunction_set_element( interface.pressure_space_.baseFunctionSet( entity ) ),
				  numSigmaBaseFunctionsElement( sigma_basefunction_set_element.numBaseFunctions() ),
				  numVelocityBaseFunctionsElement( velocity_basefunction_set_element.numBaseFunctions() ),
				  numPressureBaseFunctionsElement( pressure_basefunction_set_element.numBaseFunctions() ),
				volumeQuadratureElement( entity, ( 4 * Traits::pressureSpaceOrder ) + 1 ),
				  discrete_model(discrete_modelIn),
				  eps(1e-14),
				  viscosity(discrete_modelIn.viscosity())
			{}
		};
		struct InfoContainerInteriorFace : public InfoContainerVolume {
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
			const typename Traits::IntersectionIteratorType::Intersection& intersection;
			const typename Traits::IntersectionIteratorType::Geometry& intersectionGeometry;
			const typename Traits::FaceQuadratureType faceQuadratureElement;
			const typename Traits::FaceQuadratureType faceQuadratureNeighbour;
			const double lengthOfIntersection;
			const StabilizationCoefficients& stabil_coeff;
			const double C_11;
			const double D_11;
			typename Traits::VelocityRangeType D_12;

			InfoContainerInteriorFace (const MatrixInterfaceType& interface,
								const typename Traits::EntityType& ent,
							   const typename Traits::EntityType& nei,
							   const typename Traits::IntersectionIteratorType::Intersection& inter,
								const typename Traits::DiscreteModelType& discrete_modelIn )
				:InfoContainerVolume( interface, ent, discrete_modelIn ),
				  neighbour( nei ),
				  sigma_basefunction_set_neighbour( interface.sigma_space_.baseFunctionSet( neighbour ) ),
				  velocity_basefunction_set_neighbour( interface.velocity_space_.baseFunctionSet( neighbour ) ),
				  pressure_basefunction_set_neighbour( interface.pressure_space_.baseFunctionSet( neighbour ) ),
				  numSigmaBaseFunctionsNeighbour( sigma_basefunction_set_neighbour.numBaseFunctions() ),
				  numVelocityBaseFunctionsNeighbour( velocity_basefunction_set_neighbour.numBaseFunctions() ),
				  numPressureBaseFunctionsNeighbour( pressure_basefunction_set_neighbour.numBaseFunctions() ),
				  intersection( inter ),
				  intersectionGeometry( intersection.intersectionGlobal() ),
				  faceQuadratureElement( interface.sigma_space_.gridPart(),
																  intersection,
																  ( 4 * Traits::pressureSpaceOrder ) + 1,
																  Traits::FaceQuadratureType::INSIDE ),
				  faceQuadratureNeighbour( interface.sigma_space_.gridPart(),
																  intersection,
																  ( 4 * Traits::pressureSpaceOrder ) + 1,
																  Traits::FaceQuadratureType::OUTSIDE ),
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

		struct ApplyVolume {
			template < class IntegratorType, class InfoType >
			static void apply( IntegratorType& integrator, const InfoType& info )
			{
				integrator.applyVolume( info );
			}
		};

		struct ApplyInteriorFace {
			template < class IntegratorType, class InfoType >
			static void apply( IntegratorType& integrator, const InfoType& info )
			{
				integrator.applyInteriorFace( info );
			}
		};

		struct ApplyBoundaryFace {
			template < class IntegratorType, class InfoType >
			static void apply( IntegratorType& integrator, const InfoType& info )
			{
				integrator.applyBoundaryFace( info );
			}
		};

		template < class FunctorType, class InfoType >
		class ForEachIntegrator {
		public:
		  //! \brief Constructor
		  //! \param tuple The tuple which we want to process.
		  ForEachIntegrator(MatrixIntegratorTuple& tuple, const InfoType& info)
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
		  MatrixIntegratorTuple& tuple_;
		  const InfoType& info_;
		};

		void apply ( MatrixIntegratorTuple& matrix_integrator_tuple ) const
		{
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
				ForEachIntegrator<ApplyVolume,InfoContainerVolume>(matrix_integrator_tuple, info );

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
						ForEachIntegrator<ApplyInteriorFace,InfoContainerInteriorFace>(matrix_integrator_tuple, info );
					}
					else if ( !intersection.neighbor() && intersection.boundary() )
					{
//						applyBoundaryFace();
					}
				}
			}
		}

	};


	template < class MatrixObjectType, class Traits >
	class W
	{
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
		typedef typename Traits::SigmaJacobianRangeType
			SigmaJacobianRangeType;
		typedef typename Traits::LocalIntersectionCoordinateType
			LocalIntersectionCoordinateType;


		MatrixObjectType& matrix_object_;
		public:
			W( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						local_matrix = matrix_object_.localMatrix( info.entity, info.entity );
				const double viscosity = info.discrete_model.viscosity();
					//                                                        // we will call this one
					// (W)_{i,j} += \mu\int_{T}v_{j}\cdot(\nabla\cdot\tau_{i})dx // W's volume integral
					//                                                        // see also "W's entitity surface integral", "W's neighbour surface integral" and "W's boundary integral" below
					for ( int i = 0; i < info.numSigmaBaseFunctionsElement; ++i ) {
						for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
							double W_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
								// get x
								const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
								// get the integration factor
								const double elementVolume = info.geometry.integrationElement( x );
								// get the quadrature weight
								const double integrationWeight = info.volumeQuadratureElement.weight( quad );
								// compute v_j^t \cdot ( \nabla \cdot \tau_i^t )
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, x, v_j );
								SigmaJacobianRangeType v_j_temp_for_div
										= prepareVelocityRangeTypeForSigmaDivergence<SigmaJacobianRangeType,VelocityRangeType>( v_j );
								const double divergence_of_tau_i_times_v_j
										= info.sigma_basefunction_set_element.evaluateGradientSingle( i, info.entity,
																									 x, v_j_temp_for_div );
								W_i_j += elementVolume
									* integrationWeight
									* viscosity
									* divergence_of_tau_i_times_v_j;
							} // done sum over quadrature points
							// if small, should be zero
							if ( fabs( W_i_j ) < info.eps ) {
								W_i_j = 0.0;
							}
							else
								// add to matrix
								local_matrix.add( i, j, W_i_j );
						}
					} // done computing W's volume integral
			}

			template < class InfoContainerFaceType >
			void applyInteriorFace( const InfoContainerFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localWmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				typename MatrixObjectType::LocalMatrixType
						localWmatrixNeighbour = matrix_object_.localMatrix( info.entity, info.entity );

				//                                                                                                               // we will call this one
				// (W)_{i,j} += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's element surface integral
				//           += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's neighbour surface integral
				//                                                                                                               // see also "W's boundary integral" below
				//                                                                                                               // and "W's volume integral" above
				//                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
				for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
					// compute W's element surface integral
					for ( int i = 0; i < info.numSigmaBaseFunctionsElement; ++i ) {
						double W_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
							// get x in codim<0> and codim<1> coordinates
							const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
							const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
							// get the integration factor
							const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
							// get the quadrature weight
							const double integrationWeight = info.faceQuadratureElement.weight( quad );
							// compute \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
							const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
							SigmaRangeType tau_i( 0.0 );
							info.sigma_basefunction_set_element.evaluate( i, x, tau_i );

							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );

							VelocityJacobianRangeType v_j_dyadic_normal
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_j, outerNormal );
							typename Traits::C12 c_12( outerNormal, info.discrete_model.getStabilizationCoefficients() );
							VelocityRangeType v_j_dyadic_normal_times_C12( 0.0 );
							v_j_dyadic_normal.mv( c_12, v_j_dyadic_normal_times_C12 );

							VelocityRangeType tau_i_times_normal( 0.0 );
							tau_i.mv( outerNormal, tau_i_times_normal );

							VelocityRangeType flux_value = v_j;
							flux_value *= 0.5;
							flux_value += v_j_dyadic_normal_times_C12;
							double flux_times_tau_i_times_normal = flux_value * tau_i_times_normal;
							W_i_j -= elementVolume
									* integrationWeight
									* info.viscosity
									* flux_times_tau_i_times_normal;
						} // done sum over all quadrature points
						// if small, should be zero
						if ( fabs( W_i_j ) < info.eps ) {
							W_i_j = 0.0;
						}
						else
							// add to matrix
							localWmatrixElement.add( i, j, W_i_j );
					} // done computing W's element surface integral
					// compute W's neighbour surface integral
					for ( int i = 0; i < info.numSigmaBaseFunctionsNeighbour; ++i ) {
						double W_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
							// get x in codim<0> and codim<1> coordinates
							const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
							const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
							const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
							// get the integration factor
							const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
							// get the quadrature weight
							const double integrationWeight = info.faceQuadratureElement.weight( quad );
							// compute \hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{j}\cdot n_{T}
							const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
							SigmaRangeType tau_i( 0.0 );
							info.sigma_basefunction_set_neighbour.evaluate( i, xOutside, tau_i );

							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );

							VelocityJacobianRangeType v_j_dyadic_normal
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_j, outerNormal );
							VelocityRangeType v_j_dyadic_normal_times_C12( 0.0 );
							typename Traits::C12 c_12( outerNormal, info.discrete_model.getStabilizationCoefficients() );
							v_j_dyadic_normal.mv( c_12, v_j_dyadic_normal_times_C12 );

							VelocityRangeType tau_i_times_normal( 0.0 );
							tau_i.mv( outerNormal, tau_i_times_normal );

							VelocityRangeType flux_value = v_j;
							flux_value *= 0.5;
							flux_value -= v_j_dyadic_normal_times_C12;
							double flux_times_tau_i_times_normal = flux_value * tau_i_times_normal;
							W_i_j += elementVolume
									* integrationWeight
									* info.viscosity
									* flux_times_tau_i_times_normal;
						} // done sum over all quadrature points
						// if small, should be zero
						if ( fabs( W_i_j ) < info.eps ) {
							W_i_j = 0.0;
						}
						else
							// add to matrix
							localWmatrixNeighbour.add( i, j, W_i_j );
					} // done computing W's neighbour surface integral
				} // done computing W's surface integrals
			}
	};
} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_HH
