#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include <dune/stuff/grid.hh>
#include <dune/stuff/misc.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class Traits >
	class W
	{
		const typename Traits::DiscreteModelType&					discrete_model_;
		const typename Traits::GridPartType&						grid_part_;
		const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space_;
		const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space_;
		const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space_;

		const double eps;

		// VelocityRangeType is expected to be a FieldVector,
		// SigmaJacobianRangeType to be a Matrixmapping and
		// SigmaJacobianRangeType[i] to be a FieldVector
		//! \todo   doc me
		static typename Traits::SigmaJacobianRangeType prepareVelocityRangeTypeForSigmaDivergence(
													const typename Traits::VelocityRangeType& arg )
		{
			typename Traits::SigmaJacobianRangeType ret( 0.0 );
			assert( arg.dim() == ret[0].dim() );
			for ( int i = 0; i < int(arg.dim()) ; ++i ) {
				for ( int j = 0; j < int(arg.dim()); ++j ) {
					typename Traits::VelocityRangeType row( 0.0 );
					row[ j ] = arg[ i ];
					ret[ i * arg.dim() + j ] = row;
				}
			}
			return ret;
		}

		/**
		 *  \brief  dyadic product
		 *
		 *          Implements \f$\left(arg_{1} \otimes arg_{2}\right)_{i,j}:={arg_{1}}_{i} {arg_{2}}_{j}\f$
		 **/
		static typename Traits::SigmaRangeType dyadicProduct(
										const typename Traits::VelocityRangeType& arg1,
										const typename Traits::VelocityRangeType& arg2 )
		{
			typename Traits::SigmaRangeType ret( 0.0 );
			typedef typename Traits::SigmaRangeType::RowIterator
				MatrixRowIteratorType;
			typedef typename Traits::VelocityRangeType::ConstIterator
				ConstVectorIteratorType;
			typedef typename Traits::VelocityRangeType::Iterator
				VectorIteratorType;
			MatrixRowIteratorType rItEnd = ret.end();
			ConstVectorIteratorType arg1It = arg1.begin();
			for ( MatrixRowIteratorType rIt = ret.begin(); rIt != rItEnd; ++rIt ) {
				ConstVectorIteratorType arg2It = arg2.begin();
				VectorIteratorType vItEnd = rIt->end();
				for (   VectorIteratorType vIt = rIt->begin();
						vIt != vItEnd;
						++vIt ) {
					*vIt = *arg1It * *arg2It;
					++arg2It;
				}
				++arg1It;
			}
			return ret;
		}

		public:
			W(	const typename Traits::DiscreteModelType&					discrete_model,
				const typename Traits::GridPartType&						grid_part,
				const typename Traits::DiscreteVelocityFunctionSpaceType&	velocity_space,
				const typename Traits::DiscretePressureFunctionSpaceType&	pressure_space,
				const typename Traits::DiscreteSigmaFunctionSpaceType&		sigma_space)
			:discrete_model_(discrete_model),
			grid_part_(grid_part),
			velocity_space_(velocity_space),
			pressure_space_(pressure_space),
			sigma_space_(sigma_space),
			eps( Parameters().getParam( "eps", 1.0e-14 ) )
			{}

			template < class MatrixObjectType >
			void apply ( MatrixObjectType& matrix_object ) const
			{
				typedef typename GridType::LeafGridView
					GridView;
				typedef typename GridView::template Codim<0>::
						template Partition<Dune::All_Partition>::Iterator
					LeafIterator;
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

				const GridView gridView = grid_part_.grid().leafView();

				const LeafIterator endit = gridView.template end<0,Dune::All_Partition>();
				for ( LeafIterator entityIt = gridView.template begin<0,Dune::All_Partition>();
					 entityIt!=endit;
					 ++entityIt)
				{
					const typename Traits::EntityType& entity = *entityIt;
					const typename Traits::EntityType::Geometry& geometry = entity.geometry();

					typename MatrixObjectType::LocalMatrixType
							localWmatrixElement = matrix_object.localMatrix( entity, entity );
					const typename Traits::DiscreteSigmaFunctionSpaceType::BaseFunctionSetType
							sigma_basefunction_set_element = sigma_space_.baseFunctionSet( entity );
					const typename Traits::DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
							velocity_basefunction_set_element = velocity_space_.baseFunctionSet( entity );
					const typename Traits::DiscretePressureFunctionSpaceType::BaseFunctionSetType
							pressure_basefunction_set_element  = pressure_space_.baseFunctionSet( entity );
					const int numSigmaBaseFunctionsElement = sigma_basefunction_set_element.numBaseFunctions();
					const int numVelocityBaseFunctionsElement = velocity_basefunction_set_element.numBaseFunctions();
					const int numPressureBaseFunctionsElement = pressure_basefunction_set_element.numBaseFunctions();
					const typename Traits::VolumeQuadratureType volumeQuadratureElement( entity,
																		( 4 * Traits::pressureSpaceOrder ) + 1 );

					const double viscosity = discrete_model_.viscosity();
					//                                                        // we will call this one
					// (W)_{i,j} += \mu\int_{T}v_{j}\cdot(\nabla\cdot\tau_{i})dx // W's volume integral
					//                                                        // see also "W's entitity surface integral", "W's neighbour surface integral" and "W's boundary integral" below
					for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
						for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
							double W_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
								// get x
								const ElementCoordinateType x = volumeQuadratureElement.point( quad );
								// get the integration factor
								const double elementVolume = geometry.integrationElement( x );
								// get the quadrature weight
								const double integrationWeight = volumeQuadratureElement.weight( quad );
								// compute v_j^t \cdot ( \nabla \cdot \tau_i^t )
								VelocityRangeType v_j( 0.0 );
								velocity_basefunction_set_element.evaluate( j, x, v_j );
								const double divergence_of_tau_i_times_v_j
										= sigma_basefunction_set_element.evaluateGradientSingle( i, entity, x, prepareVelocityRangeTypeForSigmaDivergence( v_j ) );
								W_i_j += elementVolume
									* integrationWeight
									* viscosity
									* divergence_of_tau_i_times_v_j;
							} // done sum over quadrature points
							// if small, should be zero
							if ( fabs( W_i_j ) < eps ) {
								W_i_j = 0.0;
							}
							else
								// add to matrix
								localWmatrixElement.add( i, j, W_i_j );
						}
					} // done computing W's volume integral


					// walk the intersections
					const typename Traits::IntersectionIteratorType intItEnd = gridView.iend( entity );
					for (   typename Traits::IntersectionIteratorType intIt = gridView.ibegin( entity );
							intIt != intItEnd;
							++intIt )
					{
						const typename Traits::IntersectionIteratorType::Intersection& intersection = *intIt;

						// get intersection geometry
						typedef typename Traits::IntersectionIteratorType::Geometry
							IntersectionGeometryType;
						typedef typename Traits::LocalIntersectionCoordinateType
							LocalIntersectionCoordinateType;
						const IntersectionGeometryType& intersectionGeometry = intersection.intersectionGlobal();
						// get intersection quadrature, seen from inside
						const typename Traits::FaceQuadratureType faceQuadratureElement( grid_part_,
																		intersection,
																		( 4 * Traits::pressureSpaceOrder ) + 1,
																		Traits::FaceQuadratureType::INSIDE );
						const double lengthOfIntersection = Stuff::getLenghtOfIntersection( intersection );
						const StabilizationCoefficients& stabil_coeff ( discrete_model_.getStabilizationCoefficients() );
						const double C_11 = stabil_coeff.Factor("C11") * std::pow( lengthOfIntersection, stabil_coeff.Power("C11") );
						const double D_11 = stabil_coeff.Factor("D11") * std::pow( lengthOfIntersection, stabil_coeff.Power("D11") );
						//we'll leave this on 0 for the time being so it does not generate any additional  penalty terms
	//					const VelocityRangeType D_12(stabil_coeff.Factor("D12") );//TODO FIXME
						VelocityRangeType D_12( 1 );//TODO FIXME
						D_12 /= D_12.two_norm();
						D_12 *= stabil_coeff.Factor("D12");
	//					VelocityRangeType E_11(stabil_coeff.Factor("E12"));

						// if we are inside the grid
						if ( intersection.neighbor() && !intersection.boundary() )
						{
							//! DO NOT TRY TO DEREF outside() DIRECTLY
							const typename Traits::IntersectionIteratorType::EntityPointer neighbourPtr = intersection.outside();
							//there's some (copy ctor? / implicit type conversion) at play here which will make shit fail if you do
							const typename Traits::EntityType& neighbour = *neighbourPtr;
							typename MatrixObjectType::LocalMatrixType
									localWmatrixNeighbour = matrix_object.localMatrix( neighbour, entity );
							const typename Traits::DiscreteSigmaFunctionSpaceType::BaseFunctionSetType
									sigmaBaseFunctionSetNeighbour = sigma_space_.baseFunctionSet( neighbour );
							const typename Traits::DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
									velocityBaseFunctionSetNeighbour = velocity_space_.baseFunctionSet( neighbour );
							const typename Traits::DiscretePressureFunctionSpaceType::BaseFunctionSetType
									pressureBaseFunctionSetNeighbour = pressure_space_.baseFunctionSet( neighbour );
							const int numSigmaBaseFunctionsNeighbour = sigmaBaseFunctionSetNeighbour.numBaseFunctions();
							const int numVelocityBaseFunctionsNeighbour = velocityBaseFunctionSetNeighbour.numBaseFunctions();
							const int numPressureBaseFunctionsNeighbour = pressureBaseFunctionSetNeighbour.numBaseFunctions();

							// get intersection quadrature, seen from outside
							const typename Traits::FaceQuadratureType faceQuadratureNeighbour(   grid_part_,
																				intersection,
																				( 4 * Traits::pressureSpaceOrder ) + 1,
																				Traits::FaceQuadratureType::OUTSIDE );
							//                                                                                                               // we will call this one
							// (W)_{i,j} += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's element surface integral
							//           += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's neighbour surface integral
							//                                                                                                               // see also "W's boundary integral" below
							//                                                                                                               // and "W's volume integral" above
	//                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									// compute W's element surface integral
									for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
										double W_i_j = 0.0;
										// sum over all quadrature points
										for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
											// get x in codim<0> and codim<1> coordinates
											const ElementCoordinateType x = faceQuadratureElement.point( quad );
											const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
											// get the integration factor
											const double elementVolume = intersectionGeometry.integrationElement( xLocal );
											// get the quadrature weight
											const double integrationWeight = faceQuadratureElement.weight( quad );
											// compute \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
											const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
											SigmaRangeType tau_i( 0.0 );
											sigma_basefunction_set_element.evaluate( i, x, tau_i );

											VelocityRangeType v_j( 0.0 );
											velocity_basefunction_set_element.evaluate( j, x, v_j );

											VelocityJacobianRangeType v_j_dyadic_normal = dyadicProduct( v_j, outerNormal );
											typename Traits::C12 c_12( outerNormal, discrete_model_.getStabilizationCoefficients() );
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
												* viscosity
												* flux_times_tau_i_times_normal;
										} // done sum over all quadrature points
										// if small, should be zero
										if ( fabs( W_i_j ) < eps ) {
											W_i_j = 0.0;
										}
										else
											// add to matrix
											localWmatrixElement.add( i, j, W_i_j );
									} // done computing W's element surface integral
									// compute W's neighbour surface integral
									for ( int i = 0; i < numSigmaBaseFunctionsNeighbour; ++i ) {
										double W_i_j = 0.0;
										// sum over all quadrature points
										for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
											// get x in codim<0> and codim<1> coordinates
											const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
											const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
											const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
											// get the integration factor
											const double elementVolume = intersectionGeometry.integrationElement( xLocal );
											// get the quadrature weight
											const double integrationWeight = faceQuadratureElement.weight( quad );
											// compute \hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{j}\cdot n_{T}
											const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
											SigmaRangeType tau_i( 0.0 );
											sigmaBaseFunctionSetNeighbour.evaluate( i, xOutside, tau_i );

											VelocityRangeType v_j( 0.0 );
											velocity_basefunction_set_element.evaluate( j, xInside, v_j );

											VelocityJacobianRangeType v_j_dyadic_normal = dyadicProduct( v_j, outerNormal );
											VelocityRangeType v_j_dyadic_normal_times_C12( 0.0 );
											typename Traits::C12 c_12( outerNormal, discrete_model_.getStabilizationCoefficients() );
											v_j_dyadic_normal.mv( c_12, v_j_dyadic_normal_times_C12 );

											VelocityRangeType tau_i_times_normal( 0.0 );
											tau_i.mv( outerNormal, tau_i_times_normal );

											VelocityRangeType flux_value = v_j;
											flux_value *= 0.5;
											flux_value -= v_j_dyadic_normal_times_C12;
											double flux_times_tau_i_times_normal = flux_value * tau_i_times_normal;
											W_i_j += elementVolume
												* integrationWeight
												* viscosity
												* flux_times_tau_i_times_normal;
										} // done sum over all quadrature points
										// if small, should be zero
										if ( fabs( W_i_j ) < eps ) {
											W_i_j = 0.0;
										}
										else
											// add to matrix
											localWmatrixNeighbour.add( i, j, W_i_j );
									} // done computing W's neighbour surface integral
								} // done computing W's surface integrals
	//                        }

						}
						else if ( !intersection.neighbor() && intersection.boundary() )
						{

						}
					}
				}
			}
	};
} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_HH
