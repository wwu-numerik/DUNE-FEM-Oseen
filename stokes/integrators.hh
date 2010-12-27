#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include "integrators_base.hh"

namespace Dune {
namespace Stokes {
namespace Integrators {

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
