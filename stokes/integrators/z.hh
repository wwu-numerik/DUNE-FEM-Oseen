#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include <dune/stokes/integrators/base.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class Z
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
			Z( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localZmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				// (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx // Z's volume integral
				//                                                  // see also "Z's entitity surface integral", "Z's neighbour surface integral" and "Z's boundary integral" below
				for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
						double Z_i_j = 0.0;
						// sum over all quadratur points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++ quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute q_{j}\cdot(\nabla\cdot v_i)
							PressureRangeType q_j( 0.0 );
							info.pressure_basefunction_set_element.evaluate( j, x, q_j );
							const double divergence_of_v_i_times_q_j = info.velocity_basefunction_set_element.evaluateGradientSingle( i, entity, x, preparePressureRangeTypeForVelocityDivergence( q_j ) );
							Z_i_j += -1.0
								* elementVolume
								* integrationWeight
								* pressure_gradient_scaling
								* divergence_of_v_i_times_q_j;
						} // done sum over all quadrature points
						// if small, should be zero
						if ( fabs( Z_i_j ) < eps ) {
							Z_i_j = 0.0;
						}
						else
							// add to matrix
							localZmatrixElement.add( i, j, Z_i_j );
					}
				} // done computing Z's volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localZmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				typename MatrixObjectType::LocalMatrixType
						localZmatrixNeighbour = matrix_object_.localMatrix( info.entity, info.neighbour );
				// (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's element surface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{-}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's neighbour surface integral
				//                                                                                                  // see also "Z's boundary integral" below
				//                                                                                                  // and "Z's volume integral" above
//                        if ( info.discrete_model.hasPressureFlux() ) {
					for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
						// compute Z's element surface integral
						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							double Z_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_element.evaluate( j, x, q_j );
								const double v_i_times_normal = v_i * outerNormal;

								const double p_factor = ( 0.5 - ( D_12 * outerNormal ) );// (0.5 p - p D_12 ) n ) <- p+
//										const double p_factor = ( 0.5 - ( 1 ) );// (0.5 p - p D_12 ) n ) <- p+
								const double q_j_times_v_i_times_normal =  q_j * v_i_times_normal;
								Z_i_j += p_factor
									* elementVolume
									* integrationWeight
									* pressure_gradient_scaling
									* q_j_times_v_i_times_normal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( Z_i_j ) < eps ) {
								Z_i_j = 0.0;
							}
							else
								// add to matrix
								localZmatrixElement.add( i, j, Z_i_j );
						} // done computing Z's element surface integral
						// compute Z's neighbour surface integral
						for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
							double Z_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureNeighbour.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
								const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureNeighbour.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
								// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, xInside, v_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_neighbour.evaluate( j, xOutside, q_j );
								const double v_i_times_normal = v_i * outerNormal;

								const double p_factor = ( 0.5 + ( D_12 * outerNormal ) );// (0.5 p + p D_12 ) n ) <- p-
//										const double p_factor = ( 0.5 + ( 1 ) );// (0.5 p + p D_12 ) n ) <- p-
								const double q_j_times_v_i_times_normal = q_j * v_i_times_normal;
								Z_i_j += p_factor
									* elementVolume
									* integrationWeight
									* pressure_gradient_scaling
									* q_j_times_v_i_times_normal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( Z_i_j ) < eps ) {
								Z_i_j = 0.0;
							}
							else
								// add to matrix
								localZmatrixNeighbour.add( i, j, Z_i_j );
						} // done computing Z's neighbour surface integral
					} // done computing Z's surface integrals
				//}
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localZmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				// (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's boundary integral
				//                                                                                                  // see also "Z's volume integral", "Z's element surface integral" and "Z's neighbour surface integral" above
//                        if ( info.discrete_model.hasPressureFlux() ) {
					for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
						// compute the boundary integral
						for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
							double Z_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_element.evaluate( j, x, q_j );
								const double v_i_times_normal = v_i * outerNormal;
								const double q_j_times_v_times_normal = q_j * v_i_times_normal;
								Z_i_j += elementVolume
									* integrationWeight
									* pressure_gradient_scaling
									* q_j_times_v_times_normal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( Z_i_j ) < eps ) {
								Z_i_j = 0.0;
							}
							else
								// add to matrix
								localZmatrixElement.add( i, j, Z_i_j );
						}
					} // done computing Z's boundary integral
//                        }
			}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_HH
