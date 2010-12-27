#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include <dune/stokes/integrators/base.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class E
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
			E( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localEmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				// (E)_{i,j} += -\int_{T}v_{j}\cdot\nabla q_{i}dx // E's volume integral
				//                                                // see also "E's entitity surface integral", "E's neighbour surface integral" and "E's boundary integral" below
				for ( int i = 0; i < info.numPressureBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double E_i_j = 0.0;
						// sum over all quadratur points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++ quad ) {
							// get x
							ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
							double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute v_{j}\cdot(\nabla q_i)
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							PressureJacobianRangeType jacobian_of_q_i( 0.0 );
							info.pressure_basefunction_set_element.jacobian( i, x, jacobian_of_q_i );
							const VelocityRangeType gradient_of_q_i_untransposed( jacobian_of_q_i[0] );
							const JacobianInverseTransposedType jacobianInverseTransposed = info.geometry.jacobianInverseTransposed( x );
							VelocityRangeType gradient_of_q_i( 0.0 );
							jacobianInverseTransposed.mv( gradient_of_q_i_untransposed, gradient_of_q_i );
							const double gradient_of_q_i_times_v_j = gradient_of_q_i * v_j;
							E_i_j += -1.0
								* elementVolume
								* integrationWeight
								* gradient_of_q_i_times_v_j;
						} // done sum over all quadrature points
						// if small, should be zero
						if ( fabs( E_i_j ) < eps ) {
							E_i_j = 0.0;
						}
						else
							// add to matrix
							localEmatrixElement.add( i, j, E_i_j );
					}
				} // done computing E's volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localEmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				typename MatrixObjectType::LocalMatrixType
						localEmatrixNeighbour = matrix_object_.localMatrix( info.neighbour, info.entity );
				// (E)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}ds // E's element surface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}ds // E's neighbour surface integral
				//                                                                                                // see also "E's boundary integral" below
				//                                                                                                // and "E's volume integral" above
//                        if ( info.discrete_model.hasVelocityPressureFlux() ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						// compute E's element surface integral
						for ( int i = 0; i < info.numPressureBaseFunctionsElement; ++i ) {
							double E_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, x, v_j );
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_element.evaluate( i, x, q_i );
								VelocityRangeType flux_value = v_j;
								flux_value *= 0.5;
								const double v_j_times_outerNormal = v_j * outerNormal;
								VelocityRangeType jump = D_12;
								jump *= v_j_times_outerNormal;
								flux_value += jump;
								VelocityRangeType q_i_times_flux = flux_value;
								q_i_times_flux *= q_i;
								const double q_i_times_flux_times_outerNormal = q_i_times_flux * outerNormal;

								E_i_j += elementVolume
									* integrationWeight
									* q_i_times_flux_times_outerNormal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( E_i_j ) < eps ) {
								E_i_j = 0.0;
							}
							else
								// add to matrix
								localEmatrixElement.add( i, j, E_i_j );
						} // done computing E's element surface integral
						// compute E's neighbour surface integral
						for ( int i = 0; i < info.numPressureBaseFunctionsNeighbour; ++i ) {
							double E_i_j = 0.0;
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
								// compute \hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_neighbour.evaluate( i, xOutside, q_i );

								VelocityRangeType flux_value = v_j;
								flux_value *= 0.5;
								const double v_j_times_outerNormal = v_j * outerNormal;
								VelocityRangeType jump = D_12;
								jump *= v_j_times_outerNormal;
								flux_value += jump;
								VelocityRangeType q_i_times_flux = flux_value;
								q_i_times_flux *= q_i;
								const double q_i_times_flux_times_outerNormal = q_i_times_flux * outerNormal;

								E_i_j -= elementVolume
									* integrationWeight
									* q_i_times_flux_times_outerNormal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( E_i_j ) < eps ) {
								E_i_j = 0.0;
							}
							else
								// add to matrix
								localEmatrixNeighbour.add( i, j, E_i_j );
						} // done computing E's neighbour surface integral
					} // done computing E's surface integrals
//                        }
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localWmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
			}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_HH
