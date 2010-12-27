#ifndef DUNE_STOKES_INTEGRATORS_X_HH
#define DUNE_STOKES_INTEGRATORS_X_HH

#include <dune/stokes/integrators/base.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class X
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
			X( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localXmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				// (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx // X's volume integral
				//                                                  // see also "X's entitity surface integral", "X's neighbour surface integral" and "X's boundary integral" below
				for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						double X_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute \tau_{j}:\nabla v_{i}
							SigmaRangeType tau_j( 0.0 );
							info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
							const double gradient_of_v_i_times_tau_j = info.velocity_basefunction_set_element.evaluateGradientSingle( i, entity, x, tau_j );
							X_i_j += elementVolume
								* integrationWeight
								* gradient_of_v_i_times_tau_j;
						} // done sum over quadrature points
						// if small, should be zero
						if ( fabs( X_i_j ) < eps ) {
							X_i_j = 0.0;
						}
						else
							// add to matrix
							localXmatrixElement.add( i, j, X_i_j );
					}
				} // done computing X's volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localXmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				typename MatrixObjectType::LocalMatrixType
						localXmatrixNeighbour = matrix_object_.localMatrix( info.entity, info.neighbour );
				// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's element sourface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}ds // X's neighbour sourface integral
				//                                                                                                                   // see also "X's boundary integral" below
				//                                                                                                                   // and "X's volume integral" above
//                        if ( info.discrete_model.hasSigmaFlux() ) {
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						// compute X's element sourface integral
						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							double X_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x in codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								SigmaRangeType tau_j( 0.0 );
								info.sigma_basefunction_set_element.evaluate( j, x, tau_j );

								VelocityRangeType tau_j_times_normal( 0.0 );
								tau_j.mv( outerNormal, tau_j_times_normal );
								typename Traits::C12 c_12( outerNormal, info.discrete_model.getStabilizationCoefficients() );
								SigmaRangeType tau_j_times_normal_dyadic_C12
										= dyadicProduct( tau_j_times_normal, c_12 );

								SigmaRangeType flux_value = tau_j;
								flux_value *= 0.5;
								flux_value -= tau_j_times_normal_dyadic_C12;

								VelocityRangeType flux_times_normal( 0.0 );
								flux_value.mv( outerNormal, flux_times_normal );

								const double v_i_times_flux_times_normal
										= info.velocity_basefunction_set_element.evaluateSingle( i, x, flux_times_normal );
								X_i_j -= elementVolume
									* integrationWeight
									* v_i_times_flux_times_normal;

							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( X_i_j ) < eps ) {
								X_i_j = 0.0;
							}
							else
								// add to matrix
								localXmatrixElement.add( i, j, X_i_j );
						} // done computing X's element sourface integral
						// compute X's neighbour sourface integral
						for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
							double X_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureNeighbour.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
								const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureNeighbour.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								SigmaRangeType tau_j( 0.0 );
								sigmaBaseFunctionSetNeighbour.evaluate( j, xOutside, tau_j );

								VelocityRangeType tau_j_times_normal( 0.0 );
								tau_j.mv( outerNormal, tau_j_times_normal );
								typename Traits::C12 c_12( outerNormal, info.discrete_model.getStabilizationCoefficients() );
								SigmaRangeType tau_j_times_normal_dyadic_C12
										= dyadicProduct( tau_j_times_normal, c_12 );

								SigmaRangeType flux_value = tau_j;
								flux_value *= 0.5;
								flux_value += tau_j_times_normal_dyadic_C12;

								VelocityRangeType flux_times_normal( 0.0 );
								flux_value.mv( outerNormal, flux_times_normal );

								const double v_i_times_flux_times_normal
										= velocityBaseFunctionSetNeighbour.evaluateSingle( i, xInside, flux_times_normal );
								X_i_j -= elementVolume
									* integrationWeight
									* v_i_times_flux_times_normal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( X_i_j ) < eps ) {
								X_i_j = 0.0;
							}
							else
								// add to matrix
								localXmatrixNeighbour.add( i, j, X_i_j );
						} // done computing X's neighbour sourface integral
					} // done computing X's sourface integrals
//                        }
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localXmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's boundary integral
				//                                                                                                                   // see also "X's volume integral", "X's element surface integral" and "X's neighbour surface integral" above
//                        if ( info.discrete_model.hasSigmaFlux() ) {
					for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
						for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
							double X_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								SigmaRangeType tau_j( 0.0 );
								info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								VelocityRangeType tau_times_normal( 0.0 );
								tau_j.mv( outerNormal, tau_times_normal );
								const double v_i_times_tau_times_normal = v_i * tau_times_normal;
								X_i_j += -1.0
									* elementVolume
									* integrationWeight
									* v_i_times_tau_times_normal;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( X_i_j ) < eps ) {
								X_i_j = 0.0;
							}
							else
								// add to matrix
								localXmatrixElement.add( i, j, X_i_j );
						}
					} // done computing X's boundary integral
			}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_X_HH
