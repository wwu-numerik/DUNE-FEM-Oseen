#ifndef DUNE_STOKES_INTEGRATORS_R_HH
#define DUNE_STOKES_INTEGRATORS_R_HH

#include <dune/stokes/integrators/base.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class R
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
			R( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localRmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
				typename MatrixObjectType::LocalMatrixType
						localRmatrixNeighbour = matrix_object_.localMatrix( info.entity, info.neighbour );
				// (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}ds // R's element surface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}ds // R's neighbour surface integral
				//                                                                                                // see also "R's boundary integral" below
//                        if ( info.discrete_model.hasVelocityPressureFlux() ) {
					for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
						// compute R's element surface integral
						for ( int i = 0; i < info.numPressureBaseFunctionsElement; ++i ) {
							double R_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_element.evaluate( i, x, q_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_element.evaluate( j, x, q_j );
								const double q_i_times_q_j = q_i * q_j;
								R_i_j += D_11
									* elementVolume
									* integrationWeight
									* q_i_times_q_j;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( R_i_j ) < eps ) {
								R_i_j = 0.0;
							}
							else
								// add to matrix
								localRmatrixElement.add( i, j, R_i_j );
						} // done computing R's element surface integral
						// compute R's neighbour surface integral
						for ( int i = 0; i < info.numPressureBaseFunctionsNeighbour; ++i ) {
							double R_i_j = 0.0;
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
								// compute \hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_neighbour.evaluate( j, xOutside, q_j );
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_element.evaluate( i, xInside, q_i );
								const double q_i_times_q_j = q_i * q_j;
								R_i_j += -1.0
									* D_11
									* elementVolume
									* integrationWeight
									* q_i_times_q_j;
							} // done sum over all quadrature points
							// if small, should be zero
							if ( fabs( R_i_j ) < eps ) {
								R_i_j = 0.0;
							}
							else
								// add to matrix
								localRmatrixNeighbour.add( i, j, R_i_j );
						} // done computing R's neighbour surface integral
					} // done computing R's surface integrals
//                        }
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& )
			{}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_R_HH
