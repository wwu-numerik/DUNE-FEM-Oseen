#ifndef DUNE_STOKES_INTEGRATORS_Y_HH
#define DUNE_STOKES_INTEGRATORS_Y_HH

#include <dune/stokes/integrators/base.hh>
#include <dune/stuff/matrix.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class Y
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
		typedef Stuff::Matrix::LocalMatrixProxy<MatrixObjectType>
			LocalMatrixProxyType;

		MatrixObjectType& matrix_object_;
		public:
			Y( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixObjectType::LocalMatrixType
						localYmatrixElement = matrix_object_.localMatrix( info.entity, info.entity );
//                if ( info.discrete_model.isGeneralized() )
				{
				for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double Y_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute \tau_{j}:\nabla v_{i}
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, x, v_i );
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							const double v_i_times_v_j = v_i * v_j;
							Y_i_j += elementVolume
								* integrationWeight
								* info.alpha
								* v_i_times_v_j;
						} // done sum over quadrature points
						// if small, should be zero
						if ( fabs( Y_i_j ) < info.eps ) {
							Y_i_j = 0.0;
						}
						else
							// add to matrix
							localYmatrixElement.add( i, j, Y_i_j );
					}
				} // done computing Y's volume integral
				}
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localYmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localYmatrixNeighbour( matrix_object_, info.neighbour, info.entity, info.eps );
				// (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}ds // Y's element surface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}ds // Y's neighbour surface integral
				//                                                                                                         // see also "Y's boundary integral" below
//                        if ( info.discrete_model.hasSigmaFlux() ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						// compute Y's element surface integral
						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							double Y_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}
//								const VelocityRangeType /*outerNormal*/ = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, x, v_j );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								const double v_i_times_v_j = v_i * v_j;
								Y_i_j += info.C_11
									* elementVolume
									* integrationWeight
									* v_i_times_v_j;
							} // done sum over all quadrature points
							localYmatrixElement.add( i, j, Y_i_j );
						} // done computing Y's element surface integral
						// compute Y's neighbour surface integral
						for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i ) {
							double Y_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureNeighbour.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
								const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureNeighbour.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}
//								const VelocityRangeType /*outerNormal*/ = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_i );
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								const double v_i_times_v_j = v_i * v_j;
								Y_i_j += -1.0
									* info.C_11
									* elementVolume
									* integrationWeight
									* v_i_times_v_j;
							} // done sum over all quadrature points
							localYmatrixNeighbour.add( i, j, Y_i_j );
						} // done computing Y's neighbour surface integral
					} // done computing Y's surface integrals
//                        }

			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				LocalMatrixProxyType localYmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				// (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}ds // Y's boundary integral
				//                                                                                                           // see also "Y's element surface integral" and "Y's neighbour surface integral" above
//                        if ( info.discrete_model.hasSigmaFlux() ) {
					for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
						for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
							double Y_i_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
								// get the integration factor
								const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute -\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}
//                                const VelocityRangeType /*outerNormal*/ = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, x, v_j );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								const double v_i_times_v_j = v_i * v_j;
								Y_i_j += info.C_11
									* elementVolume
									* integrationWeight
									* v_i_times_v_j;
							} // done sum over all quadrature points
							localYmatrixElement.add( i, j, Y_i_j );
						}
					} // done computing Y's boundary integral
//                        }

			}
			static const std::string name;
	};

	template < class T, class R > const std::string Y<T,R>::name = "Y";

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_Y_HH
