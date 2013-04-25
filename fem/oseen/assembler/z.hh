#ifndef DUNE_OSEEN_INTEGRATORS_Z_HH
#define DUNE_OSEEN_INTEGRATORS_Z_HH

#include <dune/fem/oseen/assembler/base.hh>
#include <dune/stuff/common/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

	template < class MatrixPointerType, class Traits >
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
		typedef DSFe::LocalMatrixProxy<MatrixPointerType>
			LocalMatrixProxyType;

		MatrixPointerType& matrix_pointer_;
		public:
			Z( MatrixPointerType& matrix_object	)
				:matrix_pointer_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				typename MatrixPointerType::element_type::LocalMatrixType
						localZmatrixElement = matrix_pointer_->localMatrix( info.entity, info.entity );
				// (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx // Z's volume integral
				//                                                  // see also "Z's entitity surface integral", "Z's neighbour surface integral" and "Z's boundary integral" below
				for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++ quad ) {
					const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
					// get the integration factor
					const double elementVolume = info.geometry.integrationElement( x );
					// get the quadrature weight
					const double integrationWeight = info.volumeQuadratureElement.weight( quad );
					// compute q_{j}\cdot(\nabla\cdot v_i)
					PressureRangeType q_j( 0.0 );
					for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
						info.pressure_basefunction_set_element.evaluate( j, x, q_j );
						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							const double divergence_of_v_i_times_q_j =
									info.velocity_basefunction_set_element.evaluateGradientSingle( i, info.entity, x,
																								  preparePressureRangeTypeForVelocityDivergence<Traits>( q_j ) );
							const double Z_i_j = -1.0
								* elementVolume
								* integrationWeight
								* info.pressure_gradient_scaling
								* divergence_of_v_i_times_q_j;
							localZmatrixElement.add( i, j, Z_i_j );
						}
					}
				} // done computing Z's volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localZmatrixElement( matrix_pointer_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localZmatrixNeighbour( matrix_pointer_, info.entity, info.neighbour, info.eps );
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
								const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, x, v_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_element.evaluate( j, x, q_j );
								const double v_i_times_normal = v_i * outerNormal;

								const double p_factor = ( 0.5 - ( info.D_12 * outerNormal ) );// (0.5 p - p D_12 ) n ) <- p+
//										const double p_factor = ( 0.5 - ( 1 ) );// (0.5 p - p D_12 ) n ) <- p+
								const double q_j_times_v_i_times_normal =  q_j * v_i_times_normal;
								Z_i_j += p_factor
									* elementVolume
									* integrationWeight
									* info.pressure_gradient_scaling
									* q_j_times_v_i_times_normal;
							} // done sum over all quadrature points
							localZmatrixElement.add( i, j, Z_i_j );
						} // done computing Z's element surface integral
						// compute Z's neighbour surface integral
						for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i ) {
							double Z_i_j = 0.0;
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
								// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
								const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_i( 0.0 );
								info.velocity_basefunction_set_element.evaluate( i, xInside, v_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_neighbour.evaluate( j, xOutside, q_j );
								const double v_i_times_normal = v_i * outerNormal;

								const double p_factor = ( 0.5 + ( info.D_12 * outerNormal ) );// (0.5 p + p D_12 ) n ) <- p-
//										const double p_factor = ( 0.5 + ( 1 ) );// (0.5 p + p D_12 ) n ) <- p-
								const double q_j_times_v_i_times_normal = q_j * v_i_times_normal;
								Z_i_j += p_factor
									* elementVolume
									* integrationWeight
									* info.pressure_gradient_scaling
									* q_j_times_v_i_times_normal;
							} // done sum over all quadrature points
							localZmatrixNeighbour.add( i, j, Z_i_j );
						} // done computing Z's neighbour surface integral
					} // done computing Z's surface integrals
				//}
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				LocalMatrixProxyType localZmatrixElement( matrix_pointer_, info.entity, info.entity, info.eps );
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
								const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
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
									* info.pressure_gradient_scaling
									* q_j_times_v_times_normal;
							} // done sum over all quadrature points
							localZmatrixElement.add( i, j, Z_i_j );
						}
					} // done computing Z's boundary integral
//                        }
			}
			static const std::string name;
	};

	template < class T, class R > const std::string Z<T,R>::name = "Z";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_Z_HH

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

