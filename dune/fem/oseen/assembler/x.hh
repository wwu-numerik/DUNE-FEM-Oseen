#ifndef DUNE_OSEEN_INTEGRATORS_X_HH
#define DUNE_OSEEN_INTEGRATORS_X_HH

#include <dune/fem/oseen/assembler/base.hh>
#include <dune/stuff/common/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

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
		typedef DSFe::LocalMatrixProxy<MatrixObjectType>
			LocalMatrixProxyType;

		MatrixObjectType& matrix_object_;
		public:
			X( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				LocalMatrixProxyType localXmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				// (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx // X's volume integral
				//                                                  // see also "X's entitity surface integral", "X's neighbour surface integral" and "X's boundary integral" below
				for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						double X_i_j = 0.0;
                        for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
                            const auto x = info.volumeQuadratureElement.point( quad );
                            const double elementVolume = info.geometry.integrationElement( x );
                            const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute \tau_{j}:\nabla v_{i}
							SigmaRangeType tau_j( 0.0 );
							info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
							const double gradient_of_v_i_times_tau_j = info.velocity_basefunction_set_element.evaluateGradientSingle( i, info.entity, x, tau_j );
							X_i_j += elementVolume
								* integrationWeight
								* gradient_of_v_i_times_tau_j;
                        }
						localXmatrixElement.add( i, j, X_i_j );
					}
                }
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localXmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localXmatrixNeighbour( matrix_object_, info.entity, info.neighbour, info.eps );
				// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's element sourface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}ds // X's neighbour sourface integral
				//                                                                                                                   // see also "X's boundary integral" below
				//                                                                                                                   // and "X's volume integral" above
//                        if ( info.discrete_model.hasSigmaFlux() ) {

				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
                    const auto xInside = info.faceQuadratureElement.point( quad );
					const auto xOutside = info.faceQuadratureNeighbour.point( quad );
					const auto xLocal = info.faceQuadratureElement.localPoint( quad );
                    const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                    const double integrationWeight = info.faceQuadratureElement.weight( quad );
					// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					typename Traits::C12 c_12( outerNormal, info.discrete_model.getStabilizationCoefficients() );
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						SigmaRangeType tau_j( 0.0 );
						info.sigma_basefunction_set_element.evaluate( j, xInside, tau_j );
						// compute X's element sourface integral
						VelocityRangeType tau_j_times_normal( 0.0 );
						tau_j.mv( outerNormal, tau_j_times_normal );

						SigmaRangeType tau_j_times_normal_dyadic_C12
								= DSC::dyadicProduct<SigmaRangeType,VelocityRangeType>( tau_j_times_normal, c_12 );
						SigmaRangeType flux_value = tau_j;
						flux_value *= 0.5;
						flux_value -= tau_j_times_normal_dyadic_C12;
						VelocityRangeType flux_times_normal( 0.0 );
						flux_value.mv( outerNormal, flux_times_normal );

						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							const double v_i_times_flux_times_normal
									= info.velocity_basefunction_set_element.evaluateSingle( i, xInside, flux_times_normal );
							const double X_i_j = -elementVolume
								* integrationWeight
								* v_i_times_flux_times_normal;
							localXmatrixElement.add( i, j, X_i_j );
                        }
						// compute X's neighbour sourface integral
						info.sigma_basefunction_set_element.evaluate( j, xOutside, tau_j );
						tau_j.mv( outerNormal, tau_j_times_normal );
						tau_j_times_normal_dyadic_C12
							= DSC::dyadicProduct<SigmaRangeType,VelocityRangeType>( tau_j_times_normal, c_12 );
						flux_value = tau_j;
						flux_value *= 0.5;
						flux_value -= tau_j_times_normal_dyadic_C12;
						flux_value.mv( outerNormal, flux_times_normal );
						for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i ) {
							// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}
							const double v_i_times_flux_times_normal
								= info.velocity_basefunction_set_element.evaluateSingle( i, xInside, flux_times_normal );
							const double X_i_j = -elementVolume
								* integrationWeight
								* v_i_times_flux_times_normal;
							localXmatrixNeighbour.add( i, j, X_i_j );
                        }
                    }
				}
			}
			//                        }

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				LocalMatrixProxyType localXmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's boundary integral
				//                                                                                                                   // see also "X's volume integral", "X's element surface integral" and "X's neighbour surface integral" above
//                        if ( info.discrete_model.hasSigmaFlux() ) {
					for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
						for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
							double X_i_j = 0.0;
                            for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
                                const auto x = info.faceQuadratureElement.point( quad );
								const auto xLocal = info.faceQuadratureElement.localPoint( quad );
                                const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
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
                            }
							localXmatrixElement.add( i, j, X_i_j );
						}
                    }
			}
			static const std::string name;
	};

	template < class T, class R > const std::string X<T,R>::name = "X";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_X_HH

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

