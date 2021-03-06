#ifndef DUNE_OSEEN_INTEGRATORS_O_HH
#define DUNE_OSEEN_INTEGRATORS_O_HH

#include <dune/fem/oseen/assembler/base.hh>
#include <dune/stuff/common/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

	template < class MatrixObjectType, class Traits, class BetaFunctionType  >
	class O
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
		const BetaFunctionType& beta_;
		public:
			O( MatrixObjectType& matrix_object, const BetaFunctionType& beta)
				:matrix_object_(matrix_object),
				beta_(beta)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				applyVolume_alt2( info );
			}

			template < class InfoContainerVolumeType >
			void applyVolume_alt1( const InfoContainerVolumeType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
                            const auto x = info.volumeQuadratureElement.point( quad );
                            const double elementVolume = info.geometry.integrationElement( x );
                            const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							//calc u_h * \nabla * (v \tensor \beta )
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, x, v_i );
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							VelocityRangeType beta_eval;
							beta_lf.evaluate( x, beta_eval );

							VelocityJacobianRangeType v_i_jacobian;
							info.velocity_basefunction_set_element.jacobian( i, x, v_i_jacobian );
							VelocityJacobianRangeType v_j_jacobian;
							info.velocity_basefunction_set_element.jacobian( j, x, v_j_jacobian );

							VelocityJacobianRangeType v_i_tensor_beta
									= DSC::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, beta_eval );
							const double ret = DSC::colonProduct( v_i_tensor_beta, v_j_jacobian );

							O_i_j -= elementVolume
								* integrationWeight
								* info.convection_scaling
								* ret;
						}
						localOmatrixElement.add( i, j, O_i_j );
					}
				}
			}

			template < class InfoContainerVolumeType >
			void applyVolume_alt2( const InfoContainerVolumeType& info )
			{
                LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
                        for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
                            const auto x = info.volumeQuadratureElement.point( quad );
                            const double elementVolume = info.geometry.integrationElement( x );
                            const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							//calc u_h * \nabla * (v \tensor \beta )
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, x, v_i );
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							VelocityRangeType beta_eval;
							beta_lf.evaluate( x, beta_eval );

							VelocityJacobianRangeType v_j_jacobian;
							info.velocity_basefunction_set_element.jacobian( j, x, v_j_jacobian );
							VelocityJacobianRangeType beta_jacobian;
							beta_lf.jacobian( x, beta_jacobian );

							VelocityRangeType divergence_of_v_j_tensor_beta;

							divergence_of_v_j_tensor_beta[0] =
								beta_eval[0] * v_j_jacobian[0][0] + v_j[0] * beta_jacobian[0][0]
								   + beta_eval[0] * v_j_jacobian[1][1] + v_j[1] * beta_jacobian[0][1];

							divergence_of_v_j_tensor_beta[1] =
								beta_eval[1] * v_j_jacobian[0][0] + v_j[0] * beta_jacobian[1][0]
								   + beta_eval[1] * v_j_jacobian[1][1] + v_j[1] * beta_jacobian[1][1];


							const double u_h_times_divergence_of_beta_v_i_tensor_beta =
									v_j * divergence_of_v_j_tensor_beta;

							O_i_j -= elementVolume
									* integrationWeight
									* info.convection_scaling
									* u_h_times_divergence_of_beta_v_i_tensor_beta;
						}
						if ( fabs( O_i_j ) < info.eps ) {
							O_i_j = 0.0;
						}
						else {
							localOmatrixElement.add( i, j, O_i_j );
						}
					}
				}
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localOmatrixNeighbour( matrix_object_, info.neighbour, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType&
						beta_lf = beta_.localFunction( info.entity );

                // we call this one
				// (O)_{i,j} += \int_{ // O's element surface integral
				//           += \int_{ // O's neighbour surface integral
                // see also "O's boundary integral" below
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad )
				{
                    const auto xInside = info.faceQuadratureElement.point( quad );
                    const auto xOutside = info.faceQuadratureNeighbour.point( quad );
                    const auto xLocal = info.faceQuadratureElement.localPoint( quad );
                    const auto xLocal_neigh = info.faceQuadratureNeighbour.localPoint( quad );

					VelocityRangeType beta_eval;
					beta_lf.evaluate( xInside, beta_eval );
					VelocityRangeType beta_eval_neigh;
					beta_lf.evaluate( xOutside, beta_eval_neigh );
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					const VelocityRangeType outerNormal_neigh = info.intersection.unitOuterNormal( xLocal_neigh );
					const double beta_times_normal = beta_eval * outerNormal;
					if ( beta_times_normal > 0  )
					{
                        const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                        const double integrationWeight = info.faceQuadratureElement.weight( quad );

						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, xInside, v_i );

							for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j ) {
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								VelocityRangeType v_j_neigh( 0.0 );
								info.velocity_basefunction_set_neighbour.evaluate( j, xOutside, v_j_neigh );
								// \int_{dK} \beta * n * u_h * v ds

                                const double v_i_jump = ( (v_i * outerNormal) );
								double ret = (beta_times_normal ) * ( (v_i * v_j)*0.5 + v_i_jump );

								const double O_i_j = elementVolume
										* integrationWeight
										* info.convection_scaling
										* ret;
								localOmatrixElement.add( i, j, O_i_j );
                            }
                        }
					}
                    // compute O's neighbour surface integral
					else
					{

                        const double elementVolume = info.intersectionGeometry.integrationElement( xLocal_neigh );
                        const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
						// \int_{dK} flux_value : ( v_j \ctimes n ) ds
						for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i )
						{
                            VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_i );
							for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j )
							{
                                VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								VelocityRangeType v_j_neigh( 0.0 );
								info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_j_neigh );
								// \int_{dK} \beta * n * u_h * v ds
                                const double v_i_jump = ( (v_i * outerNormal_neigh));
								double ret = (beta_eval*outerNormal ) * ( (v_i * v_j)*0.5 - v_i_jump );

								const double O_i_j = elementVolume
										* integrationWeight
										* info.convection_scaling
										* ret;
								localOmatrixNeighbour.add( i, j, O_i_j );
                            }
						}
                    }
                }
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType&
						beta_lf = beta_.localFunction( info.entity );
				// (O)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}} STUFF n_{t}ds											// O's boundary integral
				//                                                                                                           // see also "O's element surface integral" and "Y's neighbour surface integral" above
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad )
				{
                    const auto x = info.faceQuadratureElement.point( quad );
                    const auto xLocal = info.faceQuadratureElement.localPoint( quad );
                    const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                    const double integrationWeight = info.faceQuadratureElement.weight( quad );
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					VelocityRangeType beta_eval;
					beta_lf.evaluate( x, beta_eval );
					const double beta_times_normal = beta_eval * outerNormal;

					if ( beta_times_normal < 0 )
						continue;

					for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ); ++i )
					{
						VelocityRangeType v_i( 0.0 );
						info.velocity_basefunction_set_element.evaluate( i, x, v_i );

						for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j )
						{
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );

							const double ret  = (beta_times_normal) * (v_i * v_j);
							//inner edge (self)
							const double O_i_j = elementVolume
								* integrationWeight
								* info.convection_scaling
								* ret;
							localOmatrixElement.add( i, j, O_i_j );
                        }
					}
                }
			}
			static const std::string name;
	};

	template < class T, class R, class F > const std::string O<T,R,F>::name = "O";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_O_HH

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

