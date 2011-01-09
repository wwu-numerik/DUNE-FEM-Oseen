#ifndef DUNE_STOKES_INTEGRATORS_O_HH
#define DUNE_STOKES_INTEGRATORS_O_HH

#include <dune/stokes/integrators/base.hh>
#include <dune/stuff/matrix.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

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
		typedef Stuff::Matrix::LocalMatrixProxy<MatrixObjectType>
			LocalMatrixProxyType;

		MatrixObjectType& matrix_object_;
		const BetaFunctionType& beta_;
		public:
			O( MatrixObjectType& matrix_object, const BetaFunctionType& beta)
				:matrix_object_(matrix_object),
				beta_(beta)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume_alt1( const InfoContainerVolumeType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
					// get x
					const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
					// get the integration factor
					const double elementVolume = info.geometry.integrationElement( x );
					// get the quadrature weight
					const double integrationWeight = info.volumeQuadratureElement.weight( quad );
					VelocityJacobianRangeType beta_jacobian;
					beta_lf.jacobian( x, beta_jacobian );

					VelocityRangeType beta_eval;
					beta_lf.evaluate( x, beta_eval );

					for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
						//calc u_h * \nabla * (v \tensor \beta )
						VelocityRangeType v_i( 0.0 );
						info.velocity_basefunction_set_element.evaluate( i, x, v_i );

						VelocityJacobianRangeType v_i_jacobian;
						info.velocity_basefunction_set_element.jacobian( i, x, v_i_jacobian );

						for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							VelocityJacobianRangeType v_j_jacobian;
							info.velocity_basefunction_set_element.jacobian( j, x, v_j_jacobian );
							VelocityJacobianRangeType v_i_tensor_beta
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, beta_eval );
							const double ret = Stuff::colonProduct( v_i_tensor_beta, v_j_jacobian );

							const double O_i_j = -elementVolume
								* integrationWeight
								* info.convection_scaling
								* ret;
							localOmatrixElement.add( i, j, O_i_j );
						}
					}
				}
			}

			template < class InfoContainerVolumeType >
			void applyVolume_alt2( const InfoContainerVolumeType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
					// get x
					const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
					// get the integration factor
					const double elementVolume = info.geometry.integrationElement( x );
					// get the quadrature weight
					const double integrationWeight = info.volumeQuadratureElement.weight( quad );
					VelocityJacobianRangeType beta_jacobian;
					beta_lf.jacobian( x, beta_jacobian );

					VelocityRangeType beta_eval;
					beta_lf.evaluate( x, beta_eval );

					for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
						//calc u_h * \nabla * (v \tensor \beta )
						VelocityRangeType v_i( 0.0 );
						info.velocity_basefunction_set_element.evaluate( i, x, v_i );

						VelocityJacobianRangeType v_i_jacobian;
						info.velocity_basefunction_set_element.jacobian( i, x, v_i_jacobian );

						VelocityRangeType divergence_of_v_i_tensor_beta;
						for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
							double row_result = 0;
							for ( size_t m = 0; m < beta_eval.dim(); ++m ) {
								row_result += beta_jacobian[l][m] * v_i[l] + v_i_jacobian[l][m] * beta_eval[l];
							}
							divergence_of_v_i_tensor_beta[l] = row_result;
						}
						for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
							assert( !isnan(divergence_of_v_i_tensor_beta[l]) );
						}

					   divergence_of_v_i_tensor_beta[0] = beta_eval[0] * v_i_jacobian[0][0]
							   + v_i[0] * beta_jacobian[0][0]
							   + beta_eval[0] * v_i_jacobian[1][1]
							   + v_i[1] * beta_jacobian[0][1];
					   divergence_of_v_i_tensor_beta[1] = beta_eval[1] * v_i_jacobian[0][0]
							   + v_i[0] * beta_jacobian[1][0]
							   + beta_eval[1] * v_i_jacobian[1][1]
							   + v_i[1] * beta_jacobian[1][1];

					   for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							VelocityJacobianRangeType v_j_jacobian;
							info.velocity_basefunction_set_element.jacobian( j, x, v_j_jacobian );

							const double u_h_times_divergence_of_beta_v_j_tensor_beta =
									v_j * divergence_of_v_i_tensor_beta;

							const double O_i_j_d = -elementVolume
									* integrationWeight
									* info.convection_scaling
									* u_h_times_divergence_of_beta_v_j_tensor_beta;

							localOmatrixElement.add( i, j, O_i_j_d );
						}
					}
				}
			}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				applyVolume_alt1( info );
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localOmatrixNeighbour( matrix_object_, info.neighbour, info.entity, info.eps );
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				//                                                                                                         // we call this one
				// (O)_{i,j} += \int_{ // O's element surface integral
				//           += \int_{ // O's neighbour surface integral
				//                                                                                                         // see also "O's boundary integral" below
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
						const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
						const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
					{
						const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
						const VelocityRangeType xWorld = info.geometry.global( xInside );
						const VelocityRangeType xWorld_Outside = info.geometry.global( xOutside );

						// get the integration factor
						const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = info.faceQuadratureElement.weight( quad );
						const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
						VelocityRangeType beta_eval;
						beta_lf.evaluate( xInside, beta_eval );
						const double beta_times_normal = beta_eval * outerNormal;
						const double c_s = (beta_times_normal) * 0.5;
						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {

							// compute O's element surface integral
							VelocityRangeType v_j( 0.0 );
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, xInside, v_i );

							VelocityRangeType u_h = v_i;
							VelocityJacobianRangeType mean_value
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( u_h, beta_eval );
							mean_value *= 0.5;
							VelocityJacobianRangeType u_jump
								= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, outerNormal );
							u_jump *= c_s;

							VelocityJacobianRangeType flux_value = mean_value;
							flux_value += u_jump;

							for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j ) {
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								//calc u^c_h \tensor beta * v \tensor n (self part), the flux value
								double c_s = (beta_times_normal) * 0.5;
								VelocityJacobianRangeType v_j_tensor_n
										= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_j, outerNormal );

									// \int_{dK} flux_value : ( v_j \ctimes n ) ds
								const double ret  = Stuff::colonProduct( flux_value, v_j_tensor_n );

								const double O_i_j = elementVolume
										* integrationWeight
										* info.convection_scaling
										* ret;
								localOmatrixElement.add( i, j, O_i_j );
							} // done sum over all quadrature points
						} // done computing Y's element surface integral
					}
					// compute O's neighbour surface integral
					const LocalIntersectionCoordinateType xLocal = info.faceQuadratureNeighbour.localPoint( quad );
					// get the integration factor
					const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
					// get the quadrature weight
					const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					VelocityRangeType v_i( 0.0 );
					VelocityRangeType v_j( 0.0 );
					VelocityRangeType beta_eval;
					for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i ) {
						info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_i );

						beta_lf.evaluate( xInside, beta_eval );
						const double beta_times_normal =  ( beta_eval * outerNormal );

						//calc u^c_h \tensor beta * v \tensor n (self part), the flux value
						double c_s = (beta_times_normal) * 0.5;
						VelocityRangeType u_h = v_i;
						VelocityJacobianRangeType mean_value
								= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( u_h, beta_eval );
						mean_value *= 0.5;
						VelocityJacobianRangeType u_jump
								= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, outerNormal );
						u_jump *= c_s;
						VelocityJacobianRangeType flux_value = mean_value;
						flux_value += u_jump;

						for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j ) {
							// \int_{dK} flux_value : ( v_j \ctimes n ) ds
							VelocityJacobianRangeType v_i_tensor_n
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_j, outerNormal );
							const double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );
							info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
							const double O_i_j = elementVolume
									* integrationWeight
									* info.convection_scaling
									* ret;
							localOmatrixNeighbour.add( i, j, O_i_j );
						} // done sum over all quadrature points
					} // done computing Y's neighbour surface integral
				} // done computing Y's surface integrals
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
				LocalMatrixProxyType localOmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				// (O)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}} STUFF n_{t}ds											// O's boundary integral
				//                                                                                                           // see also "O's element surface integral" and "Y's neighbour surface integral" above
				for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ); ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
							// get x codim<0> and codim<1> coordinates
							const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
							const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
							// get the integration factor
							const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
							// get the quadrature weight
							const double integrationWeight = info.faceQuadratureElement.weight( quad );
							const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
							//calc u^c_h \tensor beta * v \tensor n
							VelocityRangeType v_j( 0.0 );
							info.velocity_basefunction_set_element.evaluate( j, x, v_j );
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, x, v_i );
							VelocityRangeType beta_eval;
							beta_.localFunction( info.entity ).evaluate( x, beta_eval );
							const double beta_times_normal = beta_eval * outerNormal;

							double c_s;
							if ( beta_times_normal < 0 ) {
								c_s = beta_times_normal * 0.5;
							}
							else {
								c_s = - beta_times_normal * 0.5;
							}

							VelocityJacobianRangeType mean_value
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, beta_eval );
							mean_value *= 0.5;

							VelocityJacobianRangeType u_jump
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, outerNormal );
							u_jump *= c_s;

							VelocityJacobianRangeType flux_value = mean_value;
							flux_value += u_jump;

							VelocityJacobianRangeType v_i_tensor_n
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_j, outerNormal );
							double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );
							//inner edge (self)
							O_i_j += elementVolume
								* integrationWeight
								* info.convection_scaling
								* ret;

						} // done sum over all quadrature points
						localOmatrixElement.add( i, j, O_i_j );
					}
				} // done computing O's boundary integral
			}
			static const std::string name;
	};

	template < class T, class R, class F > const std::string O<T,R,F>::name = "O";

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_O_HH