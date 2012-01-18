#ifndef DUNE_OSEEN_INTEGRATORS_O_HH
#define DUNE_OSEEN_INTEGRATORS_O_HH

#include <dune/oseen/assembler/base.hh>
#include <dune/stuff/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

	template < class MatrixPointerType, class Traits, class BetaFunctionType  >
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
		typedef Stuff::Matrix::LocalMatrixProxy<MatrixPointerType>
			LocalMatrixProxyType;

		MatrixPointerType& matrix_pointer_;
		const BetaFunctionType& beta_;
		public:
			O( MatrixPointerType& matrix_object, const BetaFunctionType& beta)
				:matrix_pointer_(matrix_object),
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
				LocalMatrixProxyType localOmatrixElement( matrix_pointer_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
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
									= Stuff::dyadicProduct<VelocityJacobianRangeType,VelocityRangeType>( v_i, beta_eval );
							const double ret = Stuff::colonProduct( v_i_tensor_beta, v_j_jacobian );

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
//				return;
				typename MatrixPointerType::element_type::LocalMatrixType
						localOmatrixElement = matrix_pointer_->localMatrix( info.entity, info.entity );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
						beta_.localFunction( info.entity );
				for ( int i = 0; (i < info.numVelocityBaseFunctionsElement ) ; ++i ) {
					for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
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
							VelocityRangeType divergence_of_v_i_tensor_beta;

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
//				return;
				LocalMatrixProxyType localOmatrixElement( matrix_pointer_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localOmatrixNeighbour( matrix_pointer_, info.neighbour, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType&
						beta_lf = beta_.localFunction( info.entity );
//				const unsigned int inside_entity_id = beta_.space().gridPart().indexSet().index( info.entity );
//				const unsigned int outside_entity_id = beta_.space().gridPart().indexSet().index( info.neighbour );
//				if ( inside_entity_id > outside_entity_id )
//					return;
				//                                                                                                         // we call this one
				// (O)_{i,j} += \int_{ // O's element surface integral
				//           += \int_{ // O's neighbour surface integral
				//                                                                                                         // see also "O's boundary integral" below
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad )
				{
					const ElementCoordinateType xInside = info.faceQuadratureElement.point( quad );
					const ElementCoordinateType xOutside = info.faceQuadratureNeighbour.point( quad );
					const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
					const LocalIntersectionCoordinateType xLocal_neigh = info.faceQuadratureNeighbour.localPoint( quad );
					// get the integration factor

					VelocityRangeType beta_eval;
					beta_lf.evaluate( xInside, beta_eval );
					VelocityRangeType beta_eval_neigh;
					beta_lf.evaluate( xOutside, beta_eval_neigh );
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					const VelocityRangeType outerNormal_neigh = info.intersection.unitOuterNormal( xLocal_neigh );
					const double beta_times_normal = beta_eval * outerNormal;
//					const double c_star = std::abs(beta_times_normal) * 0.5;
					if ( beta_times_normal > 0  )
					{
						const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = info.faceQuadratureElement.weight( quad );

						for ( int i = 0; i < info.numVelocityBaseFunctionsElement; ++i ) {
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_element.evaluate( i, xInside, v_i );
//							VelocityRangeType u_h = v_i;

							for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j ) {
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								VelocityRangeType v_j_neigh( 0.0 );
								info.velocity_basefunction_set_neighbour.evaluate( j, xOutside, v_j_neigh );
								// \int_{dK} \beta * n * u_h * v ds

								const double v_i_jump = ( (v_i * outerNormal) );// + (v_j_neigh * outerNormal_neigh) );
//								double ret  = (beta_eval * u_h) * v_j_jump ;
								double ret = (beta_times_normal ) * ( (v_i * v_j)*0.5 + v_i_jump );
//								ret += (u_h * outerNormal) * v_j_jump * c_star;

								const double O_i_j = elementVolume
										* integrationWeight
										* info.convection_scaling
										* ret;
								localOmatrixElement.add( i, j, O_i_j );
							} // done sum over all quadrature points
						} // done computing Y's element surface integral
					}
//					continue;

						// compute O's neighbour surface integral
					else
					{

						const double elementVolume = info.intersectionGeometry.integrationElement( xLocal_neigh );
						// get the quadrature weight
						const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
						// \int_{dK} flux_value : ( v_j \ctimes n ) ds
						for ( int i = 0; i < info.numVelocityBaseFunctionsNeighbour; ++i )
						{
							// sum over all quadrature points
							VelocityRangeType v_i( 0.0 );
							info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_i );

//							VelocityRangeType u_h = v_i;
							for ( int j = 0; (j < info.numVelocityBaseFunctionsElement ); ++j )
							{
								// compute O's element surface integral
								VelocityRangeType v_j( 0.0 );
								info.velocity_basefunction_set_element.evaluate( j, xInside, v_j );
								VelocityRangeType v_j_neigh( 0.0 );
								info.velocity_basefunction_set_neighbour.evaluate( i, xOutside, v_j_neigh );
								// \int_{dK} \beta * n * u_h * v ds
								const double v_i_jump = ( (v_i * outerNormal_neigh));// + (v_j_neigh * outerNormal_neigh) );
//								double ret = (beta_eval_neigh * u_h) * v_j_jump;
								double ret = (beta_eval*outerNormal ) * ( (v_i * v_j)*0.5 - v_i_jump );
//								ret += (u_h * outerNormal_neigh) * v_j_jump * c_star;

								const double O_i_j = elementVolume
										* integrationWeight
										* info.convection_scaling
										* ret;
								localOmatrixNeighbour.add( i, j, O_i_j );
							} // done sum over all quadrature points
						}
					} // done computing Y's neighbour surface integral
				} // done computing Y's surface integrals
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{
//				return;
				LocalMatrixProxyType localOmatrixElement( matrix_pointer_, info.entity, info.entity, info.eps );
				const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType&
						beta_lf = beta_.localFunction( info.entity );
				// (O)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}} STUFF n_{t}ds											// O's boundary integral
				//                                                                                                           // see also "O's element surface integral" and "Y's neighbour surface integral" above
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad )
				{
					const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
					const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
					// get the integration factor
					const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
					// get the quadrature weight
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
						} // done sum over all quadrature points
					}
				} // done computing O's boundary integral
			}
			static const std::string name;
	};

	template < class T, class R, class F > const std::string O<T,R,F>::name = "O";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_O_HH
