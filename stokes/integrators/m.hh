#ifndef DUNE_STOKES_INTEGRATORS_M_HH
#define DUNE_STOKES_INTEGRATORS_M_HH

#include <dune/stokes/integrators/base.hh>
#include <dune/stuff/matrix.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

    template < class MatrixPointerType, class Traits >
	class M
	{
        typedef typename MatrixPointerType::element_type
            MatrixObjectType;
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
		public:
            M( MatrixPointerType& matrix_object	)
				:matrix_pointer_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				//using the proxy here would be potentially fatal because of the inversion
                LocalMatrixProxyType local_matrix ( matrix_pointer_, info.entity, info.entity, info.eps );
                ASSERT_EQ( local_matrix.rows(), info.numSigmaBaseFunctionsElement );
                ASSERT_EQ( local_matrix.cols(), info.numSigmaBaseFunctionsElement );

				// (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
				for ( int i = 0; i < info.numSigmaBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						double M_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = info.geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute \tau_{i}:\tau_{j}
							SigmaRangeType tau_i( 0.0 );
							SigmaRangeType tau_j( 0.0 );
							info.sigma_basefunction_set_element.evaluate( i, x, tau_i );
							info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
							const double tau_i_times_tau_j = Stuff::colonProduct( tau_i, tau_j );
							// compute M_i_j
							M_i_j += elementVolume
								* integrationWeight
								* tau_i_times_tau_j;
						} // done sum over quadrature points
						// if small, should be zero
						if ( fabs( M_i_j ) < info.eps ) {
							M_i_j = 0.0;
						} // else invert
						else {
							M_i_j = 1.0 / M_i_j;
							// add to matrix
                            local_matrix.add( i, j, M_i_j );
						}
					}
				} // done computing Minvs' volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& )
			{}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& )
			{}
			static const std::string name;
	};

	template < class T, class R > const std::string M<T,R>::name = "M";

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_M_HH
