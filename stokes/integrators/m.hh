#ifndef DUNE_STOKES_INTEGRATORS_HH
#define DUNE_STOKES_INTEGRATORS_HH

#include <dune/stokes/integrators/base.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {

	template < class MatrixObjectType, class Traits >
	class Dummy
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
			Dummy( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
				// (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
				for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
						double M_i_j = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = volumeQuadratureElement.weight( quad );
							// compute \tau_{i}:\tau_{j}
							SigmaRangeType tau_i( 0.0 );
							SigmaRangeType tau_j( 0.0 );
							sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
							sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
							const double tau_i_times_tau_j = Stuff::colonProduct( tau_i, tau_j );
							// compute M_i_j
							M_i_j += elementVolume
								* integrationWeight
								* tau_i_times_tau_j;
						} // done sum over quadrature points
						// if small, should be zero
						if ( fabs( M_i_j ) < eps ) {
							M_i_j = 0.0;
						} // else invert
						else {
							M_i_j = 1.0 / M_i_j;
							// add to matrix
							localMInversMatrixElement.add( i, j, M_i_j );
						}
					}
				} // done computing Minvs' volume integral
			}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& info )
			{}
	};

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_HH
