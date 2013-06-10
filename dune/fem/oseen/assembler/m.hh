#ifndef DUNE_OSEEN_INTEGRATORS_M_HH
#define DUNE_OSEEN_INTEGRATORS_M_HH

#include <dune/fem/oseen/assembler/base.hh>
#include <dune/stuff/common/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

    template < class MatrixObjectType, class Traits >
	class M
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
            M( MatrixObjectType& matrix_object	)
                :matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& info )
			{
                LocalMatrixProxyType local_matrix ( matrix_object_, info.entity, info.entity, info.eps );

				// (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
				for ( int i = 0; i < info.numSigmaBaseFunctionsElement; ++i ) {
					for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
						double M_i_j = 0.0;
                        for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++quad ) {
                            const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
                            const double elementVolume = info.geometry.integrationElement( x );
                            const double integrationWeight = info.volumeQuadratureElement.weight( quad );
							// compute \tau_{i}:\tau_{j}
							SigmaRangeType tau_i( 0.0 );
							SigmaRangeType tau_j( 0.0 );
							info.sigma_basefunction_set_element.evaluate( i, x, tau_i );
							info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
							const double tau_i_times_tau_j = DSC::colonProduct( tau_i, tau_j );
							// compute M_i_j
							M_i_j += elementVolume
								* integrationWeight
								* tau_i_times_tau_j;
                        }

                        if ( info.eps < fabs( M_i_j ) ) {
                            local_matrix.add( i, j, 1.0f / M_i_j );
                        }
					}
                }
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

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_M_HH

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

