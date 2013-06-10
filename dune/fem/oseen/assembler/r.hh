#ifndef DUNE_OSEEN_INTEGRATORS_R_HH
#define DUNE_OSEEN_INTEGRATORS_R_HH

#include <dune/fem/oseen/assembler/base.hh>
#include <dune/stuff/common/matrix.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

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
		typedef DSFe::LocalMatrixProxy<MatrixObjectType>
			LocalMatrixProxyType;

		MatrixObjectType& matrix_object_;
		public:
			R( MatrixObjectType& matrix_object	)
				:matrix_object_(matrix_object)
			{}

			template < class InfoContainerVolumeType >
			void applyVolume( const InfoContainerVolumeType& /*info*/ )
			{}

			template < class InfoContainerInteriorFaceType >
			void applyInteriorFace( const InfoContainerInteriorFaceType& info )
			{
				LocalMatrixProxyType localRmatrixElement( matrix_object_, info.entity, info.entity, info.eps );
				LocalMatrixProxyType localRmatrixNeighbour( matrix_object_, info.entity, info.neighbour, info.eps );
				// (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}ds // R's element surface integral
				//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}ds // R's neighbour surface integral
				//                                                                                                // see also "R's boundary integral" below
//                        if ( info.discrete_model.hasVelocityPressureFlux() ) {
					for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
                        for ( int i = 0; i < info.numPressureBaseFunctionsElement; ++i ) {
							double R_i_j = 0.0;
                            for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
                                const auto x = info.faceQuadratureElement.point( quad );
								const auto xLocal = info.faceQuadratureElement.localPoint( quad );
                                const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                                const double integrationWeight = info.faceQuadratureElement.weight( quad );
								// compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_element.evaluate( i, x, q_i );
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_element.evaluate( j, x, q_j );
								const double q_i_times_q_j = q_i * q_j;
								R_i_j += info.D_11
									* elementVolume
									* integrationWeight
									* q_i_times_q_j;
                            }
							localRmatrixElement.add( i, j, R_i_j );
                        }

						for ( int i = 0; i < info.numPressureBaseFunctionsNeighbour; ++i ) {
							double R_i_j = 0.0;
							for ( size_t quad = 0; quad < info.faceQuadratureNeighbour.nop(); ++quad ) {
                                const auto xInside = info.faceQuadratureElement.point( quad );
								const auto xOutside = info.faceQuadratureNeighbour.point( quad );
								const auto xLocal = info.faceQuadratureNeighbour.localPoint( quad );
                                const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                                const double integrationWeight = info.faceQuadratureNeighbour.weight( quad );
								// compute \hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}
								PressureRangeType q_j( 0.0 );
								info.pressure_basefunction_set_neighbour.evaluate( j, xOutside, q_j );
								PressureRangeType q_i( 0.0 );
								info.pressure_basefunction_set_element.evaluate( i, xInside, q_i );
								const double q_i_times_q_j = q_i * q_j;
								R_i_j += -1.0
									* info.D_11
									* elementVolume
									* integrationWeight
									* q_i_times_q_j;
                            }
							localRmatrixNeighbour.add( i, j, R_i_j );
                        }
                    }
//                        }
			}

			template < class InfoContainerFaceType >
			void applyBoundaryFace( const InfoContainerFaceType& )
			{}

			static const std::string name;
	};

	template < class T, class Y > const std::string R<T,Y>::name = "R";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_R_HH

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

