#ifndef DUNE_OSEEN_INTEGRATORS_RHS_HH
#define DUNE_OSEEN_INTEGRATORS_RHS_HH

#include <dune/fem/oseen/assembler/base.hh>

namespace Dune {
namespace Oseen {
namespace Assembler {

template < class DiscreteFunctionType, class Traits >
class H1
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


	DiscreteFunctionType& discrete_function_;
	public:
		H1( DiscreteFunctionType& df_func )
			:discrete_function_(df_func)
		{}

		template < class InfoContainerVolumeType >
		void applyVolume( const InfoContainerVolumeType& )
		{}

		template < class InfoContainerInteriorFaceType >
		void applyInteriorFace( const InfoContainerInteriorFaceType& )
		{}

		template < class InfoContainerFaceType >
		void applyBoundaryFace( const InfoContainerFaceType& info )
		{
			typename DiscreteFunctionType::LocalFunctionType
					localH1rhs = discrete_function_.localFunction( info.entity );
			//                                                                                                    // we will call this one
			// (H1)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}ds // H1's boundary integral
			if ( info.discrete_model.hasVelocitySigmaFlux() ) {
				for ( int j = 0; j < info.numSigmaBaseFunctionsElement; ++j ) {
					double H1_j = 0.0;
					// sum over all quadrature points
					for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
						// get x codim<0> and codim<1> coordinates
						const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
						const VelocityRangeType xWorld = info.geometry.global( x );
						const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
						// get the integration factor
						const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = info.faceQuadratureElement.weight( quad );
						// compute \hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}
						const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
						SigmaRangeType tau_j( 0.0 );
						info.sigma_basefunction_set_element.evaluate( j, x, tau_j );
						VelocityRangeType tau_j_times_normal( 0.0 );
						tau_j.mv( outerNormal, tau_j_times_normal );
						VelocityRangeType gD( 0.0 );
						info.discrete_model.dirichletData( info.intersection, 0.0, xWorld,  gD );
						const double gD_times_tau_j_times_normal = gD * tau_j_times_normal;
						H1_j += elementVolume
							* integrationWeight
							* info.viscosity
							* gD_times_tau_j_times_normal;
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H1_j ) < info.eps ) {
						H1_j = 0.0;
					}
					else
						// add to rhs
						localH1rhs[ j ] += H1_j;
				} // done computing H1's boundary integral
			}
		}
		static const std::string name;
};

template < class T, class R > const std::string H1<T,R>::name = "H1";

template < class DiscreteFunctionType, class Traits >
class H2
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


	DiscreteFunctionType& discrete_function_;
	public:
		H2( DiscreteFunctionType& df_func )
			:discrete_function_(df_func)
		{}

		template < class InfoContainerVolumeType >
		void applyVolume( const InfoContainerVolumeType&  info )
		{
			typename DiscreteFunctionType::LocalFunctionType
					localH2rhs = discrete_function_.localFunction( info.entity );
			//                                    // we will call this one
			// (H2)_{j} += \int_{T}f\cdot v_{j}dx // H2's volume integral
			//                                    // see also "H2's boundary integral" further down
			if ( info.discrete_model.hasForce() ) {
				for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
					double H2_j = 0.0;
					// sum over all quadratur points
					for ( size_t quad = 0; quad < info.volumeQuadratureElement.nop(); ++ quad ) {
						// get x
						const ElementCoordinateType x = info.volumeQuadratureElement.point( quad );
						const VelocityRangeType xWorld = info.geometry.global( x );
						// get the integration factor
						const double elementVolume = info.geometry.integrationElement( x );
						// get the quadrature weight
						const double integrationWeight = info.volumeQuadratureElement.weight( quad );
						// compute f\cdot v_j
						VelocityRangeType v_j( 0.0 );
						info.velocity_basefunction_set_element.evaluate( j, x, v_j );
						VelocityRangeType f( 0.0 );
#if MODEL_PROVIDES_LOCALFUNCTION
									info.discrete_model.forceF().localFunction(info.entity).evaluate( x, f );
#else
									info.discrete_model.force( 0.0, xWorld, f );
#endif
						const double f_times_v_j = f * v_j;
						H2_j += elementVolume
							* integrationWeight
							* f_times_v_j;
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H2_j ) < info.eps ) {
						H2_j = 0.0;
					}
					else
						// add to rhs
						localH2rhs[ j ] += H2_j;
				} // done computing H2's volume integral
			}
		}

		template < class InfoContainerInteriorFaceType >
		void applyInteriorFace( const InfoContainerInteriorFaceType&  )
		{}

		template < class InfoContainerFaceType >
		void applyBoundaryFace( const InfoContainerFaceType& info )
		{
			typename DiscreteFunctionType::LocalFunctionType
					localH2rhs = discrete_function_.localFunction( info.entity );
			//                                                                                                                 // we will call this one
			// (H2)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left( \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}ds         // H2's 1st boundary integral
			//                                                         -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}ds        \right) // H2's 2nd boundary integral
			//                                                                                                                 // see also "H2's volume integral" above
			if ( info.discrete_model.hasSigmaFlux() && info.discrete_model.hasPressureFlux() ) {
				for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
					double H2_j = 0.0;
					// sum over all quadrature points
					for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
						// get x codim<0> and codim<1> coordinates
						const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
						const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
//						const VelocityRangeType globalX = info.geometry.global( x );
						// get the integration factor
						const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = info.faceQuadratureElement.weight( quad );
						// prepare
						const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
						VelocityRangeType v_j( 0.0 );
						info.velocity_basefunction_set_element.evaluate( j, x, v_j );
						// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
//                                    if ( info.discrete_model.hasSigmaFlux() ) {
                            const VelocityRangeType xIntersectionGlobal = info.intersection.geometryInInside().global( xLocal );
							const VelocityRangeType xWorld = info.geometry.global( xIntersectionGlobal );
							VelocityRangeType gD( 0.0 );
							info.discrete_model.dirichletData( info.intersection, 0.0, xWorld, gD );
							SigmaRangeType gD_times_normal( 0.0 );
							gD_times_normal = DSC::dyadicProduct<SigmaRangeType,VelocityRangeType>( gD, outerNormal );
							VelocityRangeType gD_times_normal_times_normal( 0.0 );
							gD_times_normal.mv( outerNormal, gD_times_normal_times_normal );
							const double v_j_times_gD_times_normal_times_normal= v_j * gD_times_normal_times_normal;
							H2_j += info.C_11
								* elementVolume
								* integrationWeight
								* v_j_times_gD_times_normal_times_normal;
//                                    }
						// done computing H2's 1st boundary integral
						// compute -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}
//                                    if ( info.discrete_model.hasPressureFlux() ) {
							const double v_j_times_normal = v_j * outerNormal;
							const double flux_times_v_j_times_n_t = 0.0 * v_j_times_normal;
							H2_j += -1.0
								* elementVolume
								* integrationWeight
								* info.pressure_gradient_scaling
								* flux_times_v_j_times_n_t;
//                                    }
						// done computing H2's 2nd boundary integral
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H2_j ) < info.eps ) {
						H2_j = 0.0;
					}
					else
						// add to rhs
						localH2rhs[ j ] += H2_j;
				} // done computing H2's boundary integrals
			}
		}
		static const std::string name;
};

template < class T, class R > const std::string H2<T,R>::name = "H2";

template < class DiscreteFunctionType, class Traits, class BetaFunctionType >
class H2_O
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


	DiscreteFunctionType& discrete_function_;
	const BetaFunctionType& beta_;
	public:
		H2_O( DiscreteFunctionType& df_func, const BetaFunctionType& beta )
			:discrete_function_(df_func),
			beta_(beta)
		{}

		template < class InfoContainerVolumeType >
		void applyVolume( const InfoContainerVolumeType& )
		{
		}

		template < class InfoContainerInteriorFaceType >
		void applyInteriorFace( const InfoContainerInteriorFaceType& )
		{}

		template < class InfoContainerFaceType >
		void applyBoundaryFace( const InfoContainerFaceType& info )
		{
//			return;
			typename DiscreteFunctionType::LocalFunctionType
					localH2_O_rhs = discrete_function_.localFunction( info.entity );
			// (H2_O)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left(  \beta n_{T} g_D v_j ds        \right) // H2_O's boundary integral
			for ( int j = 0; j < info.numVelocityBaseFunctionsElement; ++j ) {
				double H2_O_j = 0.0;
				// sum over all quadrature points
				for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
					// get x codim<0> and codim<1> coordinates
					const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
					const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
									// get the integration factor
					const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
					// get the quadrature weight
					const double integrationWeight = info.faceQuadratureElement.weight( quad );
					// prepare
					const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
					VelocityRangeType v_j( 0.0 );
					info.velocity_basefunction_set_element.evaluate( j, x, v_j );
					// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
                    const VelocityRangeType xIntersectionGlobal = info.intersection.geometryInInside().global( xLocal );
					const VelocityRangeType xWorld = info.geometry.global( xIntersectionGlobal );
					VelocityRangeType gD( 0.0 );
					info.discrete_model.dirichletData( info.intersection, 0.0, xWorld, gD );

					VelocityRangeType beta_eval;
					beta_.localFunction(info.entity).evaluate( x, beta_eval );
					const double beta_times_normal = beta_eval * outerNormal;

					if ( beta_times_normal < 0 )
					{
						const double ret = beta_times_normal * (gD * v_j) ;
						H2_O_j -= elementVolume
								* info.convection_scaling
								* integrationWeight
								* ret;
					}
				}
				if ( fabs( H2_O_j ) < info.eps ) {
						 H2_O_j = 0.0;
				}
				else {
					localH2_O_rhs[ j ] += H2_O_j;
				}
			}
		}
		static const std::string name;
};

template < class T, class R, class F > const std::string H2_O<T,R,F>::name = "H2_O";

template < class DiscreteFunctionType, class Traits >
class H3
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


	DiscreteFunctionType& discrete_function_;
	public:
		H3( DiscreteFunctionType& df_func )
			:discrete_function_(df_func)
		{}

		template < class InfoContainerVolumeType >
		void applyVolume( const InfoContainerVolumeType& )
		{
		}

		template < class InfoContainerInteriorFaceType >
		void applyInteriorFace( const InfoContainerInteriorFaceType&  )
		{}

		template < class InfoContainerFaceType >
		void applyBoundaryFace( const InfoContainerFaceType& info )
		{
			typename DiscreteFunctionType::LocalFunctionType
					localH3rhs = discrete_function_.localFunction( info.entity );
            for ( int j = 0; j < info.numPressureBaseFunctionsElement; ++j ) {
                double H3_j = 0.0;
                // sum over all quadrature points
                for ( size_t quad = 0; quad < info.faceQuadratureElement.nop(); ++quad ) {
                    // get x codim<0> and codim<1> coordinates
                    const ElementCoordinateType x = info.faceQuadratureElement.point( quad );
                    const LocalIntersectionCoordinateType xLocal = info.faceQuadratureElement.localPoint( quad );
                    const VelocityRangeType xWorld = info.geometry.global( x );
                    // get the integration factor
                    const double elementVolume = info.intersectionGeometry.integrationElement( xLocal );
                    // get the quadrature weight
                    const double integrationWeight = info.faceQuadratureElement.weight( quad );
                    // compute -\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}
                    const VelocityRangeType outerNormal = info.intersection.unitOuterNormal( xLocal );
                    VelocityRangeType gD( 0.0 );
                    info.discrete_model.dirichletData( info.intersection, 0.0, xWorld, gD );
                    const double gD_times_normal = gD * outerNormal;
                    PressureRangeType q_j( 0.0 );
                    info.pressure_basefunction_set_element.evaluate( j, x, q_j );
                    const double q_j_times_gD_times_normal = q_j * gD_times_normal;
                    H3_j += elementVolume
                        * integrationWeight
                        * q_j_times_gD_times_normal;
                } // done sum over all quadrature points
                // if small, should be zero
                if ( fabs( H3_j ) < info.eps ) {
                    H3_j = 0.0;
                }
                else
                    // add to rhs
                    localH3rhs[ j ] += H3_j;
            } // done computing H3's boundary integral
		}
		static const std::string name;
};

template < class T, class R > const std::string H3<T,R>::name = "H3";

} // end namespace Assembler
} // end namespace Oseen
} // end namespace Dune

#endif // DUNE_OSEEN_INTEGRATORS_RHS_HH

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

