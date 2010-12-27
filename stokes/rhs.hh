#ifndef DUNE_STOKES_INTEGRATORS_RHS_HH
#define DUNE_STOKES_INTEGRATORS_RHS_HH


namespace Dune {
namespace Stokes {
namespace Integrators {

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
		void applyVolume( const InfoContainerVolumeType& info )
		{}

		template < class InfoContainerInteriorFaceType >
		void applyInteriorFace( const InfoContainerInteriorFaceType& info )
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
};

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
			//                                    // we will call this one
			// (H2)_{j} += \int_{T}f\cdot v_{j}dx // H2's volume integral
			//                                    // see also "H2's boundary integral" further down
			if ( discreteModel_.hasForce() ) {
				for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
					double H2_j = 0.0;
					// sum over all quadratur points
					for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
						// get x
						const ElementCoordinateType x = volumeQuadratureElement.point( quad );
						const VelocityRangeType xWorld = geometry.global( x );
						// get the integration factor
						const double elementVolume = geometry.integrationElement( x );
						// get the quadrature weight
						const double integrationWeight = volumeQuadratureElement.weight( quad );
						// compute f\cdot v_j
						VelocityRangeType v_j( 0.0 );
						velocityBaseFunctionSetElement.evaluate( j, x, v_j );
						VelocityRangeType f( 0.0 );
#if MODEL_PROVIDES_LOCALFUNCTION
									discreteModel_.forceF().localFunction(entity).evaluate( x, f );
#else
									discreteModel_.force( 0.0, xWorld, f );
#endif
						const double f_times_v_j = f * v_j;
						H2_j += elementVolume
							* integrationWeight
							* f_times_v_j;
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H2_j ) < eps ) {
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
			if ( discreteModel_.hasSigmaFlux() && discreteModel_.hasPressureFlux() ) {
				for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
					double H2_j = 0.0;
					// sum over all quadrature points
					for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
						// get x codim<0> and codim<1> coordinates
						const ElementCoordinateType x = faceQuadratureElement.point( quad );
						const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
						const VelocityRangeType globalX = geometry.global( x );
						// get the integration factor
						const double elementVolume = intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = faceQuadratureElement.weight( quad );
						// prepare
						const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
						VelocityRangeType v_j( 0.0 );
						velocityBaseFunctionSetElement.evaluate( j, x, v_j );
						// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
//                                    if ( discreteModel_.hasSigmaFlux() ) {
							const VelocityRangeType xIntersectionGlobal = intersection.intersectionSelfLocal().global( xLocal );
							const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
							VelocityRangeType gD( 0.0 );
							discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
							SigmaRangeType gD_times_normal( 0.0 );
							gD_times_normal = dyadicProduct( gD, outerNormal );
							VelocityRangeType gD_times_normal_times_normal( 0.0 );
							gD_times_normal.mv( outerNormal, gD_times_normal_times_normal );
							const double v_j_times_gD_times_normal_times_normal= v_j * gD_times_normal_times_normal;
							H2_j += C_11
								* elementVolume
								* integrationWeight
								* v_j_times_gD_times_normal_times_normal;
//                                    }
						// done computing H2's 1st boundary integral
						// compute -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}
//                                    if ( discreteModel_.hasPressureFlux() ) {
							const double v_j_times_normal = v_j * outerNormal;
							const double flux_times_v_j_times_n_t = 0.0 * v_j_times_normal;
							H2_j += -1.0
								* elementVolume
								* integrationWeight
								* pressure_gradient_scaling
								* flux_times_v_j_times_n_t;
//                                    }
						// done computing H2's 2nd boundary integral
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H2_j ) < eps ) {
						H2_j = 0.0;
					}
					else
						// add to rhs
						localH2rhs[ j ] += H2_j;
				} // done computing H2's boundary integrals
			}
		}
}; //end H2

template < class DiscreteFunctionType, class Traits >
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
	public:
		H2_O( DiscreteFunctionType& df_func )
			:discrete_function_(df_func)
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
			typename DiscreteFunctionType::LocalFunctionType
					localH2_O_rhs = discrete_function_.localFunction( info.entity );
			// (H2_O)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left(  \beta n_{T} g_D v_j ds        \right) // H2_O's boundary integral
			for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
				double H2_O_j = 0.0;
				// sum over all quadrature points
				for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
					// get x codim<0> and codim<1> coordinates
					const ElementCoordinateType x = faceQuadratureElement.point( quad );
					const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
									// get the integration factor
					const double elementVolume = intersectionGeometry.integrationElement( xLocal );
					// get the quadrature weight
					const double integrationWeight = faceQuadratureElement.weight( quad );
					// prepare
					const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
					VelocityRangeType v_j( 0.0 );
					velocityBaseFunctionSetElement.evaluate( j, x, v_j );
					// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
					const VelocityRangeType xIntersectionGlobal = intersection.intersectionSelfLocal().global( xLocal );
					const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
					VelocityRangeType gD( 0.0 );
					discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );

					VelocityRangeType beta_eval;
					beta_.localFunction(entity).evaluate( x, beta_eval );
					const double beta_times_normal = beta_eval * outerNormal;

					// u^c = 0.5 gD \otimes beta + Cs -gD \otimes n
					VelocityJacobianRangeType gD_tensor_beta = dyadicProduct( gD, beta_eval );
					gD_tensor_beta *= 0.5;
					double c_s;
					if ( beta_times_normal < 0 ) {
						c_s = beta_times_normal * 0.5;
					}
					else {
						c_s = - beta_times_normal * 0.5;
					}
					VelocityJacobianRangeType jump = dyadicProduct( gD, outerNormal );
					jump *= c_s;

					VelocityJacobianRangeType flux_value = gD_tensor_beta;
					flux_value += jump;

					VelocityJacobianRangeType v_j_tensor_n = dyadicProduct( v_j,  outerNormal );
					const double ret = Stuff::colonProduct( flux_value, v_j_tensor_n );
					H2_O_j += elementVolume
							* convection_scaling
							* integrationWeight
							* ret;
				}
				if ( fabs( H2_O_j ) < eps ) {
						 H2_O_j = 0.0;
				}
				else {
					localH2_O_rhs[ j ] += H2_O_j;
				}
			}
		}
}; //end H2_O

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

			// (H3)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}-\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}ds // H3's boundary integral
//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
				for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
					double H3_j = 0.0;
					// sum over all quadrature points
					for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
						// get x codim<0> and codim<1> coordinates
						const ElementCoordinateType x = faceQuadratureElement.point( quad );
						const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
						const VelocityRangeType xWorld = geometry.global( x );
						// get the integration factor
						const double elementVolume = intersectionGeometry.integrationElement( xLocal );
						// get the quadrature weight
						const double integrationWeight = faceQuadratureElement.weight( quad );
						// compute -\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}
						const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
						VelocityRangeType gD( 0.0 );
						discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
						const double gD_times_normal = gD * outerNormal;
						PressureRangeType q_j( 0.0 );
						pressureBaseFunctionSetElement.evaluate( j, x, q_j );
						const double q_j_times_gD_times_normal = q_j * gD_times_normal;
						H3_j += elementVolume
							* integrationWeight
							* q_j_times_gD_times_normal;
					} // done sum over all quadrature points
					// if small, should be zero
					if ( fabs( H3_j ) < eps ) {
						H3_j = 0.0;
					}
					else
						// add to rhs
						localH3rhs[ j ] += H3_j;
				} // done computing H3's boundary integral
//                        }
		}
}; //end H3

} // end namespace Integrators
} // end namespace Stokes
} // end namespace Dune

#endif // DUNE_STOKES_INTEGRATORS_RHS_HH
