#ifndef DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH
#define DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH

#include <dune/stokes/solver/solver_interface.hh>
#include <dune/stokes/solver/schurkomplement.hh>
namespace Dune {

	/**
		\brief Saddlepoint Solver
		The inner CG iteration is implemented by passing a custom Matrix operator to a given\n
		Dune solver. The outer iteration is a implementation of the CG algorithm as described in\n
		Kuhnibert
		Optionally the BFG scheme as described in YADDA is uesd to control the inner solver tolerance.
		/todo get references in doxygen right
	**/
	template < class StokesPassImp >
	class BiCgStabSaddlepointInverseOperator
	{
	  private:

		typedef StokesPassImp
			StokesPassType;

		typedef typename StokesPassType::Traits::DiscreteStokesFunctionWrapperType
			DiscreteStokesFunctionWrapperType;

		typedef typename StokesPassType::DomainType
			DomainType;

		typedef typename StokesPassType::RangeType
			RangeType;

		typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
			PressureDiscreteFunctionType;
		typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
			VelocityDiscreteFunctionType;


	  public:
		BiCgStabSaddlepointInverseOperator()
		{}

		/** takes raw matrices and right hand sides from pass as input, executes nested cg algorithm and outputs solution
		*/
		template <  class X_MatrixType,
					class M_invers_matrixType,
					class Y_MatrixType,
					class O_MatrixType,
					class E_MatrixType,
					class R_MatrixType,
					class Z_MatrixType,
					class W_MatrixType,
					class DiscreteSigmaFunctionType,
					class DiscreteVelocityFunctionType,
					class DiscretePressureFunctionType  >
		SaddlepointInverseOperatorInfo solve( const DomainType& /*arg*/,
					RangeType& dest,
					X_MatrixType& x_mat,
					M_invers_matrixType& m_inv_mat,
					Y_MatrixType& y_mat,
					O_MatrixType& o_mat,
					E_MatrixType& e_mat,
					R_MatrixType& r_mat,
					Z_MatrixType& z_mat,
					W_MatrixType& w_mat,
					const DiscreteSigmaFunctionType& rhs1_orig,
					const DiscreteVelocityFunctionType& rhs2_orig,
					const DiscretePressureFunctionType& rhs3 ) const
		{

			Logging::LogStream& logDebug = Logger().Dbg();
			Logging::LogStream& logInfo = Logger().Info();

			if ( Parameters().getParam( "disableSolver", false ) ) {
				logInfo.Resume();
				logInfo << "solving disabled via parameter file" << std::endl;
				return SaddlepointInverseOperatorInfo();
			}
			static const int save_interval = Parameters().getParam( "save_interval", -1 );
			static const std::string path = Parameters().getParam("fem.io.datadir", std::string("data") ) + std::string("/intermediate/");
			Stuff::testCreateDirectory( path );

			// relative min. error at which cg-solvers will abort
			const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
			// aboslute min. error at which cg-solvers will abort
			double outer_absLimit = Parameters().getParam( "absLimit", 1e-8 );
			const double inner_absLimit = Parameters().getParam( "inner_absLimit", 1e-8 );
			const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
			const int maxIter = Parameters().getParam( "maxIter", 500 );
			const bool use_velocity_reconstruct = Parameters().getParam( "use_velocity_reconstruct", true );

	#ifdef USE_BFG_CG_SCHEME
			const double tau = Parameters().getParam( "bfg-tau", 0.1 );
			const bool do_bfg = Parameters().getParam( "do-bfg", true );
	#endif
			logInfo.Resume();
			logInfo << "Begin SaddlePointInverseOperator " << std::endl;

			logDebug.Resume();
			//get some refs for more readability
			PressureDiscreteFunctionType& pressure = dest.discretePressure();
			VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();

			typedef InnerCGSolverWrapper< W_MatrixType,
									M_invers_matrixType,
									X_MatrixType,
									Y_MatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType >
				InnerCGSolverWrapperType;
	#ifdef USE_BFG_CG_SCHEME
			typedef typename InnerCGSolverWrapperType::ReturnValueType
				ReturnValueType;
			ReturnValueType a_solver_info;

			//the bfg scheme uses the outer acc. as a base
			double current_inner_accuracy = do_bfg ? tau * outer_absLimit : inner_absLimit;
			double max_inner_accuracy = current_inner_accuracy;
	#else
			double current_inner_accuracy = inner_absLimit;
	#endif
			InnerCGSolverWrapperType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1_orig.space(),rhs2_orig.space(), relLimit,
														   current_inner_accuracy, solverVerbosity > 3 );

	/*****************************************************************************************/

			int iteration = 0;
		#ifdef USE_BFG_CG_SCHEME
			int total_inner_iterations = 0;
			int min_inner_iterations = std::numeric_limits<int>::max();
			int max_inner_iterations = 0;
		#endif
			const int max_adaptions = Parameters().getParam( "max_adaptions", 2 ) ;
			int current_adaption = 0;
			double delta; //norm of residuum
			double gamma=0;
			double rho;

			// F = rhs2 - X M^{-1} * rhs1
			const double m_scale = m_inv_mat(0,0);
			DiscreteSigmaFunctionType rhs1 = rhs1_orig;
			VelocityDiscreteFunctionType v_tmp ( "v_tmp", velocity.space() );
			x_mat.apply( rhs1, v_tmp );
			v_tmp *=  m_scale;
			VelocityDiscreteFunctionType F( "f", velocity.space() );
			F.assign(rhs2_orig);
			F -= v_tmp;

			VelocityDiscreteFunctionType tmp1( "tmp1", velocity.space() );
			tmp1.clear();
			VelocityDiscreteFunctionType xi( "xi", velocity.space() );
			xi.clear();
			PressureDiscreteFunctionType tmp2( "tmp2", pressure.space() );
			VelocityDiscreteFunctionType residuum( "residuum", velocity.space() );
			VelocityDiscreteFunctionType residuum_squigly( "residuum", velocity.space() );

			PressureDiscreteFunctionType d( "d", pressure.space() );
			PressureDiscreteFunctionType h( "h", pressure.space() );

			typedef SchurkomplementOperator<    InnerCGSolverWrapperType,
												E_MatrixType,
												R_MatrixType,
												Z_MatrixType,
												M_invers_matrixType,
												DiscreteVelocityFunctionType,
												DiscretePressureFunctionType >
					Sk_Operator;
			Sk_Operator sk_op(  innerCGSolverWrapper, e_mat, r_mat, z_mat, m_inv_mat,
								velocity.space(), pressure.space() );

			// r^0 = S * p^0 - f
			residuum.assign( F );
			sk_op.apply( pressure, tmp2 );
			residuum -= tmp2;
			// r_^0 = r^0
			residuum_squigly.assign( residuum );

			// rho = (r_, r )
			rho = residuum.scalarProductDofs( residuum_squigly );
		#if 0
			while( (delta > outer_absLimit ) && (iteration++ < maxIter ) ) {
				if ( iteration > 1 ) {
					// gamma_{m+1} = delta_{m+1} / delta_m
					gamma = delta / gamma;
					// d_{m+1} = r_{m+1} + gamma_m * d_m
					d *= gamma;
					d += residuum;
				}
				if ( iteration >=  maxIter && current_adaption < max_adaptions ) {
					current_adaption++;
					iteration = 2;//do not execute first step in next iter again
					outer_absLimit /= 0.01;
	#ifdef USE_BFG_CG_SCHEME
					current_inner_accuracy /= 0.01;
					innerCGSolverWrapper.setAbsoluteLimit( current_inner_accuracy );
	#endif
					logInfo << "\n\t\t Outer CG solver reset, tolerance lowered" << std::endl;

				}

	#ifdef USE_BFG_CG_SCHEME
					if ( do_bfg ) {
						//the form from the precond. paper (does not work properly)
	//                    current_inner_accuracy = tau * std::min( 1. , absLimit / std::min ( std::pow( delta, int(iteration) ), 1.0 ) );
						//my form, works
						current_inner_accuracy = tau * std::min( 1. , outer_absLimit / std::min ( delta , 1.0 ) );
						innerCGSolverWrapper.setAbsoluteLimit( current_inner_accuracy );
						max_inner_accuracy = std::max( max_inner_accuracy, current_inner_accuracy );
						if( solverVerbosity > 1 )
							logInfo << "\t\t\t set inner limit to: " << current_inner_accuracy << "\n";
					}
	#endif
				// xi = A^{-1} ( B * d )
				tmp1.clear();
				b_mat.apply( d, tmp1 );

	#ifdef USE_BFG_CG_SCHEME
				innerCGSolverWrapper.apply( tmp1, xi, a_solver_info );
	#else
				innerCGSolverWrapper.apply( tmp1, xi );
	#endif

	#ifdef USE_BFG_CG_SCHEME
				if( solverVerbosity > 1 )
					logInfo << "\t\t inner iterations: " << a_solver_info.first << std::endl;
				total_inner_iterations += a_solver_info.first;
				min_inner_iterations = std::min( min_inner_iterations, a_solver_info.first );
				max_inner_iterations = std::max( max_inner_iterations, a_solver_info.first );
	#endif

				// h = B_t * xi  + C * d
				b_t_mat.apply( xi, h );
				tmp2.clear();
				c_mat.apply( d, tmp2 );
				h += tmp2;

				rho = delta / d.scalarProductDofs( h );

				// p_{m+1} = p_m - ( rho_m * d_m )
				pressure.addScaled( d, -rho );
				if ( !use_velocity_reconstruct ) {
					// u_{m+1} = u_m + ( rho_m * xi_m )
					velocity.addScaled( xi, +rho );
				}
				// r_{m+1} = r_m - rho_m * h_m
				residuum.addScaled( h, -rho );

				//save old delta for new gamma calc in next iter
				gamma = delta;

				// d_{m+1} = < r_{m+1} ,r_{m+1} >
				delta = residuum.scalarProductDofs( residuum );

				if( solverVerbosity > 2 )
					logInfo << "\t" << iteration << " SPcg-Iterationen  " << iteration << " Residuum:" << delta << std::endl;
				if ( save_interval > 0 && iteration % save_interval == 0 )
					dest.writeVTK( path, iteration );
			}

			if ( use_velocity_reconstruct ) {
				// u^0 = A^{-1} ( F - B * p^0 )
				schur_F.assign(rhs2);
				tmp1.clear();
				b_mat.apply( pressure, tmp1 );
				schur_F-=tmp1; // F ^= rhs2 - B * p
				innerCGSolverWrapper.apply(schur_F,velocity);
			}

			logInfo << "End SaddlePointInverseOperator " << std::endl;

			SaddlepointInverseOperatorInfo info; //left blank in case of no bfg
	#ifdef USE_BFG_CG_SCHEME
			const double avg_inner_iterations = total_inner_iterations / (double)iteration;
			if( solverVerbosity > 0 )
				logInfo << "\n #avg inner iter | #outer iter: "
						<<  avg_inner_iterations << " | " << iteration << std::endl;

			info.iterations_inner_avg = avg_inner_iterations;
			info.iterations_inner_min = min_inner_iterations;
			info.iterations_inner_max = max_inner_iterations;
			info.iterations_outer_total = iteration;
			info.max_inner_accuracy = max_inner_accuracy;
	#endif
			// ***************************
			return info;
#endif
		} //end BiCgStabSaddlepointInverseOperator::solve

	  };//end class BiCgStabBiCgStabSaddlepointInverseOperator


} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH