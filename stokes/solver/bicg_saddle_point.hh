#ifndef DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH
#define DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH

#include <dune/stokes/solver/solver_interface.hh>
#include <dune/stokes/solver/schurkomplement.hh>
namespace Dune {

	/**
		\brief Saddlepoint Solver
		The inner CG iteration is implemented by passing a custom Matrix operator to a given\n
		Dune solver. The outer iteration is a implementation of the BICGStab algorithm as described in\n
		van der Vorst: "Iterative Methods for Large Linear Systems" (2000)
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
			const std::string cg_name( "OuterCG");
			Logging::LogStream& logDebug = Logger().Dbg();
			Logging::LogStream& logInfo = Logger().Info();

			// relative min. error at which cg-solvers will abort
			const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
			// aboslute min. error at which cg-solvers will abort
			double outer_absLimit = Parameters().getParam( "absLimit", 1e-8 );
			const double inner_absLimit = Parameters().getParam( "inner_absLimit", 1e-8 );
			const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
			const int maxIter = Parameters().getParam( "maxIter", 500 );

	#ifdef USE_BFG_CG_SCHEME
			const double tau = Parameters().getParam( "bfg-tau", 0.1 );
			const bool do_bfg = Parameters().getParam( "do-bfg", true );
	#endif
			logInfo.Resume();
			logInfo << cg_name << ": Begin BICG SaddlePointInverseOperator " << std::endl;

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
														   current_inner_accuracy, solverVerbosity > 5 );

	/*****************************************************************************************/

			int iteration = 1;
		#ifdef USE_BFG_CG_SCHEME
			int total_inner_iterations = 0;
			int min_inner_iterations = std::numeric_limits<int>::max();
			int max_inner_iterations = 0;
		#endif

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
			PressureDiscreteFunctionType search_direction( "search", pressure.space() );

			PressureDiscreteFunctionType tmp2( "tmp2", pressure.space() );
			PressureDiscreteFunctionType residuum( "residuum", pressure.space() );
			PressureDiscreteFunctionType start_residuum( "residuum", pressure.space() );
			PressureDiscreteFunctionType s( "s", pressure.space() );
			PressureDiscreteFunctionType t( "s", pressure.space() );
			PressureDiscreteFunctionType v( "s", pressure.space() );
			PressureDiscreteFunctionType schur_f( "s", pressure.space() );

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

			//schur_f = - E A^{-1} F - G = -1 ( E A^{-1} F + G )
			schur_f.assign( rhs3 );
			innerCGSolverWrapper.apply( F, v_tmp );
			e_mat.apply( v_tmp, tmp2 );
			schur_f += tmp2;
			schur_f *= -1;

			// r^0 = - S * p^0 + schur_f
			sk_op.apply( pressure, residuum );
			residuum -= schur_f;
			residuum *= -1.;

			// r_^0 = r^0
			start_residuum.assign( residuum );
			search_direction.assign( residuum );

			double rho;
			double delta; //norm of residuum

			double alpha,omega,last_rho;
			if ( solverVerbosity > 3 )
				Stuff::printFunctionMinMax( logDebug, pressure );
			PressureDiscreteFunctionType residuum_T( "s", pressure.space() );
			residuum_T.assign( start_residuum );//??

			while( iteration < maxIter ) {
				rho = residuum_T.scalarProductDofs( residuum );
				if ( rho == 0.0 ) {
					if ( solverVerbosity > 3 )
						logDebug << boost::format( "%s: abort, theta = %e") % cg_name % rho << std::endl;
					break;
				}
				assert( !std::isnan(rho) );
				assert( std::isfinite(rho) );

				if ( iteration == 1 ) {
					search_direction.assign( start_residuum );
				}
				else {
					const double beta = (rho/last_rho) * (alpha/omega);
					search_direction *= beta;
					search_direction.addScaled( v, -beta*omega);
					search_direction += residuum;
				}

				sk_op.apply( search_direction, v );//v=S*p
				alpha = rho/ residuum_T.scalarProductDofs( v );
				assert( !std::isnan(alpha) );
				assert( std::isfinite(alpha) );

				s.assign( residuum );
				s.addScaled( v, -alpha );
				const double s_norm = std::sqrt( s.scalarProductDofs( s ) );
				if ( s_norm < outer_absLimit ) {
					pressure.addScaled( search_direction, alpha );
					logDebug << boost::format( "%s: iter %i\taborted: s: %e") % cg_name % iteration % s_norm << std::endl;
					break;
				}
				sk_op.apply( s, t );
				omega = t.scalarProductDofs( s ) / t.scalarProductDofs( t );

				if ( solverVerbosity > 3 )
					Stuff::printFunctionMinMax( logDebug, search_direction );
				pressure.addScaled( search_direction, alpha );
				if ( solverVerbosity > 3 )
					Stuff::printFunctionMinMax( logDebug, pressure );
				pressure.addScaled( s, omega );

				residuum.assign( s );
				residuum.addScaled( t, - omega );

				delta = std::sqrt( residuum.scalarProductDofs( residuum ) );
				if ( delta < outer_absLimit ) {
					logDebug << boost::format( "%s: aborted, iter %i\tres %e") % cg_name % iteration % delta << std::endl;
					break;
				}

				if ( solverVerbosity > 3 ) {
					logDebug << boost::format( "%s: iter %i\tres %e alpha %e \trho %e") % cg_name % iteration
								% delta % alpha %  rho << std::endl;
					Stuff::printFunctionMinMax( logDebug, pressure );
				}
				assert( omega != 0.0 );

				last_rho = rho;
				iteration++;
			} //end while
			// u = A^{-1} ( F - B * p^0 )
			v_tmp.assign(F);
			z_mat.apply( pressure, tmp1 );
			v_tmp-=tmp1; // F ^= rhs2 - B * p
			innerCGSolverWrapper.apply(v_tmp,velocity);
			logInfo << cg_name << ": End BICG SaddlePointInverseOperator " << std::endl;

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

		} //end BiCgStabSaddlepointInverseOperator::solve

	  };//end class BiCgStabBiCgStabSaddlepointInverseOperator


} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH
