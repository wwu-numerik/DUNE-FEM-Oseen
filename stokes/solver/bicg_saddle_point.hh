#ifndef DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH
#define DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH

#include <dune/stokes/solver/solver_interface.hh>
#include <dune/stokes/solver/schurkomplement.hh>
#include <dune/stuff/customprojection.hh>

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
	#else
			double current_inner_accuracy = inner_absLimit;
	#endif
			InnerCGSolverWrapperType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1_orig.space(),rhs2_orig.space(), relLimit,
														   current_inner_accuracy, solverVerbosity > 5 );

	/*****************************************************************************************/

			// F = rhs2 - X M^{-1} * rhs1
			const double m_scale = m_inv_mat(0,0);
			DiscreteSigmaFunctionType rhs1 = rhs1_orig;
			VelocityDiscreteFunctionType v_tmp ( "v_tmp", velocity.space() );
			x_mat.apply( rhs1, v_tmp );
			v_tmp *=  m_scale;
			VelocityDiscreteFunctionType F( "F", velocity.space() );
			F.assign(rhs2_orig);
			F -= v_tmp;

			VelocityDiscreteFunctionType tmp1( "tmp1", velocity.space() );
			PressureDiscreteFunctionType tmp2( "tmp2", pressure.space() );

			PressureDiscreteFunctionType schur_f( "schur_f", pressure.space() );

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
			v_tmp.clear();
			innerCGSolverWrapper.apply( F, v_tmp );
			e_mat.apply( v_tmp, tmp2 );
			schur_f += tmp2;
			schur_f *= -1;
			if ( solverVerbosity > 3 )
				Stuff::printFunctionMinMax( logDebug, schur_f );

			SOLVER_NAMESPACE::NewBicgStab< PressureDiscreteFunctionType,Sk_Operator >
					bicg( sk_op, relLimit, outer_absLimit, maxIter, solverVerbosity );
			pressure.clear();
			bicg.apply( schur_f, pressure );
//			pressure *=-1;
			//pressure mw correction
			double meanPressure_discrete = Stuff::meanValue( pressure, pressure.space() );
			typedef typename StokesPassType::Traits::DiscreteModelType::Traits::PressureFunctionSpaceType
					PressureFunctionSpaceType;
			PressureFunctionSpaceType pressureFunctionSpace;
			Stuff::ConstantFunction<PressureFunctionSpaceType> vol(pressureFunctionSpace, meanPressure_discrete );
			Dune::BetterL2Projection
				::project( 0.0, vol, tmp2 );
			pressure -= tmp2;


			// u = A^{-1} ( F - B * p^0 )
			v_tmp.assign(F);
			z_mat.apply( pressure, tmp1 );
			v_tmp-=tmp1; // F ^= rhs2 - B * p
			innerCGSolverWrapper.apply(v_tmp,velocity);
			logInfo << cg_name << ": End BICG SaddlePointInverseOperator " << std::endl;

			SaddlepointInverseOperatorInfo info; //left blank in case of no bfg
			// ***************************
			return info;

		} //end BiCgStabSaddlepointInverseOperator::solve

	  };//end class BiCgStabBiCgStabSaddlepointInverseOperator


} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_BICG_SADDLE_POINT_HH
