#ifndef DUNE_STOKES_SOLVERS_NESTED_CG_HH
#define DUNE_STOKES_SOLVERS_NESTED_CG_HH

#include <dune/stokes/solver/solver_interface.hh>
#include <dune/stokes/solver/schurkomplement.hh>

namespace Dune {

	/**   \brief Nested Conjugate Gradient Solver
		  \tparam StokesPassImp discrete function types etc. get extracted from this
		  \note not fully functional!
		  Iterative Dune solvers are used for both the inner and outer CG iterations
	 */
  template < class StokesPassImp >
  class NestedCgSaddlepointInverseOperator
  {
	private:

	  typedef StokesPassImp StokesPassType;

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
	  /** \todo Please doc me!
	   * \brief Constructor:
	  *
		**/
	  NestedCgSaddlepointInverseOperator()
	  {}

	  /** takes raw matrices from pass
	  */
	  template <  class X_MatrixType,
				  class M_invers_MatrixType,
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
				  M_invers_MatrixType& m_inv_mat,
				  Y_MatrixType& y_mat,
				  O_MatrixType& o_mat,
				  E_MatrixType& e_mat,
				  R_MatrixType& r_mat,
				  Z_MatrixType& z_mat,
				  W_MatrixType& w_mat,
				  const DiscreteSigmaFunctionType& rhs1,
				  const DiscreteVelocityFunctionType& rhs2,
				  const DiscretePressureFunctionType& rhs3 ) const
	  {
		  typedef InnerCGSolverWrapper< W_MatrixType,
								  M_invers_MatrixType,
								  X_MatrixType,
								  Y_MatrixType,
								  DiscreteSigmaFunctionType,
								  DiscreteVelocityFunctionType >
			  InnerCGSolverWrapperType;

	#ifdef USE_BFG_CG_SCHEME
		  typedef typename InnerCGSolverWrapperType::ReturnValueType
				  InnerCGSolverWrapperReturnType;
	#endif

		  typedef SchurkomplementOperator<    InnerCGSolverWrapperType,
											  E_MatrixType,
											  R_MatrixType,
											  Z_MatrixType,
											  M_invers_MatrixType,
											  DiscreteVelocityFunctionType,
											  DiscretePressureFunctionType >
				  Sk_Operator;

		  typedef SOLVER_NAMESPACE::OUTER_CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
				  Sk_Solver;
	#ifdef USE_BFG_CG_SCHEME
		  typedef typename Sk_Solver::ReturnValueType
				  SolverReturnType;
	#endif
		  //get some refs for more readability
		  PressureDiscreteFunctionType& pressure = dest.discretePressure();
		  VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();
		  Stuff::Logging::LogStream& logDebug = Logger().Dbg();
//		  Logging::LogStream& logError = Logger().Err();
		  Stuff::Logging::LogStream& logInfo = Logger().Info();

		  // relative min. error at which cg-solvers will abort
		  const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
		  // aboslute min. error at which cg-solvers will abort
		  const double absLimit = Parameters().getParam( "absLimit", 1e-3 );
		  const bool solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );

		  DiscretePressureFunctionType schur_f ( "schur_f", rhs3.space() );
		  DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
		  {
			  logInfo << "Begin NestedCgSaddlepointInverseOperator\n "
					  << " \n\tbegin calc schur_f,f_func " << std::endl;
			  logDebug.Resume();

			  DiscreteSigmaFunctionType m_tmp ( "m_tom", rhs1.space() );
			  DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
			  // f_func =  ( -X *  M_inv * rhs1 + rhs2 )
			  m_inv_mat.apply( rhs1, m_tmp );
			  x_mat.apply( m_tmp, tmp_f );
			  f_func -= tmp_f;
			  f_func += rhs2;



			  // schur_f := -1 * ( ( E * A^-1 * f_func ) - rhs3 )
			  InnerCGSolverWrapperType innerCGSolverWrapper(w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),f_func.space(),relLimit,absLimit,solverVerbosity);
			  assert( !Stuff::FunctionContainsNanOrInf( f_func ) );
	  #ifdef USE_BFG_CG_SCHEME
			  InnerCGSolverWrapperReturnType a_ret;
			  innerCGSolverWrapper.apply( f_func, tmp_f, a_ret );
	  #else
			  innerCGSolverWrapper.apply( f_func, tmp_f );
	  #endif
			  e_mat.apply( tmp_f, schur_f );
			  schur_f -= rhs3;
			  schur_f *= -1;

			  logInfo << " \n\tend calc schur_f,f_func\n "
					  << " \n\tbegin S*p=schur_f " << std::endl;
		  }

		  InnerCGSolverWrapperType innerCGSolverWrapper__ (w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),velocity.space(),relLimit,absLimit,solverVerbosity);
		  Sk_Operator sk_op(  innerCGSolverWrapper__, e_mat, r_mat, z_mat, m_inv_mat,
							  velocity.space(), pressure.space() );
		  Sk_Solver sk_solver( sk_op, relLimit, absLimit, 2000, solverVerbosity );

		  // p = S^-1 * schur_f = ( E * A^-1 * Z + rhs3 )^-1 * schur_f
  #ifdef USE_BFG_CG_SCHEME
		  SolverReturnType ret;
		  sk_solver.apply( schur_f, pressure, ret );
		  long total_inner = sk_op.getTotalInnerIterations();
		  logInfo << "\n\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;
  #else
		  sk_solver.apply( schur_f, pressure );
  #endif

		  logInfo << "\n\tend  S*p=schur_f" << std::endl;

		  DiscreteVelocityFunctionType Zp_temp ( "Zp_temp", velocity.space() );
		  // velocity = A^-1 * ( ( -1 * ( Z * pressure ) ) + f_func )
		  z_mat.apply( pressure, Zp_temp );
		  Zp_temp *= -1.0;
		  Zp_temp += f_func;
		  InnerCGSolverWrapperType innerCGSolverWrapper(w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),velocity.space(),relLimit,absLimit,solverVerbosity);
		  innerCGSolverWrapper.apply ( Zp_temp, velocity );

		  logInfo << "\nEnd NestedCgSaddlePointInverseOperator " << std::endl;

		  return SaddlepointInverseOperatorInfo();
	  }

	};


} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_NESTED_CG_HH
