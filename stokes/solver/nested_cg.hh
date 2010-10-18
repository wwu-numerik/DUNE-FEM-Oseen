#ifndef DUNE_STOKES_SOLVERS_NESTED_CG_HH
#define DUNE_STOKES_SOLVERS_NESTED_CG_HH

#include <dune/stokes/solver/solver_interface.hh>

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

	  typedef typename StokesPassType::DiscreteStokesFunctionWrapperType
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
	  template <  class XmatrixObjectType,
				  class MmatrixObjectType,
				  class YmatrixObjectType,
				  class EmatrixObjectType,
				  class RmatrixObjectType,
				  class ZmatrixObjectType,
				  class WmatrixObjectType,
				  class DiscreteSigmaFunctionType,
				  class DiscreteVelocityFunctionType,
				  class DiscretePressureFunctionType  >
	  SaddlepointInverseOperatorInfo solve( const DomainType& arg,
				  RangeType& dest,
				  XmatrixObjectType& Xmatrix,
				  MmatrixObjectType& Mmatrix,
				  YmatrixObjectType& Ymatrix,
				  YmatrixObjectType& Omatrix,
				  EmatrixObjectType& Ematrix,
				  RmatrixObjectType& Rmatrix,
				  ZmatrixObjectType& Zmatrix,
				  WmatrixObjectType& Wmatrix,
				  const DiscreteSigmaFunctionType& rhs1,
				  DiscreteVelocityFunctionType& rhs2,
				  DiscretePressureFunctionType& rhs3 ) const
	  {

		  Logging::LogStream& logDebug = Logger().Dbg();
		  Logging::LogStream& logError = Logger().Err();
		  Logging::LogStream& logInfo = Logger().Info();

		  if ( Parameters().getParam( "disableSolver", false ) ) {
			  logInfo.Resume();
			  logInfo << "solving disabled via parameter file" << std::endl;
			  return SaddlepointInverseOperatorInfo();
		  }

		  // relative min. error at which cg-solvers will abort
		  const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
		  // aboslute min. error at which cg-solvers will abort
		  const double absLimit = Parameters().getParam( "absLimit", 1e-3 );
		  const bool solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );

		  logInfo << "Begin NestedCgSaddlepointInverseOperator " << std::endl;

		  logDebug.Resume();
		  //get some refs for more readability
		  PressureDiscreteFunctionType& pressure = dest.discretePressure();
		  VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();

		  typedef typename  XmatrixObjectType::MatrixType
			  XmatrixType;
		  typedef typename  MmatrixObjectType::MatrixType
			  MmatrixType;
		  typedef typename  YmatrixObjectType::MatrixType
			  YmatrixType;
		  typedef typename  EmatrixObjectType::MatrixType
			  B_t_matrixType;                             //! renamed
		  typedef typename  RmatrixObjectType::MatrixType
			  CmatrixType;                                //! renamed
		  typedef typename  ZmatrixObjectType::MatrixType
			  BmatrixType;                                //! renamed
		  typedef typename  WmatrixObjectType::MatrixType
			  WmatrixType;

		  XmatrixType& x_mat      = Xmatrix.matrix();
		  MmatrixType& m_inv_mat  = Mmatrix.matrix();
		  YmatrixType& y_mat      = Ymatrix.matrix();
		  YmatrixType& o_mat      = Omatrix.matrix();
		  B_t_matrixType& b_t_mat = Ematrix.matrix(); //! renamed
		  CmatrixType& c_mat      = Rmatrix.matrix(); //! renamed
		  BmatrixType& b_mat      = Zmatrix.matrix(); //! renamed
		  WmatrixType& w_mat      = Wmatrix.matrix();

		  b_t_mat.scale( -1 ); //since B_t = -E

		  DiscretePressureFunctionType& g_func = rhs3;
		  g_func *= ( -1 ); //since G = -H_3

		  logInfo << " \n\tbegin calc new_f,f_func " << std::endl;
		  //Stuff::DiagonalMult( m_inv_mat, rhs1 ); //calc m_inv * H_1 "in-place"
		  DiscreteSigmaFunctionType m_tmp ( "m_tom", rhs1.space() );
		  DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
		  f_func.clear();
		  m_tmp.clear();

		  // f_func = ( ( -1 * ( X * ( M_inv * rhs1 ) ) ) + rhs2 )
		  m_inv_mat.apply( rhs1, m_tmp );
		  x_mat.apply( m_tmp, f_func );
		  f_func -= rhs2;


		  typedef InnerCGSolverWrapper< WmatrixType,
								  MmatrixType,
								  XmatrixType,
								  YmatrixType,
								  DiscreteSigmaFunctionType,
								  DiscreteVelocityFunctionType >
			  InnerCGSolverWrapperType;
  #ifdef USE_BFG_CG_SCHEME
		  typedef typename InnerCGSolverWrapperType::ReturnValueType
				  InnerCGSolverWrapperReturnType;
  #endif
		  #define innerCGSolverWrapper InnerCGSolverWrapperType(w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),f_func.space(),relLimit,absLimit,solverVerbosity)


		  DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
		  DiscretePressureFunctionType new_f ( "new_f", g_func.space() );
		  tmp_f.clear();
		  new_f.clear();

		  // new_f := ( B * A^-1 * f_func ) + g_func
  #ifdef USE_BFG_CG_SCHEME
		  InnerCGSolverWrapperReturnType a_ret;
		  innerCGSolverWrapper.apply( f_func, tmp_f, a_ret );
  #else
		  innerCGSolverWrapper.apply( f_func, tmp_f );
  #endif
  #if defined(ADAPTIVE_SOLVER) && defined(USE_BFG_CG_SCHEME)
		  if ( isnan(a_ret.second) ) {
			  logInfo << "\n\t\t NaNs detected, lowering error tolerance" << std::endl;
			  int max_adaptions = Parameters().getParam( "max_adaptions", 8 );
			  int adapt_step = 1;
			  double a_relLimit = relLimit;
			  double a_absLimit = absLimit;
			  while ( true ) {
				  if ( adapt_step > max_adaptions ) {
					  logInfo << "\n\t\t max adaption depth reached, aborting" << std::endl;
					  break;
				  }
				  a_relLimit /= 10.0;
				  a_absLimit /= 10.0;
				  logInfo << "\n\t\t\t trying with relLimit " << a_relLimit
						  << " and absLimit " << a_absLimit << std::endl;

				  InnerCGSolverWrapperType innerCGSolverWrapper_adapt( w_mat, m_inv_mat, x_mat,
																	   y_mat, o_mat, rhs1.space(),f_func.space(),
																	   relLimit, absLimit, solverVerbosity );
				  innerCGSolverWrapper_adapt.apply( f_func, tmp_f, a_ret );

				  if ( !isnan(a_ret.second) ) {
					  logInfo << "\n\t\t adaption produced NaN-free solution" << std::endl;
					  break;
				  }
				  adapt_step++;
			  }
		  }
  #endif //defined(ADAPTIVE_SOLVER) && defined(USE_BFG_CG_SCHEME)
		  b_t_mat.apply( tmp_f, new_f );
		  new_f -= g_func;

		  logInfo << " \n\tend calc new_f,f_func " << std::endl;

		  typedef SchurkomplementOperator<    InnerCGSolverWrapperType,
											  B_t_matrixType,
											  CmatrixType,
											  BmatrixType,
											  MmatrixType,
											  DiscreteVelocityFunctionType,
											  DiscretePressureFunctionType >
				  Sk_Operator;

		  typedef SOLVER_NAMESPACE::OUTER_CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
				  Sk_Solver;
  #ifdef USE_BFG_CG_SCHEME
		  typedef typename Sk_Solver::ReturnValueType
				  SolverReturnType;
  #endif

		  logInfo << " \n\tbegin S*p=new_f " << std::endl;
		  InnerCGSolverWrapperType innerCGSolverWrapper__ (w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),f_func.space(),relLimit,absLimit,solverVerbosity);
		  Sk_Operator sk_op(  innerCGSolverWrapper__, b_t_mat, c_mat, b_mat, m_inv_mat,
							  velocity.space(), pressure.space() );
		  Sk_Solver sk_solver( sk_op, relLimit, absLimit, 2000, solverVerbosity );
		  pressure.clear();

		  // p = S^-1 * new_f = ( B_t * A^-1 * B + rhs3 )^-1 * new_f
  #ifdef USE_BFG_CG_SCHEME
		  SolverReturnType ret;
		  sk_solver.apply( new_f, pressure, ret );
		  long total_inner = sk_op.getTotalInnerIterations();
		  logInfo << "\n\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;
  #else
		  sk_solver.apply( new_f, pressure );
  #endif

  #if defined(ADAPTIVE_SOLVER) && defined(USE_BFG_CG_SCHEME)
		  if ( isnan(ret.second) ) {
			  logInfo << "\n\t\t NaNs detected, lowering error tolerance" << std::endl;
			  int max_adaptions = Parameters().getParam( "max_adaptions", 8 );
			  int adapt_step = 1;
			  double a_relLimit = relLimit;
			  double a_absLimit = absLimit;
			  while ( true ) {
				  if ( adapt_step > max_adaptions ) {
					  logError << "\n\t\t max adaption depth reached, aborting" << std::endl;
					  break;
				  }
				  a_relLimit /= 10.0;
				  a_absLimit /= 10.0;
				  logInfo << "\n\t\t\t trying with relLimit " << a_relLimit
						  << " and absLimit " << a_absLimit << std::endl;

				  InnerCGSolverWrapperType innerCGSolverWrapper_adapt( w_mat, m_inv_mat, x_mat,
																	  y_mat, o_mat, rhs1.space(),f_func.space(),
																	  relLimit, absLimit, solverVerbosity );
				  Sk_Operator sk_op_adapt(  innerCGSolverWrapper_adapt, b_t_mat, c_mat, b_mat, m_inv_mat,
							  velocity.space(), pressure.space() );
				  Sk_Solver sk_solver_adapt( sk_op_adapt, a_relLimit, a_absLimit, 2000, solverVerbosity );
				  pressure.clear();
				  sk_solver_adapt.apply( new_f, pressure, ret );
				  long total_inner = sk_op_adapt.getTotalInnerIterations();
				  logInfo << "\n\t\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;
				  if ( !isnan(ret.second) ) {
					  logInfo << "\n\t\t adaption produced NaN-free solution" << std::endl;
					  break;
				  }
				  adapt_step++;
			  }
		  }
  #endif //defined(ADAPTIVE_SOLVER) && defined(USE_BFG_CG_SCHEME)
		  //
		  logInfo << "\n\tend  S*p=new_f" << std::endl;

		  pressure *= -1;//magic
		  DiscreteVelocityFunctionType Bp_temp ( "Bp_temp", velocity.space() );
		  Bp_temp.clear();
		  // velocity = A^-1 * ( ( -1 * ( B * pressure ) ) + f_func )
		  b_mat.apply( pressure, Bp_temp );
  //		Bp_temp *= -1;
		  Bp_temp += f_func;
		  innerCGSolverWrapper.apply ( Bp_temp, velocity );
		  velocity *= -1;//even more magic

		  logInfo << "\nEnd NestedCgSaddlePointInverseOperator " << std::endl;

		  return SaddlepointInverseOperatorInfo();
	  #undef innerCGSolverWrapper
	  }

	};


} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_NESTED_CG_HH
