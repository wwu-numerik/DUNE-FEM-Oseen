#ifndef DUNE_OSEEN_SOLVERS_NESTED_CG_HH
#define DUNE_OSEEN_SOLVERS_NESTED_CG_HH

#include <dune/fem/oseen/solver/solver_interface.hh>
#include <dune/fem/oseen/solver/schurkomplement.hh>

namespace Dune {

	/**   \brief Nested Conjugate Gradient Solver
		  \tparam OseenPassImp discrete function types etc. get extracted from this
		  \note not fully functional!
		  Iterative Dune solvers are used for both the inner and outer CG iterations
	 */
  template < class OseenPassImp >
  class NestedCgSaddlepointInverseOperator
  {
	private:

	  typedef OseenPassImp OseenPassType;

	  typedef typename OseenPassType::Traits::DiscreteOseenFunctionWrapperType
		  DiscreteOseenFunctionWrapperType;

	  typedef typename OseenPassType::DomainType
		  DomainType;

	  typedef typename OseenPassType::RangeType
		  RangeType;

	  typedef typename DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType
		  PressureDiscreteFunctionType;
	  typedef typename DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType
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

		  typedef typename InnerCGSolverWrapperType::ReturnValueType
				  InnerCGSolverWrapperReturnType;

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
		  typedef typename Sk_Solver::ReturnValueType
				  SolverReturnType;

		  //get some refs for more readability
		  PressureDiscreteFunctionType& pressure = dest.discretePressure();
		  VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();
          auto& logDebug = DSC_LOG_DEBUG;
//		  Logging::LogStream& logError = DSC_LOG_ERROR;
          auto& logInfo = DSC_LOG_INFO;

		  // relative min. error at which cg-solvers will abort
		  const double relLimit = DSC_CONFIG_GET( "relLimit", 1e-4 );
		  // aboslute min. error at which cg-solvers will abort
		  const double absLimit = DSC_CONFIG_GET( "absLimit", 1e-3 );
		  const bool solverVerbosity = DSC_CONFIG_GET( "solverVerbosity", 0 );

		  DiscretePressureFunctionType schur_f ( "schur_f", rhs3.space() );
		  DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
		  {
			  logInfo << "Begin NestedCgSaddlepointInverseOperator\n "
					  << " \n\tbegin calc schur_f,f_func " << std::endl;
			  logDebug.resume();

			  DiscreteSigmaFunctionType m_tmp ( "m_tom", rhs1.space() );
			  DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
			  // f_func =  ( -X *  M_inv * rhs1 + rhs2 )
			  m_inv_mat.apply( rhs1, m_tmp );
			  x_mat.apply( m_tmp, tmp_f );
			  f_func -= tmp_f;
			  f_func += rhs2;



			  // schur_f := -1 * ( ( E * A^-1 * f_func ) - rhs3 )
			  InnerCGSolverWrapperType innerCGSolverWrapper(w_mat,m_inv_mat,x_mat,y_mat,o_mat,rhs1.space(),f_func.space(),relLimit,absLimit,solverVerbosity);
			  assert( !DSFe::FunctionContainsNanOrInf( f_func ) );
			  InnerCGSolverWrapperReturnType a_ret;
			  innerCGSolverWrapper.apply( f_func, tmp_f, a_ret );

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
		  SolverReturnType ret;
		  sk_solver.apply( schur_f, pressure, ret );
		  long total_inner = sk_op.getTotalInnerIterations();
		  logInfo << "\n\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;

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

#endif // DUNE_OSEEN_SOLVERS_NESTED_CG_HH

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

