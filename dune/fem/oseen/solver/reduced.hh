#ifndef DUNE_OSEEN_SOLVERS_REDUCED_HH
#define DUNE_OSEEN_SOLVERS_REDUCED_HH

#include <dune/fem/oseen/solver/schurkomplement.hh>
#include <dune/fem/oseen/solver/cghelper.hh>

namespace Dune {


	/** the goal is to solve
		\$ ( Y + O - X M^{-1} W )u = H_2 - X M^{-1} H_1 \$
		eh, copy pasta doc fail?
		for u
	  **/
	template < class OseenLDGMethodImp >
	class ReducedInverseOperator
	{
	  private:

		typedef OseenLDGMethodImp
			OseenLDGMethodType;

		typedef typename OseenLDGMethodType::Traits::DiscreteOseenFunctionWrapperType
			DiscreteOseenFunctionWrapperType;

		typedef typename OseenLDGMethodType::DomainType
			DomainType;

		typedef typename OseenLDGMethodType::RangeType
			RangeType;

		typedef typename DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType
			PressureDiscreteFunctionType;
		typedef typename DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType
			VelocityDiscreteFunctionType;


	  public:
		ReducedInverseOperator()
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
					X_MatrixType& Xmatrix,
					M_invers_matrixType& Mmatrix,
					Y_MatrixType& Ymatrix,
					O_MatrixType& Omatrix,
					E_MatrixType& /*Ematrix*/,
					R_MatrixType& /*Rmatrix*/,
					Z_MatrixType& /*Zmatrix*/,
					W_MatrixType& Wmatrix,
					const DiscreteSigmaFunctionType& rhs1,
					const DiscreteVelocityFunctionType& rhs2,
					const DiscretePressureFunctionType& /*rhs3*/ ) const
		{

            auto& logDebug = DSC_LOG_DEBUG;
            auto& logInfo = DSC_LOG_INFO;

			if ( DSC_CONFIG_GET( "disableSolver", false ) ) {
				logInfo.resume();
				logInfo << "solving disabled via parameter file" << std::endl;
				return SaddlepointInverseOperatorInfo();
			}

			// relative min. error at which cg-solvers will abort
			const double relLimit = DSC_CONFIG_GET( "relLimit", 1e-4 );
			// aboslute min. error at which cg-solvers will abort
			const double inner_absLimit = DSC_CONFIG_GET( "inner_absLimit", 1e-8 );
			const int solverVerbosity = DSC_CONFIG_GET( "solverVerbosity", 0 );

			logInfo.resume();
			logInfo << "Begin ReducedInverseOperator " << std::endl;

			logDebug.resume();
			//get some refs for more readability
			VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();

			X_MatrixType& x_mat      = Xmatrix;
			M_invers_matrixType& m_inv_mat  = Mmatrix;
			Y_MatrixType& y_mat      = Ymatrix;
			O_MatrixType& o_mat      = Omatrix;
			W_MatrixType& w_mat      = Wmatrix;

			VelocityDiscreteFunctionType F( "f", velocity.space() );
			// F = H_2 - X M^{-1} H_1
			DiscreteSigmaFunctionType sig_tmp( "sig_tmp", rhs1.space() );
			m_inv_mat.apply( rhs1, sig_tmp );
			x_mat.apply( sig_tmp, F );
			F *= -1;
			F += rhs2;

			typedef A_InverseOperator< W_MatrixType,
									M_invers_matrixType,
									X_MatrixType,
									Y_MatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType >
				A_InverseOperatorType;
			A_InverseOperatorType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1.space(), rhs2.space(), relLimit,
														  inner_absLimit, solverVerbosity  );
			innerCGSolverWrapper.apply(F,velocity);
			logInfo << "End ReducedInverseOperator " << std::endl;

			return SaddlepointInverseOperatorInfo();
		} //end ReducedInverseOperator::solve


	  };//end class ReducedInverseOperator

} //end namespace Dune

#endif // DUNE_OSEEN_SOLVERS_REDUCED_HH

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

