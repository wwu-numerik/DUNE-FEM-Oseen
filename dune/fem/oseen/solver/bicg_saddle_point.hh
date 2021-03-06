#ifndef DUNE_OSEEN_SOLVERS_BICG_SADDLE_POINT_HH
#define DUNE_OSEEN_SOLVERS_BICG_SADDLE_POINT_HH

#include <dune/fem/oseen/solver/schurkomplement.hh>
#include <dune/stuff/fem/customprojection.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/fem/functions/analytical.hh>

namespace Dune {

	/**
		\brief Saddlepoint Solver
		The inner CG iteration is implemented by passing a custom Matrix operator to a given\n
		Dune solver. The outer iteration is a implementation of the BICGStab algorithm as described in\n
		van der Vorst: "Iterative Methods for Large Linear Systems" (2000)
	**/
    template < class OseenLDGMethodImp >
	class BiCgStabSaddlepointInverseOperator
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
            auto& logDebug = DSC_LOG_DEBUG;
            auto& logInfo = DSC_LOG_INFO;

			// relative min. error at which cg-solvers will abort
			const double relLimit = DSC_CONFIG_GET( "relLimit", 1e-4 );
			// aboslute min. error at which cg-solvers will abort
			double outer_absLimit = DSC_CONFIG_GET( "absLimit", 1e-8 );
			const double inner_absLimit = DSC_CONFIG_GET( "inner_absLimit", 1e-8 );
			const int solverVerbosity = DSC_CONFIG_GET( "solverVerbosity", 0 );
			const int maxIter = DSC_CONFIG_GET( "maxIter", 500 );

			const double tau = DSC_CONFIG_GET( "bfg-tau", 0.1 );
			const bool do_bfg = DSC_CONFIG_GET( "do-bfg", true );

			logInfo.resume();
			logInfo << cg_name << ": Begin BICG SaddlePointInverseOperator " << std::endl;

			logDebug.resume();
			//get some refs for more readability
			PressureDiscreteFunctionType& pressure = dest.discretePressure();
			VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();

            typedef A_InverseOperator< W_MatrixType,
									M_invers_matrixType,
									X_MatrixType,
									Y_MatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType >
                A_InverseOperatorType;

            typedef typename A_InverseOperatorType::ReturnValueType
				ReturnValueType;
			ReturnValueType a_solver_info;

			//the bfg scheme uses the outer acc. as a base
			double current_inner_accuracy = do_bfg ? tau * outer_absLimit : inner_absLimit;

            A_InverseOperatorType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1_orig.space(),rhs2_orig.space(), relLimit,
														   current_inner_accuracy, solverVerbosity > 5 );

	/*****************************************************************************************/

			// F = rhs2 - X M^{-1} * rhs1
            const double m_scale = m_inv_mat.matrix()(0u,0u);
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

            typedef SchurkomplementOperator<    A_InverseOperatorType,
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
				DSC::printFunctionMinMax( logDebug, schur_f );

            Dune::NewBicgStab< PressureDiscreteFunctionType,Sk_Operator >
					bicg( sk_op, relLimit, outer_absLimit, maxIter, solverVerbosity );
			pressure.clear();
			bicg.apply( schur_f, pressure );
			//pressure mw correction
            const double meanPressure_discrete = DSFe::meanValue( pressure, pressure.space() );
            typedef typename OseenLDGMethodType::Traits::DiscreteModelType::Traits::PressureFunctionSpaceType
					PressureFunctionSpaceType;
			PressureFunctionSpaceType pressureFunctionSpace;
            const DSFe::ConstantFunction<PressureFunctionSpaceType> vol(pressureFunctionSpace, meanPressure_discrete );
            DSFe::BetterL2Projection::project( 0.0, vol, tmp2 );
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

#endif // DUNE_OSEEN_SOLVERS_BICG_SADDLE_POINT_HH

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

