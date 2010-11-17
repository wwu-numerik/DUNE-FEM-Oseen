#ifndef DUNE_STOKES_SOLVERS_REDUCED_HH
#define DUNE_STOKES_SOLVERS_REDUCED_HH

#include <dune/stokes/solver/solver_interface.hh>

namespace Dune {


	/** the goal is to solve
		\$ ( Y + O - X M^{-1} W )u = H_2 + X M^{-1} H_1 \$
		for u
	  **/
	template < class StokesPassImp >
	class ReducedInverseOperator
	{
	  private:

		typedef StokesPassImp
			StokesPassType;

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
			const double inner_absLimit = Parameters().getParam( "inner_absLimit", 1e-8 );
			const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );

			logInfo.Resume();
			logInfo << "Begin ReducedInverseOperator " << std::endl;

			logDebug.Resume();
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

			typedef InnerCGSolverWrapper< W_MatrixType,
									M_invers_matrixType,
									X_MatrixType,
									Y_MatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType >
				InnerCGSolverWrapperType;
			InnerCGSolverWrapperType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1.space(), rhs2.space(), relLimit,
														  inner_absLimit, solverVerbosity  );
			SaddlepointInverseOperatorInfo info;
			innerCGSolverWrapper.apply(F,velocity);
			logInfo << "End ReducedInverseOperator " << std::endl;

			return info;
		} //end ReducedInverseOperator::solve


	  };//end class ReducedInverseOperator

} //end namespace Dune

#endif // DUNE_STOKES_SOLVERS_REDUCED_HH
