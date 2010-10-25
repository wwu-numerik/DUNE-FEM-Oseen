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
		SaddlepointInverseOperatorInfo solve( const DomainType& /*arg*/,
					RangeType& dest,
					XmatrixObjectType& Xmatrix,
					MmatrixObjectType& Mmatrix,
					YmatrixObjectType& Ymatrix,
					YmatrixObjectType& Omatrix,
					EmatrixObjectType& /*Ematrix*/,
					RmatrixObjectType& /*Rmatrix*/,
					ZmatrixObjectType& /*Zmatrix*/,
					WmatrixObjectType& Wmatrix,
					DiscreteSigmaFunctionType& rhs1,
					DiscreteVelocityFunctionType& rhs2,
					DiscretePressureFunctionType& /*rhs3*/ ) const
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

			typedef typename  XmatrixObjectType::MatrixType
				XmatrixType;
			typedef typename  MmatrixObjectType::MatrixType
				MmatrixType;
			typedef typename  YmatrixObjectType::MatrixType
				YmatrixType;
			typedef typename  WmatrixObjectType::MatrixType
				WmatrixType;

			XmatrixType& x_mat      = Xmatrix.matrix();
			MmatrixType& m_inv_mat  = Mmatrix.matrix();
			YmatrixType& y_mat      = Ymatrix.matrix();
			YmatrixType& o_mat      = Omatrix.matrix();
			WmatrixType& w_mat      = Wmatrix.matrix();

			VelocityDiscreteFunctionType F( "f", velocity.space() );
			// F = H_2 - X M^{-1} H_1
			rhs1 *=  m_inv_mat(0,0);
			x_mat.apply( rhs1, F );
			F *= -1;
			F += rhs2;

			typedef InnerCGSolverWrapper< WmatrixType,
									MmatrixType,
									XmatrixType,
									YmatrixType,
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
