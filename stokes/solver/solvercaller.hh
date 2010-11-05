#ifndef DUNE_STOKES_SOLVERCALLER_HH
#define DUNE_STOKES_SOLVERCALLER_HH

#include <dune/stokes/solver/nested_cg.hh>
#include <dune/stokes/solver/reduced.hh>
#include <dune/stokes/solver/fullsystem.hh>
#include <dune/stokes/solver/saddle_point.hh>

namespace Dune {

template<class StokesPassType>
struct SolverCaller {
	//! alternative solver implementation
	typedef NestedCgSaddlepointInverseOperator< StokesPassType >
		NestedCgSolverType;
	//! type of the used solver
	typedef SaddlepointInverseOperator< StokesPassType >
		SaddlepointSolverType;
	//! this is used for reduced (no pressure, incompress. condition) oseen pass
	typedef ReducedInverseOperator< StokesPassType >
		ReducedSolverType;
	typedef DirectKrylovSolver< StokesPassType >
		FullsytemSolverType;

	enum SolverID {
		NestedCG_Solver_ID		= 0,
		SaddlePoint_Solver_ID	= 1,
		Reduced_Solver_ID		= 2,
		FullSystem_Solver_ID	= 3
	};

	template <  class DomainType,
				class RangeType,
				class XmatrixObjectType,
				class MInversMatrixObjectType,
				class YmatrixObjectType,
				class OmatrixObjectType,
				class EmatrixObjectType,
				class RmatrixObjectType,
				class ZmatrixObjectType,
				class WmatrixObjectType,
				class DiscreteSigmaFunctionType,
				class DiscreteVelocityFunctionType,
				class DiscretePressureFunctionType  >
	static SaddlepointInverseOperatorInfo solve( const SolverID solverID,
				const DomainType& arg,
				RangeType& dest,
				const XmatrixObjectType& Xmatrix,
				const MInversMatrixObjectType& MInversMatrix,
				const YmatrixObjectType& Ymatrix,
				const OmatrixObjectType& Omatrix,
				const EmatrixObjectType& Ematrix,
				const RmatrixObjectType& Rmatrix,
				const ZmatrixObjectType& Zmatrix,
				const WmatrixObjectType& Wmatrix,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& H2rhs,
				const DiscretePressureFunctionType& H3rhs )
	{
		MatrixWrapper<XmatrixObjectType> X(Xmatrix);
		MatrixWrapper<MInversMatrixObjectType> M_invers(MInversMatrix);
		MatrixWrapper<YmatrixObjectType> Y(Ymatrix);
		MatrixWrapper<OmatrixObjectType> O(Omatrix);
		MatrixWrapper<EmatrixObjectType> E(Ematrix);
		MatrixWrapper<RmatrixObjectType> R(Rmatrix);
		MatrixWrapper<ZmatrixObjectType> Z(Zmatrix);
		MatrixWrapper<WmatrixObjectType> W(Wmatrix);

		switch ( solverID ) {
			case NestedCG_Solver_ID: return NestedCgSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
			case Reduced_Solver_ID: return ReducedSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
			case SaddlePoint_Solver_ID: return SaddlepointSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
			case FullSystem_Solver_ID: return FullsytemSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
			default: throw std::runtime_error("invalid Solver ID selected");
		}
	}
};

} //namespace Dune

#endif // DUNE_STOKES_SOLVERCALLER_HH
