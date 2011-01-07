#ifndef DUNE_STOKES_SOLVERCALLER_HH
#define DUNE_STOKES_SOLVERCALLER_HH

#include <dune/stokes/solver/nested_cg.hh>
#include <dune/stokes/solver/reduced.hh>
#include <dune/stokes/solver/fullsystem.hh>
#include <dune/stokes/solver/saddle_point.hh>
#include <dune/stokes/solver/reconstruction.hh>

namespace Dune {

template<class StokesPassType, template <class T,class S> class ReconstructionPolicyType = BruteForceReconstruction >
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
				class DiscretePressureFunctionType,
				class DataContainerType >
	static SaddlepointInverseOperatorInfo solve( RangeType& dest,
				DataContainerType* rhs_datacontainer,
				const SolverID solverID,
				const DomainType& arg,
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
				const DiscretePressureFunctionType& H3rhs,
				const DiscreteVelocityFunctionType& beta )
	{
		MatrixWrapper<MInversMatrixObjectType> M_invers(MInversMatrix, "M");
		MatrixWrapper<WmatrixObjectType> W(Wmatrix, "W");
		MatrixWrapper<XmatrixObjectType> X(Xmatrix, "X");
		MatrixWrapper<YmatrixObjectType> Y(Ymatrix, "Y");
		MatrixWrapper<OmatrixObjectType> O(Omatrix, "O");
		MatrixWrapper<ZmatrixObjectType> Z(Zmatrix, "Z");
		MatrixWrapper<EmatrixObjectType> E(Ematrix, "E");
		MatrixWrapper<RmatrixObjectType> R(Rmatrix, "R");

		SaddlepointInverseOperatorInfo result;
		switch ( solverID ) {
			case NestedCG_Solver_ID:		result = NestedCgSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case Reduced_Solver_ID:			result = ReducedSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case SaddlePoint_Solver_ID:		result = SaddlepointSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case FullSystem_Solver_ID:		result = FullsytemSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			default: throw std::runtime_error("invalid Solver ID selected");
		}
		if ( rhs_datacontainer ) {
			ReconstructionPolicyType<DataContainerType,typename StokesPassType::DiscreteModelType>
					::reconstruct(	*rhs_datacontainer, dest, beta,
									X, M_invers, Y,
									O, E, R, Z, W,
									H1rhs, H2rhs, H3rhs );
		}
		return result;
	}
};

} //namespace Dune

#endif // DUNE_STOKES_SOLVERCALLER_HH
