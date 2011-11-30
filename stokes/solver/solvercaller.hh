#ifndef DUNE_STOKES_SOLVERCALLER_HH
#define DUNE_STOKES_SOLVERCALLER_HH

#include <dune/stokes/solver/nested_cg.hh>
#include <dune/stokes/solver/reduced.hh>
#include <dune/stokes/solver/fullsystem.hh>
#include <dune/stokes/solver/saddle_point.hh>
#include <dune/stokes/solver/bicg_saddle_point.hh>
#include <dune/stokes/solver/reconstruction.hh>
#include <dune/stuff/profiler.hh>

namespace Dune {

namespace Stokes {

	namespace Solver {
		enum SolverID {
			NestedCG_Solver_ID			= 0,
			SaddlePoint_Solver_ID		= 1,
			Reduced_Solver_ID			= 2,
			FullSystem_Solver_ID		= 3,
			BiCg_Saddlepoint_Solver_ID	= 4
		};
	}

template<class StokesPassType, template <class T,class S> class ReconstructionPolicyType = BruteForceReconstruction >
struct SolverCaller {
	//! alternative solver implementation
	typedef NestedCgSaddlepointInverseOperator< StokesPassType >
		NestedCgSolverType;
	//! type of the used solver
	typedef SaddlepointInverseOperator< StokesPassType >
		SaddlepointSolverType;
	typedef BiCgStabSaddlepointInverseOperator< StokesPassType >
		BiCgSaddlepointSolverType;
	//! this is used for reduced (no pressure, incompress. condition) oseen pass
	typedef ReducedInverseOperator< StokesPassType >
		ReducedSolverType;
	typedef DirectKrylovSolver< StokesPassType >
		FullsytemSolverType;

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
				const Solver::SolverID solverID,
				const bool with_oseen_discretization,
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
		Stuff::Profiler::ScopedTiming solver_time("solver");
		MatrixWrapper<MInversMatrixObjectType> M_invers(MInversMatrix, "M");
		MatrixWrapper<WmatrixObjectType> W(Wmatrix, "W");
		MatrixWrapper<XmatrixObjectType> X(Xmatrix, "X");
		MatrixWrapper<YmatrixObjectType> Y(Ymatrix, "Y");
		MatrixWrapper<OmatrixObjectType> O(Omatrix, "O");
		MatrixWrapper<ZmatrixObjectType> Z(Zmatrix, "Z");
		MatrixWrapper<EmatrixObjectType> E(Ematrix, "E");
		MatrixWrapper<RmatrixObjectType> R(Rmatrix, "R");

		#if 0 //#ifndef NDEBUG
		{
		    // do the matlab logging stuff
		    if ( Parameters().getParam( "save_matrices", false ) ) {
			    Stuff::Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
			    Stuff::printSparseRowMatrixMatlabStyle( MInversMatrix.matrix(), "M_invers", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Wmatrix.matrix(), "W", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Omatrix.matrix(), "O", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Xmatrix.matrix(), "X", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Ymatrix.matrix(), "Y", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Zmatrix.matrix(), "Z", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Ematrix.matrix(), "E", matlabLogStream );
			    Stuff::printSparseRowMatrixMatlabStyle( Rmatrix.matrix(), "R", matlabLogStream );
			    Stuff::printDiscreteFunctionMatlabStyle( H1rhs, "H1", matlabLogStream );
			    Stuff::printDiscreteFunctionMatlabStyle( H2rhs, "H2", matlabLogStream );
			    Stuff::printDiscreteFunctionMatlabStyle( H3rhs, "H3", matlabLogStream );
			    Stuff::printDiscreteFunctionMatlabStyle( H2_O_rhs, "H_O", matlabLogStream );
			    matlabLogStream.Flush();
		    }

		//            // log local matrices
		//            Stuff::GridWalk<GridPartType> gw( gridPart_ );
		//            typedef Logging::MatlabLogStream
		//                FunctorStream;
		//            FunctorStream& functorStream = matlabLogStream;
		//            Stuff::GridWalk<DiscreteVelocityFunctionSpaceType> gw( velocitySpace_ );
		//            Stuff::LocalMatrixPrintFunctor< EmatrixType,FunctorStream> f_E ( Ematrix, functorStream, "E" );
		//            Stuff::LocalMatrixPrintFunctor< WmatrixType,FunctorStream> f_W ( Wmatrix, functorStream, "W" );
		//            Stuff::LocalMatrixPrintFunctor< XmatrixType,FunctorStream> f_X ( Xmatrix, functorStream, "X" );
		//            Stuff::LocalMatrixPrintFunctor< YmatrixType,FunctorStream> f_Y ( Ymatrix, functorStream, "Y" );
		//            Stuff::LocalMatrixPrintFunctor< ZmatrixType,FunctorStream> f_Z ( Zmatrix, functorStream, "Z" );
		//            Stuff::LocalMatrixPrintFunctor< RmatrixType,FunctorStream> f_R ( Rmatrix, functorStream, "R" );
		//                gw( f_W );
		//                gw( f_X );
		//                gw( f_Y );
		//                gw( f_Z );
		//                gw( f_E );
		//                gw( f_R );
		// do profiling


		    if ( Parameters().getParam( "outputMatrixPlots", false ) ) {
			Stuff::matrixToGnuplotFile( Ematrix.matrix(),       std::string( "mat_E.gnuplot")       );
			Stuff::matrixToGnuplotFile( Wmatrix.matrix(),       std::string( "mat_W.gnuplot")       );
			Stuff::matrixToGnuplotFile( Xmatrix.matrix(),       std::string( "mat_X.gnuplot")       );
			Stuff::matrixToGnuplotFile( Ymatrix.matrix(),       std::string( "mat_Y.gnuplot")       );
			Stuff::matrixToGnuplotFile( Zmatrix.matrix(),       std::string( "mat_Z.gnuplot")       );
			Stuff::matrixToGnuplotFile( Rmatrix.matrix(),       std::string( "mat_R.gnuplot")       );
			Stuff::matrixToGnuplotFile( MInversMatrix.matrix(), std::string( "mat_M.gnuplot")   );
		    }

		    if ( Parameters().getParam( "paranoid_checks", false ) )
		    {//paranoid checks
			    assert( !Stuff::MatrixContainsNanOrInf( Omatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Ematrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Wmatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Xmatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Ymatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Zmatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( Rmatrix.matrix() ) );
			    assert( !Stuff::MatrixContainsNanOrInf( MInversMatrix.matrix() ) );
			    assert( !Stuff::FunctionContainsNanOrInf( H1rhs ) );
			    assert( !Stuff::FunctionContainsNanOrInf( H2rhs ) );
			    assert( !Stuff::FunctionContainsNanOrInf( H3rhs ) );
			    assert( !Stuff::FunctionContainsNanOrInf( H2_O_rhs ) );
	    //				Zmatrix.matrix().scale( -1 );
	    //				assert( areTransposed( Zmatrix.matrix(), Ematrix.matrix() ));
	    //				Zmatrix.matrix().scale( -1 );
		    }
		}
		#endif

		SaddlepointInverseOperatorInfo result;

		if ( Parameters().getParam( "disableSolver", false ) ) {
			Logger().Err().Resume();
			Logger().Err() << "solving disabled via parameter file" << std::endl;
			return result;
		}

		switch ( solverID ) {
			case Solver::BiCg_Saddlepoint_Solver_ID:result = BiCgSaddlepointSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case Solver::NestedCG_Solver_ID:		result = NestedCgSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case Solver::Reduced_Solver_ID:			result = ReducedSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
			case Solver::SaddlePoint_Solver_ID:		result = SaddlepointSolverType().solve(	arg, dest,
															 X, M_invers, Y,
															 O, E, R, Z, W,
															 H1rhs, H2rhs, H3rhs );
											break;
//			case Solver::FullSystem_Solver_ID:		result = FullsytemSolverType().solve(	arg, dest,
//															 X, M_invers, Y,
//															 O, E, R, Z, W,
//															 H1rhs, H2rhs, H3rhs );
//											break;
			default: throw std::runtime_error("invalid Solver ID selected");
		}
		if ( rhs_datacontainer ) {
			rhs_datacontainer->clear();
			ReconstructionPolicyType<DataContainerType,typename StokesPassType::DiscreteModelType>
					::reconstruct(	*rhs_datacontainer, dest, beta,
									X, M_invers, Y,
									O, E, R, Z, W,
									H1rhs, H2rhs, H3rhs );
			if (!with_oseen_discretization)
				rhs_datacontainer->convection.clear();
			if (solverID == Solver::Reduced_Solver_ID)
				rhs_datacontainer->pressure_gradient.clear();
		}
		return result;
	}
};

} //namespace Stokes
} //namespace Dune

#endif // DUNE_STOKES_SOLVERCALLER_HH
