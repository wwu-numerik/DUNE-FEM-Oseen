#ifndef DUNE_OSEEN_SOLVERCALLER_HH
#define DUNE_OSEEN_SOLVERCALLER_HH

#include <dune/fem/oseen/solver/nested_cg.hh>
#include <dune/fem/oseen/solver/reduced.hh>
#include <dune/fem/oseen/solver/fullsystem.hh>
#include <dune/fem/oseen/solver/saddle_point.hh>
#include <dune/fem/oseen/solver/bicg_saddle_point.hh>
#include <dune/fem/oseen/solver/reconstruction.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>

namespace Dune {

namespace Oseen {

namespace Solver {
    enum SolverID {
        NestedCG_Solver_ID			= 0,
        SaddlePoint_Solver_ID		= 1,
        Reduced_Solver_ID			= 2,
        FullSystem_Solver_ID		= 3,
        BiCg_Saddlepoint_Solver_ID	= 4
    };
}

template<class OseenPassType, template <class S> class ReconstructionPolicyType = BruteForceReconstruction >
struct SolverCaller {
	//! alternative solver implementation
    typedef NestedCgSaddlepointInverseOperator< OseenPassType >
        NestedCgSolverType;
	//! type of the used solver
    typedef SaddlepointInverseOperator< OseenPassType >
		SaddlepointSolverType;
    typedef BiCgStabSaddlepointInverseOperator< OseenPassType >
        BiCgSaddlepointSolverType;
	//! this is used for reduced (no pressure, incompress. condition) oseen pass
    typedef ReducedInverseOperator< OseenPassType >
		ReducedSolverType;
    typedef DirectKrylovSolver< OseenPassType >
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
    static SaddlepointInverseOperatorInfo solve(
                const Solver::SolverID solverID,
                const bool with_oseen_discretization,
                DataContainerType* rhs_datacontainer,
                RangeType& dest,
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
		DSC::Profiler::ScopedTiming solver_time("solver");
		MatrixWrapper<MInversMatrixObjectType> M_invers(MInversMatrix, "M");
		MatrixWrapper<WmatrixObjectType> W(Wmatrix, "W");
		MatrixWrapper<XmatrixObjectType> X(Xmatrix, "X");
		MatrixWrapper<YmatrixObjectType> Y(Ymatrix, "Y");
		MatrixWrapper<OmatrixObjectType> O(Omatrix, "O");
		MatrixWrapper<ZmatrixObjectType> Z(Zmatrix, "Z");
		MatrixWrapper<EmatrixObjectType> E(Ematrix, "E");
		MatrixWrapper<RmatrixObjectType> R(Rmatrix, "R");

        #ifndef NDEBUG
		{
		    // do the matlab logging stuff
		    if ( DSC_CONFIG_GET( "save_matrices", false ) ) {
                auto& matlabLogStream = DSC_LOG_ERROR;
                #   define MPRINTER printSparseRowMatrixMatlabStyle
                DSC::MPRINTER( MInversMatrix->matrix(), "M_invers", matlabLogStream );
                DSC::MPRINTER( Wmatrix->matrix(), "W", matlabLogStream );
                DSC::MPRINTER( Omatrix->matrix(), "O", matlabLogStream );
                DSC::MPRINTER( Xmatrix->matrix(), "X", matlabLogStream );
                DSC::MPRINTER( Ymatrix->matrix(), "Y", matlabLogStream );
                DSC::MPRINTER( Zmatrix->matrix(), "Z", matlabLogStream );
                DSC::MPRINTER( Ematrix->matrix(), "E", matlabLogStream );
                DSC::MPRINTER( Rmatrix->matrix(), "R", matlabLogStream );
                DSC::printDiscreteFunctionMatlabStyle( H1rhs, "H1", matlabLogStream );
                DSC::printDiscreteFunctionMatlabStyle( H2rhs, "H2", matlabLogStream );
                DSC::printDiscreteFunctionMatlabStyle( H3rhs, "H3", matlabLogStream );
                matlabLogStream.flush();
                #undef MPRINTER

		    }

		    if ( DSC_CONFIG_GET( "outputMatrixPlots", false ) ) {
                DSC::matrixToGnuplotFile( Ematrix->matrix(),       std::string( "mat_E.gnuplot")       );
                DSC::matrixToGnuplotFile( Wmatrix->matrix(),       std::string( "mat_W.gnuplot")       );
                DSC::matrixToGnuplotFile( Xmatrix->matrix(),       std::string( "mat_X.gnuplot")       );
                DSC::matrixToGnuplotFile( Ymatrix->matrix(),       std::string( "mat_Y.gnuplot")       );
                DSC::matrixToGnuplotFile( Zmatrix->matrix(),       std::string( "mat_Z.gnuplot")       );
                DSC::matrixToGnuplotFile( Rmatrix->matrix(),       std::string( "mat_R.gnuplot")       );
                DSC::matrixToGnuplotFile( MInversMatrix->matrix(), std::string( "mat_M.gnuplot")   );
		    }

		    if ( DSC_CONFIG_GET( "paranoid_checks", false ) )
		    {//paranoid checks
                assert( !DSFe::MatrixContainsNanOrInf( Omatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Ematrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Wmatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Xmatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Ymatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Zmatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( Rmatrix->matrix() ) );
                assert( !DSFe::MatrixContainsNanOrInf( MInversMatrix->matrix() ) );
                assert( !DSFe::FunctionContainsNanOrInf( H1rhs ) );
                assert( !DSFe::FunctionContainsNanOrInf( H2rhs ) );
                assert( !DSFe::FunctionContainsNanOrInf( H3rhs ) );
//			    assert( !DSFe::FunctionContainsNanOrInf( H2_O_rhs ) );
	    //				Zmatrix.matrix().scale( -1 );
	    //				assert( areTransposed( Zmatrix.matrix(), Ematrix.matrix() ));
	    //				Zmatrix.matrix().scale( -1 );
		    }
		}
        #endif

		SaddlepointInverseOperatorInfo result;

		if ( DSC_CONFIG_GET( "disableSolver", false ) ) {
			DSC_LOG_ERROR.resume();
			DSC_LOG_ERROR << "solving disabled via parameter file" << std::endl;
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

            case Solver::FullSystem_Solver_ID:		result = FullsytemSolverType().solve(	arg, dest,
                                                             X, M_invers, Y,
                                                             O, E, R, Z, W,
                                                             H1rhs, H2rhs, H3rhs );
                                            break;
            default:
                throw std::runtime_error("invalid Solver ID selected");
		}
		if ( rhs_datacontainer ) {
			rhs_datacontainer->clear();
            ReconstructionPolicyType<typename OseenPassType::DiscreteModelType>
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

template<class OseenPassType >
struct SolverCallerProxy {
    template < class RangeType, class ContainerType, class... Args >
    static SaddlepointInverseOperatorInfo call( const bool do_oseen_discretization,
                                                ContainerType* container,
                                                RangeType& dest,
                                                const Args&... args )
    {
        //this lets us switch between standalone oseen and reduced oseen in  thete scheme easily
        const bool use_reduced_solver = (do_oseen_discretization && DSC_CONFIG_GET( "reduced_oseen_solver", false ))
                || DSC_CONFIG_GET( "parabolic", false );
        typedef Oseen::SolverCaller< OseenPassType>
            SolverCallerType;
        typedef Oseen::SolverCaller< OseenPassType, SmartReconstruction >
            SmartSolverCallerType;

        //Select which solver we want to use
        typename Oseen::Solver::SolverID solver_ID = do_oseen_discretization
                ? Oseen::Solver::BiCg_Saddlepoint_Solver_ID
                : Oseen::Solver::SaddlePoint_Solver_ID;

        if( !use_reduced_solver ) {
            if ( DSC_CONFIG_GET( "use_nested_cg_solver", false ) )
                solver_ID = Oseen::Solver::NestedCG_Solver_ID;
            else if ( DSC_CONFIG_GET( "use_full_solver", false ) )
                solver_ID = Oseen::Solver::FullSystem_Solver_ID;
        }
        else
            solver_ID = Oseen::Solver::Reduced_Solver_ID;

        if ( DSC_CONFIG_GET( "smart_reconstruction", false ) )
            return SmartSolverCallerType::solve(solver_ID,
                                                do_oseen_discretization,
                                                container,
                                                dest,
                                                args...);
        else
            return SolverCallerType::solve(solver_ID,
                                           do_oseen_discretization,
                                           container,
                                           dest,
                                           args...);
    }
};


} //namespace Oseen
} //namespace Dune

#endif // DUNE_OSEEN_SOLVERCALLER_HH

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
