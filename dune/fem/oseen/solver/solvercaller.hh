#ifndef DUNE_OSEEN_SOLVERCALLER_HH
#define DUNE_OSEEN_SOLVERCALLER_HH

#include <dune/fem/oseen/solver/reduced.hh>
#include <dune/fem/oseen/solver/saddle_point.hh>
#include <dune/fem/oseen/solver/bicg_saddle_point.hh>
#include <dune/fem/oseen/solver/reconstruction.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>

namespace Dune {

namespace Oseen {

namespace Solver {
    enum SolverID {
        SaddlePoint_Solver_ID		= 1,
        Reduced_Solver_ID			= 2,
        BiCg_Saddlepoint_Solver_ID	= 4
    };
}

template<class OseenPassType, template <class S> class ReconstructionPolicyType = BruteForceReconstruction >
struct SolverCaller {
	//! type of the used solver
    typedef SaddlepointInverseOperator< OseenPassType >
		SaddlepointSolverType;
    typedef BiCgStabSaddlepointInverseOperator< OseenPassType >
        BiCgSaddlepointSolverType;
	//! this is used for reduced (no pressure, incompress. condition) oseen pass
    typedef ReducedInverseOperator< OseenPassType >
		ReducedSolverType;


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
                const XmatrixObjectType& X,
                const MInversMatrixObjectType& M_invers,
                const YmatrixObjectType& Y,
                const OmatrixObjectType& O,
                const EmatrixObjectType& E,
                const RmatrixObjectType& R,
                const ZmatrixObjectType& Z,
                const WmatrixObjectType& W,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& H2rhs,
				const DiscretePressureFunctionType& H3rhs,
				const DiscreteVelocityFunctionType& beta )
	{
		DSC::Profiler::ScopedTiming solver_time("solver");

		SaddlepointInverseOperatorInfo result;

		switch ( solverID ) {
            case Solver::BiCg_Saddlepoint_Solver_ID:result = BiCgSaddlepointSolverType().solve(	arg, dest,
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
        auto solver_ID = do_oseen_discretization
                ? Oseen::Solver::BiCg_Saddlepoint_Solver_ID
                : Oseen::Solver::SaddlePoint_Solver_ID;

        if(use_reduced_solver)
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

