/**
 *  \file   postprocessing.hh
 *  \brief  postprocessing.hh
 **/

#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <dune/fem/operator/lagrangeinterpolation.hh>

#include "logging.hh"
#include "problem.hh"

template <  class ProblemImp,
            class GridPartImp,
            class DiscreteVelocityFunctionImp,
            class DiscretePressureFunctionImp >
class PostProcessor
{
    public:
        typedef ProblemImp
            ProblemType;
        typedef typename ProblemType::VelocityType
            ContinuousVelocityType;
        typedef GridPartImp
            GridPartType;

        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;



        PostProcessor( const ProblemType& problem, const GridPartType& gridPart,
                        const DiscreteVelocityFunctionSpaceType& velo_space,
                        const DiscretePressureFunctionSpaceType& press_space)
            : problem_( problem ),
            gridPart_( gridPart ),
            velocitySpace_ ( velo_space ),
            discreteExactVelocity_( "u_exact", velo_space ),
            discreteExactForce_( "f_exact", velo_space ),
            discreteExactDirichlet_( "gd_exact", velo_space ),
            discreteExactPressure_( "p_exact", press_space )
        {

        }

        ~PostProcessor()
        {
        }

        void assembleExactSolution()
        {
            Dune::LagrangeInterpolation< DiscreteVelocityFunctionType >::interpolateFunction( problem_.velocity(), discreteExactVelocity_ );
            Dune::LagrangeInterpolation< DiscreteVelocityFunctionType >::interpolateFunction( problem_.dirichletData(), discreteExactDirichlet_ );
            Dune::LagrangeInterpolation< DiscreteVelocityFunctionType >::interpolateFunction( problem_.force(), discreteExactForce_ );
            Dune::LagrangeInterpolation< DiscretePressureFunctionType >::interpolateFunction( problem_.velocity(), discreteExactPressure_ );
        }

    private:
        const ProblemType& problem_;
        //ContinuousVelocityType& continuousVelocity_;
        const GridPartType& gridPart_;
        const DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscreteVelocityFunctionType discreteExactVelocity_;
        DiscreteVelocityFunctionType discreteExactForce_;
        DiscreteVelocityFunctionType discreteExactDirichlet_;
        DiscretePressureFunctionType discreteExactPressure_;
};

#endif // end of postprocessing.hh
