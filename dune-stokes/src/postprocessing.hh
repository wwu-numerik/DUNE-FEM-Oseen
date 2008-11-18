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
            class DiscreteVelocityFunctionSpaceImp,
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
        typedef DiscreteVelocityFunctionSpaceImp
            DiscreteVelocityFunctionSpaceType;
        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;

        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;


        PostProcessor( const ProblemType& problem, const GridPartType& gridPart, const DiscreteVelocityFunctionSpaceType& velo_space )
            : problem_( problem ),
            gridPart_( gridPart ),
            discreteExactVelocity_( "u_exact", velo_space ),
            velocitySpace_ ( velo_space )
        {

        }

        ~PostProcessor()
        {
        }

        void assembleExactSolution()
        {

            //Dune::LagrangeInterpolation< Problem::VelocityType >::interpolateFunction( , exactVelocity_ );

        }

    private:
        const ProblemType& problem_;
        //ContinuousVelocityType& continuousVelocity_;
        const GridPartType& gridPart_;
        const DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscreteVelocityFunctionType discreteExactVelocity_;
};

#endif // end of postprocessing.hh
