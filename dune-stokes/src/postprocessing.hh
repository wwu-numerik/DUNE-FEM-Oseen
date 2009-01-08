/**
 *  \file   postprocessing.hh
 *  \brief  postprocessing.hh
 **/

#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <dune/fem/operator/lagrangeinterpolation.hh>

#include "logging.hh"

template <  class StokesPassImp, class ProblemImp >
class PostProcessor
{
    public:
        typedef ProblemImp
            ProblemType;

        typedef StokesPassImp
            StokesPassType;
        typedef typename StokesPassType::DiscreteStokesFunctionSpaceWrapperType
            DiscreteStokesFunctionSpaceWrapperType;
        typedef typename StokesPassType::GridPartType
            GridPartType;

        typedef typename StokesPassType::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;
        typedef typename StokesPassType::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        typedef typename StokesPassType::DiscretePressureFunctionType
            DiscretePressureFunctionType;
        typedef typename StokesPassType::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;


        PostProcessor( const StokesPassType& pass, const ProblemType& prob )
            : pass_( pass ),
            problem_( prob ),
            spaceWrapper_( pass.GetFunctionSpaceWrapper() ),
            gridPart_( spaceWrapper_.gridPart() ),
            velocitySpace_ ( spaceWrapper_.discreteVelocitySpace() ),
            discreteExactVelocity_( "u_exact", velocitySpace_ ),
            discreteExactForce_( "f_exact", velocitySpace_ ),
            discreteExactDirichlet_( "gd_exact", velocitySpace_ ),
            discreteExactPressure_( "p_exact", spaceWrapper_.discretePressureSpace() )
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
        const StokesPassType& pass_;
        const ProblemType& problem_;
        //ContinuousVelocityType& continuousVelocity_;
        const GridPartType& gridPart_;
        const DiscreteVelocityFunctionSpaceType& velocitySpace_;
        const DiscreteStokesFunctionSpaceWrapperType& spaceWrapper_;
        DiscreteVelocityFunctionType discreteExactVelocity_;
        DiscreteVelocityFunctionType discreteExactForce_;
        DiscreteVelocityFunctionType discreteExactDirichlet_;
        DiscretePressureFunctionType discreteExactPressure_;
};

#endif // end of postprocessing.hh
