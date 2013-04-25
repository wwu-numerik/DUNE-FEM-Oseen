/**
 *  \file   velocity.hh
 *
 *  \brief  contains a class Velocity with traitsclass VelocityTraits
 **/

#ifndef VELOCITY_HH
#define VELOCITY_HH

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/fem/function/common/function.hh>

#include <dune/stuff/common/logging.hh>

/**
 *  \brief  containing typedefs needed by Velocity
 *
 *  \tparam gridDim
 *          dimension of the grid (unused)
 *  \tparam VelocityFunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int gridDim, class VelocityFunctionSpaceImp >
class VelocityTraits
{
    public:
        typedef VelocityFunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType
            GradientRangeType;
        typedef typename FunctionSpaceType::HessianRangeType
            DivergenceRangeType;
};

/**
 *  \brief describes the velocity \f$u\f$ as an exact solution of a stokes
 *  problem
 *
 *  in 2 dimensions: \f$u:\mathbb{R}^{2}\mapsto\Omega =\left[-1,1\right]^{2} \f$
 *  \f[
 *      u(x) = \left(
 *          \begin{array}{c}
 *              u_{1}(x)\\
 *              u_{2}(x)
 *          \end{array}
 *      \right),
 *  \f]
 *  where
 *  \f[
 *      u_{1}(x_{1},x_{2}):= -e^{x_{1}}\Big( x_{2} cos(x_{2}) + sin(x_{2})\Big),
 *  \f]
 *  \f[
 *      u_{2}(x_{1},x_{2}) := e^{x_{1}}x_{2}sin(x_{2}).
 *  \f]
 *
 *  \tparam TraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class TraitsImp  >
class Velocity : public Dune::Fem::Function < typename TraitsImp::FunctionSpaceType , Velocity < TraitsImp > >
{
    public:
        typedef TraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::GradientRangeType
            GradientRangeType;
        typedef typename Traits::DivergenceRangeType
            DivergenceRangeType;
        typedef typename Traits::FunctionSpaceType
            FunctionSpaceType;
        typedef Dune::Fem::Function < typename TraitsImp::FunctionSpaceType , Velocity < TraitsImp > >
            BaseType;

        /**
         *  \brief constructor
         *
         *  doing nothing besides Base init
         **/
        Velocity( const FunctionSpaceType& f_space )
            : BaseType( f_space ),
            dim_( FunctionSpaceType::dimDomain )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~Velocity()
        {}

        /**
         *  \brief  evaluates the velocity
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of velocity at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            if ( dim_ == 1 ) {
                assert( !"velocity not implemented in 1D" );
            }
            else if ( dim_ == 2 ) {
                double x1 = arg[0];
                double x2 = arg[1];
#if defined(CONSTANT_PROBLEM)
                ret[0] = 0;
                ret[1] = 0;
#elif defined(GENERALIZED_STOKES_PROBLEM)
                const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
                ret[0] = tmp;
                ret[1] = -1.0 * tmp;
#else
                double exp_of_x1 = std::exp( x1 );
                double sin_of_x2 = std::sin( x2 );
                // return
                ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
                ret[1] = exp_of_x1 * x2 * sin_of_x2;
#endif
            }
            else if ( dim_ == 3 ) {
//                double x1 = arg[0];
//                double x2 = arg[1];
//                double x3 = arg[2];
#ifdef CONSTANT_PROBLEM
                ret[0] = 0;
                ret[1] = 0;
                ret[2] = 0;
#elif defined(GENERALIZED_STOKES_PROBLEM)
                assert( !"GENERALIZED_STOKES_PROBLEM not implemented in 3D" );
#elif defined(AORTA_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#else
                assert( !"velocity not implemented in 3D" );
#endif
            }
            else {
                assert( !"velocity not implemented for more than 3 dimensions" );
            }

        }

        /**
         *  \brief  evaluates the velocity
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of velocity at given point
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;

    private:
        const int dim_;
};

#endif // end of velocity.hh

/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
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

