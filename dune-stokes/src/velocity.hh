/** \file velocity.hh
    \brief contains a class Velocity with traitsclass VelocityTraits
 **/

#ifndef VELOCITY_HH
#define VELOCITY_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by Velicity
 *  \tparam int gridDim dimension of the grid
 **/
template < int gridDim >
class VelocityTraits
{
    public:
        typedef Dune::FieldVector< double, gridDim >
            DomainType;
        typedef Dune::FieldVector< double, gridDim >
            RangeType;
        typedef Dune::FieldVector< RangeType, gridDim >
            GradientRangeType;
        typedef Dune::FieldVector< double, 1 >
            DivergenceRangeType;
};

/**
 *  \brief describes the exact solution \f$u\f$ (velocity) of a stokes problem
 *
 *  in 2 dimensions: \f$u:\mathbb{R}^{2} \mapsto \Omega = \left[-1,1\right]^{2} \f$
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
 *      u_{1}(x_{1},x_{2}) := -e^{x_{1}}\Big( x_{2} cos(x_{2}) + sin(x_{2}) \Big),
 *  \f]
 *  \f[
 *      u_{2}(x_{1},x_{2}) := e^{x_{1}}x_{2}sin(x_{2}).
 *  \f]
 *
 *  \tparam int gridDim dimension of the grid
 *
 *  \todo doc test problem
 **/
template < int gridDim >
class Velocity
{
    public:
        typedef VelocityTraits< gridDim >
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::GradientRangeType
            GradientRangeType;
        typedef typename Traits::DivergenceRangeType
            DivergenceRangeType;

        /**
         *  \brief constructor
         *  doing nothing
         **/
        Velocity()
        {
        }

        /**
         *  \brief  destructor
         *  doing nothing
         **/
        ~Velocity()
        {
        }

        /**
         *  \brief evaluates the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of velocity at point arg
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of velocity at point arg
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

        /**
         *  \brief evaluates the gradient of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the gradient of the velocity at point arg
         **/
        inline void gradient( const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  evaluates the divergence of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the divergence of the velocity at point arg
         **/
        inline void divergence( const DomainType& arg, DivergenceRangeType& ret ) const;

        /**
         *  \brief  evaluates the laplacian of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the laplacian of the velocity at point arg
         **/
        inline void laplacian( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;
};

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
inline void Velocity< 2 >::evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double sin_of_x2 = std::sin( x2 );
    // return
    ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
    ret[1] = exp_of_x1 * x2 * sin_of_x2;
}

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
inline void Velocity< 2 >::gradient( const DomainType& arg, GradientRangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    assert( ret[0].dim() == 2 );
    assert( ret[1].dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double sin_of_x2 = std::sin( x2 );
    double cos_of_x2 = std::cos( x2 );
    // gradient of u_{1}
    double du1_over_dx1 = -1.0 * exp_of_x1 *
        ( x2 * cos_of_x2 + sin_of_x2 );
    double du1_over_dx2 = exp_of_x1 *
        ( x2 * sin_of_x2 - 2.0 * cos_of_x2 );
    RangeType grad_u1;
    grad_u1[0] = du1_over_dx1;
    grad_u1[1] = du1_over_dx2;
    // gradient of u_{2}
    double du2_over_dx1 = exp_of_x1 * x2 * sin_of_x2;
    double du2_over_dx2 = exp_of_x1 *
        ( x2 * cos_of_x2 + sin_of_x2 );
    RangeType grad_u2;
    grad_u2[0] = du2_over_dx1;
    grad_u2[1] = du2_over_dx2;
    // return
    ret[0] = grad_u1;
    ret[1] = grad_u2;
}

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
inline void Velocity< 2 >::divergence( const DomainType& arg, DivergenceRangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 1 );
    // return
    ret[0] = 0.0;
}

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
inline void Velocity< 2 >::laplacian( const DomainType& arg, RangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double cos_of_x2 = std::cos( x2 );
    // return
    ret[0] = -1.0 * exp_of_x1 *
        ( ( 2.0 * x2 * cos_of_x2 ) - ( 2.0 * std::sin( x2 ) ) );
    ret[1] = 2.0 * exp_of_x1 * cos_of_x2;
}

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
void Velocity< 2 >::testMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "- testing class Velocity..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "  - x: " << x[0] << std::endl;
    debugStream << "       " << x[1] << std::endl;
    RangeType u;
    evaluate( x, u );
    debugStream << "  - u(x): " << u[0] << std::endl;
    debugStream << "          " << u[1] << std::endl;
    GradientRangeType grad_u;
    gradient( x, grad_u );
    debugStream << "  - grad u(x): " << grad_u[0] << std::endl;
    debugStream << "               " << grad_u[1] << std::endl;
    DivergenceRangeType div_u;
    divergence( x, div_u );
    debugStream << "  - div u(x): " << div_u[0] << std::endl;
    RangeType laplace_u;
    laplacian( x, laplace_u );
    debugStream << "  - laplacian u(x): " << laplace_u[0] << std::endl;
    debugStream <<  "                  " << laplace_u[1] << std::endl;
    infoStream << "  ...test passed!" << std::endl;
}

#endif // end of velocity.hh
