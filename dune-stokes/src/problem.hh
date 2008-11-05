/** \file problem.hh
    \brief problem.hh
 **/

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>
//#include <assert>

#include "dune/common/fvector.hh"

#include "logging.hh"

/**
 *  \brief describes the velocity
 *  \f[u(x_{1},x_{2}):=-e^{x_{1}}\Big( x_{2} cos(x_{2}) + sin(x_{2}) \Big)\f]
 *
 *  and the gradient of the velocity
 *
 *  as the solution of the stokes test problem
 *  \todo doc test problem
 **/
template < int grid_dim >
class Velocity
{
    public:
        typedef Dune::FieldVector< double, grid_dim > DomainType;
        typedef Dune::FieldVector< double, grid_dim > RangeType;
        typedef Dune::FieldVector< RangeType, grid_dim > GradientRangeType;
        typedef Dune::FieldVector< double, 1 > DivergenceRangeType;

        /**
         *  \brief constructor
         **/
        Velocity()
            : infoStream_( Logger().Info() ),
            debugStream_( Logger().Dbg() ),
            errorStream_( Logger().Err() )
        {
        }

        /**
         *  \brief  destructor
         **/
        ~Velocity()
        {
        }

        /**
         *  \brief evaluates the velocity
         *  \arg DomainType& arg
         *  \arg RangeType& ret returns velocity at arg
         **/
        void Evaluate( const DomainType& arg, RangeType& ret );


        /**
         *  \brief evaluates the gradient of the velocity
         **/
        void Gradient( const DomainType& arg, GradientRangeType& ret );

        /**
         *  \brief  evaluates the divergence of the velocity
         **/
        void Divergence( const DomainType& arg, DivergenceRangeType& ret );

        /**
         *  \brief  evaluates the laplacian of the velocity
         **/
        void Laplacian( const DomainType& arg, RangeType& ret );

    private:
        Logging::LogStream& infoStream_;
        Logging::LogStream& debugStream_;
        Logging::LogStream& errorStream_;
};

template < >
void Velocity< 2 >::Evaluate( const DomainType& arg, RangeType& ret )
{
    assert( arg.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double sin_of_x2 = std::sin( x2 );
    // return
    ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
    ret[1] = exp_of_x1 * x2 * sin_of_x2;
}

template < >
void Velocity< 2 >::Gradient( const DomainType& arg, GradientRangeType& ret )
{
    assert( arg.dim() == 2 );
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
};

template < >
void Velocity< 2 >::Divergence( const DomainType& arg, DivergenceRangeType& ret )
{
    assert( arg.dim() == 2 );
    ret[0] = 0.0;
    ret[1] = 0.0;
}

template < >
void Velocity< 2 >::Laplacian( const DomainType& arg, RangeType& ret )
{
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double cos_of_x2 = std::cos( x2 );
    // return
    ret[0] = -1.0 * exp_of_x1 *
        ( ( 2.0 * x2 * cos_of_x2 ) - ( 2.0 * std::sin( x2 ) ) );
    ret[1] = 2.0 * exp_of_x1 * cos_of_x2;
};


#endif  // end of problem.hh

