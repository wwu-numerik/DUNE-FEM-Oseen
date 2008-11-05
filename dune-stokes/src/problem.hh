/** \file problem.hh
    \brief descibes analytical functions, that solve a stokes problem in 2d
 **/

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>
//#include <assert>

#include "dune/common/fvector.hh"

/**
 *  \brief describes the velocity
 *  \tparam int grid_dim dimension of the grid
 *
 *  \f[u(x_{1},x_{2}):=-e^{x_{1}}\Big( x_{2} cos(x_{2}) + sin(x_{2}) \Big)\f]
 *
 *  \todo doc test problem
 **/
template < int grid_dim >
class Velocity
{
    public:
        typedef Dune::FieldVector< double, grid_dim >
            DomainType;
        typedef Dune::FieldVector< double, grid_dim >
            RangeType;
        typedef Dune::FieldVector< RangeType, grid_dim >
            GradientRangeType;
        typedef Dune::FieldVector< double, 1 >
            DivergenceRangeType;

        /**
         *  \brief constructor
         **/
        Velocity()
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
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of velocity at point arg
         **/
        inline void Evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of velocity at point arg
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            Evaluate( arg, ret );
            return ret;
        };

        /**
         *  \brief evaluates the gradient of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the gradient of the velocity at point arg
         **/
        inline void Gradient( const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  evaluates the divergence of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the divergence of the velocity at point arg
         **/
        inline void Divergence( const DomainType& arg, DivergenceRangeType& ret ) const;

        /**
         *  \brief  evaluates the laplacian of the velocity
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of the laplacian of the velocity at point arg
         **/
        inline void Laplacian( const DomainType& arg, RangeType& ret ) const;

    private:
};

template < >
inline void Velocity< 2 >::Evaluate( const DomainType& arg, RangeType& ret ) const
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

template < >
inline void Velocity< 2 >::Gradient( const DomainType& arg, GradientRangeType& ret ) const
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
};

template < >
inline void Velocity< 2 >::Divergence( const DomainType& arg, DivergenceRangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // return
    ret[0] = 0.0;
    ret[1] = 0.0;
}

template < >
inline void Velocity< 2 >::Laplacian( const DomainType& arg, RangeType& ret ) const
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
};

/**
 *
 *  \brief  describes the presure
 *  \tparam int grid_dim dimension of the grid
 *
 *  \todo   doc
 **/
template < int grid_dim >
class Pressure
{
    public:
        typedef Dune::FieldVector< double, grid_dim >
            DomainType;
        typedef Dune::FieldVector< double, 1 >
            RangeType;
        typedef Dune::FieldVector< double, grid_dim >
            GradientRangeType;

        /**
         *  \brief  constructor
         **/
        Pressure()
        {
        }

        /**
         *  \brief  destructor
         **/
        ~Pressure()
        {
        }

        /**
         *  \brief evaluates the pressure
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of pressure at point arg
         **/
        inline void Evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the pressure
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of pressure at point arg
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            Evaluate( arg, ret );
            return ret;
        };

        /**
         *  \brief  evaluates the gradient of the pressure
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of gradient of the pressure at point arg
         **/
        inline void Gradient( const DomainType& arg, GradientRangeType& ret ) const;

    private:
 };

template < >
inline void Pressure< 2 >::Evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play save
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    // return
    ret[0] = 2.0 * std::exp( x1 ) * std::sin( x2 );
};

template < >
inline void Pressure< 2 >::Gradient( const DomainType& arg, GradientRangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    //return
    ret[0] = 2.0 * exp_of_x1 * std::sin( x2 );
    ret[1] = 2.0 * exp_of_x1 * std::cos( x2 );
};

/**
 *  \brief  describes the force
 *  \tparam int grid_dim dimension of the grid
 *
 *  \todo doc
 **/
template < int grid_dim >
class Force
{
    public:
        typedef Dune::FieldVector< double, grid_dim >
            DomainType;
        typedef Dune::FieldVector< double, grid_dim >
            RangeType;

        /**
         *  \brief  constructor
         *
         *  doing nothing
         **/
        Force()
        {
        }

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~Force()
        {
        }

        /**
         *  \brief  evaluates the force
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of force at point arg
         **/
        inline void Evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the force
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of force at point arg
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            Evaluate( arg, ret );
            return ret;
        };

    private:
};

template < >
inline void Force< 2 >::Evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    //return
    ret[0] = 2.0 * std::exp( x1 ) * std::cos( x2 );
    ret[1] = 0.0;
};

/**
 *  \brief  describes the dirichlet boundary data
 *  \tparam int grid_dim dimension of the grid
 *
 *  \todo doc
 **/
template < int grid_dim >
class DirichletData
{
    public:
        typedef Dune::FieldVector< double, grid_dim >
            DomainType;
        typedef Dune::FieldVector< double, grid_dim >
            RangeType;

        /**
         *  \brief  constructor
         *
         *  doing nothing
         **/
        DirichletData()
        {
        }

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
         ~DirichletData()
         {
         }

         /**
          * \brief  evaluates the dirichlet data
         *  \arg DomainType& arg point to be evaluated at
         *  \arg RangeType& ret value of dirichlet boundary data at point arg
          **/
        inline void Evaluate( DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the dirichlet data
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of dirichlet data at point arg
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            Evaluate( arg, ret );
            return ret;
        };

        private:
};

template < >
inline void DirichletData< 2 >::Evaluate( DomainType& arg, RangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double sin_of_x2 = std::sin( x2 );
    //return
    ret[0] = -1.0 * exp_of_x1 *
        ( ( x2 * std::cos( x2 ) ) + sin_of_x2 );
    ret[1] = exp_of_x1 * x2 * sin_of_x2;
};

#endif  // end of problem.hh
