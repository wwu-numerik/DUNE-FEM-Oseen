/**
 *  \file   velocity.hh
 *
 *  \brief  contains a class Velocity with traitsclass VelocityTraits
 **/

#ifndef VELOCITY_HH
#define VELOCITY_HH

#include <cmath>

#include <dune/common/fvector.hh>


#include <dune/stuff/logging.hh>

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
//        typedef Dune::FieldVector< RangeType, gridDim >
            GradientRangeType;
        typedef typename FunctionSpaceType::HessianRangeType
//        typedef Dune::FieldVector< double, 1 >
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
class Velocity : public Dune::Function < typename TraitsImp::FunctionSpaceType , Velocity < TraitsImp > >
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
        typedef Dune::Function < typename TraitsImp::FunctionSpaceType , Velocity < TraitsImp > >
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
#ifdef SIMPLE_PROBLEM
                ret[0] = 1;
                ret[1] = 0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0;
                ret[1] = 0;
#elif defined(GENRALIZED_STOKES_PROBLEM)
                const double tmp = std::cos( ( 0.5 * M_PI ) * ( x1 + x2 ) );
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
#ifdef SIMPLE_PROBLEM
                assert( !"SIMPLE_PROBLEM not implemented in 3D" );
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0;
                ret[1] = 0;
                ret[2] = 0;
#elif defined(GENRALIZED_STOKES_PROBLEM)
                assert( !"GENRALIZED_STOKES_PROBLEM not implemented in 3D" );
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
         *  \brief  evaluates the gradient of the velocity
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of the gradient of the velocity at given point
         **/
//        inline void gradient( const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  evaluates the divergence of the velocity
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of the divergence of the velocity at given point
         **/
//        inline void divergence( const DomainType& arg, DivergenceRangeType& ret ) const;

        /**
         *  \brief  evaluates the laplacian of the velocity
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of the laplacian of the velocity at given point
         **/
//        inline void laplacian( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;

    private:
        const int dim_;
};

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class TraitsImp >
//inline void Velocity< TraitsImp >::evaluate(
//        const DomainType& arg,
//        RangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
//    assert( ret.dim() == 2 );
//    double x1 = arg[0];
//    double x2 = arg[1];
//#ifdef SIMPLE_PROBLEM
//    ret[0] = 1;
//    ret[1] = 0;
//#elif defined(CONSTANT_PROBLEM)
//    ret[0] = 0;
//    ret[1] = 0;
//#else
//    double exp_of_x1 = std::exp( x1 );
//    double sin_of_x2 = std::sin( x2 );
//    // return
//    ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
//    ret[1] = exp_of_x1 * x2 * sin_of_x2;
//
//#endif
//
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class TraitsImp >
//inline void Velocity< TraitsImp >::gradient( const DomainType& arg, GradientRangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
////    assert( ret.dim() == 2 );
//    assert( ret[0].dim() == 2 );
//    assert( ret[1].dim() == 2 );
//    // some computations
//    double x1 = arg[0];
//    double x2 = arg[1];
//    double exp_of_x1 = std::exp( x1 );
//    double sin_of_x2 = std::sin( x2 );
//    double cos_of_x2 = std::cos( x2 );
//    // gradient of u_{1}
//    double du1_over_dx1 = -1.0 * exp_of_x1 *
//        ( x2 * cos_of_x2 + sin_of_x2 );
//    double du1_over_dx2 = exp_of_x1 *
//        ( x2 * sin_of_x2 - 2.0 * cos_of_x2 );
//    RangeType grad_u1;
//    grad_u1[0] = du1_over_dx1;
//    grad_u1[1] = du1_over_dx2;
//    // gradient of u_{2}
//    double du2_over_dx1 = exp_of_x1 * x2 * sin_of_x2;
//    double du2_over_dx2 = exp_of_x1 *
//        ( x2 * cos_of_x2 + sin_of_x2 );
//    RangeType grad_u2;
//    grad_u2[0] = du2_over_dx1;
//    grad_u2[1] = du2_over_dx2;
//    // return
//    ret[0] = grad_u1;
//    ret[1] = grad_u2;
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class TraitsImp >
//inline void Velocity< TraitsImp >::divergence( const DomainType& arg, DivergenceRangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
//    assert( ret.dim() == 1 );
//    // return
//    ret[0] = 0.0;
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class TraitsImp >
//inline void Velocity< TraitsImp >::laplacian( const DomainType& arg, RangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
//    assert( ret.dim() == 2 );
//    // some computations
//    double x1 = arg[0];
//    double x2 = arg[1];
//    double exp_of_x1 = std::exp( x1 );
//    double cos_of_x2 = std::cos( x2 );
//    // return
//    ret[0] = -1.0 * exp_of_x1 *
//        ( ( 2.0 * x2 * cos_of_x2 ) - ( 2.0 * std::sin( x2 ) ) );
//    ret[1] = 2.0 * exp_of_x1 * cos_of_x2;
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class TraitsImp >
//void Velocity< TraitsImp >::testMe() const
//{
//    // some logstreams
//    Logging::LogStream& infoStream = Logger().Info();
//    Logging::LogStream& debugStream = Logger().Dbg();
//    infoStream << "- testing class Velocity..." << std::endl;
//    //tests
//    DomainType x;
//    x[0] = 1.0;
//    x[1] = 1.0;
//    debugStream << "  - x: " << x[0] << std::endl;
//    debugStream << "       " << x[1] << std::endl;
//    RangeType u;
//    evaluate( x, u );
//    debugStream << "  - u(x): " << u[0] << std::endl;
//    debugStream << "          " << u[1] << std::endl;
////    GradientRangeType grad_u;
////    gradient( x, grad_u );
////    debugStream << "  - grad u(x): " << grad_u[0] << std::endl;
////    debugStream << "               " << grad_u[1] << std::endl;
////    DivergenceRangeType div_u;
////    divergence( x, div_u );
////    debugStream << "  - div u(x): " << div_u[0] << std::endl;
////    RangeType laplace_u;
////    laplacian( x, laplace_u );
////    debugStream << "  - laplacian u(x): " << laplace_u[0] << std::endl;
////    debugStream <<  "                  " << laplace_u[1] << std::endl;
//    infoStream << "  ...test passed!" << std::endl;
//}

#endif // end of velocity.hh
