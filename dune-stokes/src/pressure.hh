/**
 *  \file   pressure.hh
 *
 *  \brief  contains a class Pressure with traitsclass PressureTraits
 **/

#ifndef PRESSURE_HH
#define PRESSURE_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by Pressure
 *
 *  \tparam gridDim (unused)
 *          dimension of the grid
 *
 *  \tparam PressureFunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int gridDim, class PressureFunctionSpaceImp >
class PressureTraits
{
    public:
        typedef PressureFunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType
        //typedef Dune::FieldVector< double, gridDim >
            GradientRangeType;

};

/**
 *  \brief  describes the presure
 *
 *  \tparam PressureTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class PressureTraitsImp >
class Pressure : public Dune::Function < typename PressureTraitsImp::FunctionSpaceType, Pressure < PressureTraitsImp > >
{
    public:
        typedef PressureTraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::GradientRangeType
            GradientRangeType;
        typedef typename Traits::FunctionSpaceType
            PressureFunctionSpaceType;
        typedef Dune::Function < typename PressureTraitsImp::FunctionSpaceType , Pressure < PressureTraitsImp > >
            BaseType;
        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        Pressure( const PressureFunctionSpaceType& press_space )
            : BaseType( press_space )
        {
        }

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~Pressure()
        {
        }

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of pressure at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of pressure at given point
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

        /**
         *  \brief  evaluates the gradient of the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of gradient of the pressure at given point
         **/
//        inline void gradient(const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;
 };

/**
 *  \brief  specialization for gridDim = 2
 **/
template < class PressureTraitsImp  >
inline void Pressure< PressureTraitsImp  >::evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play save
    assert( arg.dim() == 2 );
    assert( ret.dim() == 1 );
    double x1 = arg[0];
    double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
    ret[0] = x1*x2;
#elif defined(CONSTANT_PROBLEM)
    ret[0] = 0;
#else
    ret[0] = 2.0 * std::exp( x1 ) * std::sin( x2 );
#endif
}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class PressureTraitsImp  >
//inline void Pressure< PressureTraitsImp >::gradient( const DomainType& arg, GradientRangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
////    assert( ret.dim() == 2 );
//    // some computations
//    double x1 = arg[0];
//    double x2 = arg[1];
//    double exp_of_x1 = std::exp( x1 );
//    //return
//    ret[0] = 2.0 * exp_of_x1 * std::sin( x2 );
//    ret[1] = 2.0 * exp_of_x1 * std::cos( x2 );
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
template < class PressureTraitsImp  >
void Pressure< PressureTraitsImp >::testMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "- testing class Pressure..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "  - x: " << x[0] << std::endl;
    debugStream << "       " << x[1] << std::endl;
    RangeType p;
    evaluate( x, p );
    debugStream << "  - p(x): " << p[0] << std::endl;
//    GradientRangeType grad_p;
//    gradient( x, grad_p );
//    debugStream << "  - grad p(x): " << grad_p[0] << std::endl;
//    debugStream << "               " << grad_p[1] << std::endl;
    infoStream << "  ...test passed!" << std::endl;
}

#endif // end of pressure.hh
