/**
 *  \file   force.hh
 *  \brief  contains a class Force with traitsclass ForceTraits
 **/

#ifndef FORCE_HH
#define FORCE_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by Force
 *
 *  \tparam gridDim
 *          dimension of the grid
 **/
template < int gridDim >
class ForceTraits
{
    public:
        typedef Dune::FieldVector< double, gridDim >
            DomainType;
        typedef Dune::FieldVector< double, gridDim >
            RangeType;
};

/**
 *  \brief  describes the force
 *  \tparam gridDim
 *          dimension of the grid
 *
 *  \todo   extensive docu with latex
 **/
template < int gridDim >
class Force
{
    public:
        typedef ForceTraits< gridDim >
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
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
         *  doing nothing
         **/
        ~Force()
        {
        }

        /**
         *  \brief  evaluates the force
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of force at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief  evaluates the force
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of force at given point
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
};

/**
 *  \brief  specialization for gridDim = 2
 **/
template < >
inline void Force< 2 >::evaluate( const DomainType& arg, RangeType& ret ) const
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
}

/**
 *  \brief  specialization for gridDim = 2
 **/
template < >
void Force< 2 >::testMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "- testing class Force..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "  - x: " << x[0] << std::endl;
    debugStream << "       " << x[1] << std::endl;
    RangeType f;
    evaluate( x, f );
    debugStream << "  - f(x): " << f[0] << std::endl;
    debugStream << "          " << f[1] << std::endl;
    infoStream << "  ...test passed!" << std::endl;
}

#endif // end of force.hh
