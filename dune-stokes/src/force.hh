/** \file force.hh
    \brief contains a class Force with traitsclass ForceTraits
 **/

#ifndef FORCE_HH
#define FORCE_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by Force
 *  \tparam int gridDim dimension of the grid
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
 *  \tparam int gridDim dimension of the grid
 *
 *  \todo doc
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

        /**
         *  \brief  a simple test of all class' functionalities
         *  \arg  Logging::LogStream& stream where to print
         **/
        void TestMe() const;
};

/**
 *  \brief specialization for gridDim = 2
 **/
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
 *  \brief specialization for gridDim = 2
 **/
template < >
void Force< 2 >::TestMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "\nnow testing class Force..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "\n x: " << x[0] << std::endl;
    debugStream <<   "    " << x[1] << std::endl;
    RangeType f;
    Evaluate( x, f );
    debugStream << "\n f(x): " << f[0] << std::endl;
    debugStream <<  "        " << f[1] << std::endl << std::endl;
    infoStream << "...test passed!" << std::endl;
};

#endif // end of force.hh
