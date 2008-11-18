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
 *          dimension of the grid (unused)
 *  \tparam ForceFunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int gridDim, class ForceFunctionSpaceImp  >
class ForceTraits
{
    public:
        typedef ForceFunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
};

/**
 *  \brief  describes the force
 *  \tparam ForceTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class ForceTraitsImp >
class Force : public Dune::Function < typename ForceTraitsImp::FunctionSpaceType , Force < ForceTraitsImp > >
{
    public:
        typedef ForceTraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::FunctionSpaceType
            FunctionSpaceType;
        typedef Force < Traits >
            ThisType;
        typedef Dune::Function < FunctionSpaceType, ThisType >
            BaseType;


        /**
         *  \brief  constructor
         *
         *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
         **/
        Force( const double viscosity, const FunctionSpaceType& space )
            : BaseType ( space ),
              viscosity_( viscosity )
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

    private:
        double viscosity_;
};

/**
 *  \brief  specialization for gridDim = 2
 **/
template < class ForceTraitsImp >
inline void Force< ForceTraitsImp >::evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play safe
    assert( arg.dim() == 2 );
    assert( ret.dim() == 2 );
    // some computations
    double x1 = arg[0];
    double x2 = arg[1];
    double exp_of_x1 = std::exp( x1 );
    double cos_of_x2 = std::cos( x2 );
    //return
    ret[0] = 2.0 * exp_of_x1 *
        ( ( 1.0 - viscosity_ ) * std::sin( x2 )
            + viscosity_ * cos_of_x2 );
    ret[1] = 2.0 * ( 1.0 - viscosity_ ) * exp_of_x1 * cos_of_x2;
}

/**
 *  \brief  specialization for gridDim = 2
 **/
template < class ForceTraitsImp >
void Force< ForceTraitsImp >::testMe() const
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
