/** \file pressure.hh
    \brief contains a class Pressure
 **/

#ifndef PRESSURE_HH
#define PRESSURE_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

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
         *  \return RangeType value of pressure at point arg
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

        /**
         *  \brief  a simple test of all class' functionalities
         *  \arg  Logging::LogStream& stream where to print
         **/
        void TestMe() const;

    private:
 };

template < >
inline void Pressure< 2 >::Evaluate( const DomainType& arg, RangeType& ret ) const
{
    // play save
    assert( arg.dim() == 2 );
    assert( ret.dim() == 1 );
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

template < >
void Pressure< 2 >::TestMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "\nnow testing class Pressure..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "\n x: " << x[0] << std::endl;
    debugStream <<   "    " << x[1] << std::endl;
    RangeType p;
    Evaluate( x, p );
    debugStream << "\n p(x): " << p[0] << std::endl;
    GradientRangeType grad_p;
    Gradient( x, grad_p );
    debugStream << "\n grad p(x): " << grad_p[0] << std::endl;
    debugStream <<   "            " << grad_p[1] << std::endl << std::endl;
    infoStream << "...test passed!" << std::endl;
};

#endif // end of pressure.hh
