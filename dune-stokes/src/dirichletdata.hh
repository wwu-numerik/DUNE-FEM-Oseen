/** \file dirichletdata.hh
    \brief contains a class DirichletData with traitsclass DirichletDataTraits
 **/

#ifndef DIRICHLETDATA_HH
#define DIRICHLETDATA_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by DirichletData
 *  \tparam int gridDim dimension of the grid
 **/
template < int gridDim >
class DirichletDataTraits
{
    public:
        typedef Dune::FieldVector< double, gridDim >
            DomainType;
        typedef Dune::FieldVector< double, gridDim >
            RangeType;
};

/**
 *  \brief  describes the dirichlet boundary data
 *  \tparam int gridDim dimension of the grid
 *
 *  \todo doc
 **/
template < int gridDim >
class DirichletData
{
    public:
        typedef DirichletDataTraits< gridDim >
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *  doing nothing
         **/
        DirichletData()
        {
        }

        /**
         *  \brief  destructor
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
        inline void evaluate( DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief evaluates the dirichlet data
         *  \arg DomainType& arg point to be evaluated at
         *  \return RangeType ret value of dirichlet data at point arg
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
 *  \brief specialization for gridDim = 2
 **/
template < >
inline void DirichletData< 2 >::evaluate( DomainType& arg, RangeType& ret ) const
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
}

/**
 *  \brief specialization for gridDim = 2
 **/
template < >
void DirichletData< 2 >::testMe() const
{
    // some logstreams
    Logging::LogStream& infoStream = Logger().Info();
    Logging::LogStream& debugStream = Logger().Dbg();
    infoStream << "- testing class DirichletData..." << std::endl;
    //tests
    DomainType x;
    x[0] = 1.0;
    x[1] = 1.0;
    debugStream << "  - x: " << x[0] << std::endl;
    debugStream << "       " << x[1] << std::endl;
    RangeType gd;
    evaluate( x, gd );
    debugStream << "  - gd(x): " << gd[0] << std::endl;
    debugStream << "           " << gd[1] << std::endl;
    infoStream << "  ...test passed!" << std::endl;
}

#endif  // end of dirichletdata.hh

