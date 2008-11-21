/**
 *  \file   dirichletdata.hh
 *
 *  \brief  contains a class DirichletData with traitsclass DirichletDataTraits
 **/

#ifndef DIRICHLETDATA_HH
#define DIRICHLETDATA_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"

/**
 *  \brief  containing typedefs needed by DirichletData
 *
 *  \tparam gridDim
 *          dimension of the grid
 *  \tparam FunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int griddim, class FunctionSpaceImp >
class DirichletDataTraits
{
    public:
        static const unsigned int gridDim = griddim;
        typedef FunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
};

/**
 *  \brief  describes the dirichlet boundary data
 *
 *  \tparam DirichletTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class DirichletTraitsImp >
class DirichletData : public Dune::Function < typename DirichletTraitsImp::FunctionSpaceType, DirichletData < DirichletTraitsImp > >
{
    public:
        typedef DirichletTraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::FunctionSpaceType
            DirichletFunctionSpaceType;
        typedef DirichletData < Traits >
            ThisType;
        typedef Dune::Function < DirichletFunctionSpaceType, ThisType >
            BaseType;

        static const unsigned int gridDim = Traits::gridDim;
        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        DirichletData( const DirichletFunctionSpaceType& space )
            : BaseType( space )
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
          *
          * \param  arg
          *         point to evaluate at
          * \param  ret
          *         value of dirichlet boundary data at given point
          **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const;

        /**
         *  \brief  evaluates the dirichlet data
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of dirichlet data at given point
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
template < class DirichletTraitsImp >
inline void DirichletData< DirichletTraitsImp >::evaluate( const DomainType& arg, RangeType& ret ) const
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
 *  \brief  specialization for gridDim = 2
 **/
template < class DirichletTraitsImp >
void DirichletData< DirichletTraitsImp >::testMe() const
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

