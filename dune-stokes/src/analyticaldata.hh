/**
 *  \file   analyticaldata.hh
 *  \brief  contains classes representing analytical data
 **/

#ifndef ANALYTICALDATA_HH
#define ANALYTICALDATA_HH

#include <cmath>

#include <dune/common/fvector.hh>

/**
 *  \todo   texdoc
 **/
template < class FunctionSpaceImp >
class Force : public Dune::Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
{
    public:
        typedef Force< FunctionSpaceImp >
            ThisType;
        typedef Dune::Function< FunctionSpaceImp, ThisType >
            BaseType;
        typedef typename BaseType::DomainType
            DomainType;
        typedef typename BaseType::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *
         *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
         **/
        Force( const double viscosity, const FunctionSpaceImp& space )
            : BaseType ( space ),
              viscosity_( viscosity )
        {}

        /**
         *  \brief  destructor
         *  doing nothing
         **/
        ~Force()
        {}

        /**
         *  \brief  evaluates the force
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of force at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
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
         *  \brief  evaluates the force
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of force at given point
         **/
        RangeType operator () ( const DomainType& arg )
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

    private:
        double viscosity_;
};

/**
 *  \brief  describes the dirichlet boundary data
 *
 *  \tparam DirichletTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class FunctionSpaceImp >
class DirichletData : public Dune::Function < FunctionSpaceImp, DirichletData < FunctionSpaceImp > >
{
    public:
        typedef DirichletData< FunctionSpaceImp >
            ThisType;
        typedef Dune::Function< FunctionSpaceImp, ThisType >
            BaseType;
        typedef typename BaseType::DomainType
            DomainType;
        typedef typename BaseType::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        DirichletData( const FunctionSpaceImp& space )
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
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
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
};


#endif // end of analyticaldata.hh
