/**
 *  \file   analyticaldarcydata.hh
 *
 *  \todo   doc
 **/

#ifndef ANALYTICALDARCYDATA_HH
#define ANALYTICALDARCYDATA_HH

#include <cmath>

#include <dune/common/fvector.hh>

namespace Darcy
{

/**
 *  \brief  ConstantFunction
 *
 *  \todo   doc
 **/
template < class FunctionSpaceImp >
class ConstantFunction : public Dune::Function < FunctionSpaceImp , ConstantFunction < FunctionSpaceImp > >
{
    public:
        typedef ConstantFunction< FunctionSpaceImp >
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
         *  \todo   doc
         **/
        ConstantFunction( const FunctionSpaceImp& space, const RangeType& constValue )
            : BaseType ( space ),
              dim_( FunctionSpaceImp::dimDomain ),
              constValue_( constValue )
        {}

        /**
         *  \brief  destructor
         *
         *  \todo   doc
         **/
        ~ConstantFunction()
        {}

        /**
         *  \brief  evaluate
         *
         *  \todo   doc
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            ret = constValue_;
        }

        /**
         *  \brief  evaluate
         *
         *  \todo   doc
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret, const int id ) const
        {
            evaluate( arg, ret );
        }

    private:
        const int dim_;
        const RangeType& constValue_;

};  // end of ConstantFunction

}   // end of namespace Darcy

#endif // end of analyticaldirichletdata.hh
