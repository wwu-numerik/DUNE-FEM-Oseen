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
            if ( arg.dim() == 2 ) {
//                if ( !( arg[1] < 0.0 ) && !( arg[1] > 0.0 ) ) { // micro problem for x dimension
                    if ( id == 2 ) { // faces on inner hole
                        ret[ 0 ] = 0.0;
                        ret[ 1 ] = 0.0;
                    }
                    else if ( id == 3 ) { // bottom faces
                        ret[ 0 ] = 0.0;
                        ret[ 1 ] = 0.0;
                    }
                    else if ( id == 4 ) { // right faces
                        ret[ 0 ] = 1.0;
                        ret[ 1 ] = 0.0;
                    }
                    else if ( id == 5 ) { // top faces
                        ret[ 0 ] = 0.0;
                        ret[ 1 ] = 0.0;
                    }
                    else if ( id == 6 ) { // left faces
                        ret[ 0 ] = 1.0;
                        ret[ 1 ] = 0.0;
                    }
//                }
//                else { // micro problem for y dimension
//                    if ( id == 2 ) { // faces on inner hole
//                        ret[ 0 ] = 0.0;
//                        ret[ 1 ] = 0.0;
//                    }
//                    else if ( id == 3 ) { // bottom faces
//                        ret[ 0 ] = 0.0;
//                        ret[ 1 ] = 1.0;
//                    }
//                    else if ( id == 4 ) { // right faces
//                        ret[ 0 ] = 0.0;
//                        ret[ 1 ] = 0.0;
//                    }
//                    else if ( id == 5 ) { // top faces
//                        ret[ 0 ] = 0.0;
//                        ret[ 1 ] = 1.0;
//                    }
//                    else if ( id == 6 ) { // left faces
//                        ret[ 0 ] = 0.0;
//                        ret[ 1 ] = 0.0;
//                    }
//                }
            }
            else {
                assert( !"constant dirichlet data only implemented in 2D!" );
            }
        }

    private:
        const int dim_;
        const RangeType& constValue_;

};  // end of ConstantFunction

}   // end of namespace Darcy

#endif // end of analyticaldirichletdata.hh
