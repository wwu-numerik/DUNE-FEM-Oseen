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
         *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
         **/
        Force( const double viscosity, const FunctionSpaceImp& space )
            : BaseType ( space ),
              viscosity_( viscosity ),
              dim_( FunctionSpaceImp::dimDomain )
        {}

        /**
         *  \brief  destructor
         *  doing nothing
         **/
        ~Force()
        {}

        /**
         *  \brief  evaluates the force
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of force at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            if ( dim_ == 1 ) {
                assert( !"force not implemented in 1D!" );
            }
            else if ( dim_ == 2 ) {
//                const double x1 = arg[0];
//                const double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                ret[0] = arg[1];
                ret[1] = -1.0 * arg[0];
#elif defined(MICRO_PROBLEM)
                ret[ 0 ] = 0.0;
                ret[ 1 ] = 0.0;
#else
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
#endif
            }
            else if ( dim_ == 3 ) {
//                const double x1 = arg[0];
//                const double x2 = arg[1];
//                const double x3 = arg[2];
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
                ret[2] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                assert( !"ROTATE_PROBLEM not implemented in 3D!" );
#elif defined(MICRO_PROBLEM)
                assert( !"MICRO_PROBLEM not implemented in 3D!" );
#else
                assert( !"force not implemented in 3D!" );
#endif
            }
            else {
                assert( !"force not implemented for more then 3 dimensions!" );
            }
        }

    private:
        const double viscosity_;
        const int dim_;
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
            : BaseType( space ),
              dim_( FunctionSpaceImp::dimDomain )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
         ~DirichletData()
         {}

        void evaluate( const DomainType& arg, RangeType& ret, const int id ) const
        {
            if ( dim_ == 1 ) {
                assert( !"dirichlet data not implemented in 1D!" );
            }
            else if ( dim_ == 2 ) {
                // some computations
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(ROTATE_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM)
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
#else
                const double x1 = arg[0];
                const double x2 = arg[1];
                double exp_of_x1 = std::exp( x1 );
                double sin_of_x2 = std::sin( x2 );
                double cos_of_x2 = std::cos( x2 );
                //return
                ret[0] = x2 * cos_of_x2;
                ret[0] += sin_of_x2;
                ret[0] *= -1.0 * exp_of_x1;
                ret[1] = exp_of_x1 * x2 * sin_of_x2;
#endif
            }
            else if ( dim_ == 3 ) {
//                double x1 = arg[0];
//                double x2 = arg[1];
//                double x3 = arg[2];
                // some computations
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(ROTATE_PROBLEM)
                assert( !"ROTATE_PROBLEM not implemented in 3D!" );
#elif defined(MICRO_PROBLEM)
                assert( !"MICRO_PROBLEM not implemented in 3D!" );
#else
                assert( !"dirichlet data not implemented in 3D!" );
#endif
            }
            else {
                assert( !"dirichlet data not implemented for more than 3 dimensions!" );
            }
        }

         /**
          * \brief  evaluates the dirichlet data
          * \param  arg
          *         point to evaluate at
          * \param  ret
          *         value of dirichlet boundary data at given point
          **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const {}

    private:
        const int dim_;

};


#endif // end of analyticaldata.hh
