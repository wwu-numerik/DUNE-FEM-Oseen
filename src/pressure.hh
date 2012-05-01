/**
 *  \file   pressure.hh
 *
 *  \brief  contains a class Pressure with traitsclass PressureTraits
 **/

#ifndef PRESSURE_HH
#define PRESSURE_HH

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/stuff/logging.hh>

/**
 *  \brief  containing typedefs needed by Pressure
 *
 *  \tparam gridDim (unused)
 *          dimension of the grid
 *
 *  \tparam PressureFunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int gridDim, class PressureFunctionSpaceImp >
class PressureTraits
{
    public:
        typedef PressureFunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType
        //typedef Dune::FieldVector< double, gridDim >
            GradientRangeType;

};

/**
 *  \brief  describes the presure
 *
 *  \tparam PressureTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class PressureTraitsImp >
class Pressure : public Dune::Fem::Function < typename PressureTraitsImp::FunctionSpaceType, Pressure < PressureTraitsImp > >
{
    public:
        typedef PressureTraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::GradientRangeType
            GradientRangeType;
        typedef typename Traits::FunctionSpaceType
            PressureFunctionSpaceType;
        typedef Dune::Fem::Function < typename PressureTraitsImp::FunctionSpaceType , Pressure < PressureTraitsImp > >
            BaseType;
        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        Pressure( const PressureFunctionSpaceType& press_space )
            : BaseType( press_space ),
            dim_( PressureTraitsImp::FunctionSpaceType::dimDomain )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~Pressure()
        {}

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of pressure at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            if ( dim_ == 1 ) {
                assert( !"pressure not implemented in 1D" );
            }
            else if ( dim_ == 2 ) {
                double x1 = arg[0];
                double x2 = arg[1];
#ifdef CONSTANT_PROBLEM
                ret[0] = -x2;
#elif defined(GENERALIZED_STOKES_PROBLEM)
                ret[0] = std::sin( ( M_PI_2 ) * ( x1 - x2 ) );
#else
                ret[0] = 2.0 * std::exp( x1 ) * std::sin( x2 );
#endif
            }
            else if ( dim_ == 3 ) {
//                double x1 = arg[0];
//                double x2 = arg[1];
//                double x3 = arg[2];
#ifdef SIMPLE_PROBLEM
                assert( !"SIMPLE_PROBLEM not implemented in 3D" );
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0;
                ret[1] = 0;
                ret[2] = 0;
#elif defined(GENERALIZED_STOKES_PROBLEM)
                assert( !"GENERALIZED_STOKES_PROBLEM not implemented in 3D" );
#elif defined(AORTA_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#else
                assert( !"pressure not implemented in 3D" );
#endif
            }
            else {
                assert( !"pressure not implemented for more than 3 dimensions" );
            }
        }

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of pressure at given point
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

        /**
         *  \brief  evaluates the gradient of the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of gradient of the pressure at given point
         **/
//        inline void gradient(const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;

    private:
        const int dim_;
 };

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class PressureTraitsImp  >
//inline void Pressure< PressureTraitsImp  >::evaluate( const DomainType& arg, RangeType& ret ) const
//{
//    // play save
//    assert( arg.dim() == 2 );
//    assert( ret.dim() == 1 );
//    double x1 = arg[0];
//    double x2 = arg[1];
//#ifdef SIMPLE_PROBLEM
//    ret[0] = -x2;
//#elif defined(CONSTANT_PROBLEM)
//    ret[0] = -x2;
//#else
//    ret[0] = 2.0 * std::exp( x1 ) * std::sin( x2 );
//#endif
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class PressureTraitsImp  >
//inline void Pressure< PressureTraitsImp >::gradient( const DomainType& arg, GradientRangeType& ret ) const
//{
//    // play safe
//    assert( arg.dim() == 2 );
////    assert( ret.dim() == 2 );
//    // some computations
//    double x1 = arg[0];
//    double x2 = arg[1];
//    double exp_of_x1 = std::exp( x1 );
//    //return
//    ret[0] = 2.0 * exp_of_x1 * std::sin( x2 );
//    ret[1] = 2.0 * exp_of_x1 * std::cos( x2 );
//}

/**
 *  \brief  specialization for gridDim = 2
 **/
//template < class PressureTraitsImp  >
//void Pressure< PressureTraitsImp >::testMe() const
//{
//    // some logstreams
//    Logging::LogStream& infoStream = Logger().Info();
//    Logging::LogStream& debugStream = Logger().Dbg();
//    infoStream << "- testing class Pressure..." << std::endl;
//    //tests
//    DomainType x;
//    x[0] = 1.0;
//    x[1] = 1.0;
//    debugStream << "  - x: " << x[0] << std::endl;
//    debugStream << "       " << x[1] << std::endl;
//    RangeType p;
//    evaluate( x, p );
//    debugStream << "  - p(x): " << p[0] << std::endl;
////    GradientRangeType grad_p;
////    gradient( x, grad_p );
////    debugStream << "  - grad p(x): " << grad_p[0] << std::endl;
////    debugStream << "               " << grad_p[1] << std::endl;
//    infoStream << "  ...test passed!" << std::endl;
//}

#endif // end of pressure.hh

/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

