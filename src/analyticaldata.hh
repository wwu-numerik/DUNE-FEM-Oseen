/**
 *  \file   analyticaldata.hh
 *  \brief  contains classes representing analytical data
 **/

#ifndef ANALYTICALDATA_HH
#define ANALYTICALDATA_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/function/common/function.hh>
#include <cmath>

//#include <dune/common/fvector.hh>
//#include <dune/common/fmatrix.hh>


/**
 *  \todo   texdoc
 **/
template < class FunctionSpaceImp >
class Force : public Dune::Fem::Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
{
    public:
        typedef Force< FunctionSpaceImp >
            ThisType;
        typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
            BaseType;
        typedef typename BaseType::DomainType
            DomainType;
        typedef typename BaseType::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
         **/
        Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0, const double scaling_factor = 1.0 )
            : BaseType ( space ),
              viscosity_( viscosity ),
			  alpha_( alpha ),
			  scaling_factor_( scaling_factor )
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
			dune_static_assert( dim_ > 1 && dim_ < 4, "Force_Unsuitable_WorldDim");

            if ( dim_ == 2 ) {
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                ret[0] = arg[1];
                ret[1] = -1.0 * arg[0];
#elif defined(POROSITY_PROBLEM)
                ret[ 0 ] = 0.0;
                ret[ 1 ] = 1.0;
#elif defined(POROSITY_PROBLEM_WOIDS)
                ret[ 0 ] = 0.0;
                ret[ 1 ] = 0.0;
#elif defined(GENERALIZED_STOKES_PROBLEM)
                const double x = arg[0];
                const double y = arg[1];
				const double cos_p = std::cos( M_PI_2 * (x+y) );
				const double cos_m = std::cos( M_PI_2 * (x-y) );
				const double pi_sq = M_PI_2 * M_PI ;//std::pow( M_PI_2 , 2 );
//                const double tmp = alpha_ *  + M_PI_2 * M_PI * std::cos( M_PI_2 * ( x + y ) ) + M_PI_2 * std::cos( M_PI_2 * ( x - y ) ) ;
                ret[0]  = cos_p * ( alpha_ + viscosity_ * pi_sq ) + M_PI_2 * cos_m;
                ret[1]  = cos_p * ( - alpha_ - viscosity_ * pi_sq ) - M_PI_2 * cos_m;

#elif defined(DARCY_PROBLEM)
                // im verhältnis zu [-1,1]^2
                double scaleX = DSC_CONFIG_GET( "domain_scale_x", 2.0 );
                double scaleY = DSC_CONFIG_GET( "domain_scale_y", 2.0 );
                double shiftX = DSC_CONFIG_GET( "domain_shift_x", -1.0 );
                double shiftY = DSC_CONFIG_GET( "domain_shift_y", -1.0 );
                ret[0] = ( arg[1] * scaleY ) + shiftY;
                ret[1] = -1.0 * ( ( arg[0] * scaleX ) + shiftX );
#elif defined(MICRO_PROBLEM_X)
                ret[0] = 1.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM_Y)
                ret[0] = 0.0;
                ret[1] = 1.0;
#elif defined(COCKBURN_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
#else
                DUNE_THROW(Dune::Exception, "not implemented");
#endif
            }
            else if ( dim_ == 3 ) {
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
#elif defined(POROSITY_PROBLEM)
                assert( !"POROSITY_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM_WOIDS)
                assert( !"POROSITY_PROBLEM_WOIDS not implemented in 3D!" );
#elif defined(GENERALIZED_STOKES_PROBLEM)
                assert( !"GENERALIZED_STOKES_PROBLEM not implemented in 3D!" );
#elif defined(DARCY_PROBLEM)
                assert( !"DARCY_PROBLEM not implemented in 3D!" );
#else
                DUNE_THROW(Dune::Exception, "not implemented");
#endif
            }
            ret *= scaling_factor_;
        }

    private:
        const double viscosity_;
        const double alpha_;
		const double scaling_factor_;
		static const int dim_ = FunctionSpaceImp::dimDomain;
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
class DirichletData : public Dune::Fem::Function < FunctionSpaceImp, DirichletData < FunctionSpaceImp > >
{
    public:
        typedef DirichletData< FunctionSpaceImp >
            ThisType;
        typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
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
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
         ~DirichletData()
         {}

		template < class IntersectionType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionType& intersection ) const
        {
			const int id = intersection.boundaryId();
			dune_static_assert( dim_ > 1 && dim_ < 4 , "DirichletData_Unsuitable_WorldDim");
			if ( dim_ == 2 ) {
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
#elif defined(POROSITY_PROBLEM)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(POROSITY_PROBLEM_WOIDS)
                const double x1 = arg[0];
                const double x2 = arg[1];
                if ( !( x2 > 0.0 ) ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 < 1.0 ) ) { // right faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x2 < 1.0 ) ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 > 0.0 ) ) { // left faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else {
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(GENERALIZED_STOKES_PROBLEM)
                const double x1 = arg[0];
                const double x2 = arg[1];
                const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
                ret[0] = tmp;
                ret[1] = -1.0 * tmp;
#elif defined(DARCY_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM_X)
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
#elif defined(MICRO_PROBLEM_Y)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(COCKBURN_PROBLEM)
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
#else
                DUNE_THROW(Dune::Exception, "not implemented");
#endif
            }
            else if ( dim_ == 3 ) {
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
#elif defined(POROSITY_PROBLEM)
                assert( !"POROSITY_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM_WOIDS)
                assert( !"POROSITY_PROBLEM_WOIDS not implemented in 3D!" );
#elif defined(GENERALIZED_STOKES_PROBLEM)
                assert( !"GENERALIZED_STOKES_PROBLEM not implemented in 3D!" );
#elif defined(DARCY_PROBLEM)
                assert( !"DARCY_PROBLEM not implemented in 3D!" );

#else
                DUNE_THROW(Dune::Exception, "not implemented");
#endif
            }
        }

         /**
          * \brief  evaluates the dirichlet data
          * \param  arg
          *         point to evaluate at
          * \param  ret
          *         value of dirichlet boundary data at given point
          **/
        inline void evaluate( const DomainType& /*arg*/, RangeType& /*ret*/ ) const {}

    private:
		static const int dim_ = FunctionSpaceImp::dimDomain ;
};

#endif // end of analyticaldata.hh

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

