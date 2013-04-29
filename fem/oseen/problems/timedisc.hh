#ifndef STOKES_PORBLEMS_TIMEDISC_HH
#define STOKES_PORBLEMS_TIMEDISC_HH
#include <dune/fem/function/common/function.hh>
#include <dune/stuff/common/misc.hh>
#include "common.hh"

namespace StokesProblems {
namespace TimeDisc {

ALLGOOD_SETUPCHECK;
static const std::string identifier = "TimeDisc";
static const bool hasExactSolution	= true;
static const double disc_time = 1.0;

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

        Force( const double viscosity, const FunctionSpaceImp& /*space*/, const double alpha = 0.0, const double scaling_factor = 1.0 )
            : BaseType (  ),
			  viscosity_( viscosity ),
			  alpha_( alpha ),
			  scaling_factor_( scaling_factor )
		{}

		~Force()
		{}

        inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 2  , "__CLASS__ Wrong world dim");
//			const double x			= arg[0];
//			const double y			= arg[1];
			const double v			= viscosity_;
			//laplce
			ret[0] = -2*std::pow(disc_time,3.0)*v;
			ret[1] = 0;
			//grad p
			ret[0] += disc_time;
			ret[1] += 1;
			//conv
//					  ret[0] += std::pow(disc_time,5.0)*2*x*y;
//					  ret[1] += std::pow(disc_time,5.0)*y*y;
			//dt u
//					  ret[0] += std::pow(disc_time,2.0)*3*y*y;
//					  ret[1] += 2*disc_time*x;

//					  ret *=0;
		}

	private:
		const double viscosity_;
		const double alpha_;
		const double scaling_factor_;
		static const int dim_ = FunctionSpaceImp::dimDomain;
};

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

        DirichletData( const FunctionSpaceImp& /*space*/ )
            : BaseType(  )
		{}

		 ~DirichletData()
		 {}

		template < class IntersectionType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection*/ ) const
		{
			ret[0] = std::pow(disc_time,3.0)* arg[1] * arg[1];
			ret[1] = std::pow(disc_time,2.0)* arg[0];
		}

        inline void evaluate( const DomainType& /*arg*/, RangeType& /*ret*/ ) const { assert( false ); }

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain;
};

template < class FunctionSpaceImp  >
class Velocity : public Dune::Fem::Function < FunctionSpaceImp , Velocity < FunctionSpaceImp > >
{
	public:
		typedef Velocity< FunctionSpaceImp >
			ThisType;
		typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;

        Velocity( const FunctionSpaceImp& /*f_space*/ )
            : BaseType(  )
		{}

		~Velocity()
		{}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 2  , "Wrong world dim");
			ret[0] = std::pow(disc_time,3.0)* arg[1] * arg[1];
			ret[1] = std::pow(disc_time,2.0)* arg[0];
		}

		RangeType operator () ( const DomainType& arg)
		{
			RangeType ret;
			evaluate( arg, ret );
			return ret;
		}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain;
};

template < class FunctionSpaceImp  >
class Pressure : public Dune::Fem::Function < FunctionSpaceImp , Pressure < FunctionSpaceImp > >
{
	public:
		typedef Pressure< FunctionSpaceImp >
			ThisType;
		typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;

		/**
		 *  \brief constructor
		 *
		 *  doing nothing besides Base init
		 **/
        Pressure( const FunctionSpaceImp& /*f_space*/ )
            : BaseType(  )
		{}

		/**
		 *  \brief  destructor
		 *
		 *  doing nothing
		 **/
		~Pressure()
		{}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			ret = disc_time * arg[0] + arg[1] - ( ( disc_time + 1) / 2.0 );
		}

		RangeType operator () ( const DomainType& arg)
		{
			RangeType ret;
			evaluate( arg, ret );
			return ret;
		}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain;
};

} // namespace TimeDisc {
} // namespace StokesProblems {

#endif // STOKES_PORBLEMS_TIMEDISC_HH

/** Copyright (c) 2012, Rene Milk 
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

