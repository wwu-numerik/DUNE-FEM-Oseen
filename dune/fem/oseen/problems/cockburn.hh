#ifndef STOKES_PROBLEMS_COCKBURN_HH
#define STOKES_PROBLEMS_COCKBURN_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/common/float_cmp.hh>
#include "common.hh"

namespace StokesProblems {
namespace Cockburn {

static const std::string identifier = "Simple";
static const bool hasExactSolution	= true;

struct SetupCheck {
    std::string err;
    template < class GridPart , class ...Rest >
    bool operator()( const GridPart& gridPart, const Rest&... /*rest*/ ) {
        DSG::Dimensions< typename GridPart::GridType > grid_dim( gridPart.grid() );
        bool ok =  Dune::FloatCmp::eq( grid_dim.coord_limits[0].min(), -1. )
                && Dune::FloatCmp::eq( grid_dim.coord_limits[1].min(), -1. )
                && Dune::FloatCmp::eq( grid_dim.coord_limits[0].max(), 1. )
                && Dune::FloatCmp::eq( grid_dim.coord_limits[1].max(), 1. );
        err = ( boost::format( "\n******\nSetupCheck Failed!\ngrid dimension %f,%f - %f,%f\n" )
                % grid_dim.coord_limits[0].min()
                % grid_dim.coord_limits[1].min()
                % grid_dim.coord_limits[0].max()
                % grid_dim.coord_limits[1].max() ).str();
        if (!ok)
            return false;
        const double v = DSC_CONFIG_GET( "viscosity", -10.0 );
        ok = Dune::FloatCmp::eq( v, 1.0 );
        err = ( boost::format( "viscosity %f\n" ) % v ).str();
        return ok;
    }
    std::string error() {
        return err;
    }
};


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

        Force( const double viscosity, const double alpha = 0.0, const double scaling_factor = 1.0 )
            : viscosity_( viscosity ),
			  alpha_( alpha ),
			  scaling_factor_( scaling_factor )
		{}

		inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 2, "Force_Unsuitable_WorldDim");
			ret[0] = 0.0;//arg[1];
			ret[1] = 0.0;//arg[0];
			ret *= scaling_factor_;
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

		template < class IntersectionType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection*/ ) const
		{
			dune_static_assert( dim_ == 2, "DirichletData_Unsuitable_WorldDim");
			const double x1 = arg[0];
			const double x2 = arg[1];
			const double exp_of_x1 = std::exp( x1 );
			const double sin_of_x2 = std::sin( x2 );
			const double cos_of_x2 = std::cos( x2 );
			//return
			ret[0] = x2 * cos_of_x2;
			ret[0] += sin_of_x2;
			ret[0] *= -1.0 * exp_of_x1;
			ret[1] = exp_of_x1 * x2 * sin_of_x2;
		}

		inline void evaluate( const DomainType& /*arg*/, RangeType& /*ret*/ ) const {assert(false);}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain ;
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

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			ret[0] = 2.0 * std::exp( arg[0] ) * std::sin( arg[1] );
		}

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

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			const double x1 = arg[0];
			const double x2 = arg[1];
			const double exp_of_x1 = std::exp( x1 );
			const double sin_of_x2 = std::sin( x2 );
			ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
			ret[1] = exp_of_x1 * x2 * sin_of_x2;
		}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain;
};

} // namespace Cockburn {
} // namespace StokesProblems {

#endif // STOKES_PROBLEMS_COCKBURN_HH

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

