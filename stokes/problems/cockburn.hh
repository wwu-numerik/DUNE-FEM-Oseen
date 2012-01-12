#ifndef STOKES_PROBLEMS_COCKBURN_HH
#define STOKES_PROBLEMS_COCKBURN_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

static const std::string identifier = "Simple";
static const bool hasExactSolution	= true;

struct SetupCheck {
    std::string err;
    template < class GridPart , class ...Rest >
    bool operator()( const GridPart& gridPart, const Rest&... /*rest*/ ) {
        Stuff::GridDimensions< typename GridPart::GridType > grid_dim( gridPart.grid() );
        bool ok =  Stuff::aboutEqual( grid_dim.coord_limits[0].min(), -1. )
                && Stuff::aboutEqual( grid_dim.coord_limits[1].min(), -1. )
                && Stuff::aboutEqual( grid_dim.coord_limits[0].max(), 1. )
                && Stuff::aboutEqual( grid_dim.coord_limits[1].max(), 1. );
        err = ( boost::format( "\n******\nSetupCheck Failed!\ngrid dimension %f,%f - %f,%f\n" )
                % grid_dim.coord_limits[0].min()
                % grid_dim.coord_limits[1].min()
                % grid_dim.coord_limits[0].max()
                % grid_dim.coord_limits[1].max() ).str();
        if (!ok)
            return false;
        const double v = Parameters().getParam( "viscosity", -10.0 );
        ok = Stuff::aboutEqual( v, 1.0 );
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

        DirichletData( const FunctionSpaceImp& /*space*/ )
            : BaseType(  )
		{}

		 ~DirichletData()
		 {}

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
			ret[0] = 2.0 * std::exp( arg[0] ) * std::sin( arg[1] );
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
			const double x1 = arg[0];
			const double x2 = arg[1];
			const double exp_of_x1 = std::exp( x1 );
			const double sin_of_x2 = std::sin( x2 );
			ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
			ret[1] = exp_of_x1 * x2 * sin_of_x2;
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

#endif // STOKES_PROBLEMS_COCKBURN_HH
