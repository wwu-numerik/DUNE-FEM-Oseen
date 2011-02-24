#ifndef STOKES_PROBLEMS_COCKBURN_HH
#define STOKES_PROBLEMS_COCKBURN_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/misc.hh>

static const std::string identifier = "Simple";
static const bool hasExactSolution	= true;

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

		Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0, const double scaling_factor = 1.0 )
			: BaseType ( space ),
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

		DirichletData( const FunctionSpaceImp& space )
			: BaseType( space )
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

		inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain ;
};

template < class FunctionSpaceImp  >
class Pressure : public Dune::Function < FunctionSpaceImp , Pressure < FunctionSpaceImp > >
{
	public:
		typedef Pressure< FunctionSpaceImp >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
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
		Pressure( const FunctionSpaceImp& f_space )
			: BaseType( f_space )
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
class Velocity : public Dune::Function < FunctionSpaceImp , Velocity < FunctionSpaceImp > >
{
	public:
		typedef Velocity< FunctionSpaceImp >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;

		Velocity( const FunctionSpaceImp& f_space )
			: BaseType( f_space )
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
