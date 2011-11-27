#ifndef STOKES_PROBLEMS_GENERALIZED_HH
#define STOKES_PROBLEMS_GENERALIZED_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/misc.hh>

static const std::string identifier = "Generalized";
static const bool hasExactSolution	= true;

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

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 2, "Force_Unsuitable_WorldDim");
			const double x = arg[0];
			const double y = arg[1];
			const double cos_p = std::cos( M_PI_2 * (x+y) );
			const double cos_m = std::cos( M_PI_2 * (x-y) );
			const double pi_sq = M_PI_2 * M_PI ;//std::pow( M_PI_2 , 2 );
			ret[0]  = cos_p * ( alpha_ + viscosity_ * pi_sq ) + M_PI_2 * cos_m;
			ret[1]  = cos_p * ( - alpha_ - viscosity_ * pi_sq ) - M_PI_2 * cos_m;
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
			const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
			ret[0] = tmp;
			ret[1] = -1.0 * tmp;
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert( false ); }

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
			dune_static_assert( dim_ == 2, "Velocity_Unsuitable_WorldDim");
			const double x1 = arg[0];
			const double x2 = arg[1];
			const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
			ret[0] = tmp;
			ret[1] = -1.0 * tmp;
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
			dune_static_assert( dim_ == 2, "Pressure_Unsuitable_WorldDim");
			const double x1 = arg[0];
			const double x2 = arg[1];
			ret[0] = std::sin( ( M_PI_2 ) * ( x1 - x2 ) );
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

#endif // STOKES_PROBLEMS_GENERALIZED_HH
