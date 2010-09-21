#ifndef SIMPLE_HH
#define SIMPLE_HH

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

		inline void evaluate( const DomainType& arg, RangeType& ret ) const
		{
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > __CLASS__Unsuitable_WorldDim;
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
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionType& intersection ) const
		{
			const int id = intersection.boundaryId();
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > DirichletData_Unsuitable_WorldDim;
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert( false ); }

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
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > Velocity_Unsuitable_WorldDim;
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
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > Pressure_Unsuitable_WorldDim;
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

#endif // SIMPLE_HH
