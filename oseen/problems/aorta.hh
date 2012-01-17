#ifndef PROBLEM_AORTA_HH
#define PROBLEM_AORTA_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/grid.hh>
#include "common.hh"

ALLGOOD_SETUPCHECK;

static const std::string identifier = "Simple";
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

        inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 3, "Force_Unsuitable_WorldDim");
			ret = RangeType( 0 );
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
        void evaluate( const DomainType& /*arg*/, RangeType& ret, const IntersectionType& intersection ) const
		{
				const int id = intersection.boundaryId();
				typedef Dune::FieldVector< typename IntersectionType::ctype, IntersectionType::dimension - 1 >
					LocalVectorType;

				LocalVectorType center = intersection.geometryInInside().local( intersection.geometryInInside().center() );
				RangeType normal = intersection.unitOuterNormal( center );
				ret = normal;
				double factor = Parameters().getParam( "gd_factor", 1.0 );
					switch ( id ) {
						case 1: {
							factor = 0;
							break;
						}
						case 2: {
							factor *= -1;
							break;
						}
						case 6:
						case 5:
						case 4:
						case 3: {
							factor *= 1;
							break;
						}
						default:
						assert( false );
					}
				ret *= factor;
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

        inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 3, "Velocity_Unsuitable_WorldDim");
			ret = RangeType(0);
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

        inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			dune_static_assert( dim_ == 3, "Pressure_Unsuitable_WorldDim");
			ret = RangeType(0);
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

#endif // PROBLEM_AORTA_HH
