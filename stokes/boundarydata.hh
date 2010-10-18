#ifndef BOUNDARYDATA_HH
#define BOUNDARYDATA_HH

#include <dune/stokes/boundaryinfo.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/misc.hh>

template < template < class > class DiricheltDataImp >
struct DefaultDirichletDataTraits {

	template < class FunctionSpaceImp, class GridPartImp >
	struct Implementation {
		typedef DiricheltDataImp< FunctionSpaceImp>
			AnalyticalDirichletDataType;

		template <class DiscreteStokesFunctionWrapper >
		static AnalyticalDirichletDataType getInstance( const DiscreteStokesFunctionWrapper& wrapper ) {
			return 	AnalyticalDirichletDataType( wrapper.discreteVelocitySpace() );
		}
	};
};

//! \todo RENE needs to doc me and move me to stuff
template < class FunctionSpaceImp >
class InOutFluxDirichletData : public Dune::Function < FunctionSpaceImp, InOutFluxDirichletData < FunctionSpaceImp > >
{
	public:
		enum BoundaryType {
			zeroBoundary	= 0,
			influxBoundary	= 1,
			outfluxBoundary	= 2
		};

		typedef InOutFluxDirichletData< FunctionSpaceImp >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;
		typedef std::map< int, BoundaryType >
			BoundaryIdTypeMapType;
		typedef typename BoundaryIdTypeMapType::const_iterator
			BoundaryIdTypeMapTypeConstIterator;

		/**
			*  \brief  constructor
			*
			*
			**/
		InOutFluxDirichletData( const FunctionSpaceImp& space )
			: BaseType( space ),
				dim_( FunctionSpaceImp::dimDomain )
		{
			zeroBoundaryIds_ 	= Parameters().getList( "zeroBoundaryIds" , 1 );
			influxBoundaryIds_	= Parameters().getList( "influxBoundaryIds" , 2 );
			outfluxBoundaryIds_	= Parameters().getList( "outfluxBoundaryIds" , 3 );
			setupBoundaryIdTypeMap_();
		}

		/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
		~InOutFluxDirichletData()
		{}

		template < class IntersectionIteratorType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionIteratorType& faceIter ) const
		{
			const int id = faceIter.boundaryId();

			typedef Dune::FieldVector< typename IntersectionIteratorType::ctype, IntersectionIteratorType::dimension - 1 >
				LocalVectorType;

			LocalVectorType center = Stuff::getBarycenterLocal( faceIter.intersectionSelfLocal() );
			RangeType normal = faceIter.unitOuterNormal( center );
			static const double gd_factor = Parameters().getParam( "gd_factor", 1.0 );
			ret = normal;
			BoundaryIdTypeMapTypeConstIterator id_it = boundaryIdTypeMap_.find( id );
			assert ( id_it != boundaryIdTypeMap_.end() );

			switch ( id_it->second ) {
				case zeroBoundary: {
					ret *= 0;
					return;
				}
				case influxBoundary:
					ret *= -1;
				case outfluxBoundary: {
					ret *= gd_factor;
					return;
				}
				default:
					assert( false );
					return;
			}
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert(false); }

	protected:
		void setupBoundaryIdTypeMap_() {
			Logger().Info() << "\t- Using ids ";
			for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = zeroBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \t\tfor g_d = 0 \n\t        ids ";
			for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = influxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = - n \n\t        ids ";
			for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = outfluxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = n " << std::endl;
		}

		std::vector< int > zeroBoundaryIds_;
		std::vector< int > influxBoundaryIds_;
		std::vector< int > outfluxBoundaryIds_;
		BoundaryIdTypeMapType boundaryIdTypeMap_;

		const int dim_;
};

//! \todo RENE needs to doc me and move me to stuff
template < class FunctionSpaceImp, class GridPartType >
class FirstOrderBoundaryShapeFunction : public Dune::BoundaryShapeFunctionBase< FunctionSpaceImp, GridPartType >
{
	typedef Dune::BoundaryShapeFunctionBase < FunctionSpaceImp, GridPartType >
		ParentType;

	public:
		using ParentType::p1;
		using ParentType::p2;
		using ParentType::m;
		using ParentType::scale_factor_;
		using ParentType::direction_;
		using ParentType::edge_distance_;


		FirstOrderBoundaryShapeFunction( const FunctionSpaceImp& space, typename ParentType::PointInfo pf, typename ParentType::RangeType direction, double scale_factor = 1.0 )
			: ParentType ( space, pf, direction, scale_factor )
		{}

		virtual void evaluate( const typename ParentType::DomainType& arg, typename ParentType::RangeType& ret ) const
		{
			typename ParentType::RangeType tmp = direction_;
			const double dist  = ( arg - m ).two_norm();
			tmp *= scale_factor_ * std::abs( edge_distance_ - dist ) / edge_distance_;
			ret = tmp;
		}
};

//! \todo RENE needs to doc me and move me to stuff
template < class FunctionSpaceImp, class GridPartType >
class SecondOrderBoundaryShapeFunction : public Dune::BoundaryShapeFunctionBase< FunctionSpaceImp, GridPartType >
{
	typedef Dune::BoundaryShapeFunctionBase < FunctionSpaceImp, GridPartType >
		ParentType;

	public:
		using ParentType::p1;
		using ParentType::p2;
		using ParentType::m;
		using ParentType::scale_factor_;
		using ParentType::direction_;
		using ParentType::edge_distance_;


		SecondOrderBoundaryShapeFunction( const FunctionSpaceImp& space, typename ParentType::PointInfo pf, typename ParentType::RangeType direction, double scale_factor = 1.0 )
			: ParentType ( space, pf, direction, scale_factor ),
			parabolic_stretch_( Parameters().getParam("parabolic_stretch", 1.0) )
		{
		}

		virtual void evaluate( const typename ParentType::DomainType& arg, typename ParentType::RangeType& ret ) const
		{
			typename ParentType::RangeType tmp = direction_;
			const double fac = parabolic_stretch_ * scale_factor_ * ( (arg - p1).two_norm() * (arg - p2).two_norm() );
			tmp *= fac;
			ret = tmp;
		}

	protected:
		const double parabolic_stretch_;
};

//! \todo RENE needs to doc me and move me to stuff
template < class FunctionSpaceImp, class GridPartType, template <class ,class> class BoundaryFunctionImp =  FirstOrderBoundaryShapeFunction >
class VariableDirichletData : public Dune::Function < FunctionSpaceImp, VariableDirichletData < FunctionSpaceImp, GridPartType > >
{
	public:
		enum BoundaryType {
			zeroBoundary	= 0,
			influxBoundary	= 1,
			outfluxBoundary	= 2
		};

		typedef VariableDirichletData< FunctionSpaceImp, GridPartType >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;
		typedef std::map< int, BoundaryType >
			BoundaryIdTypeMapType;
		typedef typename BoundaryIdTypeMapType::const_iterator
			BoundaryIdTypeMapTypeConstIterator;
		typedef std::map< int, RangeType>
			ID_ValueMapType;
		typedef Dune::BoundaryInfo< GridPartType >
			BoundaryInfoType;
		typedef BoundaryFunctionImp< FunctionSpaceImp, GridPartType >
			BoundaryFunctionType;
		typedef std::vector<BoundaryFunctionType*>
			BoundaryFunctionListType;
		/**
			*  \brief  constructor
			*
			*
			**/
		VariableDirichletData( const FunctionSpaceImp& space, const GridPartType& gridpart )
			: BaseType( space ),
			gridpart_( gridpart ),
			dim_( FunctionSpaceImp::dimDomain ),
			boundaryInfo_( gridpart )
		{
			zeroBoundaryIds_ 	= Parameters().getList( "zeroBoundaryIds" , 1 );
			influxBoundaryIds_	= Parameters().getList( "influxBoundaryIds" , 2 );
			outfluxBoundaryIds_	= Parameters().getList( "outfluxBoundaryIds" , 3 );
			setupBoundaryIdTypeMap_();
		}

		/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
		~VariableDirichletData()
		{}

		template < class IntersectionIteratorType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionIteratorType& faceIter ) const
		{
			const int id = faceIter.boundaryId();
//			assert( boundaryFunctionList_.size() >= id );
			const BoundaryFunctionType* boundaryFunction = boundaryFunctionList_[id];
			assert ( boundaryFunction );
			boundaryFunction->evaluate( arg, ret );
//			if ( zeroBoundaryIds_.find(id) != zeroBoundaryIds_.end() )
//				ret *=0;
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert(false); }

	protected:
		void setupBoundaryIdTypeMap_() {
			Logger().Info() << "\t- Using ids ";
			for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = zeroBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \t\tfor g_d = 0 \n\t        ids ";
			for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = influxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = - n \n\t        ids ";
			for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = outfluxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = n " << std::endl;
			const double gd_factor = Parameters().getParam( "gd_factor", 1.0 );
			boundaryFunctionList_.reserve( boundaryIdTypeMap_.size() + 1 ); //+1 since BIDs are not 0 indexed
			for ( BoundaryIdTypeMapTypeConstIterator it = boundaryIdTypeMap_.begin(); it != boundaryIdTypeMap_.end(); ++it ) {
				const int id = it->first;
				std::string paramname = std::string( "gd_" ) + Stuff::toString( id );
				std::vector< double > components = Parameters().getList( paramname, double(id) ); //senseless def val here...
				RangeType value;
				assert( components.size() == value.dim() );
				for ( int i = 0; i < value.dim(); ++i )
					value[i] = components[i];
				id_value_map_[ id ] = value;
				boundaryFunctionList_[id] = new BoundaryFunctionType(BaseType::space(), boundaryInfo_.GetPointInfo( id ),  value, gd_factor );
			}
		}

		std::vector< int > zeroBoundaryIds_;
		std::vector< int > influxBoundaryIds_;
		std::vector< int > outfluxBoundaryIds_;
		BoundaryIdTypeMapType boundaryIdTypeMap_;
		ID_ValueMapType	id_value_map_;
		const GridPartType& gridpart_;

		const int dim_;
		BoundaryInfoType boundaryInfo_;
		BoundaryFunctionListType boundaryFunctionList_;
};


#endif // BOUNDARYDATA_HH
