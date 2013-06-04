#ifndef BOUNDARYDATA_HH
#define BOUNDARYDATA_HH

#include <dune/fem/oseen/boundaryinfo.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/fem/functions/timefunction.hh>
#include <dune/stuff/grid/entity.hh>

template < template < class > class DiricheltDataImp >
struct DefaultDirichletDataTraits {

	template < class FunctionSpaceImp, class GridPartImp >
	struct Implementation {
		typedef DiricheltDataImp< FunctionSpaceImp>
			AnalyticalDirichletDataType;

		template <class DiscreteOseenFunctionWrapper >
		static AnalyticalDirichletDataType getInstance( const DiscreteOseenFunctionWrapper& wrapper ) {
			return 	AnalyticalDirichletDataType( wrapper.discreteVelocitySpace() );
		}
	};
};

//! \todo RENE needs to doc me and move me to stuff
namespace Stuff {
class BoundaryIdMapping {
    public:
	enum BoundaryType {
		zeroBoundary	= 0,
		influxBoundary	= 1,
		outfluxBoundary	= 2
	};
	typedef std::map< int, BoundaryType >
		BoundaryIdTypeMapType;
	typedef typename BoundaryIdTypeMapType::const_iterator
		BoundaryIdTypeMapTypeConstIterator;

	BoundaryIdMapping()
	{
        zeroBoundaryIds_ 	= DSC_CONFIG.getList( "zeroBoundaryIds" , int(1) );
        influxBoundaryIds_	= DSC_CONFIG.getList( "influxBoundaryIds" , 2 );
        outfluxBoundaryIds_	= DSC_CONFIG.getList( "outfluxBoundaryIds" , 3 );
	    setupBoundaryIdTypeMap_();
	}
    protected:
	void setupBoundaryIdTypeMap_() {
		DSC_LOG_INFO << "\t- Using ids ";
		for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
			boundaryIdTypeMap_[*it] = zeroBoundary;
			DSC_LOG_INFO << *it << " ";
		}
		DSC_LOG_INFO << " \t\tfor g_d = 0 \n\t        ids ";
		for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
			boundaryIdTypeMap_[*it] = influxBoundary;
			DSC_LOG_INFO << *it << " ";
		}
		DSC_LOG_INFO << " \tfor g_d = - n \n\t        ids ";
		for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
			boundaryIdTypeMap_[*it] = outfluxBoundary;
			DSC_LOG_INFO << *it << " ";
		}
		DSC_LOG_INFO << " \tfor g_d = n " << std::endl;
	}

	std::vector< int > zeroBoundaryIds_;
	std::vector< int > influxBoundaryIds_;
	std::vector< int > outfluxBoundaryIds_;
	BoundaryIdTypeMapType boundaryIdTypeMap_;
};

template < class FunctionSpaceImp >
class BoundaryFluxFunction : public Dune::Fem::Function < FunctionSpaceImp, BoundaryFluxFunction < FunctionSpaceImp > >,
	public BoundaryIdMapping
{
    public:

		typedef BoundaryFluxFunction< FunctionSpaceImp >
			ThisType;
		typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;

		BoundaryFluxFunction( const FunctionSpaceImp& space )
			: BaseType( space ),
				dim_( FunctionSpaceImp::dimDomain )
		{}

		~BoundaryFluxFunction()
		{}

		template < class IntersectionIteratorType >
        void evaluate( const DomainType& /*arg*/, RangeType& ret, const IntersectionIteratorType& faceIter ) const
		{
			const int id = faceIter.boundaryId();

			typedef Dune::FieldVector< typename IntersectionIteratorType::ctype, IntersectionIteratorType::dimension - 1 >
				LocalVectorType;

            const auto& geometry = faceIter.intersectionSelfLocal();
            LocalVectorType center = geometry.local(geometry.center());
			RangeType normal = faceIter.unitOuterNormal( center );
			static const double gd_factor = DSC_CONFIG_GET( "gd_factor", 1.0 );
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

        inline void evaluate( const DomainType& /*arg*/, RangeType& /*ret*/ ) const { assert(false); }

	protected:
		const int dim_;
};

template < class FunctionSpaceImp, class TimeProviderImp >
class InstationaryBoundaryFluxFunction :
    public DSFe::IntersectionTimeFunction < FunctionSpaceImp ,
			InstationaryBoundaryFluxFunction< FunctionSpaceImp,TimeProviderImp >,
			TimeProviderImp > ,
	public BoundaryIdMapping
{
    public:
	typedef InstationaryBoundaryFluxFunction< FunctionSpaceImp, TimeProviderImp >
		ThisType;
    typedef DSFe::IntersectionTimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
		BaseType;
	typedef typename BaseType::DomainType
		DomainType;
	typedef typename BaseType::RangeType
		RangeType;

	InstationaryBoundaryFluxFunction( const TimeProviderImp& timeprovider,
				   const FunctionSpaceImp& space,
				   const double /*constant*/ = 0.0 )
		: BaseType ( timeprovider, space )
	{}

	~InstationaryBoundaryFluxFunction()
	{}

	void evaluateTime( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
	{
	    ret = RangeType( 0 );
	}

	template < class IntersectionType >
	void evaluateTime( const double time, const DomainType& /*arg*/, RangeType& ret, const IntersectionType& intersection ) const
	{
	    const int id = intersection.boundaryId();

	    typedef Dune::FieldVector< typename IntersectionType::ctype, IntersectionType::dimension - 1 >
		    LocalVectorType;

        const auto& geometry = intersection;
        LocalVectorType center = geometry.local(geometry.center());
	    RangeType normal = intersection.unitOuterNormal( center );
	    const double gd_factor = time * DSC_CONFIG_GET( "gd_factor", 1.0 );
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
};

}// namespace Stuff

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
			parabolic_stretch_( DSC_CONFIG_GET("parabolic_stretch", 1.0) )
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
class VariableDirichletData : public Dune::Fem::Function < FunctionSpaceImp, VariableDirichletData < FunctionSpaceImp, GridPartType > >
{
	public:
		enum BoundaryType {
			zeroBoundary	= 0,
			influxBoundary	= 1,
			outfluxBoundary	= 2
		};

		typedef VariableDirichletData< FunctionSpaceImp, GridPartType >
			ThisType;
		typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
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
            zeroBoundaryIds_ 	= DSC_CONFIG.getList( "zeroBoundaryIds" , 1 );
            influxBoundaryIds_	= DSC_CONFIG.getList( "influxBoundaryIds" , 2 );
            outfluxBoundaryIds_	= DSC_CONFIG.getList( "outfluxBoundaryIds" , 3 );
            setupBoundaryIdTypeMap_();
		}

		/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
		~VariableDirichletData()
		{
		//TODO leak galore >_>
		}

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

        inline void evaluate( const DomainType& /*arg*/, RangeType& /*ret*/ ) const { assert(false); }

	protected:
		void setupBoundaryIdTypeMap_() {
			DSC_LOG_INFO << "\t- Using ids ";
			for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = zeroBoundary;
				DSC_LOG_INFO << *it << " ";
			}
			DSC_LOG_INFO << " \t\tfor g_d = 0 \n\t        ids ";
			for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = influxBoundary;
				DSC_LOG_INFO << *it << " ";
			}
			DSC_LOG_INFO << " \tfor g_d = - n \n\t        ids ";
			for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = outfluxBoundary;
				DSC_LOG_INFO << *it << " ";
			}
			DSC_LOG_INFO << " \tfor g_d = n " << std::endl;
			const double gd_factor = DSC_CONFIG_GET( "gd_factor", 1.0 );
			boundaryFunctionList_.reserve( boundaryIdTypeMap_.size() + 1 ); //+1 since BIDs are not 0 indexed
			for ( BoundaryIdTypeMapTypeConstIterator it = boundaryIdTypeMap_.begin(); it != boundaryIdTypeMap_.end(); ++it ) {
				const int id = it->first;
				std::string paramname = std::string( "gd_" ) + DSC::toString( id );
                std::vector< double > components = DSC_CONFIG.getList( paramname, double(id) ); //senseless def val here...
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

