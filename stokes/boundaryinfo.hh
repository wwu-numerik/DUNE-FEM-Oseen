#ifndef BOUNDARYINFO_HH_INCLUDED
#define BOUNDARYINFO_HH_INCLUDED

#include <dune/common/misc.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/math.hh>
#include <map>
#include <set>
#include <utility>

namespace Dune {

//! \todo RENE needs to doc me and move me to stuff
template <class GridPartType>
class BoundaryInfo
{
    public:
        typedef typename GridPartType::GridType
            GridType;
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;
		typedef typename GridPartType::template Codim< 0 >::IteratorType
            EntityIteratorType;
        typedef typename GridPartType::EntityCodim0Type
            EntityType;
        typedef typename GridType::template Codim< 1 >::Geometry
            IntersectionGeometryGlobalType;
        typedef Dune::FieldVector< typename IntersectionGeometryGlobalType::ctype, IntersectionGeometryGlobalType::coorddimension  >
            CoordType;
        typedef std::set<CoordType,int>
            CoordBoundaryIDMapType;
        typedef std::vector< CoordType >
            CoordListType;
        typedef std::pair<CoordType,int>
            GlobalListElementType;
        typedef std::vector< GlobalListElementType >
            GlobalListType;
        typedef std::map< int, std::vector< CoordType > >
            BoundaryCoordListType;
        typedef std::pair<CoordType,CoordType>
            EgdeType;
        typedef std::map<double,CoordType>
            DistancesMapType;
        typedef std::map<double,CoordType>
            VariancesMapType;
        typedef std::map<int,DistancesMapType> //idx in Globallist -> distancesMap
            CoordDistanceMapType;
        typedef std::map<int,CoordType>
            CenterCoordMapType;
        typedef std::map<int,std::pair< CoordType,CoordType > >
            OuterCoordMapType;

        struct PointInfo {
            CoordType center;
            CoordType outmost_1;
            CoordType outmost_2;
        };

        BoundaryInfo(const GridPartType& gridpart)
            : gridpart_(gridpart)
        {
            EntityIteratorType e_it = gridpart_.template begin<0>();
            for( ; e_it != gridpart_.template end<0>(); ++e_it ) {
                const EntityType& e = *e_it;
                IntersectionIteratorType intItEnd = gridpart_.iend( *e_it );
				for (   IntersectionIteratorType intIt = gridpart_.ibegin( *e_it );
						intIt != intItEnd;
						++intIt ) {
					if ( intIt.boundary() ) {
                        const int id = intIt->boundaryId();
                        const IntersectionGeometryGlobalType& globalGeo = intIt.geometry();
                        for ( int i = 0; i < globalGeo.corners(); ++i ) {
                            const CoordType& c = globalGeo[i];
                            boundaryCoordList_[id].push_back( c );
                            globalPointList_.push_back( GlobalListElementType(c,id) );
                        }
					}
				}
            }
            assert( boundaryCoordList_.size() > 0 );
            Logger().Info() << "num BIDs: " << boundaryCoordList_.size() << '\n';
            for (   typename BoundaryCoordListType::const_iterator b_it = boundaryCoordList_.begin();
                    b_it != boundaryCoordList_.end();
                    ++b_it ) {
                CoordDistanceMapType c_dists;
                const typename BoundaryCoordListType::key_type current_boundary_id = b_it->first;
                Logger().Info() << "num points for BID: " << current_boundary_id << " : " << b_it->second.size() << '\n';
                for (   typename CoordListType::const_iterator it = b_it->second.begin();
                        it != b_it->second.end();
                        ++it ) {
                    int idx_in_globallist = Stuff::getIdx( globalPointList_, GlobalListElementType( *it, current_boundary_id ) );
                    c_dists[idx_in_globallist] = getDisctances( b_it->second, *it );
                }
                assert( c_dists.size() > 0 );
                VariancesMapType m;
                for ( typename CoordDistanceMapType::const_iterator c_it = c_dists.begin(); c_it != c_dists.end(); ++c_it ) {
                    const CoordType& vertex = globalPointList_[c_it->first].first;
                    const DistancesMapType& d = c_it->second;
                    const size_t count = d.size();
                    assert ( count > 0 );
                    double variance = 0.0;
                    for (   typename DistancesMapType::const_iterator d_it = d.begin(); d_it != d.end(); ++d_it ) {
                        variance += std::abs(d_it->first)/double(count);
                    }
                    m[variance] = vertex;
                }
                assert( m.size() > 0 );
                const CoordType& center = m.begin()->second;
                int center_idx_in_globallist = Stuff::getIdx( globalPointList_, GlobalListElementType( center, current_boundary_id ) );
                assert( center_idx_in_globallist != -1 );
				Stuff::Logging::LogStream& ss = Logger().Info();
                Stuff::printFieldVector( center,   std::string("center for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                const DistancesMapType& distances_from_center = c_dists[center_idx_in_globallist];
                if( distances_from_center.size() > 1 ) {
                    const CoordType& out1 = distances_from_center.begin()->second;
                    typename DistancesMapType::const_reverse_iterator d_r_it = distances_from_center.rbegin();
					assert ( d_r_it != distances_from_center.rend() );
                    const CoordType& out2 = d_r_it->second;
                    centerCoordMap_[ current_boundary_id ] = center;
                    outerCoordMap_[ current_boundary_id ] = EgdeType( out1, out2 );
                    Stuff::printFieldVector( out1,     std::string("out1   for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                    Stuff::printFieldVector( out2,     std::string("out2   for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                }
                else
                    Logger().Info() << "no dists for BID: "<< current_boundary_id << std::endl;
            }

        }

        PointInfo GetPointInfo( int boundary_id ) const {
            typename OuterCoordMapType::const_iterator out_it = outerCoordMap_.find( boundary_id );
            assert( out_it != outerCoordMap_.end() );
            PointInfo p;
            p.outmost_1 = out_it->second.first;
            p.outmost_2 = out_it->second.second;
            typename CenterCoordMapType::const_iterator cen_it = centerCoordMap_.find( boundary_id );
            assert( cen_it != centerCoordMap_.end() );
            p.center = cen_it->second;
            return p;
        }

        /** Default destructor */
        virtual ~BoundaryInfo() {}

    protected:

        const GridPartType& gridpart_;
//        CoordBoundaryIDMapType coordBoundaryIDMap_;
        BoundaryCoordListType boundaryCoordList_;
        CenterCoordMapType centerCoordMap_;
        OuterCoordMapType outerCoordMap_;
        GlobalListType globalPointList_;

        DistancesMapType getDisctances( const CoordListType& list, const CoordType& origin ) {
            DistancesMapType ret;
            for (   typename CoordListType::const_iterator it = list.begin();
                        it != list.end();
                        ++it ) {
                if ( *it == origin )
                    continue;
                typename DistancesMapType::key_type dist = ( *it - origin ).two_norm();
                assert( dist > 0 );
                long sign = Stuff::sign( long(origin * *it) );
                if ( sign == 0 )
                    sign = origin.two_norm() > it->two_norm() ? 1 : -1;
                ret[sign*dist] = *it;
            }
            return ret;
        }
};//end BoundaryInfo

//! \todo RENE needs to doc me and move me to stuff
template < class FunctionSpaceImp, class GridPartType >
class BoundaryShapeFunctionBase : public Dune::Fem::Function < FunctionSpaceImp, BoundaryShapeFunctionBase < FunctionSpaceImp, GridPartType > >
{
	public:
		typedef BoundaryShapeFunctionBase< FunctionSpaceImp, GridPartType >
			ThisType;
		typedef Dune::Fem::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;
		typedef typename Dune::BoundaryInfo<GridPartType>::PointInfo
			PointInfo;
		BoundaryShapeFunctionBase( const FunctionSpaceImp& space, PointInfo pf, RangeType direction, double scale_factor = 1.0 )
			: BaseType( space ),
			pointInfo_( pf ),
			direction_( direction ),
			m( pointInfo_.center ),
			p1( pointInfo_.outmost_1 ),
			p2( pointInfo_.outmost_2 ),
			scale_factor_( scale_factor ),
			edge_distance_( ( (p1 - m).two_norm() + (p2 - m).two_norm() ) / 2.0 )
		{
		}

        virtual void evaluate( const DomainType& /*arg*/, RangeType& ret ) const
		{
			ret = direction_;
		}

	protected:
		const PointInfo pointInfo_;
		const RangeType direction_;
		const DomainType& m;
		const DomainType& p1;
		const DomainType& p2;
		const double scale_factor_;
		const double edge_distance_;

};

//! \todo RENE needs to doc me and move me to stuff
template < template <class,class, template <class,class> class> class AnalyticalDirichletDataImp,
			template <class,class> class BoundaryFunctionImp =  BoundaryShapeFunctionBase >
struct GeometryBasedBoundaryFunctionTraits {

	template < class FunctionSpaceImp, class GridPartImp >
	struct Implementation {
		typedef AnalyticalDirichletDataImp< FunctionSpaceImp, GridPartImp, BoundaryFunctionImp >
			AnalyticalDirichletDataType;
		typedef typename AnalyticalDirichletDataType::BoundaryFunctionType
			BoundaryFunctionType;
		template <class DiscreteStokesFunctionWrapper >
		static AnalyticalDirichletDataType getInstance( const DiscreteStokesFunctionWrapper& wrapper ) {
			return 	AnalyticalDirichletDataType( wrapper.discreteVelocitySpace(), wrapper.discreteVelocitySpace().gridPart() );
		}
	};
};

}//end namespace Dune

#endif // BOUNDARYINFO_HH_INCLUDED
