#ifndef BOUNDARYINFO_HH_INCLUDED
#define BOUNDARYINFO_HH_INCLUDED

#include <dune/common/misc.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <map>
#include <set>
#include <utility>

namespace Dune {

template <class GridPartType>
class BoundaryInfo
{
    public:
        Dune::CompileTimeChecker< 2 == GridPartType::GridType::dimension > BoundaryInfo_not_implemented_for_dim_neq_2;

        typedef typename GridPartType::GridType
            GridType;
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;
        typedef typename GridPartType::Traits::template Codim< 0 >::IteratorType
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
                        const IntersectionGeometryGlobalType& globalGeo = intIt.intersectionGlobal();
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
                if ( b_it->first != 2 )
                    continue;
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
                Logging::LogStream& ss = Logger().Info();
                Stuff::printFieldVector( center,   std::string("center for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                const DistancesMapType& distances_from_center = c_dists[center_idx_in_globallist];
                if( distances_from_center.size() > 1 ) {
                    const CoordType& out1 = distances_from_center.begin()->second;
                    typename DistancesMapType::const_reverse_iterator d_r_it = distances_from_center.rbegin();
                    //assert ( d_r_it != distances_from_center.begin() );
                    const CoordType& out2 = d_r_it->second;
                    centerCoordMap_[ current_boundary_id ] = center;
                    outerCoordMapType_[ current_boundary_id ] = EgdeType( out1, out2 );
                    Stuff::printFieldVector( out1,     std::string("out1   for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                    Stuff::printFieldVector( out2,     std::string("out2   for id: ") + Stuff::toString( current_boundary_id ), ss, "BID --- " );
                }
                else
                    Logger().Info() << "no dists for BID: "<< current_boundary_id << std::endl;
            }

        }
        /** Default destructor */
        virtual ~BoundaryInfo() {}

    protected:
        BoundaryInfo(const BoundaryInfo& other) {}
        BoundaryInfo& operator=(const BoundaryInfo& other) { return *this; }

        const GridPartType& gridpart_;
//        CoordBoundaryIDMapType coordBoundaryIDMap_;
        BoundaryCoordListType boundaryCoordList_;
        CenterCoordMapType centerCoordMap_;
        OuterCoordMapType outerCoordMapType_;
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

}//end namespace Dune

#endif // BOUNDARYINFO_HH_INCLUDED
