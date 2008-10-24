#ifndef DUNE_DEFAULTINDEXSETS_HH
#define DUNE_DEFAULTINDEXSETS_HH

//- system includes 
#include <vector>
#include <rpc/rpc.h>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/common/interfaces.hh>
#include <dune/grid/common/grid.hh>

/** @file
 @author Robert Kloefkorn
 @brief Provides default index set implementations for Level- and
 LeafIndexsets used by AlbertaGrid and ALUGrid. 
*/

namespace Dune {

/*!
  The DefaultGridIndexSet is a wrapper for the grid index which can be 
  index of entity or globalIndex. The DofMapper uses an IndexSet for
  mapping the dofs, so we can hide the real grid index behind the index
  set. Furthermore if a grid doesn't provide the consecutive index set
  then this can be calculated in the IndexSet. These two following index
  sets are just the identiy to the grid indices. 

  The DefaultGridIndexSetBase defines some methods that are needed for
  index sets that cope with adaptation, but aren't needed for the following 
  index set, so most of this methods do notin'.
*/

class DefaultEmptyIndexSet
{
  // dummy value 
  enum { myType = -1 };
protected:  
  const bool adaptive_;
public:  
  template <class IndexSetType, class EntityType,int enCodim, int codim>
  struct IndexWrapper 
  {
    static inline int index (const IndexSetType & set, const EntityType & en , int num )
    {
      return set.index(en);
    }
  };

  // if codim and enCodim are equal, then return index 
  template <class IndexSetType, class EntityType, int codim>
  struct IndexWrapper<IndexSetType,EntityType,codim,codim>
  {
    static inline int index (const IndexSetType & set, const EntityType & en , int num )
    {
      return set.index(en);
    }
  };

  // return number of vertex num
  template <class IndexSetType, class EntityType>
  struct IndexWrapper<IndexSetType,EntityType,0,3> 
  {
    static inline int index (const IndexSetType & set, const EntityType & en , int num )
    {
      return set.template subIndex<3> (en,num);
    }
  };
 
  // return number of vertex num for dim == 2 
  // return number of edge num for dim == 3
  template <class IndexSetType, class EntityType>
  struct IndexWrapper<IndexSetType,EntityType,0,2> 
  {
    static inline int index (const IndexSetType & set, const EntityType & en , int num )
    {
      return set.template subIndex<2> (en,num);
    }
  };

  // return number of vertex for dim == 1
  // return number of edge num for dim == 2 
  // return number of face num for dim == 3
  template <class IndexSetType, class EntityType>
  struct IndexWrapper<IndexSetType,EntityType,0,1> 
  {
    static inline int index (const IndexSetType & set, const EntityType & en , int num )
    {
      return set.template subIndex<1> (en,num);
    }
  };

  //! default constructor 
  DefaultEmptyIndexSet (bool adaptive) : adaptive_(adaptive) {}

  //! return false mean the no memory has to be allocated 
  //! and no compress of date has to be done 
  bool compress () { 
    return false; 
  }

  //! \brief returns true if index set can be used for adaptive 
  //! calculations 
  bool adaptive () const { return adaptive_; }
  
  //! returns true if index set gernally needs compress after adaptation 
  bool needsCompress () const { return false; }

  //! do nothing here, because fathers index should already exist 
  template <class EntityType> 
  void insertNewIndex( const EntityType & en ) {
    assert(adaptive_) ;
  }

  //! do nothing here, because fathers index should already exist 
  template <class EntityType> 
  void removeOldIndex( const EntityType & en ) {
    assert(adaptive_) ;
  }

  //! nothing to do here 
  void resize () {}

  //! no extra memory for restriction is needed
  int additionalSizeEstimate () const { return 0; }

  static int type() { return myType; }

  //! we have no old size 
  int numberOfHoles ( int codim ) const { return 0; }
  
  //! return old index, for dof manager only 
  int oldIndex (int hole, int codim ) const { return 0; }
  
  //! return new index, for dof manager only 
  int newIndex (int hole, int codim ) const { return 0; }
  
  //! write index set to xdr file 
  bool write_xdr(const std::basic_string<char> filename , int timestep) 
  {
    FILE  *file;
    XDR   xdrs;
    const char *path = "";

    std::basic_string<char> fnstr  = genFilename(path,filename, timestep);
    const char * fn = fnstr.c_str();
    file = fopen(fn, "wb");
    if (!file)
    {
        std::cerr << "\aERROR in DefaultGridIndexSet::write_xdr(..): could not open <"
                  << filename << ">!" << std::endl;
      return false;
    }

    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);

    return true;
  }
  
  //! read index set to xdr file 
  bool read_xdr(const std::basic_string<char> filename , int timestep) 
  {
    FILE   *file;
    XDR     xdrs;
    const char *path = "";

    std::basic_string<char> fnstr = genFilename(path,filename, timestep);
    const char * fn  = fnstr.c_str();
    std::cout << "Reading <" << fn << "> \n";
    file = fopen(fn, "rb");
    if(!file)
    {
      std::cerr << "\aERROR in DefaultGridIndexSet::read_xdr(..): could not open <" 
                << filename << ">!" << std::endl;
      return(false);
    }

    // read xdr 
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    return true;
  }

protected: 
  // read/write from/to xdr stream 
  bool processXdr(XDR *xdrs)
  {
    int type = myType; 
    xdr_int ( xdrs, &type);
    if(type != myType)
    {
      std::cerr << "\nERROR: DefaultGridIndex: wrong type choosen! \n\n";
      assert(type == myType);
    }
    return true;
  }
};

//! Wraps the interface methods of indexsets and adds the addiotnal needed
//! functions 
template <class IndexSetImp> 
class IndexSetWrapper : public DefaultEmptyIndexSet 
{
public:
  //! The types of the iterator 
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename IndexSetImp::template Codim<cd>::
        template Partition<pitype>::Iterator  Iterator;
    };
  };

  //! store const reference to set 
  IndexSetWrapper(const IndexSetImp & set, bool adaptive = false) 
    : DefaultEmptyIndexSet(adaptive)
    , set_(set) 
  {}
  
  //! store const reference to set 
  IndexSetWrapper(const IndexSetWrapper<IndexSetImp> & s) 
    : DefaultEmptyIndexSet(s.adaptive_), set_(s.set_) {}
 
  //! return size of set for codim  
  int size ( GeometryType type ) const   
  {
    return set_.size(type);
  }

  //! return size of grid entities per level and codim 
  int size ( int codim ) const   
  {
    return set_.size(codim);
  }

  //! return index of en 
  template <class EntityType> 
  int index (const EntityType & en) const
  {
    return set_.index(en); 
  }

  //! return sub index of given entities sub entity with codim and number 
  template <int codim,class EntityType> 
  int subIndex (const EntityType & en, int num) const
  {
    return set_.template subIndex<codim> (en,num);
  }

  //! wrap geomTypes method of set 
  const std::vector< GeometryType > & geomTypes (int codim) const 
  {
    return set_.geomTypes(codim); 
  }

  //! returns true if this set provides an index for given entity
  template<class EntityType>
  bool contains (const EntityType& en) const
  {
    return set_.contains(en); 
  }

  //! return index 
  template <int codim, class EntityType> 
  int index (const EntityType & en, int num) const
  {
    enum { enCodim = EntityType::codimension };
    return IndexWrapper<IndexSetImp,EntityType,enCodim,codim>::index(set_,en,num);
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin () const
  { 
    return set_.template begin<cd,pitype> ();
  }

  /** @brief Iterator to one past the last entity of given codim for partition type
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end () const
  { 
    return set_.template end<cd,pitype> ();
  }

private: 
  const IndexSetImp & set_; 
};

//! default base class for index set implementations for FR numerics
template <class GridType>
class DefaultGridIndexSetBase : public DefaultEmptyIndexSet
{
  // dummy value 
  enum { myType = -1 };
public:  
  enum { ncodim = GridType::dimension + 1 };

  //! Conschdrugdor 
  DefaultGridIndexSetBase (const GridType & grid ) 
    // here false, because methods have to be overloaded
    : DefaultEmptyIndexSet(false) 
    , grid_ (grid) {}
protected:
  // the corresponding grid 
  const GridType & grid_;
};

//! Wraps LevelIndexSet for use with LagrangeFunctionSpace 
template <class GridType>
class WrappedLevelIndexSet 
: public IndexSetWrapper< typename GridType :: Traits:: LevelIndexSet > 
{
  // my type, to be revised 
  enum { myType = 1 };

  typedef typename GridType :: Traits :: LevelIndexSet LevelIndexSetType;
  typedef WrappedLevelIndexSet < GridType >  ThisType;
public:
  
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };

  //! Constructor getting grid and level for Index Set 
  WrappedLevelIndexSet (const GridType & grid , const int level ) 
    : IndexSetWrapper<  LevelIndexSetType > (grid.levelIndexSet(level)) {}
 
  //! return type of index set (for input/output)
  static int type() { return myType; }
  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid,grid.maxLevel());
    return set;
  }

};

//! compile time chooser for hierarchic or leaf index set
template <class GridImp> 
class HierarchicIndexSetSelector
{
  
  // true if GridImp has HierarchicIndexSet 
  enum { hasHierarchicIndexSet = Conversion<GridImp,HasHierarchicIndexSet>::exists };

  template <class GridType, bool hasHSet> 
  struct HSetChooser
  {
    typedef typename GridType::Traits::LeafIndexSet IndexSetType;
    static const IndexSetType & hierarchicIndexSet(const GridType & grid)  
    { 
      return grid.leafIndexSet();
    }
  };
  
  template <class GridType> 
  struct HSetChooser<GridType,true>
  {
    typedef typename GridImp:: HierarchicIndexSet IndexSetType;
    static const IndexSetType & hierarchicIndexSet(const GridType & grid)  
    { 
      return grid.hierarchicIndexSet();
    }
  };
  
public: 
  //! \brief type of HierarchicIndexSet, default is LeafIndexSet 
  typedef typename HSetChooser<GridImp,hasHierarchicIndexSet>::IndexSetType HierarchicIndexSet; 
 
  //! \brief return reference to hierarchic index set 
  static const HierarchicIndexSet & hierarchicIndexSet(const GridImp & grid) 
  { 
    return HSetChooser<GridImp,hasHierarchicIndexSet>::hierarchicIndexSet(grid);
  }

  //! return true if index set can be used for adapitve calculations 
  static bool adaptive () { return hasHierarchicIndexSet; }
};

//! Wraps HierarchicIndex Sets of AlbertaGrid and ALUGrid 
template <class GridType>
class WrappedHierarchicIndexSet
: public IndexSetWrapper< typename HierarchicIndexSetSelector<GridType> :: HierarchicIndexSet >
{
  // my type, to be revised 
  enum { myType = 0 };

  //! type of hset selector 
  typedef HierarchicIndexSetSelector<GridType> SelectorType;
  
  // my index set type 
  typedef typename SelectorType :: HierarchicIndexSet HSetType;

  typedef WrappedHierarchicIndexSet<GridType> ThisType;
public:
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };

  //! constructor 
  WrappedHierarchicIndexSet ( const GridType & grid , const int level =-1 ) 
    : IndexSetWrapper< HSetType > ( SelectorType :: hierarchicIndexSet(grid) 
                                  , SelectorType :: adaptive() ) {}
     
  //! return type (for Grape In/Output)
  static int type() { return myType; }

  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid);
    return set;
  }

};

//! Wraps LeafIndexSet of Dune Grids for use with LagrangeFunctionSpace 
template <class GridType>
class WrappedLeafIndexSet
:  public IndexSetWrapper<typename GridType :: Traits :: LeafIndexSet> 
{
  // my type, to be revised 
  enum { myType = 5 };

  // my index set type 
  typedef typename GridType :: Traits :: LeafIndexSet IndexSetType;
  typedef WrappedLeafIndexSet<GridType> ThisType;
public:
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };
  //! constructor 
  WrappedLeafIndexSet ( const GridType & grid , const int level =-1 ) 
    : IndexSetWrapper < IndexSetType > (grid.leafIndexSet()) {}

  //! return type (for Grape In/Output)
  static int type() { return myType; }
  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid,grid.maxLevel());
    return set;
  }
};

//*********************************************************************
/*! \brief 
 * DefaultLevelIndexSet generates a level index set for a grid out of a
 * grids hierarchic index set by storing for each entity in the grid 
 * a number in an array. 
 */
//*********************************************************************

template <class DefaultLevelIndexSetType, int codim>
struct CheckLevelForCodim 
{
  static void recursive(DefaultLevelIndexSetType & d)
  {
    d.template checkLevelIndexForCodim<codim> (); 
    CheckLevelForCodim<DefaultLevelIndexSetType,codim-1>::recursive(d);
  }
};

template <class DefaultLevelIndexSetType>
struct CheckLevelForCodim<DefaultLevelIndexSetType,0>
{
  static void recursive(DefaultLevelIndexSetType & d)
  {
    d.template checkLevelIndexForCodim<0> ();
  }
};


template <class GridImp>
struct DefaultLevelIteratorTypes
{
  //! The types
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename GridImp::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator Iterator;
    };
  };
};

/*! \brief 
  DefaultLevelIndexSet creates a LevelIndexSet for a Grid by using its
  HierarchicIndexSet 
 */ 
template <class GridImp>
class DefaultLevelIndexSet :
  public IndexSetDefaultImplementation< GridImp, 
                   DefaultLevelIndexSet <GridImp>,
                   DefaultLevelIteratorTypes<GridImp> > 

{
  typedef GridImp GridType;
  enum { dim = GridType :: dimension };

  //! type of used arrays to store indices 
  typedef std::vector<int> IndexArrayType;
  
public:
  enum { ncodim = GridType::dimension + 1 };
private:
  typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

  typedef DefaultLevelIndexSet<GridType> ThisType;

  template <class EntityType, int codim> 
  struct InsertEntity
  {
    template <class HierarchicIndexSet>
    static void insert(const EntityType & en, 
                       const HierarchicIndexSet & hset, 
                       IndexArrayType (&index)[ncodim], 
                       int (&num)[ncodim])
    {
      IndexArrayType & idx = index[codim];
      for(int i=0; i<en.template count<codim>(); ++i) 
      {
        const int id = hset.template subIndex<codim>(en,i);
        if( idx[id] < 0) 
        {
          idx[id] = num[codim];
          ++num[codim];
        }
      }
      InsertEntity<EntityType,codim-1>::insert(en,hset,index,num);
    }
  };
  
  template <class EntityType> 
  struct InsertEntity<EntityType,0>
  {
    template <class HierarchicIndexSet>
    static void insert(const EntityType & en, 
                       const HierarchicIndexSet & hset, 
                       IndexArrayType (&index)[ncodim], 
                       int (&num)[ncodim])
    {
      enum { codim = 0 };
      IndexArrayType & idx = index[codim];
      const int id = hset.index(en); 
      if( idx[id] < 0 ) 
      {
        idx[id] = num[codim];
        ++num[codim];
      }
    }
  };
  
public:
  //! create LevelIndex by using the HierarchicIndexSet of a grid 
  //! for the given level 
  DefaultLevelIndexSet(const GridType & grid , int level ) :
    grid_(grid) , level_(level) , hIndexSet_ ( grid.hierarchicIndexSet() )
    , size_ ( ncodim )
  {
    calcNewIndex ();
  }

  //! return LevelIndex of given entity 
  template<class EntityType>
  int index (const EntityType & en) const 
  {
    enum { cd = EntityType :: codimension };
    // this must not be true for vertices 
    // therefore only check other codims 
#ifndef NDEBUG 
    const int codim = cd;
    assert( (codim == dim) ? (1) : (level_ == en.level() )); 
    assert( levelIndex_[codim][ hIndexSet_.index(en) ] >= 0 );
#endif
    return levelIndex_[cd][ hIndexSet_.index(en) ];
  }
 
  //! return LevelIndex of given entity 
  template<int cd>
  int index (const typename GridImp::template Codim<cd>::Entity& en) const 
  {
    // this must not be true for vertices 
    // therefore only check other codims 
#ifndef NDEBUG 
    const int codim = cd;
    assert( (codim == dim) ? (1) : (level_ == en.level() )); 
    assert( levelIndex_[codim][ hIndexSet_.index(en) ] >= 0 );
#endif
    return levelIndex_[cd][ hIndexSet_.index(en) ];
  }
 
  //! return subIndex (LevelIndex) for a given Entity of codim = 0 and a
  //! given SubEntity codim and number of SubEntity
  template <int cd>
  int subIndex (const typename GridType::template Codim<0>::Entity & en, int i) const
  {
    // this must not be true for vertices 
    // therefore only check other codims 
#ifndef NDEBUG
    const int codim = cd;
    assert( (codim != 0) ? (1) : (level_ == en.level() )); 
    assert(levelIndex_[codim][ hIndexSet_.template subIndex<cd>(en,i) ] >= 0 );
#endif
    return levelIndex_[cd][ hIndexSet_.template subIndex<cd>(en,i) ];
  }

  //! returns true if this set provides an index for given entity
  template<class EntityType>
  bool contains (const EntityType& en) const
  {
    enum { cd = EntityType :: codimension };
    return (levelIndex_[cd][ hIndexSet_.index(en) ] >= 0);
  }

  //! return size of IndexSet for a given level and codim 
  int size ( int codim ) const
  {
    assert( codim >= 0 && codim <= GridType::dimension );
    return size_[codim];
  }

  //! return size of IndexSet for a given level and codim 
  //! this method is to be revised 
  int size ( GeometryType type ) const
  {
    if( typeNotValid(type) ) return 0;
    return size_[GridType::dimension-type.dim()];
  }

  //! do calculation of the index set, has to be called when grid was
  //! changed or if index set is created 
  void calcNewIndex ()
  {
    // resize arrays to new size 
    for(int cd=0; cd<ncodim; ++cd) 
    {
      resizeVectors(levelIndex_[cd], hIndexSet_.size(cd));
    }

    // walk grid and store index 
    typedef typename DefaultLevelIteratorTypes<GridImp>:: template Codim<0>::
        template Partition<All_Partition> :: Iterator IteratorType;

    // we start with zero for all codims 
    int num[ncodim];
    for(int cd=0; cd<ncodim; ++cd) num[cd] = 0;

    IteratorType endit  = this->template end  < 0, All_Partition > ();
    for(IteratorType it = this->template begin< 0, All_Partition > (); 
        it != endit; ++it)
    {
      assert( it->level() == level_ );
      insertEntity(*it,num);
    }

    // remember the number of entity on level and cd = 0
    for(int cd=0; cd<ncodim; ++cd) 
    {
      size_[cd] = num[cd];
      assert( size_[cd] == grid_.size(level_,cd) );
    }

#ifndef NDEBUG 
    CheckLevelForCodim<ThisType,dim>::recursive(*this);
#endif
  }

  // calculate index for the codim 
  template <int cd>
  void checkLevelIndexForCodim ()
  {
    IndexArrayType & levIndex = levelIndex_[cd];
    // resize memory if necessary 
    // walk grid and store index 
    typedef typename DefaultLevelIteratorTypes<GridImp>:: template Codim<cd>::
      template Partition<All_Partition> :: Iterator LevelIterator;

    LevelIterator endit  = this->template end  < cd , All_Partition > ();
    for(LevelIterator it = this->template begin< cd , All_Partition > (); it != endit; ++it)
    { 
      int no = hIndexSet_.index(*it);
      assert( levIndex[no] >= 0 );
    }
  }

  //! deliver all geometry types used in this grid
  const std::vector<GeometryType>& geomTypes (int codim) const
  {
    return grid_.geomTypes(codim);
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLevelIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
  { 
    return this->grid_.template lbegin<cd,pitype> (level_);
  }

  /** @brief Iterator to one past the last entity of given codim for partition type
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLevelIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
  { 
    return this->grid_.template lend<cd,pitype> (level_);
  }

  //! returns true if this set provides an index for given entity
  bool containsIndex (int cd , int idx) const
  {
    assert( cd >= 0 );
    assert( cd < ncodim );
    
    assert( idx >= 0);
    assert( idx < (int) levelIndex_[cd].size());
    return (levelIndex_[cd][ idx ] >= 0);
  }

private:
  // return whether set has this type stored or not 
  bool typeNotValid (const GeometryType & type) const 
  {
    int codim = GridType :: dimension - type.dim();
    const std::vector<GeometryType> & geomT = geomTypes(codim); 
    for(size_t i=0; i<geomT.size(); ++i) if(geomT[i] == type) return false;
    return true;  
  }
  // calculate index for the codim 
  template <class EntityType>
  void insertEntity(EntityType & en, int (&num)[ncodim])
  {
    InsertEntity<EntityType,dim>::insert(en,hIndexSet_,levelIndex_,num);
  }

  // resize vectors of index set 
  void resizeVectors(IndexArrayType &a, int newNumberOfEntries)
  {
    if(newNumberOfEntries > 0)
    {
      a.resize(newNumberOfEntries);
    }
    for(size_t i=0; i<a.size(); i++) a[i] = -1;
  }

  // method prints indices of given codim, for debugging
  void print (int codim) const 
  {
    for(size_t i=0; i<levelIndex_[codim].size(); i++) 
    {
      std::cout << "levelind[" << i << "] = " << levelIndex_[codim][i] << "\n";
    }
  }

  // grid this level set belongs to 
  const GridType & grid_;

  // the level for which this index set is created 
  const int level_;
  
  // the grids HierarchicIndexSet
  const HierarchicIndexSetType & hIndexSet_;

  // number of entitys of each level an codim 
  IndexArrayType size_;

  //*********************************************************
  // Methods for mapping the hierarchic Index to index on Level
  IndexArrayType levelIndex_[ncodim];
  
};


//! LeafIterator tpyes for all codims and partition types
template <class GridImp>
struct DefaultLeafIteratorTypes
{
  //! The types of the iterator 
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename GridImp::Traits::template Codim<cd>::
        template Partition<pitype>::LeafIterator  Iterator;
    };
  };
};


//! Default LeafIndexSet 
template <class GridImp>
class DefaultLeafIndexSet :
  public IndexSetDefaultImplementation< GridImp, 
                   DefaultLeafIndexSet <GridImp>,
                   DefaultLeafIteratorTypes<GridImp> > 

{
  typedef GridImp GridType;
  enum { dim = GridType :: dimension };

  //! type of used arrays to store indices 
  typedef std::vector<int> IndexArrayType;
public:
  enum { ncodim = dim + 1 };
  typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

private:
  typedef DefaultLeafIndexSet<GridType> ThisType;

  template <class EntityType, int codim> 
  struct InsertEntity
  {
    template <class HierarchicIndexSet>
    static void insert(const EntityType & en, 
                       const HierarchicIndexSet & hset, 
                       IndexArrayType (&index)[ncodim], 
                       int (&num)[ncodim])
    {
      IndexArrayType & idx = index[codim];
      for(int i=0; i<en.template count<codim>(); ++i) 
      {
        const int id = hset.template subIndex<codim>(en,i);
        if( idx[id] < 0) 
        {
          idx[id] = num[codim];
          ++num[codim];
        }
      }
      InsertEntity<EntityType,codim-1>::insert(en,hset,index,num);
    }
  };
  
  template <class EntityType> 
  struct InsertEntity<EntityType,0>
  {
    template <class HierarchicIndexSet>
    static void insert(const EntityType & en, 
                       const HierarchicIndexSet & hset, 
                       IndexArrayType (&index)[ncodim], 
                       int (&num)[ncodim])
    {
      enum { codim = 0 };
      IndexArrayType & idx = index[codim];
      const int id = hset.index(en); 
      if( idx[id] < 0 ) 
      {
        idx[id] = num[codim];
        ++num[codim];
      }
    }
  };
  
public:
  //! create LevelIndex by using the HierarchicIndexSet of a grid 
  //! for the given level 
  DefaultLeafIndexSet(const GridType & grid) 
    : grid_(grid) 
    , hIndexSet_ ( grid.hierarchicIndexSet() )
    , size_ ( ncodim )
  {
    calcNewIndex ();
  }

  //! return LeafIndex of given entity 
  template<class EntityType>
  int index (const EntityType & en) const 
  {
    enum { cd = EntityType :: codimension };
    // this must not be true for vertices 
    // therefore only check other codims 
    assert( index_[cd][ hIndexSet_.index(en) ] >= 0 );
    return index_[cd][ hIndexSet_.index(en) ];
  }
 
  template<int cd>
  int index (const typename GridImp::template Codim<cd>::Entity& en) const 
  {
    // this must not be true for vertices 
    // therefore only check other codims 
    assert( index_[cd][ hIndexSet_.index(en) ] >= 0 );
    return index_[cd][ hIndexSet_.index(en) ];
  }
 
  //! return subIndex (LevelIndex) for a given Entity of codim = 0 and a
  //! given SubEntity codim and number of SubEntity
  template <int cd>
  int subIndex (const typename GridType::template Codim<0>::Entity & en, int i) const
  {
    // this must not be true for vertices 
    // therefore only check other codims 
    assert(index_[cd][ hIndexSet_.template subIndex<cd>(en,i) ] >= 0 );
    return index_[cd][ hIndexSet_.template subIndex<cd>(en,i) ];
  }

  //! returns true if this set provides an index for given entity
  template<class EntityType>
  bool contains (const EntityType& en) const
  {
    enum { cd = EntityType :: codimension };
    return (index_[cd][ hIndexSet_.index(en) ] >= 0 );
  }

  //! return size of IndexSet for a given codim 
  int size ( int codim ) const
  {
    assert( codim >= 0 && codim <= GridType::dimension );
    return size_[codim];
  }

  //! return size of IndexSet for a codim 
  //! this method is to be revised 
  int size ( GeometryType type ) const
  {
    if( typeNotValid(type) ) return 0;
    return size_[GridType::dimension-type.dim()];
  }

  //! do calculation of the index set, has to be called when grid was
  //! changed or if index set is created 
  void calcNewIndex ()
  {
    // resize arrays to new size 
    for(int cd=0; cd<ncodim; ++cd) 
    {
      resizeVectors(index_[cd], hIndexSet_.size(cd));
    }

    // walk grid and store index 
    typedef typename DefaultLeafIteratorTypes<GridImp>:: template Codim<0>::
        template Partition<All_Partition> :: Iterator IteratorType;

    // we start with zero for all codims 
    int num[ncodim];
    for(int cd=0; cd<ncodim; ++cd) num[cd] = 0;

    int elems=0;

    IteratorType endit  = this->template end  < 0, All_Partition > ();
    for(IteratorType it = this->template begin< 0, All_Partition > (); 
        it != endit; ++it)
    {
      ++elems;
      insertEntity(*it,num);
    }
    
    // remember the number of entity on level and cd = 0
    for(int cd=0; cd<ncodim; ++cd) 
    {
      size_[cd] = num[cd];
      //if( size_[cd] != grid_.size(cd) ) 
      //  std::cout << size_[cd] << " calc | grid " << grid_.size(cd)  << "\n";
      //assert( size_[cd] == grid_.size(cd) );
    }
  }

  //! deliver all geometry types used in this grid
  const std::vector<GeometryType>& geomTypes (int codim) const
  {
    return grid_.geomTypes(codim);
  }


  /** @brief Iterator to first entity of given codimension and partition type.
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
  { 
    return this->grid_.template leafbegin<cd,pitype> ();
  }

  /** @brief Iterator to one past the last entity of given codim for partition type
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
  { 
    return this->grid_.template leafend<cd,pitype> ();
  }

private:
  // return whether set has this type stored or not 
  bool typeNotValid (const GeometryType & type) const 
  {
    int codim = GridType :: dimension - type.dim();
    const std::vector<GeometryType> & geomT = geomTypes(codim); 
    for(size_t i=0; i<geomT.size(); ++i) if(geomT[i] == type) return false;
    return true;  
  }
  // calculate index for the codim 
  template <class EntityType>
  void insertEntity(EntityType & en, int (&num)[ncodim])
  {
    InsertEntity<EntityType,dim>::insert(en,hIndexSet_,index_,num);
  }

  // resize vectors of index set 
  void resizeVectors(IndexArrayType &a, int newNumberOfEntries)
  {
    if( newNumberOfEntries > 0 )
    {
      a.resize(newNumberOfEntries);
    }
    for(size_t i=0; i<a.size(); ++i) a[i] = -1;
  }

  // method prints indices of given codim, for debugging
  void print (int codim) const 
  {
    for(size_t i=0; i<index_[codim].size(); i++) 
    {
      std::cout << "levelind[" << i << "] = " << index_[codim][i] << "\n";
    }
  }

  // grid this level set belongs to 
  const GridType & grid_;

  // the grids HierarchicIndexSet
  const HierarchicIndexSetType & hIndexSet_;

  // number of entitys of each level an codim 
  IndexArrayType size_;

  //*********************************************************
  // Methods for mapping the hierarchic Index to index on Level
  IndexArrayType index_[ncodim];
};


//! class for combining 2 index sets together for adaptation process
template <class A, class B >
class CombinedAdaptProlongRestrict
{
  //! space A and B 
  const A & _a;
  const B & _b;
public: 
  //! constructor storing the two references 
  CombinedAdaptProlongRestrict ( const A & a, const B & b ) : _a ( a ) , _b ( b )
  {}
  
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    _a.restrictLocal(father,son,initialize);
    _b.restrictLocal(father,son,initialize);
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    _a.prolongLocal(father,son,initialize);
    _b.prolongLocal(father,son,initialize);
  }
};

} // end namespace Dune 
#endif