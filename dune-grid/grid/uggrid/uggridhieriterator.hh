#ifndef DUNE_UGHIERITERATOR_HH
#define DUNE_UGHIERITERATOR_HH

/** \file
 * \brief The UGGridHierarchicIterator class
 */

#include <dune/common/stack.hh>

namespace Dune {

//**********************************************************************
//
// --UGGridHierarchicIterator
// --HierarchicIterator
  /** \brief Iterator over the descendants of an entity.
   * \ingroup UGGrid
  Mesh entities of codimension 0 ("elements") allow to visit all entities of
  codimension 0 obtained through nested, hierarchic refinement of the entity.
  Iteration over this set of entities is provided by the HierarchicIterator,
  starting from a given entity.
 */

template<class GridImp>
class UGGridHierarchicIterator :
        public Dune::UGGridEntityPointer <0,GridImp>,
        public HierarchicIteratorDefaultImplementation <GridImp,UGGridHierarchicIterator>
{

    friend class UGGridEntity<0,GridImp::dimension,GridImp>;

public:
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! the default Constructor
    UGGridHierarchicIterator(int maxLevel) 
        : maxlevel_(maxLevel) 
    {}

    void increment();

  //! max level to go down 
  int maxlevel_;

    Stack<typename UG_NS<GridImp::dimension>::Element*> elementStack_;

};

    // Include class method definitions
#include "uggridhieriterator.cc"

}  // end namespace Dune

#endif
