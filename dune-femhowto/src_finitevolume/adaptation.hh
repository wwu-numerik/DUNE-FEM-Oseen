#ifndef ADAPTATION_HH 
#define ADAPTATION_HH 

#include <map>

// include capabilities 
#include <dune/grid/common/capabilities.hh>

//#include "datacollector.hh"
//#include "fvdatahandle.hh"

#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/dgspace/dgadaptmanager.hh>  //needed for "typedef RestrictProlongDiscontinuousSpace< DiscreteFunctionType, 0 >  RestrictProlongType;"
using namespace Dune;

//forward declaration
template <class DiscreteFunctionImp, bool ad>
struct Adaptation;


//specialized class for the case adaptation==false
template <class DiscreteFunctionImp> 
struct Adaptation<DiscreteFunctionImp,false>
{
  typedef DiscreteFunctionImp DiscreteFunctionType;
  
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  // do nothing here 
  template<class VectorType>
  static void adapt(GridPartType& gridPart,
                    VectorType& solution,
                    const int minLevel, const int maxLevel)
  {
  } 
  static void adaptFem(GridPartType& gridPart, DiscreteFunctionType& solution,
		       const int minLevel, const int maxLevel)
  {
  } 
  
  static void markEntities(GridPartType& gridPart,const DiscreteFunctionType& solution,
                           const int minLevel, 
                           const int maxLevel,
                           const double refinetol, 
                           const double coarsentol)		  
  {
  }
};


// class doing adaptation of grid and data restriction/prolongation
template <class DiscreteFunctionImp, bool ad> 
class Adaptation
{
  typedef DiscreteFunctionImp DiscreteFunctionType;
  
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  // type of grid 
  typedef typename GridPartType :: GridType GridType;

  // first we extract the dimensions of the grid  
  enum { dim = GridType :: dimension };
  enum { dimworld = GridType ::dimensionworld };

  // type used for coordinates in the grid
  typedef typename GridType::ctype ctype;

  // type of used index set  
  typedef typename GridPartType :: IndexSetType IndexSetType;

  // local id set type 
  typedef typename GridType :: Traits :: LocalIdSet IdSetType; 

  // type of local id 
  typedef typename IdSetType::IdType IdType;

  // type of level iterator 
  typedef typename GridType:: template Codim<0> :: LevelIterator LevelIterator;

  // iterator type 
  typedef typename GridPartType :: template Codim<0>:: IteratorType IteratorType;

  // intersection iterator type
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

  // entity pointer type
  typedef typename GridType :: template Codim<0>::EntityPointer EntityPointer;
  // type of entity 
  typedef typename GridType :: template Codim<0>::Entity   Entity;
  // tpye of geometry 
  typedef typename GridType :: template Codim<0>::Geometry Geometry;

  // type of hierarchic iterator 
  typedef typename Entity :: HierarchicIterator HierarchicIterator;

  // typedef  RestrictProlongPieceWiseConstantData<DiscreteFunctionType> RestrictProlongType;
 // typedef RestrictProlongDiscontinuousSpace< DiscreteFunctionType, 0 >  RestrictProlongType; //Polorder = 0
 // typedef AdaptationManager < GridType, RestrictProlongType > ADOperatorType;

  //-------------------------------------------------------------
  //
  //  constructor
  //
  //-------------------------------------------------------------
  
  // Adaption(GridPartType& gridPart,DiscreteFunctionType& solution,int minLevel, int maxLevel) :
  //gridPart_(gridPart
  //......bad idea!!!!!!!!


  //-------------------------------------------------------------
  //
  //  private methods 
  //
  //-------------------------------------------------------------
  template<class VectorType>
  static double calculateIndicator(const DiscreteFunctionType& solution,  
                                   VectorType& indicator)
  { 

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    
    //get DiscreteFunctionSpace
    const  DiscreteFunctionSpaceType& dfspace = solution.space();

    //get gridPart
    const GridPartType& gridPart = dfspace.gridPart();
    
    const IndexSetType& indexSet = gridPart.indexSet();

    double globalmax = -1E100;
    double globalmin =  1E100;

    // compute update vector and optimum dt in one grid traversal
    const IteratorType endit = dfspace.end();     
    for (IteratorType it     = dfspace.begin(); 
         it!=endit; ++it)
    {
      const Entity& entity = *it;
      
      LocalFunctionType enLocalFunction = solution.localFunction( entity );
      const unsigned int numDofs = enLocalFunction.numDofs();
      // std::cout<<"Dofanzahl :"<<numDofs<<std::endl;
      double enMeanvalue = 0.0;
      for( unsigned int i = 0; i < numDofs; ++i )
	{
	  //RangeType psi;
	  //baseFunctionSet.evaluate( i, quadrature[ qp ], psi );
	  //localFunction[ i ] += weight * (phi * psi);
	  enMeanvalue += enLocalFunction[ i ];
	}
      enMeanvalue /= (double) numDofs;
      
      const int enIdx = indexSet.index( entity ); 
 
      // global min/max
      globalmax = std::max(globalmax,enMeanvalue);
      globalmin = std::min(globalmin,enMeanvalue);

      // iterate over neighbors 
      const IntersectionIteratorType nbend = gridPart.iend( entity ); 
      for (IntersectionIteratorType nb = gridPart.ibegin( entity ); 
           nb != nbend; ++nb)
      {
        if (nb.neighbor())
        {
          // access neighbor
          EntityPointer outside = nb.outside();
          const Entity& neighbor = *outside; 
          const int nbIdx = indexSet.index( neighbor );

          // handle face from one side only
          if ( entity.level() > neighbor.level() ||
              (entity.level() == neighbor.level() && enIdx < nbIdx) )
          {
	    LocalFunctionType nbLocalFunction = solution.localFunction( neighbor );
	    const unsigned int numDofs = nbLocalFunction.numDofs();
	    double nbMeanvalue = 0.0;
	    for( unsigned int i = 0; i < numDofs; ++i )
	      {
		//RangeType psi;
		//baseFunctionSet.evaluate( i, quadrature[ qp ], psi );
		//localFunction[ i ] += weight * (phi * psi);
		nbMeanvalue += nbLocalFunction[ i ];
	      }
	    nbMeanvalue /= (double) numDofs; 
	    
	    double localdelta = std::abs(enMeanvalue - nbMeanvalue);
            indicator[enIdx] = std::max(indicator[enIdx],localdelta);
            indicator[nbIdx] = std::max(indicator[nbIdx],localdelta);
          }
        }
      }
    }
    return globalmax-globalmin;
  }
  /***************************************************************************
  template<class VectorType>
  static double calculateIndicator(const GridPartType& gridPart, 
                                   const VectorType& solution,  
                                   VectorType& indicator)
  {
    const IndexSetType& indexSet = gridPart.indexSet();

    double globalmax = -1E100;
    double globalmin =  1E100;

    // compute update vector and optimum dt in one grid traversal
    const IteratorType endit = gridPart.template end<0>();     
    for (IteratorType it     = gridPart.template begin<0>(); 
         it!=endit; ++it)
    {
      const Entity& entity = *it;
      const int enIdx = indexSet.index( entity ); 
 
      // global min/max
      globalmax = std::max(globalmax,solution[enIdx]);
      globalmin = std::min(globalmin,solution[enIdx]);

      // iterate over neighbors 
      const IntersectionIteratorType nbend = gridPart.iend( entity ); 
      for (IntersectionIteratorType nb = gridPart.ibegin( entity ); 
           nb != nbend; ++nb)
      {
        if (nb.neighbor())
        {
          // access neighbor
          EntityPointer outside = nb.outside();
          const Entity& neighbor = *outside; 
          const int nbIdx = indexSet.index( neighbor );

          // handle face from one side only
          if ( entity.level() > neighbor.level() ||
              (entity.level() == neighbor.level() && enIdx < nbIdx) )
          {
            double localdelta = std::abs(solution[enIdx] - solution[nbIdx]);
            indicator[enIdx] = std::max(indicator[enIdx],localdelta);
            indicator[nbIdx] = std::max(indicator[nbIdx],localdelta);
          }
        }
      }
    }
    return globalmax-globalmin;
  }
  *******************************************************************************/
 public:
  
  static void markEntities(GridPartType& gridPart,const DiscreteFunctionType& solution,
                           const int minLevel, 
                           const int maxLevel,
                           const double refinetol, 
                           const double coarsentol)
  {
   
    
    //get DiscreteFunctionSpace
//    const  DiscreteFunctionSpaceType& dfspace = solution.space();

    //get gridPart
//    const GridPartType& gridPart = dfspace.gridPart();

    const IndexSetType& indexSet = gridPart.indexSet();
    GridType& grid = gridPart.grid();
    
    // create indicator 
//    std::vector<double> indicator( indexSet.size(), -1e100);
    std::vector<double> indicator( indexSet.size(0), -1e100);
    //    std::cout<< "solution size = " << solution.size() << std::endl;
    //    std::cout<< "Grid size = " << grid.size(0) << std::endl;                  // is exactly the same

    // calculate indicator 
    double globalDelta = calculateIndicator(solution,indicator);
   

    // mark cells for refinement/coarsening
    const double localRefineTol = refinetol * globalDelta;
    const double localCoarsenTol = coarsentol * globalDelta;

    
//for test purpose 
    int refCounter = 0;
    int coarsCounter = 0;

    // compute update vector and optimum dt in one grid traversal
    const IteratorType endit = gridPart.template end<0>();     
    for (IteratorType it     = gridPart.template begin<0>(); 
         it!=endit; ++it)
    {
      // entity 
      const Entity& entity = *it;

      // get local error indicator 
      double localIndicator = indicator[indexSet.index(entity)]; 

      if (localIndicator > localRefineTol 
          && (entity.level() < maxLevel) || ! entity.isRegular())
      {
        // mark for refinement 
	++refCounter;
        grid.mark(1,it);
        // mark all neighbors                     warum????????
        const IntersectionIteratorType nbend = gridPart.iend(entity); 
        for (IntersectionIteratorType nb = gridPart.ibegin(entity); 
             nb != nbend; ++nb)
        {
          if (nb.neighbor())
          {
            EntityPointer outside = nb.outside();
            const Entity& neighbor = *outside; 
            if (neighbor.level() < maxLevel || !neighbor.isRegular())
            {
              // mark for refinement 
              grid.mark(1, outside);
            }
          }
        }
      }
      if (localIndicator < localCoarsenTol && entity.level() > minLevel)
      {
        // mark for coarsening 
	++coarsCounter;
        grid.mark(-1, it);
      }
    }
//    std::cout<<std::endl<<"    "<<refCounter<<" Entities marked for refinement"<<std::endl
//	     <<"    "<<coarsCounter<<" Entities marked for coarsening"<<std::endl;
  }

  /***************************************************************************************
  template<class VectorType>
  static void markEntities(GridPartType& gridPart, 
                           VectorType& solution,
                           const int minLevel, 
                           const int maxLevel,
                           const double refinetol, 
                           const double coarsentol)
  {
    // create indicator 
    std::vector<double> indicator( solution.size(), -1e100);

    // calculate indicator 
    double globalDelta = calculateIndicator(gridPart,solution,indicator);

    // mark cells for refinement/coarsening
    const double localRefineTol = refinetol * globalDelta;
    const double localCoarsenTol = coarsentol * globalDelta;

    const IndexSetType& indexSet = gridPart.indexSet();
    GridType& grid = gridPart.grid();

    // compute update vector and optimum dt in one grid traversal
    const IteratorType endit = gridPart.template end<0>();     
    for (IteratorType it     = gridPart.template begin<0>(); 
         it!=endit; ++it)
    {
      // entity 
      const Entity& entity = *it;

      // get local error indicator 
      double localIndicator = indicator[indexSet.index(entity)]; 

      if (localIndicator > localRefineTol 
          && (entity.level() < maxLevel) || ! entity.isRegular())
      {
        // mark for refinement 
        grid.mark(1,it);

        // mark all neighbors                     warum????????
        const IntersectionIteratorType nbend = gridPart.iend(entity); 
        for (IntersectionIteratorType nb = gridPart.ibegin(entity); 
             nb != nbend; ++nb)
        {
          if (nb.neighbor())
          {
            EntityPointer outside = nb.outside();
            const Entity& neighbor = *outside; 
            if (neighbor.level() < maxLevel || !neighbor.isRegular())
            {
              // mark for refinement 
              grid.mark(1, outside);
            }
          }
        }
      }
      if (localIndicator < localCoarsenTol && entity.level() > minLevel)
      {
        // mark for coarsening 
        grid.mark(-1, it);

      }
    }
   
  }


  // do restriction of element and all children 
  template <class VectorType, class DataMapType>  
  static bool hierarchicRestrict(const Entity& entity, 
                                 VectorType& solution,
                                 const IndexSetType& indexSet,
                                 const IdSetType& idSet,
                                 DataMapType& dataMap)
  {
    if(! entity.isLeaf() )
    {
      // true means we are going to restrict data 
      bool doRestrict = true;

      // if the children have children then we have to go deeper 
      const int childLevel = entity.level() + 1;

      // check all children first 
      {
        const HierarchicIterator endit  = entity.hend  ( childLevel );
        for(HierarchicIterator it = entity.hbegin( childLevel ); it != endit; ++it)
        {
          doRestrict &= hierarchicRestrict( *it , solution, indexSet, idSet, dataMap);
        }
      }

      // if doRestrict is still true, restrict data 
      if(doRestrict)
      {
        // true for first child, otherwise false 
        int count = 0;
        double fatherValue = 0.0;

        // restrict for all children 
        const HierarchicIterator endit  = entity.hend  ( childLevel );
        for(HierarchicIterator it = entity.hbegin( childLevel ); it != endit; ++it)
        {
          assert( it->isLeaf() );
          fatherValue += solution[indexSet.index(*it)]; 
          ++count;
        }

        // take mean value 
        fatherValue /= (double) count;
        // store in data map 
        dataMap[idSet.id(entity)] = fatherValue;
      }
    }
    else 
    {
      // save leaf data 
      dataMap[idSet.id(entity)] = solution[indexSet.index(entity)]; 
    } 

    // returns true if entity meight be removed during next adaptation
    return entity.mightBeCoarsened();
  }


  template <class VectorType, class DataMapType>  
  static void hierarchicProlong(const Entity& entity, 
                                VectorType& solution,
                                const IndexSetType& indexSet,
                                const IdSetType& idSet,
                                DataMapType& dataMap,
                                const int maxLevel)
  {
    // check entity  
    if( entity.isLeaf() && ! entity.wasRefined() )
    {
      solution[indexSet.index(entity)] = dataMap[idSet.id(entity)]; 
    }

    // hierarchically prolong data 
    const HierarchicIterator endit = entity.hend ( maxLevel );
    for(HierarchicIterator it = entity.hbegin( maxLevel );
        it != endit; ++it)
    {
      // get entity 
      const Entity& son = *it;

      // true if entity was created during last adaptation 
      const bool wasRefined = son.wasRefined();

      if( son.isLeaf () ) 
      {
        // prolong data from father
        if( wasRefined )
        {
          // get father 
          EntityPointer father = son.father(); 
          // copy entries 
          solution[indexSet.index(son)] = dataMap[idSet.id(*father)]; 
        }
        else 
        {
          // if element has not changed, just copy data 
          solution[indexSet.index(son)] = dataMap[idSet.id(son)]; 
        }
      }
      else 
      {
        if( wasRefined )
        {
          // get father 
          EntityPointer father = son.father(); 
          // copy entries 
          dataMap[idSet.id(son)] = dataMap[idSet.id(*father)]; 
        }
      }
    }
  } 


  //-------------------------------------------------------------
  //
  //  public adapt() method
  //
  //-------------------------------------------------------------
public:
  //adaptFem()-Method
  static void adaptFem(GridPartType &gridPart,DiscreteFunctionType &solution, 
                    const int minLevel, const int maxLevel)
  {
    // tolerance value for refinement strategy
    const double refineTol  = 0.1;
    // tolerance value for coarsening strategy
    const double coarsenTol = 0.1 * refineTol;

    //get DisceteFunctionSpace
//    DiscreteFunctionSpaceType& dfspace = solution.space();

    //get gridPart
//    GridPartType& gridPart = dfspace.gridPart();

    std::cout<<"  Marking entities...";

    // mark entities for coarsening and refinement 
    markEntitiesFem(gridPart,solution,minLevel,maxLevel,refineTol,coarsenTol);
    std::cout<<"  done"<<std::endl;
    // get grid 
    GridType&grid = gridPart.grid();
    
    // type of data map 
//    typedef std::map<IdType,double> DataMapType;

    // map to store all data during adaptation 
//    DataMapType dataMap;

    RestrictProlongType rp(solution);
    ADOperatorType adop(grid,rp);

    #if 1
    adop.adapt();
    #else
    // check if elements might be removed in next adaptation cycle 
    bool mightCoarsen = grid.preAdapt();
    // if elements might be removed 
    if(mightCoarsen)
    {
      // restrict data and save leaf level 
      const LevelIterator endit = grid.template lend<0>(0);
      for (LevelIterator it = grid.template lbegin<0>(0);
           it != endit; ++it)
      {
        hierarchicRestrict(*it, solution, gridPart.indexSet(), 
                           grid.localIdSet(), dataMap);
      }
    }
    else 
    {
	    // only copy leaf level 			//wozu????
      const IteratorType endit = gridPart.template end<0>();
      for (IteratorType it = gridPart.template begin<0>();
           it != endit; ++it)
      {
        hierarchicRestrict(*it, solution, gridPart.indexSet(), 
                           grid.localIdSet(), dataMap);
      }
    }
    
    //assert( dataMap.size() >= solution.size() );

    // adapt grid, returns if new elements were created 
    bool refined = grid.adapt();
    // re-balance grid (parallel stuff)
    {
      #if 1
      // create data collector 
      DataCollectorImp<GridType,DataMapType> dataCollector(grid,dataMap);
      #else
      typedef Dune::DataCollector<GridType,DataMapType> DataCollectorType;
      DataCollectorType dataCollector(grid,dataMap);
      #endif
      // call load balance 
      grid.loadBalance( dataCollector );
    }


    // resize solution vector if elements might have been removed 
    // or were created 
    if( refined || mightCoarsen ) 
    {
      solution.resize( gridPart.indexSet().size(0) );
    }

    // interpolate all new cells 
    {
      const LevelIterator endit = grid.template lend<0>(0); 
      for (LevelIterator it = grid.template lbegin<0>(0);
           it != endit; ++it)
      {
        hierarchicProlong(*it, solution, gridPart.indexSet(), 
                          grid.localIdSet(), dataMap, grid.maxLevel());
      }
    }

    // cleanup adaptation markers 
    grid.postAdapt();
    #endif
  }

  //adapt() method
  template<class VectorType>
  static void adapt(GridPartType& gridPart, 
                    VectorType& solution, 
                    const int minLevel, const int maxLevel)
  {
    // tolerance value for refinement strategy
    const double refineTol  = 0.1;
    // tolerance value for coarsening strategy
    const double coarsenTol = 0.1 * refineTol;

    // mark entities for coarsening and refinement 
    markEntities(gridPart,solution,minLevel,maxLevel,refineTol,coarsenTol);

    // get grid 
    GridType& grid = gridPart.grid();

    // type of data map 
    typedef std::map<IdType,double> DataMapType;

    // map to store all data during adaptation 
    DataMapType dataMap;

    // check if elements might be removed in next adaptation cycle 
    bool mightCoarsen = grid.preAdapt();

    // if elements might be removed 
    if(mightCoarsen)
    {
      // restrict data and save leaf level 
      const LevelIterator endit = grid.template lend<0>(0);
      for (LevelIterator it = grid.template lbegin<0>(0);
           it != endit; ++it)
      {
        hierarchicRestrict(*it, solution, gridPart.indexSet(), 
                           grid.localIdSet(), dataMap);
      }
    }
    else 
    {
	    // only copy leaf level 			//wozu????
      const IteratorType endit = gridPart.template end<0>();
      for (IteratorType it = gridPart.template begin<0>();
           it != endit; ++it)
      {
        hierarchicRestrict(*it, solution, gridPart.indexSet(), 
                           grid.localIdSet(), dataMap);
      }
    }

    //assert( dataMap.size() >= solution.size() );

    // adapt grid, returns if new elements were created 
    bool refined = grid.adapt();

    // re-balance grid (parallel stuff)
    {
      #if 1
      // create data collector 
      DataCollectorImp<GridType,DataMapType> dataCollector(grid,dataMap);
      #else
      typedef Dune::DataCollector<GridType,DataMapType> DataCollectorType;
      DataCollectorType dataCollector(grid,dataMap);
      #endif
      // call load balance 
      grid.loadBalance( dataCollector );
    }


    // resize solution vector if elements might have been removed 
    // or were created 
    if( refined || mightCoarsen ) 
    {
      solution.resize( gridPart.indexSet().size(0) );
    }

    // interpolate all new cells 
    {
      const LevelIterator endit = grid.template lend<0>(0); 
      for (LevelIterator it = grid.template lbegin<0>(0);
           it != endit; ++it)
      {
        hierarchicProlong(*it, solution, gridPart.indexSet(), 
                          grid.localIdSet(), dataMap, grid.maxLevel());
      }
    }

    // cleanup adaptation markers 
    grid.postAdapt();
  }
  ******************************************************************************************/
};

#endif 
