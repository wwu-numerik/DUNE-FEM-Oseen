#ifndef FINITEVOLUMESCHEME_HH
#define FINITEVOLUMESCHEME_HH

// include dune reference elements
#include <dune/grid/common/referenceelements.hh>

// include field vector
#include <dune/common/fvector.hh>

// geometry types
#include <dune/common/geometrytype.hh>

//Quadrature type
#include <dune/fem/quadrature/cachequad.hh>

// local includes
#include "transportproblem.hh"
//#include "fvdatahandle.hh"

using namespace Dune;

template <class DiscreteFunctionImp> 
struct FiniteVolumeScheme
{
//----------------------------------------------------------
//
//  Typedefs
//
//----------------------------------------------------------

  typedef DiscreteFunctionImp DiscreteFunctionType;
  
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
 
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  
  typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
  // type of grid
  typedef typename GridPartType :: GridType GridType;

  //problem type
  typedef ProblemImp<GridType> Problem;

  // first we extract the dimensions of the grid  
  enum { dim = GridType :: dimension };
  enum { dimworld = GridType ::dimensionworld };

  // type used for coordinates in the grid
  typedef typename GridType::ctype ctype;

  // index set type 
  typedef typename GridPartType :: IndexSetType IndexSetType; 
  
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
  // iterator type 
  //typedef typename GridPartType :: template Codim<0>:: IteratorType IteratorType;

  // intersection iterator type
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

  // entity pointer type
  typedef typename GridType :: template Codim<0>::EntityPointer EntityPointer;
  typedef typename GridType :: template Codim<0>::Entity   Entity;
  typedef typename GridType :: template Codim<0>::Geometry Geometry;

  // type of domain vector space element
  typedef FieldVector<ctype,dim> DomainType; 
  //typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  
  // type of face domain vector space 
  typedef FieldVector<ctype,dim-1> FaceDomainType; 

  typedef ReferenceElements<ctype,dim-1> FaceReferenceElementContainerType;
  typedef ReferenceElements<ctype,dim> ReferenceElementContainerType;


  //----------------------------------------------------------
  //
  //  Methods
  //
  //----------------------------------------------------------

 

  
  // ************** initialize()
  // initialize the solution discrete function with initial value
  static void initialize ( DiscreteFunctionType &solution)
  {
    calcExactSolution (solution, 0.0);
  }


  // ************** calcExactSolution()
  // calculate the exact solution on a discrete function at time t.
  //template<class DiscreteFunctionType>
  static void  calcExactSolution ( DiscreteFunctionType &solution, double t )
  {
    solution.clear();
    const DiscreteFunctionSpaceType &dfSpace = solution.space();
    
    // grid traversal
    const IteratorType endit = dfSpace.end();
    for (IteratorType it = dfSpace.begin();
         it!=endit; ++it)
    {
      // get entity
      const Entity& entity = *it;

      const GeometryType geomType=entity.type();
      double vol =  ReferenceElementContainerType::general(geomType).volume();

      // cell geometry type
      const Geometry& geo = entity.geometry(); 
     

      /**************************************************************
      // ****** mean value over corners 
      int corners = geo.corners();
      double value = 0.0;
      for(int i=0; i<corners; ++i)
      {
	// evaluate initial data on geometry corners 
        value += Problem::exact( geo[i], t ); 
      }
      value/=corners;
      ***********************************************************/      


      LocalFunctionType localFunction = solution.localFunction( entity );
      const unsigned int numDofs = localFunction.numDofs();
     

      
      const BaseFunctionSetType &baseFunctionSet = localFunction.baseFunctionSet();
      
      QuadratureType quadrature( entity, 2*dfSpace.order() );
//      QuadratureType quadrature( entity, 2 );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
	{
	  RangeType phi;
	  phi = Problem::exact( geo.global( quadrature.point( qp ) ),t );
	  
	  const RangeFieldType weight = quadrature.weight( qp );
	  const double factor = weight/vol;

	  for( unsigned int i = 0; i < numDofs; ++i )
	    {
	      RangeType psi;
	      baseFunctionSet.evaluate( i, quadrature[ qp ], psi );
	      localFunction[ i ] += factor * (phi * psi);
//	      std::cout << "psi = "<<psi<<std::endl;
	      //localFunction[ i ] += (phi * psi)/corners;
	    }
	}
      
      //for( unsigned int i = 0; i < numDofs; ++i )
      //	localFunction[i] = value;
    }
    
    /**** Parallel stuff
    {
      // create data handle 
      FVDataHandle<GridType,IndexSetType,VectorType> dataHandle(indexSet,solution);

      // communicate data 
      gridPart.communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
    }
    ********/

  }


  // ************** printFunctionValues()
  // initialize the solution discrete function with initial value
  static void printFunctionValues (const DiscreteFunctionType &solution )
  {
    const DiscreteFunctionSpaceType &dfSpace = solution.space();

    // grid traversal
    const IteratorType endit = dfSpace.end();
    for (IteratorType it = dfSpace.begin();
         it!=endit; ++it)
    {
      // get entity
      const Entity& entity = *it;

      
      // cell geometry type
      const Geometry& geo = entity.geometry(); 

      // ****** mean value over corners 
      double value = 0.0;
      for(int i=0; i<geo.corners(); ++i)
      {
        // evaluate initial data on geometry corners 
	value += Problem::initial( geo[i] ); 
      }
      LocalFunctionType localFunction = solution.localFunction( entity );
      const unsigned int numDofs = localFunction.numDofs();
//      const BaseFunctionSetType &baseFunctionSet = localFunction.baseFunctionSet();
      
      //QuadratureType quadrature( entity, 2*dfSpace.order() );
      QuadratureType quadrature( entity, 2 );                 
      const unsigned int numQuadraturePoints = quadrature.nop();
      
      for( unsigned int i = 0; i < numDofs; ++i )
	{
	 
	  std::cout << "DofNr. " <<i<<" localFunctionValue "<< localFunction[ i ];
	    //	  	  <<" BaseFunctionValue "<<psi<<std::endl;
	  for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
	    {
	      // RangeType psi;    //psi is 1 if numDof is 1 as well.
	      //baseFunctionSet.evaluate( i, quadrature[ qp ], psi );
	      //std::cout <<" BaseFunctionValue "<<psi<<std::endl;

	      RangeType phi; 
	      const RangeFieldType weight = quadrature.weight( qp );
	      phi = Problem::initial( geo.global( quadrature.point( qp ) ) );
	      std::cout <<"Quadraturpunkt Nr. " <<qp<< ", Weight : "<< weight << " Wert: "<< phi  <<std::endl;
	      //std::cout << "Initialwert "<< phi <<" an Stelle "<<geo.global( quadrature.point( qp ) )<<", Quadraturpunkt Nr. " <<qp<<std::endl;
	      
	      
	    }
	}
      //take mean value 
	 //solution[enIdx] = value/((ctype) geo.corners());
    }
    
    /**** Parallel stuff
    {
      // create data handle 
      FVDataHandle<GridType,IndexSetType,VectorType> dataHandle(indexSet,solution);

      // communicate data 
      gridPart.communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
    }
    ********/

  }

 // ************** finiteVolumeSchemeFem()
  template<class NumericalFluxType>
  static double finiteVolumeScheme (const NumericalFluxType& flux,
				       const DiscreteFunctionType& solution,
				       DiscreteFunctionType& update,
				       const double time,
				       const double cfl)
  {
    
    //set update to zero
    update.clear();

    
    //get discretefunctionspace
    const DiscreteFunctionSpaceType &dfSpace = solution.space();
    
    const GridPartType& gridPart = dfSpace.gridPart();
    
    // get index set 
    const IndexSetType& indexSet = gridPart.indexSet();
    
    // time step size
    double dt = 1e308;

    // evaluate velocity at face center
    DomainType velocity;

    // ******* loop over all leaf elements *********
    // compute update function and optimum dt in one grid traversal
   
    
    // grid traversal
    const IteratorType endit = dfSpace.end();
    for (IteratorType it = dfSpace.begin();
         it!=endit; ++it)
      {
	double dtEstimate = 0.0;                  // estimate for local dt
	const Entity& entity = *it;               // get entity
	const Geometry& geo = entity.geometry();  // cell geometry type
	const double enVolume = geo.volume();     // cell volume
	const double enVolume_1 = 1.0/enVolume;   // 1 over cell volume
	const int enIdx = indexSet.index(entity);
	// **** run through all intersections with neighbors and boundary	
	const IntersectionIteratorType nbend = gridPart.iend(entity); 
	for (IntersectionIteratorType nb = gridPart.ibegin(entity); 
	     nb != nbend; ++nb)
	  {
	    // get geometry type of face
	    const GeometryType faceGeomType = nb.intersectionGlobal().type();
	    // center in face's reference element
	    const FaceDomainType& 
	      xLocal = FaceReferenceElementContainerType::general(faceGeomType).position(0,0);
	    // get normal vector scaled with volume
	    DomainType integrationOuterNormal = nb.integrationOuterNormal(xLocal);
	    // multiply integrationOuterNormal with volume of reference element
	    integrationOuterNormal *= FaceReferenceElementContainerType::general(faceGeomType).volume();
	    // center of face in global coordinates
	    const DomainType xGlobal = nb.intersectionGlobal().global(xLocal);
	    
	    // evaluate velocity at face center
	    Problem::velocity(time,xGlobal,velocity);

	    LocalFunctionType solutionEnLocal = solution.localFunction(entity);  //get local function of the solution on the entity
	    LocalFunctionType updateEnLocal = update.localFunction(entity);   //get local function of the update function on the entity

	    // handle interior face
	    if (nb.neighbor()) 
	      {
		// access neighbor
		EntityPointer outside = nb.outside();
		const Entity& neighbor = *outside;
		const int nbIdx = indexSet.index(neighbor);
		
		// compute flux from one side only
		// this should become easier with the new IntersectionIterator functionality!
		if ( entity.level() >  neighbor.level() || 
		     ( entity.level() == neighbor.level() && enIdx < nbIdx ) ||
		     neighbor.partitionType() != InteriorEntity 
		     )
		  {
		    double gLeft,gRight; 
		    const double nbVolume_1 = 1.0 / neighbor.geometry().volume();
		    LocalFunctionType solutionNbLocal = solution.localFunction(neighbor);   //get local function of the solution function on the neighbor entity 
		    LocalFunctionType updateNbLocal = update.localFunction(neighbor);       //get local function of the update function on the neighbor entity
		    // apply numerical flux 
		    dtEstimate += flux.numericalFlux(integrationOuterNormal,velocity,  
						     solutionEnLocal[0],solutionNbLocal[0],
						     gLeft,gRight); 
		    // calc update of entity
		    updateEnLocal[0] -= gLeft * enVolume_1;
		    // calc update of neighbor
		    updateNbLocal[0] += gRight * nbVolume_1;
		  }
	      }
	    
	    // handle boundary face
	    if (nb.boundary())
	      {
		double bndFlux;
		
		// apply boundary flux                                                    //FEM in all three lines
		dtEstimate += flux.boundaryFlux(time,xGlobal,
						integrationOuterNormal,velocity,
						solutionEnLocal[0],bndFlux); 
		// calc update of entity 
		updateEnLocal[0] -= bndFlux * enVolume_1;
	      }
	    
	    
	  }// ***********end intersectiontraversal *****************
	       
	// compute dt restriction
	dt = std::min(dt, enVolume/dtEstimate);
	
      }// *********** end grid traversal **********
        
    #if 0
    //  communicate data (parallel stuff)
    {
      // create data handle 
      FVDataHandle<GridType,IndexSetType,VectorType> dataHandle(indexSet,update);             //FEM parallel stuff
      // communicate data 
      gridPart.communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
    }
    #endif
    // scale dt with CFL number 
    dt *= cfl;
    return dt;
  } // end method finiteVolumeSchemeFem()


}; // end struct FiniteVolumeScheme
#endif
