#ifndef BOUNDARYTREAT_HH
#define BOUNDARYTREAT_HH

#include <dune/fem/quadrature/cachequad.hh>

namespace Dune {

// Assembler for right hand side
// Note: this is not an L2-Projection
// discreteFunction is an output parameter (kind of return value)
template< int polOrd,
          class FunctionType,
          class DiscreteFunctionType>
static void assembleRHS( const FunctionType &function,
            DiscreteFunctionType &discreteFunction )
{
  // some standard typedefs
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType
    LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: GridType GridType;

  enum { dimension = GridType :: dimension };

  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
  typedef typename EntityType :: Geometry GeometryType;

  // get discrete function space from discrete function
  const DiscreteFunctionSpaceType &discreteFunctionSpace
    = discreteFunction.space();
  discreteFunction.clear();

  // variables for return values
  RangeType val, phi;

  // loop over all elements
  const IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
  {
    // get entity
    const EntityType& entity = *it;
  
    // get geometry
    const GeometryType &geometry = entity.geometry();
  
    // get local function of discrete function
    LocalFunctionType localFunction = discreteFunction.localFunction( entity );
  
    // get base function set from discrete function space
    const BaseFunctionSetType baseFunctionSet =
      discreteFunctionSpace.baseFunctionSet(entity);
  
    // create quadrature for element integration
    CachingQuadrature< GridPartType, 0 > quadrature( entity, polOrd );
  
  
    const int numDofs = localFunction.numDofs();
    const int numQuadraturePoints = quadrature.nop();
    // loop over all quadrature points
    for( int qP = 0; qP < numQuadraturePoints; ++qP )
    {
      // calculate integration weight
      const double intel
        = geometry.integrationElement( quadrature.point( qP ) )
          * quadrature.weight( qP );
  
      // evaluate function
      function.evaluate( geometry.global( quadrature.point( qP ) ), val );
  
      // evaluate value to add to local function
      for( int i = 0; i < numDofs; ++i )
      {
        // evaluate base function i on quadrature point qP
        baseFunctionSet.evaluate( i, quadrature[ qP ], phi );
        // add to local function
        localFunction[ i ] += intel * (val * phi);
      }
    }
  }
}



//! set the dirichlet points to the (given) boundary values
template< class EntityType,
          class FunctionType,
          class DiscreteFunctionType >
static void boundaryTreatment( const EntityType &entity,
                               const FunctionType& boundaryValue,
                               DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType :: FunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType LagrangePointSetType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  enum { faceCodim = 1 };
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef typename LagrangePointSetType :: template Codim<faceCodim>
                                :: SubEntityIteratorType FaceDofIteratorType;


  // get function space 
  const DiscreteFunctionSpaceType &discreteFunctionSpace = rhs.space();
  // get grid part 
  const GridPartType &gridPart = discreteFunctionSpace.gridPart();

  // saves boundary value
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType; 
  RangeType bnd;

  // run over all intersections
  const IntersectionIteratorType endit = gridPart.iend( entity );
  for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != endit; ++it )
  {
    // check if intersection is on the boundary
    if( it.boundary() )
    {
      // get local dofs 
      LocalFunctionType rhsLocal = rhs.localFunction( entity );

      // get lagrange point set 
      const LagrangePointSetType
         &lagrangePointSet = discreteFunctionSpace.lagrangePointSet( entity );

      typedef typename EntityType :: Geometry Geometry;
      const Geometry& geo = entity.geometry();

      // get face number 
      const int face = it.numberInSelf();

      // get dof iterator on face 
      const FaceDofIteratorType
          faceEndIt = lagrangePointSet.template endSubEntity<faceCodim>( face );

      // set all dof entries on this face to the (given) boundary value
      for( FaceDofIteratorType faceIt
              = lagrangePointSet.template beginSubEntity<faceCodim>( face ); faceIt != faceEndIt; ++faceIt )
      {
        const int localDof = *faceIt;

        // evaluate the boundary value
        boundaryValue.evaluate( geo.global( lagrangePointSet.point(localDof)), bnd );

        // set the boundary value
        rhsLocal[ localDof ] = bnd[0];
      }
    }
  }
}


} // end name space Dune

#endif
