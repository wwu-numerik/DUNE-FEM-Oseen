#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

//- Dune includes
#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/operator/matrix/ontheflymatrix.hh>

#include "matrixoperator.hh"


namespace Dune
{
  // **** forward declaration 
  template< class DiscreteFunctionImp >
  class LaplaceOperator;



  // **** traits class
  template <class DiscreteFunctionImp > 
  struct LaplaceOperatorTraits 
  {
    // implementation type 
    typedef LaplaceOperator<DiscreteFunctionImp> MatrixOperatorType;

    // type of discrete function 
    typedef DiscreteFunctionImp DiscreteFunctionType; 

    // field type 
    typedef typename DiscreteFunctionImp :: RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    // type of system matrix
#if USE_DUNE_ISTL 
    typedef ISTLMatrixObject<DiscreteFunctionSpaceType,
                             DiscreteFunctionSpaceType> MatrixType;
#else 
    typedef SparseRowMatrixObject< DiscreteFunctionSpaceType ,
                                   DiscreteFunctionSpaceType >  MatrixType;
    //typedef OnTheFlyMatrixObject< DiscreteFunctionSpaceType ,
    //                              DiscreteFunctionSpaceType >  MatrixType;
    //typedef BlockMatrixObject< DiscreteFunctionSpaceType ,
    //                           DiscreteFunctionSpaceType >  MatrixType;
#endif
  };



  //**** The Laplace operator
  template< class DiscreteFunctionImp >
  class LaplaceOperator : public MatrixOperator< LaplaceOperatorTraits<DiscreteFunctionImp> >
  {
  public:
    // type of discrete function
    typedef DiscreteFunctionImp DiscreteFunctionType;

    // type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType ;

    // some more typedefs
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType  JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType   BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType  GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType  GridType;

    // polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };
    // The grid's dimension
    enum { dimension = GridType :: dimension };

    // type of quadrature to be used
    //typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef ElementQuadrature< GridPartType, 0 > QuadratureType;

    // type of jacobian inverse type
    typedef FieldMatrix< double, dimension, dimension > JacobianInverseType;

  private:
    // type of this class 
    typedef LaplaceOperator< DiscreteFunctionType > ThisType;
    // type of base class
    typedef MatrixOperator< LaplaceOperatorTraits<DiscreteFunctionImp> > BaseType;

  private:
    mutable JacobianRangeType grad;
    mutable JacobianRangeType othGrad;

    enum { maxBaseFunctions = 100 };
    mutable JacobianRangeType mygrad[ maxBaseFunctions ];


  public:
    // constructor
    LaplaceOperator( const DiscreteFunctionSpaceType &discreteFunctionSpace , bool bnd = true )
    : BaseType( discreteFunctionSpace, bnd )
    {
    }


    //******** assemble the local matrix **********
    // In the Finite Element Method, most base functions are zero on a given
    // entity. To build the matrix on one entity, we only need to store a very
    // limited amount of matrix entries, the socalled local matrix.
    template< class  EntityType, class LocalMatrixType >
    void assembleLocalMatrix ( const EntityType &entity,
                               LocalMatrixType &matrix ) const
    {
      // type of geometry 
      typedef typename EntityType :: Geometry GeometryType;
      const GeometryType &geometry = entity.geometry();

      // get discrete function space 
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();

      // get base function set 
      const BaseFunctionSetType baseSet = dfSpace.baseFunctionSet( entity );

      // get number of base functions 
      const int numBaseFunctions = baseSet.numBaseFunctions();

      assert( numBaseFunctions <= maxBaseFunctions );

      // create quadrature 
      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );

      const int numQuadraturePoints = quad.nop();
      // loop over all quadrature points 
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
      {
        // get Jacobian inverse 
        const JacobianInverseType&inv = 
            geometry.jacobianInverseTransposed( quad.point( pt ) );

        // get integration element 
        const double weight = geometry.integrationElement( quad.point( pt ) )
                              * quad.weight( pt );

        // evaluate gradients of base functions 
        for( int row = 0; row < numBaseFunctions; ++row ) 
        {
          // evaluate gradient of base function on reference element 
          baseSet.jacobian( row , quad[ pt ], grad ); 

          // set grad to zero (umv does only add)
          mygrad[ row ] = 0;
          // multiply with transpose of jacobian inverse 
          // mygrad += inv * grad 
          inv.umv( grad[ 0 ], mygrad[ row ][ 0 ] );
        }

        // add products to matrix  
        for( int row = 0; row < numBaseFunctions; ++row ) 
        {
          // add diagonal 
          {
            RangeFieldType val = (mygrad[ row ][ 0 ] * mygrad[ row ][ 0 ]) * weight;
            matrix.add( row, row , val);
          }

          // add other entries (use symetric structure)
          for ( int col = 0; col < row; ++col ) 
          {
            RangeFieldType val = (mygrad[ row ][ 0 ] * mygrad[ col ][ 0 ]) * weight;
            matrix.add( row, col , val);
            matrix.add( col, row , val);
          }
        }
      }
    }
  };

} // end namespace 
#endif
