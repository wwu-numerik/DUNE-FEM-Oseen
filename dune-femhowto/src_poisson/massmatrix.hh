#ifndef DUNE_MASSMATRIX_HH
#define DUNE_MASSMATRIX_HH

//- Dune includes
#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include "matrixoperator.hh"

namespace Dune
{

  // forward declaration 
  template< class DiscreteFunctionImp >
  class MassOperator;



  template <class DiscreteFunctionImp >
  struct MassOperatorTraits
  {
    // implementation type 
    typedef MassOperator<DiscreteFunctionImp> MatrixOperatorType;

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
    //typedef BlockMatrixObject< DiscreteFunctionSpaceType ,
    //                           DiscreteFunctionSpaceType >  MatrixType;
#endif
  };



  // **** The MassMatrix operator
  template< class DiscreteFunctionImp>
  class MassOperator : public MatrixOperator< MassOperatorTraits<DiscreteFunctionImp> >
  {
  public:
    // type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;
    // type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
    // some more types
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    // polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };
    // The grid's dimension
    enum { dimension = GridType :: dimension };

    // type of quadrature to be used
    //typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef ElementQuadrature< GridPartType, 0 > QuadratureType;

  private:
    // type of base class 
    typedef MatrixOperator< MassOperatorTraits<DiscreteFunctionImp> > BaseType;

    // maximal allowed number of base functions
    enum { maxBaseFunctions = 100 };
    mutable RangeType phi_[ maxBaseFunctions ];

  public:
    // constructor
    MassOperator( const DiscreteFunctionSpaceType &discreteFunctionSpace , bool bnd = true)
    : BaseType( discreteFunctionSpace , bnd )
    {
    }


    //*********** assemble the local matrix ******************
    // In the Finite Element Method, most base functions are zero on a given
    // entity. To build the matrix on one entity, we only need to store a very
    // limited amount of matrix entries, the socalled local matrix.
    template<class  EntityType, class LocalMatrixType >
    void assembleLocalMatrix ( const EntityType &entity,
                              LocalMatrixType &matrix ) const
    {
      // type of geometry 
      typedef typename EntityType :: Geometry GeometryType;
      // get geometry 
      const GeometryType &geometry = entity.geometry();

      // get discrete function space 
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();

      // get base function set 
      const BaseFunctionSetType baseSet
        = dfSpace.baseFunctionSet( entity );

      // get number of base functions 
      const int numBaseFunctions = baseSet.numBaseFunctions();

      assert( numBaseFunctions <= maxBaseFunctions );

      // create quadrature 
      QuadratureType quad( entity, 2 * (polynomialOrder) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
      {
        // calculate integration weight 
        const double weight = quad.weight( pt ) * 
                              geometry.integrationElement( quad.point( pt ) );

        // evaluate all base functions 
        for( int row = 0; row < numBaseFunctions; ++row ) 
        {
          baseSet.evaluate( row , quad[ pt ], phi_[ row ] ); 
        }

        // create matrix entries 
        for( int row = 0; row < numBaseFunctions; ++row ) 
        {
          // calculate diagonal 
          RangeFieldType val = (phi_[row] * phi_[row]) * weight;
          matrix.add( row, row, val);

          // other entries (use symetric structure)
          for ( int col = 0; col < row; ++col ) 
          {
            RangeFieldType val = (phi_[row] * phi_[col]) * weight;
            matrix.add( row, col, val);
            matrix.add( col, row, val);
          }
        }
      }
    }
  };

} // end namespace 
#endif
