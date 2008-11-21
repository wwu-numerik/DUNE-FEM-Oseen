#ifndef DUNE_ELLIPTICFEOP_HH
#define DUNE_ELLIPTICFEOP_HH

//- Dune includes
#include <dune/fem/quadrature/quadrature.hh>

//- local includes
#include "feop.hh"

namespace Dune 
{

  /** \brief Mass-matrix and stiff-matrix of a simple linear elliptic problem: \epsilon a(x) u(x) - \div A(x) \grad u(x) = f(x)
   */
  template< class DiscreteFunctionImp, class TensorImp, class MassImp > // Imp = Implementation 
  class EllipticFEOp
  : public FEOp< DiscreteFunctionImp,
                 SparseRowMatrix< typename DiscreteFunctionImp :: RangeFieldType >,
                 EllipticFEOp< DiscreteFunctionImp, TensorImp, MassImp > >
  {
  public:
    //! type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: FunctionSpaceType
      DiscreteFunctionSpaceType;
        
    //! type of an element of the jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;
    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;
    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };
        
    //! type of tensor
    typedef TensorImp TensorType;
    
    //! type of mass-term 
    typedef MassImp MassType;

    //! type of system matrix
    typedef SparseRowMatrix< RangeFieldType > MatrixType;

    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
            
  private:
    typedef EllipticFEOp< DiscreteFunctionType, TensorType, MassType > ThisType;
    typedef FEOp< DiscreteFunctionType, MatrixType, ThisType > BaseType;

  public:
    //! Operation mode for Finite Element Operator
    typedef typename BaseType :: OpMode OpMode;
 
  private:
    mutable JacobianRangeType grad;
    mutable JacobianRangeType othGrad;

    enum { maxnumOfBaseFct = 100 };
    mutable JacobianRangeType mygrad[ maxnumOfBaseFct ];
    mutable RangeType mymass[ maxnumOfBaseFct ];
        
  private:

    TensorType *stiffTensor_;
    MassType *massTerm_;
     
    const RangeFieldType epsilon_; 
       
  public:
    //! constructor - no Tensor (A(x)), no MassTerm (a(x)) and no \epsilon-term
    EllipticFEOp( const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( NULL ),
      massTerm_( NULL ),
      epsilon_(0)
    {
    }
  
    //! constructor - no Tensor (A(x) and no MassTerm (a(x))
    EllipticFEOp( const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode, const RangeFieldType &epsilon)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( NULL ),
      massTerm_( NULL ),
      epsilon_(epsilon)
    {
    }
    
    //! constructor - Tensor (A(x) but no MassTerm (a(x))
    EllipticFEOp( TensorType &stiff,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode, const RangeFieldType &epsilon)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),  
      massTerm_( NULL ),
      epsilon_(epsilon)
    { 
    }
    
    //! constructor - no Tensor (A(x) but MassTerm (a(x))
    EllipticFEOp( MassType &mass,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode, const RangeFieldType &epsilon)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( NULL ),  
      massTerm_( &mass ),
      epsilon_(epsilon)
    { 
    }
    
    //! constructor - both: Tensor (A(x) and MassTerm (a(x))
    EllipticFEOp( TensorType &stiff,
                 MassType &mass,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode, const RangeFieldType &epsilon)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),  
      massTerm_( &mass ),
      epsilon_(epsilon)
    { 
    }
        
    //! Returns the actual matrix if it is assembled
    const MatrixType* getMatrix () const
    {
      assert( this->matrix_ );
      return this->matrix_;
    }
        
    //! Creates a new empty matrix
    MatrixType* newEmptyMatrix () const 
    {
      return new MatrixType( this->functionSpace_.size(),
                             this->functionSpace_.size(), 
                             15 * dimension );
    }
        
    //! Prepares the local operator before calling apply()
    void prepareGlobal ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest )
    {
      this->arg_  = &arg;
      this->dest_ = &dest;
      this->dest_.clear();
    }
        
    //! return the matrix entr that belong to basis function i and j 
    template< class EntityType >
    double getLocalMatrixEntry( const EntityType &entity,
                                const int i,
                                const int j ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseSet
        = discreteFunctionSpace.baseFunctionSet( entity );

      double val = 0;
      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) { 
        baseSet.jacobian( i, quad, pt, grad );

        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimension, dimension > &inv
          = geometry.jacobianInverseTransposed( quad.point( pt ) );
        const double vol = geometry.integrationElement( quad.point( pt ) );
            
        // multiply with transpose of jacobian inverse 
        grad[ 0 ] = FMatrixHelp :: mult( inv, grad[ 0 ] );
                    
        if( i != j )
	{
          baseSet.jacobian( j, quad, pt, othGrad );
                            
          // multiply with transpose of jacobian inverse 
          
	  othGrad[ 0 ] = FMatrixHelp :: multTransposed( inv, othGrad[ 0 ] );
          val += (grad[ 0 ] * othGrad[ 0 ] ) * quad.weight( pt ) * vol;
        } else
          val += ( grad[ 0 ] * grad[ 0 ] ) * quad.weight( pt ) * vol;
	  
      }
      return val;
    }
 
    template< class  EntityType, class LocalMatrixType > //Corresponding matrix to only one element of the grid
    void getLocalStiffMatrix( const EntityType &entity,
                              const int matrixSize,
                              LocalMatrixType &stiffMatrix ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;
	
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseSet
        = discreteFunctionSpace.baseFunctionSet( entity );
            
      assert( matrixSize <= maxnumOfBaseFct ); 
      for( int i = 0; i < matrixSize; ++i ) 
        for( int j = 0; j <= i; ++j ) 
          stiffMatrix[ j ][ i ] = 0;

      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) {
        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimension, dimension > &inv
                    = geometry.jacobianInverseTransposed( quad.point( pt ) );
        const double volume = geometry.integrationElement( quad.point( pt ) );

        for( int i = 0; i < matrixSize; ++i ) { 
          baseSet.jacobian( i, quad, pt, mygrad[ i ] ); 
      
          // multiply with transpose of jacobian inverse 
          mygrad[ i ][ 0 ] = FMatrixHelp :: mult( inv, mygrad[ i ][ 0 ] );
        }
        
        RangeFieldType weight = quad.weight( pt ) * volume;
        
	// if there is a stiffTensor A(x) (if stiffTensor != NULL), do the following:
	
	if( stiffTensor_ )
        {
          RangeType y;
          stiffTensor_->evaluate( geometry.global( quad.point( pt ) ), y );
          weight *= y[ 0 ];
        }
	
        for( int i = 0; i < matrixSize; ++i ) 
          for ( int j = 0; j <= i; ++j )
	  {
            stiffMatrix[ j ][ i ] += (mygrad[ i ][ 0 ] * mygrad[ j ][ 0 ]) * weight;
	  }
      }
      
      // symmetrize matrix
      for( int i = 0; i < matrixSize; ++i ) 
        for( int j = matrixSize; j > i; --j ) 
          stiffMatrix[ j ][ i ] = stiffMatrix[ i ][ j ];
    }

    
    template< class  EntityType, class LocalMatrixType >
    void getLocalMassMatrix( const EntityType &entity,
                             const int matrixSize,
                             LocalMatrixType &massMatrix ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseSet
        = discreteFunctionSpace.baseFunctionSet( entity );
            
      assert( matrixSize <= maxnumOfBaseFct );
      for( int i = 0; i < matrixSize; ++i )
        for( int j = 0; j <= i; ++j ) 
          massMatrix[ j ][ i ] = 0;

      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) {
     
        const double volume = geometry.integrationElement( quad.point( pt ) );

        for( int i = 0; i < matrixSize; ++i )
	{
	  baseSet.evaluate( i, quad, pt, mymass[ i ] );
        }
        
        RangeFieldType weight = quad.weight( pt ) * volume;
	
	// if there is an epsilon (if epsilon_ != NULL), do the following:
	if ( epsilon_ )
	{
	   // if there is a massTerm_ a(x) (if massTerm_ != NULL), do the following:
           if( massTerm_ )
	   {
             RangeType y;
             massTerm_->evaluate( geometry.global( quad.point( pt ) ), y ); 
             weight *= y[ 0 ];
           }
	   else 
	   {
             weight *= epsilon_;
	   }
	}
	   
        for( int i = 0; i < matrixSize; ++i ) 
          for ( int j = 0; j <= i; ++j )
	  {
	    massMatrix[ j ][ i ] += (mymass[ i ][ 0 ] * mymass[ j ][ 0 ]) * weight;
	  }
      }
      
      // symmetrize matrix
      for( int i = 0; i < matrixSize; ++i ) 
        for( int j = matrixSize; j > i; --j ) 
          massMatrix[ j ][ i ] = massMatrix[ i ][ j ];
    }

    
    template< class  EntityType, class LocalMatrixType >
    void getLocalMatrix( const EntityType &entity,
                                 const int matrixSize,
                                 LocalMatrixType &matrix ) const
    {
      LocalMatrixType stiffMatrix, massMatrix;
    
      // you only need to care about the mass matrix if there is a mass-term (a(x)) or an \epsilon. So if there is, do the following:
      if ( massTerm_ || epsilon_) 
	   {
            getLocalMassMatrix( entity, matrixSize, massMatrix );
           }
      
      getLocalStiffMatrix( entity, matrixSize, stiffMatrix );

      // if there is a mass-term (a(x)) or an \epsilon, do the following:
      if ( massTerm_ || epsilon_) 
	   {
	    for( int i = 0; i < matrixSize; ++i ) 
               for ( int j = 0; j < matrixSize; ++j )
	       {
	          matrix[ j ][ i ] = massMatrix[ j ][ i ] + stiffMatrix[ j ][ i ];
	       } 
           }
       // if there is no mass-term (a(x)) and no \epsilon, do the following:
       else
          {
	   for( int i = 0; i < matrixSize; ++i ) 
               for ( int j = 0; j < matrixSize; ++j )
	       {
	          matrix[ j ][ i ] = stiffMatrix[ j ][ i ];
	       } 
	  }      
    }


  }; // end class


  // Assembler for right rand side
  template< class DiscreteFunctionImp >
  class RightHandSideAssembler
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;
    
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
  
  public:
    // discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    static void assemble( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;
      
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();
  
      discreteFunction.clear(); //discreteFunction auf Null setzen

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        //it* Pointer auf ein Element der Entity
        const GeometryType &geometry = (*it).geometry(); //Referenz auf Geometrie
      
        LocalFunctionType localFunction = discreteFunction.localFunction( *it ); 
        const BaseFunctionSetType &baseFunctionSet //BaseFunctions leben immer auf Refernzelement!!!
          = discreteFunctionSpace.baseFunctionSet( *it ); 

        CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrd ); //0 --> codim 0
        const int numDofs = localFunction.numDofs(); //Dofs = Freiheitsgrade (also die Unbekannten)
        for( int i = 0; i < numDofs; ++i )
        {
          RangeType y, z; //return values
        
          const int numQuadraturePoints = quadrature.nop();
          for( int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
          {
            const double det
              = geometry.integrationElement( quadrature.point( quadraturePoint ) );
            function.evaluate( geometry.global( quadrature.point( quadraturePoint ) ), y );
            baseFunctionSet.evaluate( i, quadrature, quadraturePoint, z ); //i = i'te Basisfunktion;
            localFunction[ i ] += det * quadrature.weight( quadraturePoint ) * (y * z);
          }
        }
      }
    }
  };

} // end namespace 

#endif

