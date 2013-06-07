#ifndef DUNE_OSEEN_SPMATRIX_HH
#define DUNE_OSEEN_SPMATRIX_HH

//- system includes 
#include <vector>
#include <set>
#include <algorithm>

//- local includes 
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/operator/common/localmatrix.hh> 
#include <dune/fem/operator/common/localmatrixwrapper.hh> 
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/misc/functor.hh>

#ifdef ENABLE_UMFPACK 
#include <umfpack.h>
#endif

namespace Dune
{
  
    //! PortedSparseRowMatrixObject
  template <class RowFunctionImp, class ColFunctionImp, class TraitsImp>
  class PortedSparseRowMatrixObject : public OEMMatrix
  {
      typedef typename RowFunctionImp::DiscreteFunctionSpaceType
          DomainSpace;
      typedef typename ColFunctionImp::DiscreteFunctionSpaceType
          RangeSpace;
  public:
      typedef RowFunctionImp RowDiscreteFunctionType;
      typedef ColFunctionImp ColumnDiscreteFunctionType;
    //! type of traits 
    typedef TraitsImp Traits;

    //! type of stencil class 
    typedef typename Traits :: StencilType StencilType;
    
    typedef DomainSpace DomainSpaceType;
    typedef RangeSpace RangeSpaceType;

  private:  
    typedef PortedSparseRowMatrixObject< RowFunctionImp, ColFunctionImp, Traits > ThisType;

  protected:
    typedef typename DomainSpaceType :: GridType GridType;

    typedef typename DomainSpace :: EntityType  ColumnEntityType ;
    typedef typename RangeSpace :: EntityType   RowEntityType ;

    template< class MatrixObject >
    struct LocalMatrixTraits;

    template< class MatrixObject >
    class LocalMatrix;
    
  public:  
    typedef SparseRowMatrix< double > MatrixType;
    typedef MatrixType PreconditionMatrixType;

  public:
    //! type of local matrix 
    typedef LocalMatrix<ThisType> ObjectType;
    typedef ThisType LocalMatrixFactoryType;
    typedef Fem :: ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
    //! type of local matrix 
    typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  protected:
    const DomainSpaceType &domainSpace_;
    const RangeSpaceType &rangeSpace_;
    
    int sequence_;

    mutable MatrixType matrix_;
    bool preconditioning_;

    mutable LocalMatrixStackType localMatrixStack_;

  public:
    //! setup matrix handler 
    inline PortedSparseRowMatrixObject( const DomainSpaceType &domainSpace,
                                  const RangeSpaceType &rangeSpace,
                                  const std::string &paramfile = "" )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      sequence_( -1 ),
      matrix_(),
      preconditioning_( false ),
      localMatrixStack_( *this )
    {
      int precon = 0;
      if( paramfile != "" )
      {
        readParameter( paramfile, "Preconditioning", precon );
      }
      else 
      {
        precon = Parameter :: getValue("Preconditioning", precon );
      }
      preconditioning_ = (precon > 0) ? true : false;
    }

    //! return reference to stability matrix
    inline MatrixType &matrix () const
    {
      return matrix_;
    }

    //! interface method from LocalMatrixFactory
    inline ObjectType *newObject () const
    {
      return new ObjectType( *this, domainSpace_, rangeSpace_ );
    }

    //! return local matrix 
    inline LocalMatrixType localMatrix( const RowEntityType &rowEntity,
                                        const ColumnEntityType &colEntity ) const
    {
      return LocalMatrixType( localMatrixStack_, rowEntity, colEntity );
    }

    //! resize all matrices and clear them 
    inline void clear ()
    {
      matrix_.clear();
    }

    //! return true if precoditioning matrix is provided 
    bool hasPreconditionMatrix () const
    {
      return preconditioning_;
    }

    //! return reference to preconditioner 
    const PreconditionMatrixType &preconditionMatrix () const 
    { 
      return matrix_;
    }

    //! reserve memory corresponnding to size of spaces
    inline void reserve(bool /*verbose*/ = false )
    {
      if( sequence_ != domainSpace_.sequence() )
      {
#ifndef DUNE_FEM_DONT_CHECKITERATORS_OF_SPACE
        // if empty grid do nothing (can appear in parallel runs)
        if( (domainSpace_.begin() != domainSpace_.end())
            && (rangeSpace_.begin() != rangeSpace_.end()) )
#endif
        {        
          // upper estimate for number of non-zeros 
          const int nonZeros = std::max( StencilType :: nonZerosEstimate( rangeSpace_ ), matrix_.numNonZeros() );
          matrix_.reserve( domainSpace_.size(), rangeSpace_.size(), nonZeros, 0.0 );
        }
        sequence_ = domainSpace_.sequence();
      }
    }

    //! solve A dest = arg using the UMFPACK direct solver 
    template< class DomainFunction, class RangeFunction >
    void solveUMF( const DomainFunction &arg, RangeFunction &dest ) const
    {
      DUNE_THROW(NotImplemented,"solveUMF only implemented for AdaptiveDiscreteFunctions");
    }

    //! solve A dest = arg using the UMFPACK direct solver 
    void solveUMF ( const AdaptiveDiscreteFunction< DomainSpaceType > &arg, 
                    AdaptiveDiscreteFunction< RangeSpaceType> &dest ) const
    {
      matrix_.solveUMF( arg.leakPointer(), dest.leakPointer() );
    }

    //! apply matrix to discrete function
    template< class DomainFunction, class RangeFunction >
    void apply ( const DomainFunction &arg, RangeFunction &dest ) const
    {
      // do matrix vector multiplication 
      matrix_.apply( arg, dest );

      // communicate data 
      dest.communicate();
    }

    //! apply matrix to discrete function
    void apply ( const AdaptiveDiscreteFunction< DomainSpaceType > &arg, 
                 AdaptiveDiscreteFunction< RangeSpaceType> &dest ) const
    {
      // do matrix vector multiplication 
      matrix_.multOEM( arg.leakPointer(), dest.leakPointer() );

      // communicate data 
      dest.communicate();
    }

    //! apply transposed matrix to discrete function
    template< class DomainFunction, class RangeFunction >
    void apply_t ( const RangeFunction &arg, DomainFunction &dest ) const
    {
      // do matrix vector multiplication 
      matrix_.apply_t( arg, dest );

      // communicate data 
      dest.communicate();
    }

    //! apply transposed matrix to discrete function
    void apply_t ( const AdaptiveDiscreteFunction< RangeSpaceType > &arg, 
                   AdaptiveDiscreteFunction< DomainSpaceType> &dest ) const
    {
      // do matrix vector multiplication 
      matrix_.multOEM_t( arg.leakPointer(), dest.leakPointer() );

      // communicate data 
      dest.communicate();
    }


    //! mult method of matrix object used by oem solver
    double ddotOEM( const double *v, const double *w ) const
    {
      typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
      DomainFunctionType V( "ddot V", domainSpace_, v );
      DomainFunctionType W( "ddot W", domainSpace_, w );
      return V.scalarProductDofs( W );
    }

    //! mult method of matrix object used by oem solver
    void multOEM( const double *arg, double *dest ) const
    {
      typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
      typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunctionType;

      DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
      RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );
      apply( farg, fdest );
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort() 
    {
      matrix_.resort();
    }

    void createPreconditionMatrix()
    { 
    }

    /** \brief delete all row belonging to a hanging node and rebuild them */
    template <class HangingNodesType> 
    void changeHangingNodes(const HangingNodesType& hangingNodes) 
    {
      {
        typedef typename HangingNodesType :: IteratorType IteratorType;
        const IteratorType end = hangingNodes.end();
        for( IteratorType it = hangingNodes.begin(); it != end; ++it)
        {
          insertHangingRow( hangingNodes, (*it).first , (*it).second );
        }
      }

      /*
      {
        typedef typename HangingNodesType :: SlaveNodesType SlaveNodesType; 
        typedef typename SlaveNodesType :: const_iterator iterator;
        iterator end =  hangingNodes.slaveNodes().end();
        for( iterator it =  hangingNodes.slaveNodes().begin(); 
             it != end; ++it )
        {
          matrix().unitRow( *it );
          matrix().set( *it, *it, 0.0 );
        }
      }
      */
    }
protected:
    /** \brief insert row to be a row for a hanging node */
    template <class HangingNodesType, class ColumnVectorType>
    void insertHangingRow( const HangingNodesType& hangingNodes,
                           const int row, const ColumnVectorType& colVec)
    {
      const size_t cols = colVec.size();

      // distribute row to associated rows 
      const int nonZeros = matrix().numNonZeros( row );
      for( int c = 0; c < nonZeros; ++c)
      {
        std::pair< double, int > val =  matrix().realValue( row, c );
        for( size_t j = 0; j < cols; ++ j)
        {
          const double value = colVec[j].second * val.first;
          assert( ! hangingNodes.isHangingNode( colVec[j].first ) );
          matrix().add( colVec[j].first , val.second, value );
        }
      }

      // replace hanging row 
      matrix().unitRow( row );
      for( size_t j = 0; j < cols; ++ j)
      {
        matrix().add( row, colVec[j].first , -colVec[j].second );
      }
    }
  };



  template< class DomainFunction, class RangeFunction, class TraitsImp >
  template< class MatrixObject >
  struct PortedSparseRowMatrixObject< DomainFunction, RangeFunction, TraitsImp >::LocalMatrixTraits
  {
      typedef typename DomainFunction::DiscreteFunctionSpaceType
        DomainSpace;
    typedef DomainSpace DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType
        RangeSpace;
    typedef RangeSpace RangeSpaceType;
    
    typedef typename PortedSparseRowMatrixObject< DomainFunction, RangeFunction, TraitsImp >
      :: template LocalMatrix< MatrixObject > LocalMatrixType;

    typedef typename RangeSpaceType :: RangeFieldType RangeFieldType;
    typedef RangeFieldType LittleBlockType;
  };



  //! LocalMatrix 
  template< class DomainFunction, class RangeFunction, class TraitsImp >
  template< class MatrixObject >
  class PortedSparseRowMatrixObject< DomainFunction, RangeFunction, TraitsImp > :: LocalMatrix
  : public LocalMatrixDefault< LocalMatrixTraits< MatrixObject > >
  {
      typedef typename DomainFunction::DiscreteFunctionSpaceType
        DomainSpace;
      typedef typename RangeFunction::DiscreteFunctionSpaceType
        RangeSpace;
  public:
    //! type of matrix object 
    typedef MatrixObject MatrixObjectType;

    //! type of the traits
    typedef LocalMatrixTraits< MatrixObjectType > Traits;

  private:
    typedef LocalMatrixDefault< Traits > BaseType;

  public:
    //! type of matrix 
    typedef typename MatrixObjectType :: MatrixType MatrixType;

    //! type of entries of little blocks 
    typedef typename Traits :: RangeFieldType RangeFieldType;

    //! type of the DoFs
    typedef RangeFieldType DofType;

    //! type of little blocks 
    typedef typename Traits :: LittleBlockType LittleBlockType;

  protected:
    MatrixType &matrix_; 
    
    //! global row numbers 
    std :: vector< int > row_;
    //! global col numbers  
    std :: vector< int > col_;

    using BaseType :: domainSpace_;
    using BaseType :: rangeSpace_;
    
  public:  
    //! constructor taking entity and spaces for using mapToGlobal
    //! class RowSpaceType, class ColSpaceType> 
    inline LocalMatrix( const MatrixObjectType &matrixObject,
                        const DomainSpaceType &domainSpace,
                        const RangeSpaceType &rangeSpace )
    : BaseType( domainSpace, rangeSpace),
      matrix_( matrixObject.matrix() )
    {
    }
    
  private: 
    // prohibit copying 
    LocalMatrix( const LocalMatrix & );

  public:
    void init( const RowEntityType &rowEntity, const ColumnEntityType &colEntity )
    {
      // initialize base functions sets 
      BaseType :: init ( rowEntity , colEntity );
        
      row_.resize( domainSpace_.baseFunctionSet( rowEntity ).numBaseFunctions() );
      col_.resize( rangeSpace_.baseFunctionSet( colEntity ).numBaseFunctions() );

      // Martin: shouldn't domainSpace and rangeSpace be flipped, here?
      typedef typename DomainSpaceType::MapperType::DofMapIteratorType DomainMapIterator;
      const DomainMapIterator dmend = domainSpace_.mapper().end( rowEntity );
      for( DomainMapIterator dmit = domainSpace_.mapper().begin( rowEntity ); dmit != dmend; ++dmit )
      {
        assert( dmit.global() == domainSpace_.mapToGlobal( rowEntity, dmit.local() ) );
        row_[ dmit.local() ] = dmit.global();
      }

      typedef typename RangeSpaceType::MapperType::DofMapIteratorType RangeMapIterator;
      const RangeMapIterator rmend = rangeSpace_.mapper().end( colEntity );
      for( RangeMapIterator rmit = rangeSpace_.mapper().begin( colEntity ); rmit != rmend; ++rmit )
      {
        assert( rmit.global() == rangeSpace_.mapToGlobal( colEntity, rmit.local() ) );
        col_[ rmit.local() ] = rmit.global();
      }

#if 0
      const size_t rows = row_.size();
      for( size_t i = 0; i < rows; ++i )
        row_[ i ] = domainSpace_.mapToGlobal( rowEntity, i );
      
      const size_t cols = col_.size();
      for( size_t i = 0; i < cols; ++i )
        col_[ i ] = rangeSpace_.mapToGlobal( colEntity, i );
#endif
      // rows are determined by the range space
      row_.resize( domainSpace_.mapper().numDofs( colEntity ) );
      domainSpace_.mapper().mapEach( colEntity, Fem::AssignFunctor< std::vector< int > >( row_ ) );

      // columns are determind by the domain space
      col_.resize( rangeSpace_.mapper().numDofs( rowEntity ) );
      rangeSpace_.mapper().mapEach( rowEntity, Fem::AssignFunctor< std::vector< int > >( col_ ) );
    }

    //! return number of rows 
    int rows () const
    {
      return row_.size();
    }

    //! return number of columns 
    int columns () const
    {
      return col_.size();
    }

    //! add value to matrix entry
    void add( int localRow, int localCol, const DofType value )
    {
	  assert( value == value );
      assert( (localRow >= 0) );
	  assert( (localRow < rows()) );
      assert( (localCol >= 0) );
      assert( (localCol < columns()) );

      matrix_.add( row_[ localRow ], col_[ localCol ], value );
    }

    //! get matrix entry 
    DofType get( int localRow, int localCol ) const
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );

      return matrix_( row_[ localRow ], col_[ localCol ] );
    }

    //! set matrix entry to value 
    void set( int localRow, int localCol, const DofType value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );

      matrix_.set( row_[ localRow ], col_[ localCol ], value );
    }

    //! set matrix row to zero except diagonla entry 
    void unitRow( const int localRow )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      matrix_.unitRow( row_[ localRow ] );
    }

    //! set matrix row to zero
    void clearRow( const int localRow )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      matrix_.clearRow( row_[localRow]);
    }

    //! set matrix column to zero
    void clearCol ( const int localCol )
    {
      assert( (localCol >= 0) && (localCol < columns()) );
      matrix_.clearCol( col_[localCol] );
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      const int row = rows();
      for( int i = 0; i < row; ++i )
        matrix_.clearRow( row_[ i ] );
    }

    //! scale local matrix with a certain value 
    void scale ( const DofType& value ) 
    {
      const int row = rows();
      for( int i = 0; i < row; ++i )
        matrix_.scaleRow( row_[ i ] , value );
    }

    //! resort all global rows of matrix to have ascending numbering 
    void resort ()
    {
      const int row = rows();
      for( int i = 0; i < row; ++i )
        matrix_.resortRow( row_[ i ] );
    }
  };



  // SparseRowMatrixOperator
  // -----------------------

  template< class DomainFunction, class RangeFunction, class TraitsImp >
  class PortedSparseRowMatrixOperator
  : public PortedSparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType, TraitsImp >,
    public Operator< typename DomainFunction::RangeFieldType, typename RangeFunction::RangeFieldType, DomainFunction, RangeFunction >
  {
    typedef PortedSparseRowMatrixOperator< DomainFunction, RangeFunction, TraitsImp > This;
    typedef PortedSparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType, TraitsImp > Base;

  public:
    typedef typename Base::DomainSpaceType DomainSpaceType;
    typedef typename Base::RangeSpaceType RangeSpaceType;

    using Base::apply;

    PortedSparseRowMatrixOperator ( const std::string &name,
                              const DomainSpaceType &domainSpace,
                              const RangeSpaceType &rangeSpace,
                              const std::string &paramfile = "" )
    : Base( domainSpace, rangeSpace, paramfile )
    {}

    virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const
    {
      Base::apply( arg, dest );
    }

    const Base &systemMatrix () const
    {
      return *this;
    }
  };

} // end namespace Dune 


#endif  // DUNE_OSEEN_SPMATRIX_HH
