#ifndef DUNE_MATRIXOPERATOR_HH
#define DUNE_MATRIXOPERATOR_HH

//- Dune includes
#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/ontheflymatrix.hh>

namespace Dune
{

  //**** Base class for Laplace and Massmatrix operator
  template< class OperatorTraits >
  class MatrixOperator
    : public Operator< typename OperatorTraits :: RangeFieldType,
                       typename OperatorTraits :: RangeFieldType,
                       typename OperatorTraits :: DiscreteFunctionType,
                       typename OperatorTraits :: DiscreteFunctionType > 
  {

    typedef MatrixOperator<OperatorTraits> ThisType;
  public:
    // type of discrete functions
    typedef typename OperatorTraits :: DiscreteFunctionType  DiscreteFunctionType;
    typedef typename OperatorTraits :: MatrixOperatorType  MatrixOperatorType;

    // type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

    // some more typedefs
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType  JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType   RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType  RangeType;
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType  GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType  GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;

    // polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };
    // The grid's dimension
    enum { dimension = GridType :: dimension };

    // type of system matrix
    typedef typename OperatorTraits :: MatrixType MatrixObjectImp;

    typedef typename MatrixObjectImp :: LocalMatrixType LocalMatrixType;
    typedef typename MatrixObjectImp :: MatrixType MatrixType;
    typedef typename MatrixObjectImp :: PreconditionMatrixType PreconditionMatrixType;

    typedef ThisType MatrixObjectType ;


  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

    // pointer to the system matrix
    mutable MatrixObjectImp matrixObj_;

    // flag indicating whether the system matrix has been assembled
    mutable bool matrix_assembled_;

    // true if boundary correction should be done
    const bool boundaryCorrect_;


  protected:
    //**** constructor is protected
    MatrixOperator( const DiscreteFunctionSpaceType &discreteFunctionSpace, bool bnd )
      : discreteFunctionSpace_( discreteFunctionSpace ),
        matrixObj_( discreteFunctionSpace, discreteFunctionSpace , "" ),
        matrix_assembled_( false ),
        boundaryCorrect_(bnd)
    {
    }

  public:
    virtual ~MatrixOperator ()
    {
    }


    //! \brief apply the operator
    virtual void operator() ( const DiscreteFunctionType &u, 
                              DiscreteFunctionType &w ) const 
    {
      // -- on-the-fly version (without explicitly saving the system matrix)
      /*
      {
        typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
        typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType; 

        const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
        w.clear();
  
        const IteratorType end = dfSpace.end();
        for(IteratorType it = dfSpace.begin();
            it != end; ++it )
        {
          const EntityType& entity = *it;
          // get local matrix object 
          LocalMatrixType localMatrix = matrixObj_.localMatrix(entity,entity);
          localMatrix.clear();

          LocalFunctionType uLf = u.localFunction( entity ); 
          LocalFunctionType wLf = w.localFunction( entity ); 

          // assemble local matrix 
          assembleLocalMatrix( entity, localMatrix );

          // multiply 
          localMatrix.multiplyAdd(uLf, wLf);

          if( boundaryCorrect_ && entity.hasBoundaryIntersections() )
          {
            boundaryCorrectOnEntity( entity, uLf, wLf );
          }
        }
      }
      */
      {
        systemMatrix().multOEM(u.leakPointer(), w.leakPointer());
      }
    }


    //! \brief apply the operator
    void multOEM ( const double* arg, 
                   double* dest) const 
    {
      DiscreteFunctionType Arg("ARG", discreteFunctionSpace_ , arg);
      DiscreteFunctionType Dest("DEST", discreteFunctionSpace_ , dest);
      this->operator () (Arg,Dest);
    }


    /*! \brief obtain a reference to the system matrix 
     *
     *  The assembled matrix is returned. If the system matrix has not been
     *  assembled, yet, the assembly is performed.
     *
     *  \returns a reference to the system matrix
     */
    /*
    const MatrixObjectType& systemMatrix () const 
    {
      if( !matrix_assembled_ )
      {
        assemble();
      }
      return const_cast<MatrixObjectType&> (*this);
      //return matrixObj_;
    }
    */
    const MatrixObjectImp& systemMatrix () const
    {
      if( !matrix_assembled_ )
      {
        assemble();
      }

      return matrixObj_;
    }


    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType & preconditionMatrix () const {
      systemMatrix();
      return matrixObj_.pcMatrix();
    }


    //! print the system matrix into a stream
    void print ( std :: ostream out = std :: cout ) const 
    {
      systemMatrix().print( out );
    }


    const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }


    /*! 
     *   assemble: perform grid-walkthrough and assemble global matrix
     * 
     *   If the matrix storage is 
     *   not allocated, new storage is allocated by newEmptyMatrix.
     *   the begin and end iterators are determined and the assembling
     *   of the global matrix initiated by call of assembleOnGrid and 
     *   bndCorrectOnGrid. The assemled flag is set. 
     */
    void assemble () const
    {
      // reserve memory and assemble structure 
      matrixObj_.reserve();
      // clear matrix
      matrixObj_.clear();
      // including boundary correcting 
      assembleOnGrid();
      // matrix is assembled 
      matrix_assembled_ = true;
    }


  protected:
    /*! perform grid walkthrough and assemble matrix
     *
     *  For each element, the local element matrix is determined into the
     *  given local matrix storage and distributed into the global matrix.
     *  Distribution is performed by an add(row,col,val) method on the 
     *  global matrix class.
     */
    void assembleOnGrid () const
    {
      Timer timer; 
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType& entity = *it;
        assembleOnEntity( entity );

        if( boundaryCorrect_ && entity.hasBoundaryIntersections() )
        {
          // get local matrix object 
          LocalMatrixType localMatrix = matrixObj_.localMatrix(entity, entity);

          boundaryCorrectOnEntity( entity , localMatrix );
        }
      }

      std::cout << "Matrix assembly took " << timer.elapsed() << " sec!" << std::endl;
    }



    /*! perform matrix assemble for one entity
     *
     *  \param[in] entity entity for current local update
     */
    void assembleOnEntity ( const EntityType &entity ) const
    {
      // get local matrix object 
      LocalMatrixType localMatrix = matrixObj_.localMatrix(entity,entity);

      // assemble local matrix 
      assembleLocalMatrix( entity, localMatrix );

      // resort for ascending order 
      localMatrix.resort();
    }



    /*! assemble the local matrix
     *
     *  In the Finite Element Method, most base functions are zero on a given
     *  entity. To build the matrix on one entity, we only need to store a very
     *  limited amount of matrix entries, the socalled local matrix.
     */
    template< class LocalMatrixType >
    void assembleLocalMatrix ( const EntityType &entity,
                               LocalMatrixType &matrix ) const
    {
      // Barton-Nackman trick: avoid using virtual functions:
      // Calls the corresponding method in the derived class.
      // This method has to be implemented in the derived class!
      asImp().assembleLocalMatrix(entity, matrix);
    }



    /*! treatment of Dirichlet-DoFs for one entity
     *
     *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
     *
     *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
     *
     *   \param[in]  entity  entity to perform Dirichlet treatment on
     *   \param[out] localMatrix local matrix to correct 
     */
    void boundaryCorrectOnEntity ( const EntityType &entity ,
                                   LocalMatrixType& localMatrix ) const
    {
      typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType  LagrangePointSetType;
      enum { faceCodim = 1 };
      typedef typename GridPartType :: IntersectionIteratorType  IntersectionIteratorType;
      typedef typename LagrangePointSetType :: template Codim<faceCodim> :: SubEntityIteratorType
          FaceDofIteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();

      const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );
 
      const IntersectionIteratorType endit = gridPart.iend( entity );
      for(IntersectionIteratorType it = gridPart.ibegin( entity ); it != endit ; ++it )
      {
        if( it.boundary() )
        {
          const int face = it.numberInSelf();

          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.template endSubEntity< faceCodim >( face );
          for(FaceDofIteratorType faceIt = lagrangePointSet.template beginSubEntity<faceCodim>( face );
              faceIt != faceEndIt; ++faceIt )
          {
            localMatrix.unitRow( *faceIt );
          }
        }
      }
    }


    // This version is only needed by the "on-the-fly" version of operator()
    /*! treatment of Dirichlet-DoFs for one entity
     *
     *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
     *
     *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
     *
     *   \param[in]  entity  entity to perform Dirichlet treatment on
     *   \param[in]  uLf left hand side 
     *   \param[out] wLf right hand side to copy 
     */
    template <class LocalFunctionType>
    void boundaryCorrectOnEntity ( const EntityType &entity ,
                                   const LocalFunctionType& uLf,
                                   LocalFunctionType& wLf) const
    {
      typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType  LagrangePointSetType;
      enum { faceCodim = 1 };
      typedef typename GridPartType :: IntersectionIteratorType  IntersectionIteratorType;
      typedef typename LagrangePointSetType :: template Codim<faceCodim> :: SubEntityIteratorType
          FaceDofIteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();

      const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );
 
      const IntersectionIteratorType endit = gridPart.iend( entity );
      for(IntersectionIteratorType it = gridPart.ibegin( entity ); it != endit ; ++it )
      {
        if( it.boundary() )
        {
          const int face = it.numberInSelf();

          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.template endSubEntity< faceCodim >( face );
          for(FaceDofIteratorType faceIt = lagrangePointSet.template beginSubEntity< faceCodim >( face );
              faceIt != faceEndIt; ++faceIt ) 
          {
            int faceIdx = *faceIt;
            wLf[faceIdx] = uLf[faceIdx];
          }
        }
      }
    }


   private:
     // Barton-Nackman Trick
     MatrixOperatorType& asImp() { return static_cast<MatrixOperatorType&> (*this); }
     const MatrixOperatorType& asImp() const { return static_cast<const MatrixOperatorType&> (*this); }
  };

} // end namespace 
#endif
