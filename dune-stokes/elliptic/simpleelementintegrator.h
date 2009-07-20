#ifndef SIMPLEELEMENTINTEGRATOR_H
#define SIMPLEELEMENTINTEGRATOR_H

#include <dune/fem/operator/elementintegrators.hh>

  /** \class SimpleElementMatrixIntegrator
   *
   *  \brief The SimpleElementMatrixIntegrator class implements some
   *         default functionality.
   *
   *  The class can be used for deriving own ElementMatrixIntegrator classes.
   *  Mainly it provides model access and storage functionality. It does not
   *  implement the addElementMatrix method, so this still has to be implemented
   *  in derived classes. But the class provides building blocks for use in
   *  future addElementMatrix methods. In particular methods
   *  addMassElementMatrix, addDiffusiveFluxElementMatrix,
   *  addConvectiveFluxElementMatrix.
   *
   *  So a derivation of an own class is very simple, e.g. done in
   *  dune/fem/examples/elliptic/elliptic.cc
   */
  template <class TraitsImp, class ModelImp >
  class SimpleElementMatrixIntegrator
  : public ElementMatrixIntegratorInterface< TraitsImp,
                                             typename ModelImp :: LinearEllipticModelInterfaceType,
                                             SimpleElementMatrixIntegrator<TraitsImp,ModelImp> >
  {
  public:
    typedef SimpleElementMatrixIntegrator
      < TraitsImp, ModelImp >
      ThisType;
    typedef ElementMatrixIntegratorInterface
      < TraitsImp,
        typename ModelImp :: LinearEllipticModelInterfaceType,
        SimpleElementMatrixIntegrator<TraitsImp,ModelImp> >
      BaseType;

  public:
    typedef typename BaseType :: TraitsType TraitsType;
    typedef typename BaseType :: ModelType ModelType;
    typedef typename TraitsType::ElementMatrixType    ElementMatrixType;
    typedef typename TraitsType::EntityType           EntityType;
    typedef typename TraitsType::EntityPointerType    EntityPointerType;

    //! derived typedefs for ommitting long type-dereferencing in the methods
    typedef typename TraitsType :: ElementQuadratureType ElementQuadratureType;
    typedef typename TraitsType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename TraitsType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename TraitsType::RangeType RangeType;
    typedef typename TraitsType::DomainType DomainType;
    typedef typename TraitsType::JacobianRangeType JacobianRangeType;
    typedef typename TraitsType::IntersectionIteratorType
                     IntersectionIteratorType;
    typedef typename TraitsType::IntersectionQuadratureType
                     IntersectionQuadratureType;
    typedef typename TraitsType::GridPartType GridPartType;
    //! used grid type
    typedef typename GridPartType :: GridType GridType;

    //! dimension of world
    enum { dimworld = GridType :: dimensionworld };

    /**  constructor: initialize ElementMatrixIntegrator with given model
     *
     *   Default implementation storing a reference to the
     *   model with the data functions. The constructor addditionally
     *   allocates storage for temporary basis-function gradients, which are
     *   required in case of diffusive contributions.
     *
     *   perhaps extension by a further parameter (global matrix) might be useful,
     *   if a local-matrix interface is used, which directly accesses the
     *   global matrix?
     *
     *   \param[in]  model    model providing the (continuous) data
     *   \param[in]  dfSpace  discrete function space
     *   \param[in]  verbose  optional verbosity flag
     */
    SimpleElementMatrixIntegrator ( ModelType& model,
                                     const DiscreteFunctionSpaceType &dfSpace,
                                     int verbose = 0 )
    : model_( model ),
      discreteFunctionSpace_( dfSpace ),
      verbose_( verbose )
    {
      // determine number of basis functions on first entity
      if (verbose_)
          std::cout << "entered constructor of SimpleElementMatrixIntegrator" << std :: endl;

      if (verbose_)
          std::cout << "got discrete functionspace" << std :: endl;

      EntityPointerType ep = dfSpace.begin();
      EntityType& entity = *ep;

      if (verbose_)
          std::cout << "got first entity" << std :: endl;

//            const typename EntityType::Geometry& geo = entity.geometry();

      if (verbose_)
          std::cout << "successful got geometry" << std :: endl;

      const BaseFunctionSetType baseFunctionSet
        = dfSpace.baseFunctionSet( entity );
      numBaseFunctions_ = baseFunctionSet.numBaseFunctions();
      gradPhiPtr_ = new JacobianRangeType[ numBaseFunctions_ ];

      if (verbose_)
          std::cout << "allocated temporary storage for gradients\n";
    }

    /*!  access function for model
     *
     *   Default implementation is return of the stored reference.
     *
     *   \return reference to the locally stored model reference
     */
    inline const ModelType& model () const
    {
      return model_;
    }

    inline const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    /*!
     *   destructor: free of temporary memory for basis-fct-gradients
     *               in mygrad
     */
    ~SimpleElementMatrixIntegrator ()
    {
      if( verbose_ )
        std::cout << "entered destructor of SimpleElementMatrixIntegrator"
                  << std :: endl;
      delete[] gradPhiPtr_;
    }

    inline void addElementMatrix ( const EntityType &entity,
                                   ElementMatrixType& matrix,
                                   double coefficient = 1 ) // const
    {
      addDiffusiveFluxElementMatrix( entity, matrix, coefficient );

      if( entity.hasBoundaryIntersections() )
      {
        if( ModelType :: Properties :: hasRobinValues )
          addRobinElementMatrix( entity, matrix, coefficient );
      }
    }

    /** \brief accumulate diffusive contributions
     *
     *  The method is used for the diffusive flux of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity $E$:
     *  \f[
     *     L_ij :=  \int_E   [a     grad(phi_j) ]^T  grad(phi_i)
     *  \f]
     *  The model class is assumed to have a diffusiveFlux() method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param[in]  entity       entity over which the intrgration is performed
     *  \param      matrix       local matrix to update
     *  \param[in]  coefficient  optional weighting coefficient (defaults to 1)
     */
    template< class EntityImp, class ElementMatrixImp >
    void addDiffusiveFluxElementMatrix( const EntityImp &entity,
                                        ElementMatrixImp &matrix,
                                        double coefficient = 1 ) const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;

      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &discreteFunctionSpace = this->discreteFunctionSpace();

      const GeometryType &geometry = entity.geometry();

      // get local basis
      const BaseFunctionSetType baseSet
        =  discreteFunctionSpace.baseFunctionSet( entity );
      int numBaseFunctions = baseSet.numBaseFunctions();

      // assert that allocated space for gradPhiPtr is sufficient!!
      assert( numBaseFunctions <= numBaseFunctions_ );

      // assert that matrix allocation is sufficient
      assert( matrix.rows() >= numBaseFunctions );
      assert( matrix.cols() >= numBaseFunctions );

      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree );
      const int numQuadraturePoints = quadrature.nop();
      for ( int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const CoordinateType &x = quadrature.point( pt );

        const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
        const double volume = geometry.integrationElement( x );

        for( int i = 0; i < numBaseFunctions; ++i )
        {
          JacobianRangeType &gradPhi = gradPhiPtr_[ i ];
          baseSet.jacobian( i, quadrature[ pt ], gradPhi );
          // multiply with transposed of the jacobian inverse
          gradPhi[ 0 ] = FMatrixHelp :: mult( inv, gradPhi[ 0 ] );
        }

        // evaluate diffusiveFlux for all gradients of basis functions
        const double factor = coefficient * quadrature.weight( pt ) * volume;
        for( int j = 0; j < numBaseFunctions; ++j )
        {
          JacobianRangeType psi;
          model.diffusiveFlux( entity, quadrature[ pt ], gradPhiPtr_[ j ], psi );
          for( int i = 0; i < numBaseFunctions; ++i )
            matrix.add( i, j, factor * (psi[ 0 ] * gradPhiPtr_[ i ][ 0 ]) );
        }
      }
    } // end addDiffusiveFluxElementMatrix


    /*! addRobinElementMatrix: accumulate Robin boundary contributions
     *
     *  The method is used for the robin boundary of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[
     *     L_{ij} :=  +  \int_{\Gamma_R} alpha      phi_i        phi_j
     *  \f]
     *  The model class is assumed to have an alpha() and a
     *  boundaryType() member method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param[in]  entity       entity over which the intrgration is performed
     *  \param      matrix       local matrix to update
     *  \param[in]  coefficient  optional weighting coefficient (defaults to 1)
     */
    template< class EntityType, class ElementMatrixType >
    void addRobinElementMatrix( const EntityType &entity,
                                ElementMatrixType &matrix,
                                double coefficient = 1 ) //  const
    {
      assert( ModelType :: Properties :: hasRobinValues );
      const ModelType &model = this->model();

      enum { quadratureDegree = TraitsType :: quadDegree };

      // for all intersections check whether boundary
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();

      const GridPartType &gridPart = dfSpace.gridPart();

      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        // check for boundary
        if( !it.boundary() )
          continue;

        // check for robin boundary values
        if( model.boundaryType( it ) != ModelType :: Robin )
          continue;

        const BaseFunctionSetType baseFunctionSet
          = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();

        // integrate over intersection
        IntersectionQuadratureType quadrature
          ( gridPart, it, quadratureDegree, IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop();
        for( int pt = 0; pt < numQuadraturePoints; ++pt )
        {
          const double volume
            = it.intersectionGlobal().integrationElement( quadrature.localPoint( pt ) );
          const double alpha = model.robinAlpha( it, quadrature, pt );
          const double factor = coefficient * alpha * quadrature.weight( pt ) * volume;

          for( int i = 0; i < numBaseFunctions; ++i )
          {
            RangeType phi_i;
            baseFunctionSet.evaluate( i, quadrature[ pt ], phi_i );
            for( int j = 0; j < numBaseFunctions; ++j )
			      {
              RangeType phi_j;
              baseFunctionSet.evaluate( j, quadrature[ pt ], phi_j );
              matrix.add( i, j, factor * (phi_i[ 0 ] * phi_j[ 0 ]) );
            }
          }
        } // end loop over quadraturepoints
      } // end loop over intersections
    } // end method addRobinElementMatrix

  protected:
    //! Reference to the data model
    const ModelType& model_;
    //! the discrete function space
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    //! verbosity flag
    int verbose_;
    //! number of basis functions
    int numBaseFunctions_;
    //! storage for basis-function gradients
    JacobianRangeType *gradPhiPtr_;
  }; // end of SimpleElementMatrixIntegrator class



#endif // SIMPLEELEMENTINTEGRATOR_H
