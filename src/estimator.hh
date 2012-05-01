#ifndef ERROR_ESTIMATOR_HH 
#define ERROR_ESTIMATOR_HH 

//- Dune-fem includes 
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh> 
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

namespace Dune
{

  // Estimator
  // ---------
  template< class DiscreteFunction >
  class Estimator
  {
    typedef Estimator< DiscreteFunction > ThisType;

  public:
    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType; 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType; 
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType; 
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType; 
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType; 
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType :: Intersection IntersectionType;

    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
    typedef typename GridType :: template Codim< 0 > :: EntityPointer
      ElementPointerType;

    static const int dimension = GridType :: dimension;

    typedef CachingQuadrature< GridPartType, 0 > ElementQuadratureType;
    typedef CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
    
    typedef std :: vector< double > ErrorIndicatorType;

  private:
    const DiscreteFunctionType &uh_;
    const DiscreteFunctionSpaceType &dfSpace_;
    GridPartType &gridPart_;
    const IndexSetType &indexSet_;
    GridType &grid_;
    ErrorIndicatorType indicator_;

  public:
    explicit Estimator ( const DiscreteFunctionType &uh )
    : uh_( uh ),
      dfSpace_( uh.space() ),
      gridPart_( dfSpace_.gridPart() ),
      indexSet_( gridPart_.indexSet() ),
      grid_( gridPart_.grid() ),
      indicator_( indexSet_.size( 0 ) )
    {}

    void clear ()
    {
      typedef typename ErrorIndicatorType :: iterator IteratorType;
      const IteratorType end = indicator_.end();
      for( IteratorType it = indicator_.begin(); it != end; ++it )
        *it = 0.0;
    }

    template< class RHSFunctionType >
    double estimate ( const RHSFunctionType &rhs )
    {
      clear();
      const IteratorType end = dfSpace_.end();
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
        estimateLocal( rhs, *it );

      double error = 0.0;
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
        error += indicator_[ indexSet_.index( *it ) ];
      std :: cout << "Estimated H1-Error: " << sqrt( error ) << std :: endl;
      return sqrt( error );
    }

    template< class RHSFunctionType >
    bool estimateAndMark ( const RHSFunctionType &rhs, const double tolerance )
    {
      const double error = estimate( rhs );
      return (error < tolerance ? false : mark( 0.9 * tolerance ));
    }

    //! calculates the error estimator and also marks the 
    //! grid for adaptation 
    template< class RHSFunctionType >
    static bool estimateAndMark ( const DiscreteFunctionType &uh,
                                  const RHSFunctionType &rhs,
                                  const double tolerance )
    {
      ThisType estimator( uh );
      return estimator.estimateAndMark( rhs, tolerance );
    }

    //! mark all elements due to given tolerance 
	bool mark ( const double /*tolerance*/ ) const
    {
      bool marked = false;
      // loop over all elements 
      const IteratorType end = dfSpace_.end();
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      {
        const ElementType &entity = *it;
        // check local error indicator 
        //if( indicator_[ indexSet_.index( entity ) ] > localTol2 )
        {
          // mark entity for refinement 
          grid_.mark( 1, entity );
          marked = true;
        }
      }
      return marked;
    }

  private:
    //! caclulate error on element 
    template< class RHSFunctionType >
    void estimateLocal ( const RHSFunctionType &rhs, const ElementType &entity )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();

      const double volume = geometry.volume();
      const double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ));
      const int index = indexSet_.index( entity );
      const LocalFunctionType uLocal = uh_.localFunction( entity );

      ElementQuadratureType quad( entity, 2*(dfSpace_.order() + 1) );
      const int numQuadraturePoints = quad.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        RangeType y;
        rhs.evaluate( geometry.global( quad.point( qp ) ), y );
        const double weight = quad.weight( qp ) * geometry.integrationElement( quad.point( qp ) );
        indicator_[ index ] += h2 * weight * (y * y);
      }

      IntersectionIteratorType end = gridPart_.iend( entity );
      for( IntersectionIteratorType it = gridPart_.ibegin( entity ); it != end; ++it )
      {
        const IntersectionType &intersection = *it;
        if( intersection.neighbor() )
          estimateIntersection( intersection, entity, uLocal );
      }
    }

    //! caclulate error on boundary intersections 
    void estimateIntersection ( const IntersectionType &intersection,
                                const ElementType &inside,
                                const LocalFunctionType &uInside )
    {
      const ElementPointerType pOutside = intersection.outside();
      const ElementType &outside = *pOutside;

      const int quadOrder = 2*(dfSpace_.order()-1);

      const int insideIndex = indexSet_.index( inside );
      const int outsideIndex = indexSet_.index( outside );

      const bool isOutsideInterior = (outside.partitionType() == InteriorEntity);
      if( !isOutsideInterior || (insideIndex < outsideIndex) )
      {
        const LocalFunctionType uOutside = uh_.localFunction( outside );

        double error;

#ifdef OLD_DUNE_GRID_VERSION
        typedef TwistUtility< GridType > TwistUtilityType;
        if( TwistUtilityType::conforming( gridPart_.grid() , intersection ) )
#else 
        if( intersection.conforming() )
#endif
        {
          FaceQuadratureType quadInside( gridPart_, intersection, quadOrder, FaceQuadratureType :: INSIDE );
          FaceQuadratureType quadOutside( gridPart_, intersection, quadOrder, FaceQuadratureType :: OUTSIDE );
          error = estimateIntersection( intersection, quadInside, quadOutside, uInside, uOutside );
        }
        else 
        {
          typedef typename FaceQuadratureType :: NonConformingQuadratureType
            NonConformingQuadratureType;
          NonConformingQuadratureType quadInside( gridPart_, intersection, quadOrder, FaceQuadratureType :: INSIDE );
          NonConformingQuadratureType quadOutside( gridPart_, intersection, quadOrder, FaceQuadratureType :: OUTSIDE );
          error = estimateIntersection( intersection, quadInside, quadOutside, uInside, uOutside );
        }

        if( error > 0.0 )
        {
          const double volume
            = 0.5 * (inside.geometry().volume() + outside.geometry().volume());
          const double h = (dimension == 1 ? volume : std :: pow( volume, 1.0 / (double)dimension ));
          indicator_[ insideIndex ] += h * error;
          if( isOutsideInterior )
            indicator_[ outsideIndex ] += h * error;
        }
      }
    }

    //! caclulate error on element intersections 
    template< class Quadrature >
    double estimateIntersection ( const IntersectionType &intersection,
                                  const Quadrature &quadInside,
                                  const Quadrature &quadOutside,
                                  const LocalFunctionType &uInside,
                                  const LocalFunctionType &uOutside ) const
    {
      double error = 0.0;
      const int numQuadraturePoints = quadInside.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        DomainType integrationNormal
          = intersection.integrationOuterNormal( quadInside.localPoint( qp ) );
        const double integrationElement = integrationNormal.two_norm();
        
        // evaluate | (d u_l * n_l) + (d u_r * n_r) | = | (d u_l - d u_r) * n_l |
        JacobianRangeType jacobianInside;
        JacobianRangeType jacobianOutside;
        uInside.jacobian( quadInside[ qp ], jacobianInside );
        uOutside.jacobian( quadOutside[ qp ], jacobianOutside );

        RangeType jump;
        jacobianInside -= jacobianOutside;
        jacobianInside.mv( integrationNormal, jump );
        
        error += quadInside.weight( qp ) * (jump * jump) / integrationElement;
      }
      return error;
    }
  };

}

#endif

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

