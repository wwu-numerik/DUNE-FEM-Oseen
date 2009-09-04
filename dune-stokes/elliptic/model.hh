/**************************************************************************
**       Title: examples of elliptic models
**    $RCSfile$
**   $Revision: 2745 $$Name$
**       $Date: 2007-12-05 21:08:12 +0100 (Wed, 05 Dec 2007) $
**   Copyright: GPL $Author: nolte $
** Description: implementation of default elliptic models for
**              use in solvers such as FEOp. A model is a class providing
**              all data functions of a problem. The model is assumed to
**              be instantiated once and passed to solvers. It is depending
**              on a traits class providing typedefs.
**
**************************************************************************/
#ifndef DUNE_ELLIPTICMODEL_HH
#define DUNE_ELLIPTICMODEL_HH

#include <sstream>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/stokes/saddlepoint_inverse_operator.hh>


#ifndef NLOG // if we want logging, should be removed in the end
    #include <dune/stuff/printing.hh>
    #include <dune/stuff/misc.hh>
    #include <dune/stuff/logging.hh>
#endif

#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

#include <dune/stuff/profiler.hh>
#include <dune/stuff/parametercontainer.hh>
#include <vector>
#include <string>

#include <iostream>
#include <cmath>
#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/stokes/discretestokesfunctionspacewrapper.hh>
#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include "../src/analyticaldata.hh"

#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/postprocessing.hh>
#include <dune/stuff/profiler.hh>

namespace Dune
{

    struct PoissonModelProperties
    {
        enum { hasDirichletValues = true };
        enum { hasNeumannValues = true };
        enum { hasRobinValues = true };
        enum { hasGeneralizedNeumannValues = true };
        enum { hasConvectiveFlux = false };
        enum { hasMass = false };
        enum { hasSource = true };
    };

    /** \class PoissonModel
     *  \brief The PoissonModel class provides a default model for an elliptic
     *  problem to be handled by FEOp
     *
     *  The model problem simply is
     *  \f[ - \mathrm{div} \nabla u = n \pi^2 \prod_{i=1}^n \sin( \pi x_i ). \f]
     *  Using homogeneous Dirichlet boundary values, the exact solution on the unit
     *  square is \f$u( x ) = \prod_{i=1}^n \sin( \pi x_i ) \f$
     *
     *  All types are extracted from the TraisImp, which defaults to
     *  DefaultElementMatrixTraits.
     */
    template< class FunctionSpaceImp >
    class PoissonModel
                : public LinearEllipticModelDefault
                < FunctionSpaceImp, PoissonModel< FunctionSpaceImp >, PoissonModelProperties >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

        typedef PoissonModelProperties Properties;

    private:
        typedef PoissonModel< FunctionSpaceType > ThisType;
        typedef LinearEllipticModelDefault
        < FunctionSpaceType, ThisType, Properties >
        BaseType;

    public:
        typedef typename BaseType :: BoundaryType BoundaryType;

        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        using BaseType :: diffusiveFlux;
        using BaseType :: source;

    public:
        template< class IntersectionType >
        inline BoundaryType boundaryType( const IntersectionType &intersection ) const
        {
            return BaseType :: Dirichlet;
        }

        template< class IntersectionType, class QuadratureType >
        inline void dirichletValues( const IntersectionType &intersection,
                                     const QuadratureType &quadrature,
                                     int p,
                                     RangeType &ret ) const
        {
            typedef typename IntersectionType::Entity EntityType;

            const int dimension = DomainType::dimension;

            const DomainType &x = intersection.inside()->geometry().global( quadrature.point( p ) );

            ret[ 0 ] = 1.0;
            for ( int i = 0; i < dimension; ++i )
                ret[ 0 ] *= sin( M_PI * x[ i ] );
        }

        template< class IntersectionType, class QuadratureType >
        inline void neumannValues( const IntersectionType &intersection,
                                   const QuadratureType &quadrature,
                                   int p,
                                   RangeType &ret ) const
        {
            std :: cout << "Neumann boundary values are not implemented." << std :: endl;
            assert( false );

            ret[ 0 ] = 0.0;
        }

        //! determine robin value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void robinValues( const IntersectionType &intersection,
                                 const QuadratureType &quadrature,
                                 int p,
                                 RangeType &ret ) const
        {
            std :: cout << "Robin boundary values are not implemented." << std :: endl;
            assert( false );

            //const DomainType &x = entity.geometry().global( quadrature.point( p ) );

            ret[ 0 ] = 0.0;
        }

        template< class EntityType, class PointType >
        inline void source ( const EntityType &entity,
                             const PointType &x,
                             RangeType &ret ) const
        {
            const int dimension = DomainType :: dimension;

            const DomainType &global = entity.geometry().global( coordinate( x ) );

            ret[ 0 ] = (dimension * M_PI * M_PI);
            for ( int i = 0; i < dimension; ++i )
                ret[ 0 ] *= sin( M_PI * global[ i ] );
        }

        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradient,
                                    JacobianRangeType &flux ) const
        {
            flux = gradient;
        }

        template< class IntersectionType, class QuadratureType >
        inline double robinAlpha( const IntersectionType &intersection,
                                  const QuadratureType &quadrature,
                                  int p ) const
        {
            return 1.0;
        }
    };  // end of PoissonModel class



    /*======================================================================*/
    /*!
     *  \class PoissonExactSolution
     *  \brief The class provides the exact solution for the model given by
     *         the PoissonModel class
     *
     *  The function represents u = x(1-x)y(1-y), which is the solution of
     *  the model problem  - div grad u = 2(x(1-x)+ y(1-y)) on
     *  the unit square with homogeneous Dirichlet boundary values.
     *  The function can be used for EOC calculation
     */
    /*======================================================================*/
    template< class FunctionSpaceImp >
    class PoissonExactSolution
                : public Function< FunctionSpaceImp, PoissonExactSolution< FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef PoissonExactSolution< FunctionSpaceType > ThisType;
        typedef Function< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        inline PoissonExactSolution ( const FunctionSpaceType &functionSpace )
                : BaseType( functionSpace )
        {
        }

        /*
        //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
        void evaluate (const DomainType & x , RangeType & ret) const
        {
          ret = 1.0;
          for(int i=0; i<DomainType::dimension; i++)
              ret *= ( x[i] - SQR(x[i]) );
          //    ret[0] += x[0]+x[1];
          // ret += 1.0;   // add dirichlet-values!!
        }
        */

        //! u( x ) = sin( pi x_1 ) * ... * sin( pi x_n )
        inline void evaluate ( const DomainType &x, RangeType &ret ) const
        {
            enum { dimension = DomainType :: dimension };

            ret[ 0 ] = 1.0;
            for ( int i = 0; i < dimension; ++i )
                ret[ 0 ] *= sin( M_PI * x[ i ] );
        }

        inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
        {
            evaluate( x , ret );
        }
    };



    /*======================================================================*/
    /*!
     *  \class Elliptic2dModel
     *  \brief The Elliptic2dModel class provides a complete model
     *         for an elliptic problem to be handled by FEOp
     *
     *  The model problem is defined on the 2D unit square with
     *        - div ( a grad u - b u) + cu = f
     *            u = g_D on Dirichlet-boundary (upper and lower edge)
     *      (a grad u - bu ) n = g_N on Neuman boundary (left edge)
     *      (a grad u - bu ) n + alpha u= g_R on Robin boundary (right edge)
     *
     *  The data functions are parametrized by nonnegative scalar q,r,s for
     *  switching on/off certain contributions and regulating the stability
     *     alpha = 1
     *     a = [(1+q), -q; -q,  (1+q)]
     *     b = [1; 1] * s * y
     *     c = xy * r
     *     f = 2q + r* x^2y^2 + r* x^2y + s * (x+y+y^2 + 2xy)
     *     g_D = xy+x
     *     g_N = -(1+q)(y+1)
     *     g_R = 2 - s*y^2 + (2+q-s)* y
     *
     *  The solution is simply u = xy+x
     *
     *  All types are extracted from the DefaultElementMatrixTraits class
     *  additionally, the model contains member variables and methods.
     */
    /*======================================================================*/
    template< class FunctionSpaceImp >
    class Elliptic2dModel
                : public LinearEllipticModelDefault< FunctionSpaceImp, Elliptic2dModel< FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef Elliptic2dModel< FunctionSpaceType > ThisType;
        typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

        static const double q = 1.0; // q>0 => non-unit diffusivity
        static const double r = 1.0; // r>0 => mass term activated
        static const double s = 1.0; // 0.001; // s>0 => convective term activated

    public:
        typedef typename BaseType :: BoundaryType BoundaryType;

    public:
        using BaseType :: diffusiveFlux;
        using BaseType :: convectiveFlux;
        using BaseType :: mass;
        using BaseType :: source;

    public:
        //! constructor with functionspace argument such that the space and the
        //! grid is available
        Elliptic2dModel()
        {
            // currently implementation is fitted to 2 dimensions
            assert( dimworld == 2 );
        }

        //! return boundary type of a boundary point p used in a quadrature
        template< class IntersectionType >
        inline BoundaryType boundaryType ( const IntersectionType &intersection ) const
        {
            const int boundaryId = intersection.boundaryId();

            switch ( boundaryId )
            {
            case 1:
                return BaseType :: Dirichlet;

            case 2:
                return BaseType :: Neumann;

            case 3:
                return BaseType :: Robin;

            default:
                std :: ostringstream stream;
                stream << "Unknown boundary id: " << boundaryId << ".";
                DUNE_THROW( RangeError, stream.str() );
            }
        }

        //! determine dirichlet value in a boundary point used in a quadrature
        template <class IntersectionType, class QuadratureType >
        inline void dirichletValues ( const IntersectionType &intersection,
                                      const QuadratureType& quad, int p,
                                      RangeType& ret ) const
        {
            const DomainType &x = intersection.inside()->geometry().global( quad.point( p ) );
            ret = x[ 0 ] * (1 + x[ 1 ]);
        }

        //! determine neumann value in a boundary point used in a quadrature
        template <class IntersectionType, class QuadratureType>
        inline void neumannValues ( const IntersectionType &intersection,
                                    const QuadratureType& quad, int p,
                                    RangeType& ret) const
        {
            const DomainType &x = intersection.inside()->geometry().global( quad.point( p ) );
            ret = -(1 + q) * (x[ 1 ] + 1);
        }

        //! determine robin value in a boundary point used in a quadrature
        template <class IntersectionType, class QuadratureType>
        inline void robinValues ( const IntersectionType &intersection,
                                  const QuadratureType& quad, int p,
                                  RangeType& ret) const
        {
            const DomainType &x = intersection.inside()->geometry().global( quad.point( p ) );
            ret = 2 - s * SQR( x[ 1 ] ) + (2 + q - s) * x[ 1 ];
            //ret = (2 + q - s * x[ 1 ]) * (1 + x[ 1 ]) - q;
        }

        template< class EntityType, class PointType >
        inline void mass ( const EntityType &entity,
                           const PointType &x,
                           RangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret = r * global[ 0 ] * global[ 1 ];
        }

        template< class EntityType, class PointType >
        inline void source( const EntityType &entity,
                            const PointType &x,
                            RangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret = 2 * q
                  + s * (global[ 0 ] + global[ 1 ]) * (1 + global[ 1 ])
                  + s * global[ 0 ] * global[ 1 ]
                  + r * SQR( global[ 0 ] ) * global[ 1 ] * (1 + global[ 1 ]);
        }

        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradphi,
                                    JacobianRangeType &ret ) const
        {
            ret[ 0 ][ 0 ] = (1 + q) * gradphi[ 0 ][ 0 ] - q * gradphi[ 0 ][ 1 ];
            ret[ 0 ][ 1 ] = (1 + q) * gradphi[ 0 ][ 1 ] - q * gradphi[ 0 ][ 0 ];
        }

        template< class EntityType, class PointType >
        inline void convectiveFlux( const EntityType &entity,
                                    const PointType &x,
                                    const RangeType &phi,
                                    JacobianRangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret[ 0 ][ 0 ] = -global[ 1 ] * s * phi[ 0 ];
            ret[ 0 ][ 1 ] = -global[ 1 ] * s * phi[ 0 ];
        }

        //! the coefficient for robin boundary condition
        template< class IntersectionType, class QuadratureType >
        inline RangeFieldType robinAlpha ( const IntersectionType &intersection,
                                           const QuadratureType &quadrature,
                                           int pt ) const
        {
            return 1;
        }
    };  // end of Elliptic2dModel class



    /*======================================================================*/
    /*!
     *  \class Elliptic2dExactSolution
     *  \brief The class provides the exact solution for the model given by
     *         the Elliptic2dModel class
     *
     *  The function represents u = x y + x, which is the solution of
     *  the model problem on the unit square with inhomogeneous Dirichlet
     *  Neumann and Dirichlet boundary values.
     *  Function can be used for EOC calculation
     */
    /*======================================================================*/
    template< class FunctionSpaceImp >
    class Elliptic2dExactSolution
                : public Function< FunctionSpaceImp, Elliptic2dExactSolution< FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef Elliptic2dExactSolution< FunctionSpaceType > ThisType;
        typedef Function< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        inline Elliptic2dExactSolution ( const FunctionSpaceType &functionSpace )
                : BaseType( functionSpace )
        {
        }

        void evaluate ( const DomainType &x, RangeType &ret ) const
        {
            ret = x[0]*x[1] + x[0];
        }

        void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
        {
            evaluate ( x , ret );
        }
    }; // end class Elliptic2dExactSolution



    /*======================================================================*/
    /*!
     *  \class Elliptic3dModel
     *  \brief The Elliptic3dModel class provides a complete model
     *         for an elliptic problem to be handled by FEOp
     *
     *  The model problem is defined on the 3D unit cube with
     *        - div ( a grad u - b u) + cu = f
     *      (a grad u - bu ) n = g_N on Neuman boundary (x=0 face)
     *      (a grad u - bu ) n + alpha u= g_R on Robin boundary (x=1 face)
     *            u = g_D on Dirichlet-boundary (remaining 4 faces)
     *
     *  The data functions are
     *     alpha = 1
     *     a = [3 -1 -1; -1 3 -1; -1 -1 3]
     *     b = [1; 1; 1] y
     *     c = xy
     *     f = 2z+3y+3x+y^2z+2xyz+xy^2+x^2y^2z+x^2y
     *     g_D = xyz+x
     *     g_N = -3yz-3
     *     g_R = 4yz-2y-z+4-y^2z
     *
     *  The solution is simply u = xyz+x
     *
     *  All types are extracted from the DefaultElementMatrixTraits class
     *  additionally, the model contains member variables and methods.
     */
    /*======================================================================*/
    template< class FunctionSpaceImp >
    class Elliptic3dModel
                : public LinearEllipticModelDefault< FunctionSpaceImp, Elliptic3dModel< FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef Elliptic3dModel< FunctionSpaceType > ThisType;
        typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename BaseType :: BoundaryType BoundaryType;

        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        using BaseType :: diffusiveFlux;
        using BaseType :: convectiveFlux;
        using BaseType :: mass;
        using BaseType :: source;

    public:
        //! constructor with functionspace argument such that the space and the
        //! grid is available
        inline Elliptic3dModel ()
        {
            // currently implementation is fitted to 3 dimensions
            assert(dimworld==3);
        }

        //! return boundary type of a boundary point p used in a quadrature
        template< class IntersectionType >
        inline BoundaryType boundaryType( const IntersectionType &intersection ) const
        {
            const int boundaryId = intersection.boundaryId();

            switch ( boundaryId )
            {
            case 1:
                return BaseType :: Dirichlet;

            case 2:
                return BaseType :: Neumann;

            case 3:
                return BaseType :: Robin;

            default:
                DUNE_THROW( RangeError, "Unknown boundary id." );
            }
        }

        //! determine dirichlet value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void dirichletValues( const IntersectionType &intersection,
                                     const QuadratureType& quad, int p,
                                     RangeType& ret) const
        {
            const DomainType& glob = intersection.inside()->geometry().global(quad.point(p));
            ret[0] = glob[0] * ( 1.0 + glob[1]*glob[2]);
        }

        //! determine neumann value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void neumannValues( const IntersectionType &intersection,
                                   const QuadratureType& quad, int p,
                                   RangeType& ret) const
        {
            const DomainType& glob = intersection.inside()->geometry().global(quad.point(p));
            ret[0] = -3.0 * glob[1]*glob[2] - 3.0;
        }

        //! determine robin value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void robinValues( const IntersectionType &intersection,
                                 const QuadratureType& quad, int p,
                                 RangeType& ret) const
        {
            const DomainType& glob = intersection.inside()->geometry().global(quad.point(p));
            ret[0] = 4 * glob[1] * glob[2] - 2* glob[1] - glob[2] +
                     4.0 - SQR(glob[1])*glob[2];
        }

        //! function c from top doc
        template< class EntityType, class PointType >
        inline void mass ( const EntityType &entity,
                           const PointType &x,
                           RangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret = global[ 0 ] * global[ 1 ];
        }

        //! function f from top doc
        template< class EntityType, class PointType >
        inline void source ( const EntityType &entity,
                             const PointType &x,
                             RangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret = 2 * global[2] + 3* global[1] + 3 * global[0] +
                  SQR(global[1])*global[2] + 2* global[0] * global[1]* global[2] +
                  global[0] * SQR(global[1]) + SQR(global[0]*global[1]) * global[2] +
                  SQR(global[0]) * global[1];
        }

        //! function a from top doc
        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradient,
                                    JacobianRangeType &flux ) const
        {
            const DomainType &grad = gradient[ 0 ];
            flux[ 0 ][ 0 ] = 3 * grad[ 0 ] - grad[ 1 ]- grad[ 2 ];
            flux[ 0 ][ 1 ] = -grad[ 0 ] + 3 * grad[ 1 ] - grad[ 2 ];
            flux[ 0 ] [2 ] = -grad[ 0 ] - grad[ 1 ] + 3 * grad[ 2 ];
        }

        //! function b from top doc
        template< class EntityType, class PointType >
        inline void convectiveFlux( const EntityType &entity,
                                    const PointType &x,
                                    const RangeType &phi,
                                    JacobianRangeType &ret ) const
        {
            const DomainType global = entity.geometry().global( coordinate( x ) );
            ret[ 0 ][ 0 ] = -global[ 1 ] * phi[ 0 ];
            ret[ 0 ][ 1 ] = -global[ 1 ] * phi[ 0 ];
            ret[ 0 ][ 2 ] = -global[ 1 ] * phi[ 0 ];
        }

        //! the coefficient for robin boundary condition
        template< class IntersectionType, class QuadratureType >
        inline RangeFieldType robinAlpha ( const IntersectionType &intersection,
                                           const QuadratureType &quadrature,
                                           int pt ) const
        {
            return 1;
        }
    };  // end of Elliptic3dModel class

    /*======================================================================*/
    /*!
     *  \class Elliptic3dExactSolution
     *  \brief The class provides the exact solution for the model given by
     *         the Elliptic2dModel class
     *
     *  The function represents u = x y z + x, which is the solution of
     *  the model problem on the unit cube with inhomogeneous Dirichlet
     *  Neumann and Dirichlet boundary values.
     *  Function can be used for EOC calculation
     */
    /*======================================================================*/
    template< class FunctionSpaceImp >
    class Elliptic3dExactSolution
                : public Function< FunctionSpaceImp, Elliptic3dExactSolution< FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef Elliptic3dExactSolution< FunctionSpaceType > ThisType;
        typedef Function< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        inline Elliptic3dExactSolution ( const FunctionSpaceType &functionSpace )
                : BaseType( functionSpace )
        {
        }

        inline void evaluate (const DomainType & x , RangeType & ret) const
        {
            ret = x[0]*x[1]*x[2] + x[0];
        }

        inline void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
        {
            evaluate ( x , ret );
        }
    }; // end class Elliptic3dExactSolution

    struct AortaModelProperties
    {
        enum { hasDirichletValues = true };
        enum { hasNeumannValues = true };
        enum { hasRobinValues = false };
        enum { hasGeneralizedNeumannValues = false };
        enum { hasConvectiveFlux = false };
        enum { hasMass = false };
        enum { hasSource = false };
    };

    template< class FunctionSpaceImp >
    class AortaModel
                : public LinearEllipticModelDefault< FunctionSpaceImp, AortaModel < FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp FunctionSpaceType;

    private:
        typedef AortaModel< FunctionSpaceType > ThisType;
        typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType > BaseType;

    public:
        typedef typename BaseType :: BoundaryType BoundaryType;

        typedef typename FunctionSpaceType :: DomainType DomainType;
        typedef typename FunctionSpaceType :: RangeType RangeType;
        typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    public:
        using BaseType :: diffusiveFlux;
        using BaseType :: convectiveFlux;
        using BaseType :: mass;
        using BaseType :: source;

    public:
        //! constructor with functionspace argument such that the space and the
        //! grid is available
        inline AortaModel ()
        {
            // currently implementation is fitted to 3 dimensions
            assert(dimworld==3);
        }

        //! return boundary type of a boundary point p used in a quadrature
        template< class IntersectionType >
        inline BoundaryType boundaryType( const IntersectionType &intersection ) const
        {
            const int boundaryId = intersection.boundaryId();
            switch ( boundaryId )
            {
            case 2:
                return BaseType :: Dirichlet;
            default:
                return BaseType :: Neumann;
            }
        }

        //! determine dirichlet value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void dirichletValues( const IntersectionType &intersection,
                                     const QuadratureType& quad, int p,
                                     RangeType& ret) const
        {
            double fac = 10.0;
            ret = RangeType( fac );
        }

        //! determine neumann value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void neumannValues( const IntersectionType &intersection,
                                   const QuadratureType& quad, int p,
                                   RangeType& ret) const
        {
            const int boundaryId = intersection.boundaryId();
            double fac = 100.0;
            switch ( boundaryId )
            {
            case 1:
                ret[0] = 0;
                break;

            case 2:
                ret[0] = fac;
                break;

            case 3:
            case 4:
            case 5:
            case 6:
                ret[0] = -1 * fac;
                break;

            default:
                DUNE_THROW( RangeError, "Unknown boundary id." );
            }
        }

        //! determine robin value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void robinValues( const IntersectionType &intersection,
                                 const QuadratureType& quad, int p,
                                 RangeType& ret) const
        {
            assert( false );
        }

        //! function c from top doc
        template< class EntityType, class PointType >
        inline void mass ( const EntityType &entity,
                           const PointType &x,
                           RangeType &ret ) const
        {
            const DomainType &global = entity.geometry().global( coordinate( x ) );
            ret = global[ 0 ] * global[ 1 ];
        }

        //! function f from top doc
        template< class EntityType, class PointType >
        inline void source ( const EntityType &entity,
                             const PointType &x,
                             RangeType &ret ) const
        {
            ret = RangeType( 0.0 );
        }

        //! function a from top doc
        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradient,
                                    JacobianRangeType &flux ) const
        {
            const DomainType &grad = gradient[ 0 ];
            flux[ 0 ][ 0 ] = 3 * grad[ 0 ] - grad[ 1 ]- grad[ 2 ];
            flux[ 0 ][ 1 ] = -grad[ 0 ] + 3 * grad[ 1 ] - grad[ 2 ];
            flux[ 0 ] [2 ] = -grad[ 0 ] - grad[ 1 ] + 3 * grad[ 2 ];
        }

        //! function b from top doc
        template< class EntityType, class PointType >
        inline void convectiveFlux( const EntityType &entity,
                                    const PointType &x,
                                    const RangeType &phi,
                                    JacobianRangeType &ret ) const
        {
            const DomainType global = entity.geometry().global( coordinate( x ) );
            ret[ 0 ][ 0 ] = -global[ 1 ] * phi[ 0 ];
            ret[ 0 ][ 1 ] = -global[ 1 ] * phi[ 0 ];
            ret[ 0 ][ 2 ] = -global[ 1 ] * phi[ 0 ];
        }

        //! the coefficient for robin boundary condition
        template< class IntersectionType, class QuadratureType >
        inline RangeFieldType robinAlpha ( const IntersectionType &intersection,
                                           const QuadratureType &quadrature,
                                           int pt ) const
        {
            assert( false );
        }
    };  // end of AortaModel class

    struct DarcyModelProperties
    {
        enum { hasDirichletValues = false };
        enum { hasNeumannValues = true };
        enum { hasRobinValues = false };
        enum { hasGeneralizedNeumannValues = false };
        enum { hasConvectiveFlux = false };
        enum { hasMass = false };
        enum { hasSource = false };
    };

    template< class FunctionSpaceImp >
    class DarcyModel
                : public LinearEllipticModelDefault< FunctionSpaceImp, DarcyModel < FunctionSpaceImp > >
    {
    public:
        typedef FunctionSpaceImp
            FunctionSpaceType;

    private:
        typedef DarcyModel< FunctionSpaceType >
            ThisType;

        typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType >
            BaseType;

        typedef Dune::FieldMatrix< double, dimworld, dimworld >
            PermeabilityTensorType;

    public:
        typedef typename BaseType :: BoundaryType
            BoundaryType;

        typedef typename FunctionSpaceType :: DomainType
            DomainType;

        typedef typename FunctionSpaceType :: RangeType
            RangeType;

        typedef typename FunctionSpaceType :: JacobianRangeType
        JacobianRangeType;

        typedef typename FunctionSpaceType :: DomainFieldType
            DomainFieldType;

        typedef typename FunctionSpaceType :: RangeFieldType
            RangeFieldType;

    public:
        using BaseType::diffusiveFlux;
        using BaseType::convectiveFlux;
        using BaseType::mass;
        using BaseType::source;

    public:
        //! constructor with functionspace argument such that the space and the
        //! grid is available
        inline DarcyModel()
        {
            assert( dimworld == 2 );
            permeabilityTensor_ = 0.0;

            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();

            infoStream << "\nComputing permeability tensor for the darcy equation..." << std::endl;

            /* ********************************************************************** *
             * initialize the grid                                                    *
             * ********************************************************************** */
            const int gridDim = GridType::dimensionworld;
            Dune::GridPtr< GridType > gridPtr( "micro_2d.dgf" );
            const int refine_level = Parameters().getParam( "micro_refine", 0 ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
            gridPtr->globalRefine( refine_level );

            typedef Dune::AdaptiveLeafGridPart< GridType >
                MicroGridPartType;
            MicroGridPartType microGridPart( *gridPtr );
            infoStream << "\tInitialised the grid." << std::endl;

            /* ********************************************************************** *
             * initialize problem                                                     *
             * ********************************************************************** */
            const int microPolOrder = MICRO_POLORDER;
            const double microViscosity = Parameters().getParam( "micro_viscosity", 1.0 );

            // model traits
            typedef Dune::DiscreteStokesModelDefaultTraits<
                            MicroGridPartType,
                            MicroForce,
                            MicroDirichletData,
                            gridDim,
                            microPolOrder,
                            microPolOrder,
                            microPolOrder >
                MicroStokesModelTraitsImp;
            typedef Dune::DiscreteStokesModelDefault< MicroStokesModelTraitsImp >
                MicroStokesModelImpType;

//            // treat as interface
//            typedef Dune::DiscreteStokesModelInterface< MicroStokesModelTraitsImp >
//                MicroStokesModelType;

            // function wrapper for the solutions
            typedef MicroStokesModelTraitsImp::DiscreteStokesFunctionSpaceWrapperType
                MicroDiscreteStokesFunctionSpaceWrapperType;

            MicroDiscreteStokesFunctionSpaceWrapperType
                microDiscreteStokesFunctionSpaceWrapper( microGridPart );

            typedef MicroStokesModelTraitsImp::DiscreteStokesFunctionWrapperType
                MicroDiscreteStokesFunctionWrapperType;

            MicroDiscreteStokesFunctionWrapperType
                microSolutions( "micro_",
                                microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                dummy( "dummy_", microDiscreteStokesFunctionSpaceWrapper );

            typedef MicroStokesModelTraitsImp::AnalyticalForceType
                MicroAnalyticalForceType;
            MicroAnalyticalForceType microAnalyticalForce( microViscosity , microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );

            typedef MicroStokesModelTraitsImp::AnalyticalDirichletDataType
                MicroAnalyticalDirichletDataType;
            MicroAnalyticalDirichletDataType microAnalyticalDirichletData( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );

            MicroStokesModelImpType microStokesModel(   Dune::StabilizationCoefficients::StabilizationCoefficients( 1, 1, 1, 1, 1, 0, 1, 0 ),
                                                        microAnalyticalForce,
                                                        microAnalyticalDirichletData,
                                                        microViscosity );

            infoStream << "\tInitialised problem." << std::endl;
            /* ********************************************************************** *
             * initialize passes                                                      *
             * ********************************************************************** */

            typedef Dune::StartPass< MicroDiscreteStokesFunctionWrapperType, -1 >
                MicroStartPassType;
            MicroStartPassType microStartPass;

            typedef Dune::StokesPass< MicroStokesModelImpType, MicroStartPassType, 0 >
                MicroStokesPassType;
            MicroStokesPassType microStokesPass(    microStartPass,
                                                    microStokesModel,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            microSolutions.discretePressure().clear();
            microSolutions.discreteVelocity().clear();
            dummy.discretePressure().clear();
            dummy.discreteVelocity().clear();

            microStokesPass.apply( dummy, microSolutions );

            infoStream << "\tMicroPass done." << std::endl;
//
//    /* ********************************************************************** *
//     * Problem postprocessing
//     * ********************************************************************** */
//    infoStream << "\n- postprocesing" << std::endl;
//
//
//    profiler().StartTiming( "Problem/Postprocessing" );
//
//#ifndef COCKBURN_PROBLEM //bool tpl-param toggles ana-soltion output in post-proc
//    typedef Problem< gridDim, DiscreteStokesFunctionWrapperType, false >
//        ProblemType;
//#else
//    typedef Problem< gridDim, DiscreteStokesFunctionWrapperType, true >
//        ProblemType;
//#endif
//    ProblemType problem( viscosity , computedSolutions );
//
//    typedef PostProcessor< StokesPassType, ProblemType >
//        PostProcessorType;
//
//    PostProcessorType postProcessor( discreteStokesFunctionSpaceWrapper, problem );
//
//    postProcessor.save( *gridPtr, computedSolutions, refine_level );
//    info.L2Errors = postProcessor.getError();
//    typedef Dune::StabilizationCoefficients::ValueType
//        Pair;
//    info.c11 = Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
//    info.c12 = Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
//    info.d11 = Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
//    info.d12 = Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
//    info.bfg = Parameters().getParam( "do-bfg", true );
//    info.gridname = gridPart.grid().name();
//    info.refine_level = refine_level;
//
//    info.polorder_pressure = StokesModelTraitsImp::pressureSpaceOrder;
//    info.polorder_sigma = StokesModelTraitsImp::sigmaSpaceOrder;
//    info.polorder_velocity = StokesModelTraitsImp::velocitySpaceOrder;
//
//    info.solver_accuracy = Parameters().getParam( "absLimit", 1e-4 );
//    info.inner_solver_accuracy = Parameters().getParam( "inner_absLimit", 1e-4 );
//    info.bfg_tau = Parameters().getParam( "bfg-tau", 0.1 );
//
//    profiler().StopTiming( "Problem/Postprocessing" );
//    profiler().StopTiming( "SingleRun" );
//
//    return info;

        }

        //! return boundary type of a boundary point p used in a quadrature
        template< class IntersectionType >
        inline BoundaryType boundaryType( const IntersectionType &intersection ) const
        {
//            const int boundaryId = intersection.boundaryId();
//            switch ( boundaryId )
//            {
//            case 2:
//                return BaseType :: Dirichlet;
//            default:
            return BaseType :: Neumann;
//            }
        }

        //! determine dirichlet value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void dirichletValues( const IntersectionType &intersection,
                                     const QuadratureType& quad, int p,
                                     RangeType& ret) const
        {
            assert( !"There should be no Dirichlet Values" );
//            double fac = 10.0;
//            ret = RangeType( fac );
        }

        //! determine neumann value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void neumannValues( const IntersectionType &intersection,
                                   const QuadratureType& quad, int p,
                                   RangeType& ret) const
        {
//            const int boundaryId = intersection.boundaryId();
//            double fac = 100.0;
//            switch ( boundaryId )
//            {
//            case 1:
//                ret[0] = 0;
//                break;
//
//            case 2:
//                ret[0] = fac;
//                break;
//
//            case 3:
//            case 4:
//            case 5:
//            case 6:
//                ret[0] = -1 * fac;
//                break;
//
//            default:
//                DUNE_THROW( RangeError, "Unknown boundary id." );
//            }
            ret = 0.0;
        }

        //! determine robin value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        inline void robinValues( const IntersectionType &intersection,
                                 const QuadratureType& quad, int p,
                                 RangeType& ret) const
        {
            assert( !"There should be no Robin Values!" );
        }

        //! function c from top doc
        template< class EntityType, class PointType >
        inline void mass ( const EntityType &entity,
                           const PointType &x,
                           RangeType &ret ) const
        {
//            const DomainType &global = entity.geometry().global( coordinate( x ) );
//            ret = global[ 0 ] * global[ 1 ];
            ret = 0.0;
        }

        //! function f from top doc
        template< class EntityType, class PointType >
        inline void source ( const EntityType &entity,
                             const PointType &x,
                             RangeType &ret ) const
        {
//            ret = RangeType( 0.0 );
            ret = 1.0;
        }

        //! function a from top doc
        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradient,
                                    JacobianRangeType &flux ) const
        {
            const DomainType& grad = gradient[ 0 ];
            const DomainType& ret = flux[ 0 ];
            flux = 1.0;
        }

        //! function b from top doc
        template< class EntityType, class PointType >
        inline void convectiveFlux( const EntityType &entity,
                                    const PointType &x,
                                    const RangeType &phi,
                                    JacobianRangeType &ret ) const
        {
//            const DomainType global = entity.geometry().global( coordinate( x ) );
//            ret[ 0 ][ 0 ] = -global[ 1 ] * phi[ 0 ];
//            ret[ 0 ][ 1 ] = -global[ 1 ] * phi[ 0 ];
//            ret[ 0 ][ 2 ] = -global[ 1 ] * phi[ 0 ];
            ret = 0.0;
        }

        //! the coefficient for robin boundary condition
        template< class IntersectionType, class QuadratureType >
        inline RangeFieldType robinAlpha ( const IntersectionType &intersection,
                                           const QuadratureType &quadrature,
                                           int pt ) const
        {
            assert( false );
        }

    private:

        PermeabilityTensorType permeabilityTensor_;


    };  // end of DarcyModel class

} // end namespace Dune

#endif
