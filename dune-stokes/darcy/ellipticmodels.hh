/**
 *  \brief  ellipticmodels.hh
 *
 *  \todo   doc
 **/

#ifndef ELLIPTICMODELS_HH
#define ELLIPTICMODELS_HH

#ifndef ENABLE_ALUGRID
    #define ENABLE_ALUGRID
#endif
#ifndef HAVE_ALUGRID
    #define HAVE_ALUGRID
#endif
#ifndef ALUGRID_SIMPLEX
    #define ALUGRID_SIMPLEX
#endif
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/gridpart/periodicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/model/linearellipticmodel.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

#include "analyticaldarcydata.hh"

namespace Darcy
{

template< class FunctionSpaceImp, class MacroGridPartImp >
class DarcyModel
            : public Dune::LinearEllipticModelDefault< FunctionSpaceImp, DarcyModel < FunctionSpaceImp, MacroGridPartImp > >
{
    public:
        typedef FunctionSpaceImp
            FunctionSpaceType;

    private:
        typedef DarcyModel< FunctionSpaceType, MacroGridPartImp >
            ThisType;

        typedef Dune::LinearEllipticModelDefault< FunctionSpaceType, ThisType >
            BaseType;

        typedef Dune::FieldMatrix< double, dimworld, dimworld >
            PermeabilityTensorType;

        typedef MacroGridPartImp
            MacroGridPartType;

        typedef typename MacroGridPartType::GridType
            MacroGridType;

    public:
        typedef typename BaseType::BoundaryType
            BoundaryType;

        typedef typename FunctionSpaceType::DomainType
            DomainType;

        typedef typename FunctionSpaceType::RangeType
            RangeType;

        typedef typename FunctionSpaceType::JacobianRangeType
            JacobianRangeType;

        typedef typename FunctionSpaceType::DomainFieldType
            DomainFieldType;

        typedef typename FunctionSpaceType::RangeFieldType
            RangeFieldType;

        typedef Dune::ALUSimplexGrid< gridDim, gridDim >
            MicroGridType;

        typedef Dune::AdaptiveLeafGridPart< MicroGridType >
            MicroGridPartType;

        using BaseType::diffusiveFlux;
        using BaseType::convectiveFlux;
        using BaseType::mass;
        using BaseType::source;

        //! constructor with functionspace argument such that the space and the
        //! grid is available
        inline DarcyModel(  const int verbosity = 0 )
            : verbosity_( verbosity )
        {

            static const int gridDim = dimworld;

            if ( gridDim != 2 ) {
                assert( !"DarcyModel implemented in 2D only!" );
            }

            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();

            if ( verbosity_ == 0) {
                infoStream.Suspend();
                debugStream.Suspend();
            }
            else if ( verbosity_ == 1 ) {
                debugStream.Suspend();
            }

            /*
             * grid
             */

            debugStream << "\tInitialising micro grid..." << std::endl;

            std::string microGridFile( "micro_grid_2d.dgf" );
            if ( gridDim == 2 ) {
                microGridFile = Dune::Parameter::getValue( "micro_grid_2d", std::string("micro_grid_2d.dgf") );
            }
            else if ( gridDim == 3 ) {
                assert( !"Darcy only implemented in 2D!" );
                microGridFile = Dune::Parameter::getValue( "micro_grid_3d", std::string("micro_grid_3d.dgf") );
            }
            else {
                assert( !"Darcy only implemented in 2D and 3D!" );
            }

            Dune::GridPtr< MicroGridType > microGridPointer( microGridFile );
            const int refine_level = Dune::Parameter::getValue( "micro_refine", 0 ) * Dune::DGFGridInfo< MicroGridType >::refineStepsForHalf();
            microGridPointer->globalRefine( refine_level );

            MicroGridPartType microGridPart( *microGridPointer );

            infoStream << "\tInitialised micro grid (with " << microGridPart.grid().size( 0 ) << " elements)." << std::endl;

            /*
             * micro problem
             */

            debugStream << "\tInitialising micro problem..." << std::endl;

            const int microPolOrder = MICRO_POLORDER;
            const double microViscosity = Dune::Parameter::getValue( "micro_viscosity", 1.0 );

            typedef Dune::DiscreteStokesModelDefaultTraits<
                            MicroGridPartType,
                            Darcy::ConstantFunction,
                            Darcy::ConstantFunction,
                            gridDim,
                            microPolOrder,
                            microPolOrder,
                            microPolOrder >
                MicroStokesModelTraitsImp;

            typedef Dune::DiscreteStokesModelDefault< MicroStokesModelTraitsImp >
                MicroStokesModelImpType;

            // treat as interface
            typedef Dune::DiscreteStokesModelInterface< MicroStokesModelTraitsImp >
                MicroStokesModelType;

            // function wrapper for the solutions
            typedef typename MicroStokesModelType::DiscreteStokesFunctionSpaceWrapperType
                MicroDiscreteStokesFunctionSpaceWrapperType;

            MicroDiscreteStokesFunctionSpaceWrapperType
                microDiscreteStokesFunctionSpaceWrapper( microGridPart );

            typedef typename MicroStokesModelType::DiscreteStokesFunctionWrapperType
                MicroDiscreteStokesFunctionWrapperType;

            MicroDiscreteStokesFunctionWrapperType
                microSolutionsX( "microX", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                microSolutionsY( "microY", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                dummy( "dummy", microDiscreteStokesFunctionSpaceWrapper );

            typedef typename MicroStokesModelTraitsImp::AnalyticalForceType
                MicroAnalyticalForceType;

            typedef typename MicroStokesModelTraitsImp::AnalyticalDirichletDataType
                MicroAnalyticalDirichletDataType;

            // micro solution for x dimnesion
            typename MicroAnalyticalForceType::RangeType unitVectorX( 0.0 );
            unitVectorX[ 0 ] = 1.0;

            MicroAnalyticalForceType microAnalyticalForceX( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorX );

            MicroAnalyticalDirichletDataType microAnalyticalDirichletDataX( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorX );

            // model x
            MicroStokesModelImpType microStokesModelX(  Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
                                                        microAnalyticalForceX,
                                                        microAnalyticalDirichletDataX,
                                                        microViscosity );

            // micro solution for y dimnesion
            typename MicroAnalyticalForceType::RangeType unitVectorY( 0.0 );
            unitVectorX[ 1 ] = 1.0;

            MicroAnalyticalForceType microAnalyticalForceY( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorY );

            MicroAnalyticalDirichletDataType microAnalyticalDirichletDataY( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorY );

            // model y
            MicroStokesModelImpType microStokesModelY(  Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
                                                        microAnalyticalForceY,
                                                        microAnalyticalDirichletDataY,
                                                        microViscosity );

            /*
             * micro pass
             */

            typedef Dune::StartPass< MicroDiscreteStokesFunctionWrapperType, -1 >
                MicroStartPassType;
            MicroStartPassType microStartPass;

            typedef Dune::StokesPass< MicroStokesModelImpType, MicroStartPassType, 0 >
                MicroStokesPassType;
            MicroStokesPassType microStokesPassX(   microStartPass,
                                                    microStokesModelX,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            MicroStokesPassType microStokesPassY(   microStartPass,
                                                    microStokesModelY,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            infoStream << "\tInitialised micro problem." << std::endl;

            debugStream << "\tSolving micro system..." << std::endl;

            dummy.discretePressure().clear();
            dummy.discreteVelocity().clear();
            microSolutionsX.discretePressure().clear();
            microSolutionsX.discreteVelocity().clear();

            microStokesPassX.apply( dummy, microSolutionsX );

            dummy.discretePressure().clear();
            dummy.discreteVelocity().clear();
            microSolutionsY.discretePressure().clear();
            microSolutionsY.discreteVelocity().clear();

            microStokesPassY.apply( dummy, microSolutionsY );

            infoStream << "\tMicro system solved." << std::endl;

            debugStream << "\tWriting micro output..." << std::endl;

            typedef Dune::VTKIO< MicroGridPartType >
                MicroVTKWriterType;

            MicroVTKWriterType microVtkWriter( microGridPart );

            microVtkWriter.addVertexData( microSolutionsX.discreteVelocity() );
            microVtkWriter.write( "data/microVelocityX" );
            microVtkWriter.clear();

            microVtkWriter.addVertexData( microSolutionsX.discretePressure() );
            microVtkWriter.write( "data/microPressureX" );
            microVtkWriter.clear();

            microVtkWriter.addVertexData( microSolutionsY.discreteVelocity() );
            microVtkWriter.write( "data/microVelocityY" );
            microVtkWriter.clear();

            microVtkWriter.addVertexData( microSolutionsY.discretePressure() );
            microVtkWriter.write( "data/microPressureY" );
            microVtkWriter.clear();

            infoStream << "\tMicro Output written." << std::endl;

            debugStream << "\tComputing permeability tensor...";

            PermeabilityTensorType permeabilityTensor_ = computePermeabilityTensor( microSolutionsX.discreteVelocity(),
                                                                                    microSolutionsY.discreteVelocity() );

            Stuff::printFieldMatrix( permeabilityTensor_, "permeabilityTensor_", debugStream, "\t" );
            debugStream << std::endl;

            infoStream << "\tPermeability tensor computed." << std::endl;

            if ( verbosity_ == 0) {
                infoStream.Resume();
                debugStream.Resume();
            }
            else if ( verbosity_ == 1 ) {
                debugStream.Resume();
            }
        }

    private:

        template< class DiscreteFunctionImp >
        PermeabilityTensorType computePermeabilityTensor( DiscreteFunctionImp& discreteFunctionX, DiscreteFunctionImp& discreteFunctionY )
        {
            PermeabilityTensorType permeabilityTensor( 0.0 );

            double gradientMicroX_times_gradientMicroX( 0.0 );
            double gradientMicroX_times_gradientMicroY( 0.0 );
            double gradientMicroY_times_gradientMicroX( 0.0 );
            double gradientMicroY_times_gradientMicroY( 0.0 );

            typedef typename MicroGridPartType::template Codim< 0 >::IteratorType
                EntityIteratorType;

            typedef typename MicroGridType::template Codim< 0 >::Entity
                EntityType;

            typedef typename EntityType::Geometry
                EntityGeometryType;

            typedef Dune::CachingQuadrature< MicroGridPartType, 0 >
                VolumeQuadratureType;

            typedef DiscreteFunctionImp
                DiscreteFunctionType;

            typedef typename DiscreteFunctionType::LocalFunctionType
                LocalFunctionType;

            typedef typename LocalFunctionType::JacobianRangeType
                JacobianRangeType;

            typedef typename Dune::FieldMatrix< typename EntityGeometryType::ctype,
                                                EntityGeometryType::coorddimension,
                                                EntityGeometryType::mydimension >
                JacobianInverseTransposedType;

            bool doOnce = true;

//            double min( 0.0 );
//            double max( 0.0 );


            EntityIteratorType entityIteratorEnd = discreteFunctionX.space().end();
            for (   EntityIteratorType entityIterator = discreteFunctionX.space().begin();
                    entityIterator != entityIteratorEnd;
                    ++entityIterator ) {

                // get entity and geometry
                const EntityType& entity = *entityIterator;
                const EntityGeometryType& geometry = entity.geometry();

                // get the local functions
                LocalFunctionType localFunctionX = discreteFunctionX.localFunction( entity );
                LocalFunctionType localFunctionY = discreteFunctionY.localFunction( entity );

                // quadrature
                VolumeQuadratureType quadrature( entity, ( 2 * discreteFunctionX.space().order() ) + 1 );

                for ( int quad = 0; quad < quadrature.nop(); ++quad ) {

                    // quadrature point
                    const DomainType x = quadrature.point( quad );

                    // get the integration factor
                    const double elementVolume = geometry.integrationElement( x );
                    // get the quadrature weight
                    const double integrationWeight = quadrature.weight( quad );

                    // get gradients
                    const JacobianInverseTransposedType jacobianInverseTransposed = geometry.jacobianInverseTransposed( x );

                    JacobianRangeType gradient_untransposed_x( 0.0 );
                    localFunctionX.jacobian( x, gradient_untransposed_x );
                    JacobianRangeType gradient_x = Stuff::rowWiseMatrixMultiplication( jacobianInverseTransposed, gradient_untransposed_x );

                    JacobianRangeType gradient_untransposed_y( 0.0 );
                    localFunctionY.jacobian( x, gradient_untransposed_y );
                    JacobianRangeType gradient_y = Stuff::rowWiseMatrixMultiplication( jacobianInverseTransposed, gradient_untransposed_y );

//                    for ( int i = 0; i < 2; ++i ) {
//                        for ( int j = 0; j < 2; ++j )
//                        {
//                            min = std::min( gradient_x[i][j], min );
//                            max = std::max( gradient_x[i][j], max );
//                        }
//                    }


                    const DomainType xWorld = geometry.global( x );
                    if ( !( ( xWorld[0] < 0.0 ) || ( xWorld[0] > 1.0 ) ) ) { // only in unit square
                        if ( !( ( xWorld[1] < 0.0 ) || ( xWorld[1] > 1.0 ) ) ) {
                            gradientMicroX_times_gradientMicroX += elementVolume * integrationWeight * Stuff::colonProduct( gradient_x, gradient_x );
                            gradientMicroX_times_gradientMicroY += elementVolume * integrationWeight * Stuff::colonProduct( gradient_x, gradient_y );
                            gradientMicroY_times_gradientMicroX += elementVolume * integrationWeight * Stuff::colonProduct( gradient_y, gradient_x );
                            gradientMicroY_times_gradientMicroY += elementVolume * integrationWeight * Stuff::colonProduct( gradient_y, gradient_y );
                        }
                    } // end of unit suqare test
                } // end of element quadrature
            } // end of gridwalk

            permeabilityTensor[0][0] = gradientMicroX_times_gradientMicroX;
            permeabilityTensor[0][1] = gradientMicroX_times_gradientMicroY;
            permeabilityTensor[1][0] = gradientMicroY_times_gradientMicroX;
            permeabilityTensor[1][1] = gradientMicroY_times_gradientMicroY;

//            std::cout << "\t\t\tmin: " << min << std::endl;
//            std::cout << "\t\t\tmax: " << max << std::endl;

            return permeabilityTensor;
        }

    public:

        //! return boundary type of a boundary point p used in a quadrature
        template< class IntersectionType >
        BoundaryType boundaryType( const IntersectionType &intersection ) const
        {
            return BaseType::Neumann;
        }

        //! determine dirichlet value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        void dirichletValues( const IntersectionType &intersection,
                                     const QuadratureType& quad, int p,
                                     RangeType& ret) const
        {
            assert( !"There should be no Dirichlet Values" );
        }

        //! determine neumann value in a boundary point used in a quadrature
        template< class IntersectionType, class QuadratureType >
        void neumannValues( const IntersectionType &intersection,
                                   const QuadratureType& quad, int p,
                                   RangeType& ret) const
        {
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
            ret = 0.0;
        }

        //! function f from top doc
        template< class EntityType, class PointType >
        inline void source ( const EntityType &entity,
                             const PointType &x,
                             RangeType &ret ) const
        {
            ret = 0.0;
        }

        //! function a from top doc
        template< class EntityType, class PointType >
        inline void diffusiveFlux ( const EntityType &entity,
                                    const PointType &x,
                                    const JacobianRangeType &gradient,
                                    JacobianRangeType &flux ) const
        {
            const DomainType& grad = gradient[ 0 ];
            DomainType& ret = flux[ 0 ];
            ret = grad;
        }

        //! function b from top doc
        template< class EntityType, class PointType >
        inline void convectiveFlux( const EntityType &entity,
                                    const PointType &x,
                                    const RangeType &phi,
                                    JacobianRangeType &ret ) const
        {
            ret = 0.0;
        }

        //! the coefficient for robin boundary condition
        template< class IntersectionType, class QuadratureType >
        inline RangeFieldType robinAlpha ( const IntersectionType &intersection,
                                           const QuadratureType &quadrature,
                                           int pt ) const
        {
            assert( !"There should be no robin alpha!" );
        }

    private:

        PermeabilityTensorType permeabilityTensor_;
        const int verbosity_;

};  // end of DarcyModel

}   // end of namespace Darcy

#endif
