/**
 *  \brief  ellipticmodels.hh
 *
 *  \todo   doc
 **/

#ifndef ELLIPTICMODELS_HH
#define ELLIPTICMODELS_HH

#if ! defined(MICRO_GRIDTYPE)
    #warning ("MICRO_GRIDTYPE undefined, defaulting to ALUSimplexGrid!")
    #define MICRO_GRIDTYPE=ALUSimplexGrid
    #include <dune/grid/alugrid.hh>
#elif MICRO_GRIDTYPE==ALUSimplexGrid
    #include <dune/grid/alugrid.hh>
#else
    #warning ("only ALUSimplexGrid allowed for MICRO_GRIDTYPE, defaulting to ALUSimplexGrid!")
    #define MICRO_GRIDTYPE=ALUSimplexGrid
    #include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/gridpart/periodicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/model/linearellipticmodel.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/grid.hh>

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

        using BaseType::diffusiveFlux;
        using BaseType::convectiveFlux;
        using BaseType::mass;
        using BaseType::source;

        //! constructor with functionspace argument such that the space and the
        //! grid is available
        inline DarcyModel(  const int verbosity = 0 )
            : verbosity_( verbosity )
        {

            static const int gridDim = MacroGridType::dimensionworld;

            if ( gridDim != 2 ) {
                assert( !"DarcyModel implemented in 2D only!" );
            }

            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();

            if ( verbosity_ == 0 ) {
                infoStream.Suspend();
                debugStream.Suspend();
            }

            debugStream << "\tComputing permeability tensor..." << std::endl;

            /*
             * grid
             */

            debugStream << "\tInitialising micro grid..." << std::endl;

            typedef Dune::MICRO_GRIDTYPE< gridDim, gridDim >
                MicroGridType;

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

            typedef Dune::PeriodicLeafGridPart< MicroGridType >
//            typedef Dune::AdaptiveLeafGridPart< MicroGridType >
                MicroGridPartType;
            MicroGridPartType microGridPart( *microGridPointer );

            infoStream << "\tInitialised micro grid." << std::endl;

            /*
             * micro problem
             */

            infoStream << "\tInitialising micro problem..." << std::endl;

            const int microPolOrder = MICRO_POLORDER;
            const double microViscosity = 1.0;

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
                microSolutions( "micro_",
                                microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                dummy( "dummy_", microDiscreteStokesFunctionSpaceWrapper );

            // analytcal data
            typedef typename MicroStokesModelTraitsImp::AnalyticalForceType
                MicroAnalyticalForceType;

            typename MicroAnalyticalForceType::RangeType unitVector( 0.0 );
            unitVector[ 0 ] = 1.0;

            MicroAnalyticalForceType microAnalyticalForce( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVector );

            typedef typename MicroStokesModelTraitsImp::AnalyticalDirichletDataType
                MicroAnalyticalDirichletDataType;

            typename MicroAnalyticalDirichletDataType::RangeType zero( 0.0 );

            MicroAnalyticalDirichletDataType microAnalyticalDirichletData( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), zero );

            // model
            MicroStokesModelImpType microStokesModel(   Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
                                                        microAnalyticalForce,
                                                        microAnalyticalDirichletData,
                                                        microViscosity );

            infoStream << "\tInitialised micro problem." << std::endl;

            /*
             * micro pass
             */

            debugStream << "\tInitialising micro system..." << std::endl;

            Stuff::getGridInformation( microGridPart, microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), debugStream );

//            typedef Dune::StartPass< MicroDiscreteStokesFunctionWrapperType, -1 >
//                MicroStartPassType;
//            MicroStartPassType microStartPass;
//
//            typedef Dune::StokesPass< MicroStokesModelImpType, MicroStartPassType, 0 >
//                MicroStokesPassType;
//            MicroStokesPassType microStokesPass(    microStartPass,
//                                                    microStokesModel,
//                                                    microGridPart,
//                                                    microDiscreteStokesFunctionSpaceWrapper );
//
//            microSolutions.discretePressure().clear();
//            microSolutions.discreteVelocity().clear();
//            dummy.discretePressure().clear();
//            dummy.discreteVelocity().clear();
//
//            microStokesPass.apply( dummy, microSolutions );

//            infoStream << "\tMicro system solved." << std::endl;
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

            infoStream << "Computed permeability tensor." << std::endl;

            if ( verbosity_ == 0 ) {
                infoStream.Resume();
                debugStream.Resume();
            }
        }

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
