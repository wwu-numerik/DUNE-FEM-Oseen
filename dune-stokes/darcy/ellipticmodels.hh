/**
 *  \brief  ellipticmodels.hh
 *
 *  \todo   doc
 **/

#ifndef ELLIPTICMODELS_HH
#define ELLIPTICMODELS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/model/linearellipticmodel.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>

namespace Darcy
{

template< class FunctionSpaceImp >
class DarcyModel
            : public Dune::LinearEllipticModelDefault< FunctionSpaceImp, DarcyModel < FunctionSpaceImp > >
{
public:
    typedef FunctionSpaceImp
        FunctionSpaceType;

private:
    typedef DarcyModel< FunctionSpaceType >
        ThisType;

    typedef Dune::LinearEllipticModelDefault< FunctionSpaceType, ThisType >
        BaseType;

    typedef Dune::FieldMatrix< double, dimworld, dimworld >
        PermeabilityTensorType;

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
    inline DarcyModel()
    {
//            assert( dimworld == 2 );
//            permeabilityTensor_ = 0.0;
//
//            Logging::LogStream& infoStream = Logger().Info();
//            Logging::LogStream& debugStream = Logger().Dbg();
//
//            infoStream << "\nComputing permeability tensor for the darcy equation..." << std::endl;

        /* ********************************************************************** *
         * initialize the grid                                                    *
         * ********************************************************************** */
//            const int gridDim = GridType::dimensionworld;
//            Dune::GridPtr< GridType > gridPtr( "micro_2d.dgf" );
//            const int refine_level = Parameters().getParam( "micro_refine", 0 ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
//            gridPtr->globalRefine( refine_level );
//
//            typedef Dune::AdaptiveLeafGridPart< GridType >
//                MicroGridPartType;
//            MicroGridPartType microGridPart( *gridPtr );
//            infoStream << "\tInitialised the grid." << std::endl;

        /* ********************************************************************** *
         * initialize problem                                                     *
         * ********************************************************************** */
//            const int microPolOrder = MICRO_POLORDER;
//            const double microViscosity = Parameters().getParam( "micro_viscosity", 1.0 );

        // model traits
//            typedef Dune::DiscreteStokesModelDefaultTraits<
//                            MicroGridPartType,
//                            MicroForce,
//                            MicroDirichletData,
//                            gridDim,
//                            microPolOrder,
//                            microPolOrder,
//                            microPolOrder >
//                MicroStokesModelTraitsImp;
//            typedef Dune::DiscreteStokesModelDefault< MicroStokesModelTraitsImp >
//                MicroStokesModelImpType;

//            // treat as interface
//            typedef Dune::DiscreteStokesModelInterface< MicroStokesModelTraitsImp >
//                MicroStokesModelType;

        // function wrapper for the solutions
//            typedef MicroStokesModelTraitsImp::DiscreteStokesFunctionSpaceWrapperType
//                MicroDiscreteStokesFunctionSpaceWrapperType;
//
//            MicroDiscreteStokesFunctionSpaceWrapperType
//                microDiscreteStokesFunctionSpaceWrapper( microGridPart );
//
//            typedef MicroStokesModelTraitsImp::DiscreteStokesFunctionWrapperType
//                MicroDiscreteStokesFunctionWrapperType;
//
//            MicroDiscreteStokesFunctionWrapperType
//                microSolutions( "micro_",
//                                microDiscreteStokesFunctionSpaceWrapper );
//
//            MicroDiscreteStokesFunctionWrapperType
//                dummy( "dummy_", microDiscreteStokesFunctionSpaceWrapper );
//
//            typedef MicroStokesModelTraitsImp::AnalyticalForceType
//                MicroAnalyticalForceType;
//            MicroAnalyticalForceType microAnalyticalForce( microViscosity , microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );
//
//            typedef MicroStokesModelTraitsImp::AnalyticalDirichletDataType
//                MicroAnalyticalDirichletDataType;
//            MicroAnalyticalDirichletDataType microAnalyticalDirichletData( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );
//
//            MicroStokesModelImpType microStokesModel(   Dune::StabilizationCoefficients::StabilizationCoefficients( 1, 1, 1, 1, 1, 0, 1, 0 ),
//                                                        microAnalyticalForce,
//                                                        microAnalyticalDirichletData,
//                                                        microViscosity );
//
//            infoStream << "\tInitialised problem." << std::endl;
        /* ********************************************************************** *
         * initialize passes                                                      *
         * ********************************************************************** */

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
//
//            infoStream << "\tMicroPass done." << std::endl;
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


};  // end of DarcyModel

}   // end of namespace Darcy

#endif
