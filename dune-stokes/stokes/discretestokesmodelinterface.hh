/**
 *  \file   discretestokesmodelinterface.hh
 *
 *  \brief  containing a class DiscreteStokesModelInterface
 **/
#ifndef DUNE_DISCRESTOKESTEMODELINTERFACE_HH
#define DUNE_DISCRESTOKESTEMODELINTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

namespace Dune
{

/**
 *  \brief  Interface class for stokes problem definition in the LDG context.
 *  \todo   doc with tex
 **/
template < class DiscreteStokesModelTraits >
class DiscreteStokesModelInterface
{
    public:


        typedef DiscreteStokesModelInterface< DiscreteStokesModelTraits >
            ThisType;

        //! Traits class defined by the user
        typedef DiscreteStokesModelTraits
            Traits;

        //! Implementation type for Barton-Nackman trick
        typedef typename Traits::DiscreteModelType
            DiscreteModelType;

        //! volume quadrature type used in pass
        typedef typename Traits::VolumeQuadratureType
            VolumeQuadratureType;

        //! face quadrature type used in pass
        typedef typename Traits::FaceQuadratureType
            FaceQuadratureType;

        //! Velocity function space
        typedef typename Traits::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! Sigma function space
        typedef typename Traits::DiscreteSigmaFunctionSpaceType
            DiscreteSigmaFunctionSpaceType;

        //! Pressure function space
        typedef typename Traits::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! Coordinate type (world coordinates)
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            DomainType;

        //! Vector type of the velocity's discrete function space's range
        typedef typename DiscreteVelocityFunctionSpaceType::RangeType
            VelocityRangeType;

        //! vector type of sigmas' discrete functions space's range
        typedef typename DiscreteSigmaFunctionSpaceType::RangeType
            SigmaRangeType;

        //! Vector type of the pressure's discrete function space's range
        typedef typename DiscretePressureFunctionSpaceType::RangeType
            PressureRangeType;

        //! Local function of type Velocity
        typedef typename Traits::LocalVelocityFunctionType
            LocalVelocityFunctionType;

        //! local function of type sigma
        typedef typename Traits::LocalSigmaFunctionType
            LocalSigmaFunctionType;

        //! Local function of type pressure
        typedef typename Traits::LocalPressureFunctionType
            LocalPressureFunctionType;

        //! Type of GridPart
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;

        //! Type of the grid
        typedef typename GridPartType::GridType
            GridType;

        //! Intersection iterator of the grid
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;

        //! Element (codim 0 entity) of the grid
        typedef typename GridType::template Codim<0>::Entity
            EntityType;

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{\sigma}\f$
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> velocitySigmaFlux and
         *              velocitySigmaBoundaryFlux as well
         **/
        bool hasVelocitySigmaFlux() const
        {
            return asImp().hasVelocitySigmaFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{p}\f$
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> velocityPressureFlux and
         *              velocityPressureBoundaryFlux as well.
         **/
        bool hasVelocityPressureFlux() const
        {
            return asImp().hasVelocityPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{p}\f$
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> pressureFlux and
         *              pressureBoundaryFlux as well.
         **/
        bool hasPressureFlux() const
        {
            return asImp().hasPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{\sigma}\f$
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> sigmaFlux and
         *              sigmaBoundaryFlux as well.
         **/
        bool hasSigmaFlux() const
        {
            return asImp().hasSigmaFlux;
        }

        /**
         *  \brief  Returns true if problem has a force contribution \f$f\f$
         *  \attention  If you let this method return true, make sure to
         *              implement force as well.
         **/
        bool hasForce() const
        {
            return asImp().hasForce();
        }

        /**
         *  \brief
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \tparam LocalVelocityFunctionType
         *          type of local function (of type velocity)
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at
         *  \param  uInner
         *          local function (of type velocity) on given entity
         *  \param  uOuter
         *          local function (of type velocity) on neighbour of
         *          given entity
         *  \todo   latex doc
         **/
        template < class FaceDomainType >
        void velocitySigmaFlux(         const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const LocalVelocityFunctionType& uInner,
                                        const LocalVelocityFunctionType& uOuter,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& uContribOuter,
                                        VelocityRangeType& emptyContribInner,
                                        VelocityRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaFlux(  it,
                                            time,
                                            x,
                                            uInner,
                                            uOuter,
                                            uContribInner,
                                            uContribOuter,
                                            emptyContribInner,
                                            emptyContribOuter ) );
        }

        /**
         *  \brief
         *  \todo   doc like velocitySigmaFlux
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const LocalVelocityFunctionType& uInner,
                                        const LocalVelocityFunctionType& uOuter,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& uContribOuter,
                                        VelocityRangeType& emptyContribInner,
                                        VelocityRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaBoundaryFlux(
                                        it,
                                        time,
                                        x,
                                        uInner,
                                        uOuter,
                                        uContribInner,
                                        uContribOuter,
                                        emptyContribInner,
                                        emptyContribOuter ) );
        }

        /**
         *  \brief
         *  \todo   latex doc
         **/
         template < class FaceDomainType >
        void velocityPressureFlux(      const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const LocalVelocityFunctionType& uInner,
                                        const LocalVelocityFunctionType& uOuter,
                                        const LocalPressureFunctionType& pInner,
                                        const LocalPressureFunctionType& pOuter,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& uContribOuter,
                                        VelocityRangeType& pContribInner,
                                        VelocityRangeType& pContribOuter,
                                        VelocityRangeType& emptyContribInner,
                                        VelocityRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureFlux(   it,
                                                time,
                                                x,
                                                uInner,
                                                uOuter,
                                                pInner,
                                                pOuter,
                                                uContribInner,
                                                uContribOuter,
                                                pContribInner,
                                                pContribOuter,
                                                emptyContribInner,
                                                emptyContribOuter) );
        }

        /**
         *  \brief
         *  \todo   doc like velocityPressureFlux
         **/
         template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const LocalVelocityFunctionType& uInner,
                                    const LocalVelocityFunctionType& uOuter,
                                    const LocalPressureFunctionType& pInner,
                                    const LocalPressureFunctionType& pOuter,
                                    VelocityRangeType& uContribInner,
                                    VelocityRangeType& uContribOuter,
                                    VelocityRangeType& pContribInner,
                                    VelocityRangeType& pContribOuter,
                                    VelocityRangeType& emptyContribInner,
                                    VelocityRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureBoundaryFlux(
                                    it,
                                    time,
                                    x,
                                    uInner,
                                    uOuter,
                                    pInner,
                                    pOuter,
                                    uContribInner,
                                    uContribOuter,
                                    pContribInner,
                                    pContribOuter,
                                    emptyContribInner,
                                    emptyContribOuter ) );
        }

        /**
         *  \brief
         *  \todo latex doc
         **/
        template < class FaceDomainType >
        void pressureFlux(  const IntersectionIteratorType& it,
                            const double time,
                            const FaceDomainType& x,
                            const LocalPressureFunctionType& pInner,
                            const LocalPressureFunctionType& pOuter,
                            PressureRangeType& pContribInner,
                            PressureRangeType& pContribOuter,
                            PressureRangeType& emptyContribInner,
                            PressureRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureFlux(   it,
                                        time,
                                        x,
                                        pInner,
                                        pOuter,
                                        pContribInner,
                                        pContribOuter,
                                        emptyContribInner,
                                        emptyContribOuter ) );
        }

        /**
         *  \brief
         *  \todo   doc like pressureFlux
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const LocalPressureFunctionType& pInner,
                                    const LocalPressureFunctionType& pOuter,
                                    PressureRangeType& pContribInner,
                                    PressureRangeType& pContribOuter,
                                    PressureRangeType& emptyContribInner,
                                    PressureRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureBoundaryFlux(   it,
                                                time,
                                                x,
                                                pInner,
                                                pOuter,
                                                pContribInner,
                                                pContribOuter,
                                                emptyContribInner,
                                                emptyContribOuter) );
        }

        /**
         *  \brief
         *  \todo   latex doc
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const LocalVelocityFunctionType& uInner,
                        const LocalVelocityFunctionType& uOuter,
                        const LocalSigmaFunctionType& sigmaInner,
                        const LocalSigmaFunctionType& sigmaOuter,
                        SigmaRangeType& sigmaContribInner,
                        SigmaRangeType& sigmaContribOuter,
                        SigmaRangeType& uContribInner,
                        SigmaRangeType& uContribOuter,
                        SigmaRangeType& emptyContribInner,
                        SigmaRangeType& emptyContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaFlux(  it,
                                    time,
                                    x,
                                    uInner,
                                    uOuter,
                                    sigmaInner,
                                    sigmaOuter,
                                    sigmaContribInner,
                                    sigmaContribOuter,
                                    uContribInner,
                                    uContribOuter,
                                    emptyContribInner,
                                    emptyContribOuter ) );
        }

        /**
         *  \brief
         *  \todo   doc like sigmaFlux
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const LocalVelocityFunctionType& uInner,
                                const LocalVelocityFunctionType& uOuter,
                                const LocalSigmaFunctionType& sigmaInner,
                                const LocalSigmaFunctionType& sigmaOuter,
                                SigmaRangeType& sigmaContribInner,
                                SigmaRangeType& sigmaContribOuter,
                                SigmaRangeType& uContribInner,
                                SigmaRangeType& uContribOuter,
                                SigmaRangeType& emptyContribInner,
                                SigmaRangeType& emptyContribOuter)
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaBoundaryFlux(  it,
                                            time,
                                            x,
                                            uInner,
                                            uOuter,
                                            sigmaInner,
                                            sigmaOuter,
                                            sigmaContribInner,
                                            sigmaContribOuter,
                                            uContribInner,
                                            uContribOuter,
                                            emptyContribInner,
                                            emptyContribOuter) );
        }

        /**
         *  \brief
         *  \todo   latex doc
         **/
        template < class FaceDomainType >
        void force( const IntersectionIteratorType& it,
                    const double time,
                    const FaceDomainType& x,
                    VelocityRangeType& forceContribInner,
                    VelocityRangeType& forceContribOuter )
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().force(
                ) );
        }

    protected:
        //! for Barton-Nackmann trick
        DiscreteModelType& asImp()
        {
            return static_cast<DiscreteModelType&>(*this);
        }

        //! for Barton-Nackmann trick
        const DiscreteModelType& asImp() const
        {
            return static_cast<const DiscreteModelType&>(*this);
        }

};

// forward declaration
template < class DiscreteStokesModelDefaultTraits >
class DiscreteStokesModelDefault;

/**
 *  \brief  traits class for DiscreteStokesModelDefault
 *  \todo   doc
 **/
template < class GridPartImp, int gridDim, int polOrder >
class DiscreteStokesModelDefaultTraits
{
    public:
        //! for Barton-Nackmann trick
        typedef DiscreteStokesModelDefault < DiscreteStokesModelDefaultTraits >
            DiscreteModelType;

        //! we use caching quadratures for the entities
        typedef Dune::CachingQuadrature< GridPartImp, 0 >
            VolumeQuadratureType;

        //! and for the faces
        typedef Dune::CachingQuadrature< GridPartImp, 1 >
            FaceQuadratureType;

        //! function space type for the velocity
        typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
            VelocityFunctionSpaceType;

        //! discrete function type space for the velocity
        typedef Dune::DiscontinuousGalerkinSpace<   VelocityFunctionSpaceType,
                                                    GridPartImp,
                                                    polOrder >
            DiscreteVelocityFunctionSpaceType;

        //! discrete function type for the velocity
        typedef Dune::AdaptiveDiscreteFunction< DiscreteVelocityFunctionSpaceType >
            DiscreteVelocityFunctionType;

        //! local function type (like a local base function) for velocity
        typedef typename DiscreteVelocityFunctionType::LocalFunctionType
            LocalVelocityFunctionType;

        //! function space type for sigma
        typedef Dune::MatrixFunctionSpace<  double,
                                            double,
                                            gridDim,
                                            gridDim,
                                            gridDim >
            SigmaFunctionSpaceType;

        //! discrete function space type for sigma
        typedef Dune::DiscontinuousGalerkinSpace<   SigmaFunctionSpaceType,
                                                    GridPartImp,
                                                    polOrder >
            DiscreteSigmaFunctionSpaceType;

        //! discrete function type for sigma
        typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
            DiscreteSigmaFunctionType;

          //! local function type (like a local base function) for sigma
        typedef typename DiscreteSigmaFunctionType::LocalFunctionType
            LocalSigmaFunctionType;

      //! function space type for the pressure
        typedef Dune::FunctionSpace< double, double, gridDim, 1 >
            PressureFunctionSpaceType;

        //! discrete function space type for the pressure
        typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
                                                    GridPartImp,
                                                    polOrder >
            DiscretePressureFunctionSpaceType;

        //! discrete function type for the pressure
        typedef Dune::AdaptiveDiscreteFunction< DiscretePressureFunctionSpaceType >
            DiscretePressureFunctionType;

        //! local function type (like a local base function) for sigma
        typedef typename DiscretePressureFunctionType::LocalFunctionType
            LocalPressureFunctionType;


};


/**
 *  \brief  definition of an ldg method for a stokes problem
 *  \todo   texdoc
 **/
template < class DiscreteStokesModelDefaultTraitsImp >
class DiscreteStokesModelDefault : public DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
{
};

}; // end of namespace Dune

#endif // end of discretestokesmodelinterface.hh
