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
 *
 *  \todo   doc with tex
 **/
template < class DiscreteStokesModelTraits >
class DiscreteStokesModelInterface
{
    public:


    typedef DiscreteStokesModelInterface<DiscreteStokesModelTraits>
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
        typedef typename Traits::DiscreteSigmaFucntionSpaceType
            DiscreteSigmaFucntionSpaceType;

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
        typedef typename DiscreteSigmaFucntionSpaceType::RangeType
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
         *
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> VelocitySigmaFlux and
         *              VelocitySigmaBoundaryFlux as well
         **/
        bool hasVelocitySigmaFlux() const
        {
            return asImp().hasVelocitySigmaFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{p}\f$
         *
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> VelocityPressureFlux and
         *              VelocityPressureBoundaryFlux as well.
         **/
        bool hasVelocityPressureFlux() const
        {
            return asImp().hasVelocityPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{p}\f$
         *
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> PressureFlux and
         *              PressureBoundaryFlux as well.
         **/
        bool hasPressureFlux() const
        {
            return asImp().hasPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{\sigma}\f$
         *
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> SigmaFlux and
         *              SigmaBoundaryFlux as well.
         **/
        bool hasSigmaFlux() const
        {
            return asImp().hasSigmaFlux;
        }

        /**
         *  \brief  Returns true if problem has a force contribution \f$f\f$
         *
         *  \attention  If you let this method return true, make sure to
         *              implement Force as well.
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
        template < class FaceDomainType, class LocalVelocityFunctionType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
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
                                            emptyContribOuter )
            );
            asImp().velocitySigmaFlux(  it,
                                        time,
                                        x,
                                        uInner,
                                        uOuter,
                                        uContribInner,
                                        uContribOuter,
                                        emptyContribInner,
                                        emptyContribOuter );
        }

        /**
         *  \brief
         **/
        template < class ArgumentTuple, class FaceDomainType >
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
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocitySigmaBoundaryFlux() );
            asImp().velocitySigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void velocityPressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocityPressureFlux() );
            asImp().velocityPressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void velocityPressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocityPressureBoundaryFlux() );
            asImp().velocityPressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void pressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().pressureFlux() );
            asImp().pressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void pressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().pressureBoundaryFlux() );
            asImp().pressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void sigmaFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().sigmaFlux() );
            asImp().sigmaFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void sigmaBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().sigmaBoundaryFlux() );
            asImp().sigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        void force()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().force() );
            asImp().force();
        }

    protected:
        DiscreteModelType& asImp()
        {
            return static_cast<DiscreteModelType&>(*this);
        }
        const DiscreteModelType& asImp() const
        {
            return static_cast<const DiscreteModelType&>(*this);
        }

};

}; // end of namespace Dune

#endif // end of discretestokesmodelinterface.hh
