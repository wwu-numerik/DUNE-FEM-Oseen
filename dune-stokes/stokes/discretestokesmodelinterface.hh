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
        /**
         *  \brief  Traits class defined by the user
         **/
        typedef DiscreteStokesModelTraits
            Traits;
        /**
         *  \brief  Implementation type for Barton-Nackman trick
         **/
        typedef typename Traits::DiscreteModelType
            DiscreteModelType;

        //! volume quadrature type used in pass
        typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
        //! face quadrature type used in pass
        typedef typename Traits::FaceQuadratureType FaceQuadratureType;
        /**
         *  \brief  Velocity function space
         **/
        typedef typename Traits::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;
        /**
         *  \brief  Pressure function space
         **/
        typedef typename Traits::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;
        /**
         *  \brief  Coordinate type (world coordinates)
         **/
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            DomainType;
        /**
         *  \brief  Vector type of the velocity's discrete function space's range
         **/
        typedef typename DiscreteVelocityFunctionSpaceType::RangeType
            VelocityRangeType;
        /**
         *  \brief  Vector type of the pressure's discrete function space's range
         **/
        typedef typename DiscretePressureFunctionSpaceType::RangeType
            PressureRangeType;
        /**
         *  \brief  Type of GridPart
         **/
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;
        /**
         *  \brief  Type of the grid
         **/
        typedef typename GridPartType::GridType
            GridType;
        /**
         *  \brief  Intersection iterator of the grid
         **/
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;
        /**
         *  \brief  Element (codim 0 entity) of the grid
         **/
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
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at
         *  \param  uInner
         *
         **/
        template < class ArgumentTuple, class FaceDomainType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const ArgumentTuple& uInner,
                                const ArgumentTuple& uOuter,
                                VelocityRangeType& uContribInner,
                                VelocityRangeType& uContribOuter,
                                VelocityRangeType& emptyContribInner,
                                VelocityRangeType& emptyContribOuter )
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocitySigmaFlux() );
            return asImp().velocitySigmaFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        template < class ArgumentTuple, class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const ArgumentTuple& uInner,
                                        const ArgumentTuple& uOuter,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& uContribOuter,
                                        VelocityRangeType& emptyContribInner,
                                        VelocityRangeType& emptyContribOuter )
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocitySigmaBoundaryFlux() );
            return asImp().velocitySigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void velocityPressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocityPressureFlux() );
            return asImp().velocityPressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void velocityPressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().velocityPressureBoundaryFlux() );
            return asImp().velocityPressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void pressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().pressureFlux() );
            return asImp().pressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void pressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().pressureBoundaryFlux() );
            return asImp().pressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void sigmaFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().sigmaFlux() );
            return asImp().sigmaFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void sigmaBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().sigmaBoundaryFlux() );
            return asImp().sigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        void force()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().force() );
            return asImp().force();
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

#endif // end of dixcretestokesmodelinterface.hh
