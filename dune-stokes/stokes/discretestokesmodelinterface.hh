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

        //! (volume) quadrature type used in pass
        typedef typename Traits::QuadratureType QuadratureType;
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
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        void VelocitySigmaFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().VelocitySigmaFlux() );
            return asImp().VelocitySigmaFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        void VelocitySigmaBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().VelocitySigmaBoundaryFlux() );
            return asImp().VelocitySigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void VelocityPressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().VelocityPressureFlux() );
            return asImp().VelocityPressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocityPressureFlux contribution.
         **/
        void VelocityPressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().VelocityPressureBoundaryFlux() );
            return asImp().VelocityPressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void PressureFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().PressureFlux() );
            return asImp().PressureFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a PressureFlux contribution.
         **/
        void PressureBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().PressureBoundaryFlux() );
            return asImp().PressureBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void SigmaFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().SigmaFlux() );
            return asImp().SigmaFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a SigmaFlux contribution.
         **/
        void SigmaBoundaryFlux()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().SigmaBoundaryFlux() );
            return asImp().SigmaBoundaryFlux();
        }

        /**
         *  \brief  Empty implementation that fails if problem claims to have
         *          a VelocitySigmaFlux contribution.
         **/
        void Force()
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().Force() );
            return asImp().Force();
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
