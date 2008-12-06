/**
 *  \file   discretestokesmodelinterface.hh
 *  \brief  containing a class DiscreteStokesModelInterface
 *          and a class DiscreteStokesModelDefault with traits class
 *          DiscreteStokesModelDefaultTraits
 *  \todo   check compatibility with certain doxygen versions
 **/
#ifndef DUNE_DISCRESTOKESTEMODELINTERFACE_HH
#define DUNE_DISCRESTOKESTEMODELINTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include "../src/stuff.hh" // should be removed in the end

namespace Dune
{

/**
 *  \brief  Interface class for stokes problem definition in the LDG context.
 *
 *          A discrete model implementation of the user should be derived from
 *          this interface class to be compatible with the StokesPass.
 *          In the LDG context a weak discrete formulation of a stokes problem
 *          enforces
 *          \f$
 *              \forall T \in \mathcal{T}, \quad \forall \tau \in \Sigma,
 *              v \in V, q \in Q
 *          \f$
 *          \f{eqnarray*}
 *              \int\limits_{T}{\sigma:\tau dx} &=&
 *                  \int\limits_{\partial T}{
 *                      \hat{u}_{\sigma} \cdot \tau \cdot n_{T} ds}
 *                  -\int\limits_{T}{u\cdot\left(\nabla\cdot\tau\right)dx},\\
 *              \mu\int\limits_{T}{\sigma :\nabla v dx}
 *                  &-&\mu\int\limits_{\partial T}{
 *                      v\cdot\hat{\sigma}\cdot n_{T} ds}
 *                  -\int\limits_{T}{p\cdot\left(\nabla\cdot v\right)dx}\\
 *              &+&\int\limits_{\partial T}{\hat{p}\cdot v\cdot n_{T}ds}
 *                  = \int\limits_{T}{f\cdot v},\\
 *              \int\limits_{\partial T}{\hat{u}_{p}\cdot n_{T}q ds}
 *                  &-&\int\limits_{T}{u\cdot\nabla q dx}=0,
 *          \f}
 *          where \f$\mathcal{T}\f$ is a triangulation and \f$\Sigma\f$,
 *          \f$V\f$, \f$Q\f$ are discrete function spaces.\n
 *          (For a detailed description see B. Cockburn, G. Kanschat,
 *          D. Schötzau, C. Schwab: <EM>Local Discontinuous Galerkin Methods
 *          for the Stokes System</EM> (2000).\n
 *          The fluxes \f$\hat{u}_{\sigma}\f$, \f$\hat{\sigma}\f$,
 *          \f$\hat{p}\f$, \f$\hat{u}_{p}\f$ in the corresponding surface
 *          integrals are implemented in the methods velocitySigmaFlux(),
 *          sigmaFlux(), pressureFlux(), velocityPressureFlux().
 *          If the face in consideration is on the boundary of \f$\Omega\f$, the
 *          computation is done by velocitySigmaBoundaryFlux(),
 *          sigmaBoundaryFlux(), pressureBoundaryFlux() and
 *          velocityPressureBoundaryFlux().\n
 *          The fluxes are designed to take values of functions on the
 *          intersection, once seen from the inside (from the entity in
 *          consideration) and once from the outside (the entities neighbour over
 *          the given intersection). Accordingly the fluxes return all
 *          contributions (to coefficients and right hand side) seen from
 *          both entities, thus saving  computational effort on half of the
 *          entities.
 *
 *  \tparam DiscreteStokesModelTraits
 *          traits class defined by the user, should provide all types needed
 *          by this interface
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

        //! discrete function type for the velocity
        typedef typename Traits::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! discrete function space type for the velocity
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! function space type for the velocity
        typedef typename DiscreteVelocityFunctionSpaceType::FunctionSpaceType
            VelocityFunctionSpaceType;

        //! discrete function type for sigma
        typedef typename Traits::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

        //! discrete function space type for sigma
        typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
            DiscreteSigmaFunctionSpaceType;

        //! function space type for sigma
        typedef typename DiscreteSigmaFunctionSpaceType::FunctionSpaceType
            SigmaFunctionSpaceType;

        //! discrete function type for the pressure
        typedef typename Traits::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! discrete function space type for the pressure
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! function space type for the pressure
        typedef typename DiscretePressureFunctionSpaceType::FunctionSpaceType
            PressureFunctionSpaceType;

        //! function type for analytical force
        typedef typename Traits::AnalyticalForceType
            AnalyticalForceType;

        typedef typename Traits::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

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

        /** \name Dummy types needed to comply to LocalPass
         *  \{
         */
        //! dummy return value of the pass
        typedef typename Traits::DestinationType
            DestinationType;

        //! dummy discrete function space belonging to DestinationType
        typedef typename Traits::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;
        /**
         *  \}
         **/


        /**
         *  \brief  constructor
         *
         *  doing nothing
         **/
        DiscreteStokesModelInterface()
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~DiscreteStokesModelInterface()
        {}

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> velocitySigmaFlux() and
         *              velocitySigmaBoundaryFlux() as well
         **/
        bool hasVelocitySigmaFlux() const
        {
            return asImp().hasVelocitySigmaFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{p}\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> velocityPressureFlux() and
         *              velocityPressureBoundaryFlux() as well.
         **/
        bool hasVelocityPressureFlux() const
        {
            return asImp().hasVelocityPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{p}\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> pressureFlux() and
         *              pressureBoundaryFlux() as well.
         **/
        bool hasPressureFlux() const
        {
            return asImp().hasPressureFlux();
        }

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement <b>both</b> sigmaFlux() and
         *              sigmaBoundaryFlux() as well.
         **/
        bool hasSigmaFlux() const
        {
            return asImp().hasSigmaFlux;
        }

        /**
         *  \brief  Returns true if problem has a force contribution \f$f\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement force() as well.
         **/
        bool hasForce() const
        {
            return asImp().hasForce();
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}
         *          \f$ for faces, that are inside \f$\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{u}_{\sigma}(u) = \hat{u}_{\sigma}^{U}
         *              + \hat{u}_{\sigma}^{RHS}
         *          \f$, where \f$\hat{u}_{\sigma}^{U}\f$ is this fluxes
         *          contribution to the coefficients of \f$U\f$ and
         *          \f$\hat{u}_{\sigma}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param  uContribInner
         *          \f$\hat{u}_{\sigma}^{U}\f$ (seen from the inside)
         *  \param  uContribOuter
         *          \f$\hat{u}_{\sigma}^{U}\f$ (seen from the outside)
         *  \param  rhsContribInner
         *          \f$\hat{u}_{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param  rhsContribOuter
         *          \f$\hat{u}_{\sigma}^{RHS}\f$ (seen from the outside)
         **/
        template < class FaceDomainType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& uInner,
                                const VelocityRangeType& uOuter,
                                VelocityRangeType& uContribInner,
                                VelocityRangeType& uContribOuter,
                                VelocityRangeType& rhsContribInner,
                                VelocityRangeType& rhsContribOuter ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaFlux(  it,
                                            time,
                                            x,
                                            uInner,
                                            uOuter,
                                            uContribInner,
                                            uContribOuter,
                                            rhsContribInner,
                                            rhsContribOuter ) );
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}
         *          \f$ for faces on \f$\partial\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{u}_{\sigma}(u) = \hat{u}_{\sigma}^{U}
         *              + \hat{u}_{\sigma}^{RHS}
         *          \f$, where \f$\hat{u}_{\sigma}^{U}\f$ is this fluxes
         *          contribution to the coefficients of \f$U\f$ and
         *          \f$\hat{u}_{\sigma}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param  uContribInner
         *          \f$\hat{u}_{\sigma}^{U}\f$
         *  \param  rhsContribInner
         *          \f$\hat{u}_{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const VelocityRangeType& uInner,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& rhsContribInner ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaBoundaryFlux(  it,
                                                    time,
                                                    x,
                                                    uInner,
                                                    uContribInner,
                                                    rhsContribInner ) );
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{p}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{u}_{p}(u,p):\Omega\rightarrow R^{d}
         *          \f$ for faces, that are inside \f$\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{u}_{p}(u,p) = \hat{u}_{p}^{U}
         *              + \hat{u}_{p}^{P} + \hat{u}_{p}^{RHS}
         *          \f$, where \f$\hat{u}_{p}^{U}\f$ is this fluxes
         *          contribution to the coefficients of \f$U\f$,
         *          \f$\hat{u}_{p}^{P}\f$ the contribution to the
         *          coefficients of \f$P\f$ and
         *          \f$\hat{u}_{p}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the inside)
         *  \param  uContribOuter
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the outside)
         *  \param  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the inside)
         *  \param  pContribOuter
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the outside)
         *  \param  rhsContribInner
         *          \f$\hat{u}_{p}^{RHS}\f$ (seen from the inside)
         *  \param  rhsContribOuter
         *          \f$\hat{u}_{p}^{RHS}\f$ (seen from the outside)
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const VelocityRangeType& uInner,
                                    const VelocityRangeType& uOuter,
                                    const PressureRangeType& pInner,
                                    const PressureRangeType& pOuter,
                                    VelocityRangeType& uContribInner,
                                    VelocityRangeType& uContribOuter,
                                    VelocityRangeType& pContribInner,
                                    VelocityRangeType& pContribOuter,
                                    VelocityRangeType& rhsContribInner,
                                    VelocityRangeType& rhsContribOuter ) const
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
                                                rhsContribInner,
                                                rhsContribOuter) );
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{p}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{u}_{p}(u,p):\Omega\rightarrow R^{d}
         *          \f$ for faces on \f$\partial\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{u}_{p}(u,p) = \hat{u}_{p}^{U}
         *              + \hat{u}_{p}^{P} + \hat{u}_{p}^{RHS}
         *          \f$, where \f$\hat{u}_{p}^{U}\f$ is this fluxes
         *          contribution to the coefficients of \f$U\f$,
         *          \f$\hat{u}_{p}^{P}\f$ the contribution to the
         *          coefficients of \f$P\f$ and
         *          \f$\hat{u}_{p}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param  pInner
         *          value of \f$p\f$ in \f$x\f$
         *  \param  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$
         *  \param  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$
         *  \param  rhsContribInner
         *          \f$\hat{u}_{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const VelocityRangeType& uInner,
                                    const PressureRangeType& pInner,
                                    VelocityRangeType& uContribInner,
                                    VelocityRangeType& pContribInner,
                                    VelocityRangeType& rhsContribInner ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureBoundaryFlux(   it,
                                                        time,
                                                        x,
                                                        uInner,
                                                        pInner,
                                                        uContribInner,
                                                        pContribInner,
                                                        rhsContribInner ) );
        }

        /**
         *  \brief  implementation of \f$\hat{p}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{p}(p):\Omega\rightarrow R
         *          \f$ for faces, that are inside \f$\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{p}(p) = \hat{p}^{P}
         *              + \hat{p}^{RHS}
         *          \f$, where \f$\hat{p}^{P}\f$ is this fluxes
         *          contribution to the coefficients of \f$P\f$ and
         *          \f$\hat{p}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param  pContribInner
         *          \f$\hat{p}^{P}\f$ (seen from the inside)
         *  \param  pContribOuter
         *          \f$\hat{p}^{P}\f$ (seen from the outside)
         *  \param  rhsContribInner
         *          \f$\hat{p}^{RHS}\f$ (seen from the inside)
         *  \param  rhsContribOuter
         *          \f$\hat{p}^{RHS}\f$ (seen from the outside)
         **/
        template < class FaceDomainType >
        void pressureFlux(  const IntersectionIteratorType& it,
                            const double time,
                            const FaceDomainType& x,
                            const PressureRangeType& pInner,
                            const PressureRangeType& pOuter,
                            PressureRangeType& pContribInner,
                            PressureRangeType& pContribOuter,
                            PressureRangeType& rhsContribInner,
                            PressureRangeType& rhsContribOuter ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureFlux(   it,
                                        time,
                                        x,
                                        pInner,
                                        pOuter,
                                        pContribInner,
                                        pContribOuter,
                                        rhsContribInner,
                                        rhsContribOuter ) );
        }

        /**
         *  \brief  implementation of \f$\hat{p}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{p}(p):\Omega\rightarrow R
         *          \f$ for faces on \f$\partial\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{p}(p) = \hat{p}^{P}
         *              + \hat{p}^{RHS}
         *          \f$, where \f$\hat{p}^{P}\f$ is this fluxes
         *          contribution to the coefficients of \f$P\f$ and
         *          \f$\hat{p}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param  pContribInner
         *          \f$\hat{p}^{P}\f$
         *  \param  rhsContribInner
         *          \f$\hat{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& pInner,
                                    PressureRangeType& pContribInner,
                                    PressureRangeType& rhsContribInner ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureBoundaryFlux(   it,
                                                time,
                                                x,
                                                pInner,
                                                pContribInner,
                                                rhsContribInner) );
        }

        /**
         *  \brief  implementation of \f$\hat{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}
         *          \f$ for faces, that are inside \f$\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{\sigma}(u,\sigma) = \hat{\sigma}^{\sigma}
         *              + \hat{\sigma}^{U} + \hat{\sigma}^{RHS}
         *          \f$, where \f$\hat{\sigma}^{\sigma}\f$ is this fluxes
         *          contribution to the coefficients of \f$\sigma\f$,
         *          \f$\hat{\sigma}^{U}\f$ the contribution to the
         *          coefficients of \f$U\f$ and
         *          \f$\hat{\sigma}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param  sigmaOuter
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the outside)
         *  \param  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the inside)
         *  \param  sigmaContribOuter
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the outside)
         *  \param  uContribInner
         *          \f$\hat{\sigma}^{U}\f$ (seen from the inside)
         *  \param  uContribOuter
         *          \f$\hat{\sigma}^{U}\f$ (seen from the outside)
         *  \param  rhsContribInner
         *          \f$\hat{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param  rhsContribOuter
         *          \f$\hat{\sigma}^{RHS}\f$ (seen from the outside)
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const VelocityRangeType& uInner,
                        const VelocityRangeType& uOuter,
                        const SigmaRangeType& sigmaInner,
                        const SigmaRangeType& sigmaOuter,
                        SigmaRangeType& sigmaContribInner,
                        SigmaRangeType& sigmaContribOuter,
                        SigmaRangeType& uContribInner,
                        SigmaRangeType& uContribOuter,
                        SigmaRangeType& rhsContribInner,
                        SigmaRangeType& rhsContribOuter ) const
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
                                    rhsContribInner,
                                    rhsContribOuter ) );
        }

        /**
         *  \brief  implementation of \f$\hat{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *
         *          Implements
         *          \f$
         *              \hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}
         *          \f$ for faces on \f$\partial\Omega\f$.\n
         *          <b>Assumption:</b> the flux can be written as
         *          \f$
         *              \hat{\sigma}(u,\sigma) = \hat{\sigma}^{\sigma}
         *              + \hat{\sigma}^{U} + \hat{\sigma}^{RHS}
         *          \f$, where \f$\hat{\sigma}^{\sigma}\f$ is this fluxes
         *          contribution to the coefficients of \f$\sigma\f$,
         *          \f$\hat{\sigma}^{U}\f$ the contribution to the
         *          coefficients of \f$U\f$ and
         *          \f$\hat{\sigma}^{RHS}\f$ the contribution
         *          to the right hand side.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param  it
         *          faceiterator
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at (on the face)
         *  \param  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$
         *  \param  uContribInner
         *          \f$\hat{\sigma}^{U}\f$
         *  \param  rhsContribInner
         *          \f$\hat{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& uInner,
                                const SigmaRangeType& sigmaInner,
                                SigmaRangeType& sigmaContribInner,
                                SigmaRangeType& uContribInner,
                                SigmaRangeType& rhsContribInner ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaBoundaryFlux(  it,
                                            time,
                                            x,
                                            uInner,
                                            sigmaInner,
                                            sigmaContribInner,
                                            uContribInner,
                                            rhsContribInner ) );
        }

        /**
         *  \brief  implementation of \f$f\f$.
         *          Calls the implementation of the derived class.
         *
         *  \tparam DomainType
         *          domain type in entity
         *  \param  time
         *          global time
         *  \param  x
         *          point to evaluate at
         *  \param  forceContrib
         *          value of \f$f\f$ in \f$x\f$
         **/
        template < class DomainType >
        void force( const double time,
                    const DomainType& x,
                    VelocityRangeType& forceContrib ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().force(  time,
                                x,
                                forceContrib ) );
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
template < class GridPartImp, class AnalyticalForceImp, class AnalyticalDirichletDataImp, int gridDim, int polOrder >
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

        //!
        typedef Dune::Function< VelocityFunctionSpaceType, AnalyticalForceImp >
            AnalyticalForceType;

        typedef Dune::Function< VelocityFunctionSpaceType, AnalyticalDirichletDataImp >
            AnalyticalDirichletDataType;


        /** \name Dummy types needed to comply with LocalPass
         *  \{
         */
        //! dummy return value of the pass
        typedef DiscreteVelocityFunctionType
            DestinationType;

        //! dummy discrete function space belonging to DestinationType
        typedef DiscreteVelocityFunctionSpaceType
            DiscreteFunctionSpaceType;
        /**
         *  \}
         **/
};


/**
 *  \brief  definition of an ldg method for a stokes problem
 *  \todo   texdoc
 **/
template < class DiscreteStokesModelDefaultTraitsImp >
class DiscreteStokesModelDefault : public DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
{
    public:

        typedef DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
            BaseType;

        typedef typename BaseType::IntersectionIteratorType
            IntersectionIteratorType;

        typedef typename BaseType::VelocityRangeType
            VelocityRangeType;

        typedef typename BaseType::SigmaRangeType
            SigmaRangeType;

        typedef typename BaseType::PressureRangeType
            PressureRangeType;

        typedef typename BaseType::AnalyticalForceType
            AnalyticalForceType;

        typedef typename BaseType::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

        /**
         *  \brief  constructor
         *
         *  set \f$C_{11}\in R\f$, \f$C_{12}\in R\f$, \f$D_{11}\in R^{d}\f$,
         *  \f$D_{12}\in R^{d}\f$
         **/
        DiscreteStokesModelDefault( const AnalyticalForceType& force,
                                    const AnalyticalDirichletDataType& dirichletData )
            : C_12_( 1.0 ),
            D_12_( 1.0 ),
            force_( force ),
            dirichletData_( dirichletData )
        {
            C_11_ = 1.0;
            D_11_ = 1.0;
        }

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~DiscreteStokesModelDefault()
        {}

        /**
         *  \brief
         *  \todo   doc me
         **/
        bool hasVelocitySigmaFlux() const
        {
            return true;
        }

        /**
         *  \brief
         *  \todo   doc me
         **/
        bool hasVelocityPressureFlux() const
        {
            return true;
        }

        /**
         *  \brief
         *  \todo   doc me
         **/
        bool hasPressureFlux() const
        {
            return true;
        }

        /**
         *  \brief
         *  \todo   doc me
         **/
        bool hasSigmaFlux() const
        {
            return true;
        }

        /**
         *  \brief
         *  \todo   doc me
         **/
        bool hasForce() const
        {
            return true;
        }

        /**
         *  \brief
         *  \todo   doc me
         *  \attention  assumption: \f$n_{inner}=-1*n_{outer}\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& uInner,
                                const VelocityRangeType& uOuter,
                                VelocityRangeType& uContribInner,
                                VelocityRangeType& uContribOuter,
                                VelocityRangeType& rhsContribInner,
                                VelocityRangeType& rhsContribOuter ) const
        {
            // some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType innerNormal = outerNormal;
            innerNormal *= -1.0;

            // contribution to u vector ( from inside entity )
            SigmaRangeType innerJump = uTypeMatrixJump( uInner,
                                                        uOuter,
                                                        outerNormal );
            innerJump.mv( C_12_, uContribInner );
            uContribInner += meanValue( uInner, uOuter );

            // contribution to u vector ( from outside entity )
            SigmaRangeType outerJump = uTypeMatrixJump( uOuter,
                                                        uInner,
                                                        innerNormal );
            outerJump.mv( C_12_, uContribOuter );
            uContribOuter += meanValue( uOuter, uInner );

            // contribution to rhs ( from inside entity )
            rhsContribInner = 0.0;

            // contribution to rhs  ( from outside entity )
            rhsContribOuter = 0.0;
        }

        /**
         *  \brief
         *  \todo   doc like velocitySigmaFlux
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const VelocityRangeType& uInner,
                                        VelocityRangeType& uContribInner,
                                        VelocityRangeType& rhsContribInner ) const
        {
            // some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType global = it->intersectionSelfLocal().global( x );

            // contribution to u vector ( from inside entity )
            uContribInner = 0.0;

            // contribution to rhs ( from inside entity )
            dirichletData_.evaluate( global,  rhsContribInner );
        }

        /**
         *  \brief
         *  \todo   doc me
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const VelocityRangeType& uInner,
                                    const VelocityRangeType& uOuter,
                                    const PressureRangeType& pInner,
                                    const PressureRangeType& pOuter,
                                    VelocityRangeType& uContribInner,
                                    VelocityRangeType& uContribOuter,
                                    VelocityRangeType& pContribInner,
                                    VelocityRangeType& pContribOuter,
                                    VelocityRangeType& rhsContribInner,
                                    VelocityRangeType& rhsContribOuter ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType innerNormal = outerNormal;
            innerNormal *= -1.0;

            // contribution to u vector ( from inside entity )
            uContribInner = D_12_;
            uContribInner *= uTypeJump( uInner,
                                        uOuter,
                                        outerNormal );
            uContribInner += meanValue( uInner,
                                        uOuter );

            // contribution to u vector ( from outside entity )
            uContribOuter = D_12_;
            uContribOuter = uTypeJump(  uOuter,
                                        uInner,
                                        innerNormal );
            uContribInner += meanValue( uOuter,
                                        uInner );

            // contribution to p vector ( from inside entity )
            pContribInner = pTypeJump(  pInner,
                                        pOuter,
                                        outerNormal );
            pContribInner *= D_11_;

            // contribution to p vector ( from outside entity )
            pContribOuter = pTypeJump(  pOuter,
                                        pInner,
                                        innerNormal );
            pContribInner *= D_11_;

            // contribution to rhs ( from inside entity )
            rhsContribInner = 0.0;

            // contribution to rhs ( from outside entity )
            rhsContribInner = 0.0;
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
                                    const VelocityRangeType& uInner,
                                    const PressureRangeType& pInner,
                                    VelocityRangeType& uContribInner,
                                    VelocityRangeType& pContribInner,
                                    VelocityRangeType& rhsContribInner ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType global = it->intersectionSelfLocal().global( x );

            // contribution to u vector ( from inside entity )
            uContribInner = 0.0;

            // contribution to p vector ( from inside entity )
            pContribInner = 0.0;

            // contribution to rhs ( from inside entity )
            dirichletData_.evaluate( global, rhsContribInner );
        }

        /**
         *  \brief
         *  \todo latex doc
         **/
        template < class FaceDomainType >
        void pressureFlux(  const IntersectionIteratorType& it,
                            const double time,
                            const FaceDomainType& x,
                            const PressureRangeType& pInner,
                            const PressureRangeType& pOuter,
                            PressureRangeType& pContribInner,
                            PressureRangeType& pContribOuter,
                            PressureRangeType& rhsContribInner,
                            PressureRangeType& rhsContribOuter ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType innerNormal = outerNormal;
            innerNormal *= -1.0;

            // contribution to p vector ( from inside entity )
            pContribInner = -1.0
                * ( D_12_ * pTypeJump( pInner, pOuter, outerNormal ) );
            pContribInner += meanValue( pInner, pOuter );

            // contribution to p vector ( from outside entity )
            pContribInner = -1.0
                * ( D_12_ * pTypeJump( pOuter, pInner, innerNormal ) );
            pContribInner += meanValue( pOuter, pInner );

            // contribution to rhs ( from inside entity )
            rhsContribInner = 0.0;

            // contribution to rhs ( from outside entity )
            rhsContribOuter = 0.0;
        }

        /**
         *  \brief
         *  \todo   doc like pressureFlux
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& pInner,
                                    PressureRangeType& pContribInner,
                                    PressureRangeType& rhsContribInner ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType global = it->intersectionSelfLocal().global( x );

            // contribution to p vector ( from inside entity )
            pContribInner = pInner;

            // contribution to rhs ( from inside entity )
            rhsContribInner = 0.0;
        }

        /**
         *  \brief
         *  \todo   latex doc
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const VelocityRangeType& uInner,
                        const VelocityRangeType& uOuter,
                        const SigmaRangeType& sigmaInner,
                        const SigmaRangeType& sigmaOuter,
                        SigmaRangeType& sigmaContribInner,
                        SigmaRangeType& sigmaContribOuter,
                        SigmaRangeType& uContribInner,
                        SigmaRangeType& uContribOuter,
                        SigmaRangeType& rhsContribInner,
                        SigmaRangeType& rhsContribOuter ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType innerNormal = outerNormal;
            innerNormal *= -1.0;

            // contribution to sigma vector ( from inside entity )
            sigmaContribInner = dyadicProduct(
                sigmaTypeJump(
                    sigmaInner, sigmaOuter, outerNormal ),
                C_12_ );
            sigmaContribInner += meanValue( sigmaInner, sigmaOuter );

            // contribution to sigma vector ( from outside entity )
            sigmaContribInner = dyadicProduct(
                sigmaTypeJump(
                    sigmaOuter, sigmaInner, innerNormal ),
                C_12_ );
            sigmaContribInner += meanValue( sigmaOuter, sigmaInner );

            // contribution to u vector ( from inside entity )
            uContribInner = uTypeMatrixJump( uInner, uOuter, outerNormal );
            uContribInner *= -1.0 * C_11_;

            // contribution to u vector ( from outside entity )
            uContribInner = uTypeMatrixJump( uOuter, uInner, innerNormal );
            uContribInner *= -1.0 * C_11_;

            // contribution to rhs ( from inside entity )
            rhsContribInner = 0.0;

            // contribution to rhs ( from outside entity )
            rhsContribOuter = 0.0;
        }

        /**
         *  \brief
         *  \todo   doc like sigmaFlux
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& uInner,
                                const SigmaRangeType& sigmaInner,
                                SigmaRangeType& sigmaContribInner,
                                SigmaRangeType& uContribInner,
                                SigmaRangeType& rhsContribInner ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType global = it->intersectionSelfLocal().global( x );

            // contribution to sigma vector ( from inside entity )
            sigmaContribInner = sigmaInner;

            // contribution to u vector ( from inside entity )
            uContribInner = dyadicProduct( uInner, outerNormal );
            uContribInner *= ( -1.0 * C_11_ );

            // contribution to rhs ( from inside entity )
            VelocityRangeType gD( 0.0 );
            dirichletData_.evaluate( global, gD );
            rhsContribInner = dyadicProduct( gD, outerNormal );
            rhsContribInner *= C_11_;
        }

        /**
         *  \brief
         *  \todo   latex doc
         **/
        template < class DomainType >
        void force( const double time,
                    const DomainType& x,
                    VelocityRangeType& forceContrib ) const
        {
            force_.evaluate( x, forceContrib );
        }

    private:

        double C_11_, D_11_;
        const VelocityRangeType C_12_, D_12_;
        const AnalyticalForceType& force_;
        const AnalyticalDirichletDataType& dirichletData_;

        /**
         *  \brief  jump for pressure-type functions
         *
         *  \f$\left[\left[\p\right]\right]:=\left(p^{+} + p^{-}\right)n^{+}\in R^{d}\f$,
         *  where \f$n^{+}\f4 is the unit outer normal,
         *  \f$p^{+}\f$ is the value of p on the inside and
         *  \f$p^{-}\f$ the value of p at the outside
         *  \attention  assumption: \f$n_{inner}=-1*n_{outer}\f$
         **/
        template < class NormalType >
        VelocityRangeType pTypeJump(    const PressureRangeType& pInner,
                                        const PressureRangeType& pOuter,
                                        const NormalType& outerNormal ) const
        {
            VelocityRangeType ret = outerNormal;
            ret *= ( pInner - pOuter );
            return ret;
        }

        /**
         *  \brief  jump for velocity-type functions
         *  \todo doc like pTypeJump
         **/
        double uTypeJump(   const VelocityRangeType& uInner,
                            const VelocityRangeType& uOuter,
                            const VelocityRangeType& outerNormal ) const
        {
            return ( uInner - uOuter ) * outerNormal;
        }

        /**
         *  \brief  matrix valued jump for velocity-type functions
         *  \todo   doc like pTypeJump
         **/
        SigmaRangeType uTypeMatrixJump( const VelocityRangeType& uInner,
                                        const VelocityRangeType& uOuter,
                                        const VelocityRangeType& outerNormal ) const
        {
            SigmaRangeType ret( 0.0 );
            VelocityRangeType uDiff = uInner - uOuter;
            ret = dyadicProduct( uDiff, outerNormal );
            return ret;
        }

        /**
         *  \brief  dyadic product
         *  \todo   doc
         **/
        SigmaRangeType dyadicProduct(   const VelocityRangeType& arg1,
                                        const VelocityRangeType& arg2 ) const
        {
            SigmaRangeType ret( 0.0 );
            typedef typename SigmaRangeType::RowIterator
                MatrixRowIteratorType;
            typedef typename VelocityRangeType::ConstIterator
                ConstVectorIteratorType;
            typedef typename VelocityRangeType::Iterator
                VectorIteratorType;
            MatrixRowIteratorType rItEnd = ret.end();
            ConstVectorIteratorType arg1It = arg1.begin();
            for ( MatrixRowIteratorType rIt = ret.begin(); rIt != rItEnd; ++rIt ) {
                ConstVectorIteratorType arg2It = arg2.begin();
                VectorIteratorType vItEnd = rIt->end();
                for (   VectorIteratorType vIt = rIt->begin();
                        vIt != vItEnd;
                        ++vIt ) {
                    *vIt = *arg1It * *arg2It;
                    ++arg2It;
                }
                ++arg1It;
            }
            return ret;
        }

        /**
         *  \brief  jump for sigma-type functions
         *  \todo   doc
         **/
        VelocityRangeType sigmaTypeJump(    const SigmaRangeType& sInner,
                                            const SigmaRangeType& sOuter,
                                            const VelocityRangeType& outerNormal ) const
        {
            VelocityRangeType ret( 0.0 );
            SigmaRangeType sDiff = sInner;
            sDiff -= sOuter;
            sDiff.mv( outerNormal, ret );
            return ret;
        }

        /**
         *  \brief  mean value of two functions (of same type)
         *  \todo   texdoc example
         **/
        template < class DiscreteFunctionImp >
        DiscreteFunctionImp meanValue( const DiscreteFunctionImp& funcInner,
                                    const DiscreteFunctionImp& funcOuter ) const
        {
            DiscreteFunctionImp ret( 0.0 );
            ret += funcInner;
            ret += funcOuter;
            ret *= 0.5;
            return ret;
        }


};

}; // end of namespace Dune

#endif // end of discretestokesmodelinterface.hh
