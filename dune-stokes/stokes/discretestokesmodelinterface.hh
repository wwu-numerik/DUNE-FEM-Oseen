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

#include <dune/stokes/discretefunctionspacepair.hh>
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
 *          both entities, thus enabling the StokesPass to save computational
 *          effort on half of the entities.
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

        //! matrix type of sigmas' discrete functions space's range
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{\sigma}^{U}\f$ (seen from the inside)
         *  \param[out]  uContribOuter
         *          \f$\hat{u}_{\sigma}^{U}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{u}_{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{\sigma}^{U}\f$
         *  \param[out]  rhsContribInner
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the inside)
         *  \param[out]  uContribOuter
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the outside)
         *  \param[out]  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the inside)
         *  \param[out]  pContribOuter
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{u}_{p}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$
         *  \param[out]  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$
         *  \param[out]  rhsContribInner
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  pContribInner
         *          \f$\hat{p}^{P}\f$ (seen from the inside)
         *  \param[out]  pContribOuter
         *          \f$\hat{p}^{P}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{p}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[out]  pContribInner
         *          \f$\hat{p}^{P}\f$
         *  \param[out]  rhsContribInner
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[in]  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  sigmaOuter
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the inside)
         *  \param[out]  sigmaContribOuter
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the outside)
         *  \param[out]  uContribInner
         *          \f$\hat{\sigma}^{U}\f$ (seen from the inside)
         *  \param[out]  uContribOuter
         *          \f$\hat{\sigma}^{U}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param[out]  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$
         *  \param[out]  uContribInner
         *          \f$\hat{\sigma}^{U}\f$
         *  \param[out]  rhsContribInner
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
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at
         *  \param[out]  forceContrib
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

        //! pair of velocity and pressure space for the pass
        typedef Dune::DiscreteFunctionSpacePair< Dune::DiscreteFunctionSpacePairTraits< DiscreteVelocityFunctionSpaceType, DiscretePressureFunctionSpaceType > >
            DiscreteFunctionSpacePair;

        //! discrete function type for the pressure
        typedef Dune::AdaptiveDiscreteFunction< DiscretePressureFunctionSpaceType >
            DiscretePressureFunctionType;

        //! function type for the analytical force
        typedef Dune::Function< VelocityFunctionSpaceType, AnalyticalForceImp >
            AnalyticalForceType;

        //! function type for the analytical boundary values
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
 *  \brief  A default implementation of a discrete stokes model.
 *
 *          Implements the fluxes needed for the ldg method
 *          (see DiscreteStokesModelInterface).
 *          The fluxes \f$\hat{u}_{\sigma}\f$, \f$\hat{\sigma}\f$,
 *          \f$\hat{p}\f$ and \f$\hat{u}_{p}\f$ are implemented as proposed in
 *          B. Cockburn, G. Kanschat, D. Schötzau, C. Schwab: <EM>Local
 *          Discontinuous Galerkin Methodsfor the Stokes System</EM> (2000).\n\n
 *          To use this model, a user has to implement the analytical force
 *          \f$f\f$ and the dirichlet data \f$g_{D}\f$ as a Dune::Function
 *          (only the method evaluate( arg, ret ) is needed) and specify the
 *          types of this functions as template arguments for the traits class
 *          DiscreteStokesModelDefaultTraits.\n\n
 *          <b>Notation:</b> Given simplices \f$T_{+}\f$ and
 *          \f$T_{-}\f$ and a face \f$\varepsilon\f$ between them, the values
 *          of a function \f$u\f$ on the face \f$\varepsilon\f$ are denoted by \f$u^{+}\f$,
 *          if seen from \f$T_{+}\f$ and \f$u^{-}\f$, if seen from \f$T_{-}\f$.
 *          The outer normals of \f$T_{+,-}\f$ in a given point on
 *          the face \f$\varepsilon\f$ are denoted by \f$n_{+,-}\f$,
 *          accordingly.\n
 *          We define\n
 *          the <b>mean values</b>\n
 *          - \f$\{\{p\}\}\in R\f$ for a \f$p\in R\f$ as
 *              \f[
 *                  \{\{p\}\}:=\frac{1}{2}\left(p^{+}+p^{-}\right),
 *              \f]
 *          - \f$\{\{u\}\}\in R^{d}\f$ for a \f$u\in R^{d}\f$ as
 *              \f[
 *                  \{\{u\}\}:=\frac{1}{2}\left(u^{+}+u^{-}\right),
 *              \f]
 *          - \f$\{\{\sigma\}\}\in R^{d\times d}\f$ for a \f$\sigma\in R^{d\times d}\f$ as
 *              \f[
 *                  \{\{\sigma\}\}:=\frac{1}{2}\left(\sigma^{+}+\sigma^{-}\right)
 *              \f]
 *
 *          and the <b>jumps</b>\n
 *          - \f$\left[\left[p\right]\right]\in R^{d}\f$ for a \f$p\in R\f$ as
 *              \f[
 *                  \left[\left[p\right]\right]:=p^{+}n^{+}+p^{-}n^{-},
 *              \f]
 *          - \f$\left[\left[u\right]\right]\in R\f$ for a \f$u\in R^{d}\f$ as
 *              \f[
 *                  \left[\left[u\right]\right]:=u^{+}\cdot n^{+}+u^{-}\cdot n^{-},
 *              \f]
 *          - \f$\underline{\left[\left[u\right]\right]}\in R^{d\times d}\f$ for a \f$u\in R^{d}\f$ as
 *              \f[
 *                  \underline{\left[\left[u\right]\right]}:=u^{+}\otimes n^{+}+u^{-}\otimes n^{-},
 *              \f]
 *          - \f$\left[\left[\sigma\right]\right]\in R^{d}\f$ for a \f$\sigma\in R^{d\times d}\f$ as
 *              \f[
 *                  \left[\left[\sigma\right]\right]:=\sigma^{+}\cdot n^{+}+\sigma^{-}\cdot n^{-}.
 *              \f]
 *
 *          We also denote by \f$\mathcal{E}_{D}\f$ the set of those faces
 *          \f$\varepsilon\f$ that lie on the boundary \f$\partial\Omega\f$ and
 *          by \f$\mathcal{E}_{I}\f$ those which are inside \f$\Omega\f$.\n
 *          For a detailed definition of this notation see B. Cockburn, G. Kanschat, D. Schötzau, C. Schwab: <EM>Local
 *          Discontinuous Galerkin Methodsfor the Stokes System</EM> (2000), again.\n
 *          <b>Attention:</b> For reasons of simplicity the assumtion \f$n^{-}=-1\cdot n^{+}\f$ is used.
 *          This may be not true for nonconforming grids.\n\n
 *          With this notation at hand the fluxes can de described as
 *          - \f$\hat{u}_{\sigma}:\Omega\rightarrow R^{d}\f$ for an inner face (see velocitySigmaFlux() )
 *              \f[
 *                  \hat{u}_{\sigma}(u):=\{\{u\}\}+\underline{\left[\left[u\right]\right]}\cdot C_{12}\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{u}_{\sigma}:\Omega\rightarrow R^{d}\f$ for a boundary face (see velocitySigmaBoundaryFlux() )
 *              \f[
 *                  \hat{u}_{\sigma}(u):=g_{D}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{\sigma}:\Omega\rightarrow R^{d\times d}\f$ for an inner face (see sigmaFlux() )
 *              \f[
 *                  \hat{\sigma}(u,\sigma):=\{\{\sigma\}\}-C_{11}\underline{\left[\left[u\right]\right]}+\left[\left[\sigma\right]\right]\otimes C_{12}\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{\sigma}:\Omega\rightarrow R^{d\times d}\f$ for a boundary face (see sigmaBoundaryFlux() )
 *              \f[
 *                  \hat{\sigma}(u,\sigma):=\sigma^{+}-C_{11}\left(u^{+}-g_{D}\right)\otimes n^{+}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{p}:\Omega\rightarrow R\f$ for an inner face (see pressureFlux() )
 *              \f[
 *                  \hat{p}(p):=\{\{p\}\}-D_{12}\left[\left[p\right]\right]\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{p}:\Omega\rightarrow R\f$ for a boundary face (see pressureBoundaryFlux() )
 *              \f[
 *                  \hat{p}(p):=p^{+}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{u}_{p}:\Omega\rightarrow R^{d}\f$ for an inner face (see velocityPressureFlux() )
 *              \f[
 *                  \hat{u}_{p}(u,p):=\{\{u\}\}+D_{11}\left[\left[p\right]\right]+D_{12}\left[\left[u\right]\right]\quad\quad\varepsilon\in\mathcal{E}_{I}
 *              \f]
 *
 *          and
 *          - \f$\hat{u}_{p}:\Omega\rightarrow R^{d}\f$ for a boundary face (see velocityPressureBoundaryFlux() )
 *              \f[
 *                  \hat{u}_{p}(u,p):=g_{D}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *
 *          where \f$C_{11},\;\;D_{11}\in R\f$ are the stability coefficients
 *          and \f$C_{12},\;\;D_{12}\in R^{d}\f$ are the coefficients
 *          concerning efficiency and accuracy.
 **/
template < class DiscreteStokesModelDefaultTraitsImp >
class DiscreteStokesModelDefault : public DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
{
    public:

        //! interface class
        typedef DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
            BaseType;

        //! iterator over intersections
        typedef typename BaseType::IntersectionIteratorType
            IntersectionIteratorType;

        //! Vector type of the velocity's discrete function space's range
        typedef typename BaseType::VelocityRangeType
            VelocityRangeType;

        //! Matrix type of the sigma's discrete function space's range
        typedef typename BaseType::SigmaRangeType
            SigmaRangeType;

        //! Vector type of the pressure's discrete function space's range
        typedef typename BaseType::PressureRangeType
            PressureRangeType;

        //! type of analytical force (usually Dune::Function)
        typedef typename BaseType::AnalyticalForceType
            AnalyticalForceType;

        //! type of analytical dirichlet data (usually Dune::Function)
        typedef typename BaseType::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

        /**
         *  \brief  constructor
         *
         *  sets the coefficients
         *  \param[in]  C_11
         *          \f$C_{11}\in R\f$
         *  \param[in]  C_12
         *          \f$C_{12}\in R^{d}\f$
         *  \param[in]  D_11
         *          \f$D_{11}\in R\f$
         *  \param[in]  D_12
         *          \f$D_{12}\in R^{d}\f$
         *  \param[in]  force
         *          analytical force
         *  \param[in]  dirichletData
         *          analytical dirichlet data
         **/
        DiscreteStokesModelDefault( const double C_11,
                                    const VelocityRangeType& C_12,
                                    const double D_11,
                                    const VelocityRangeType& D_12,
                                    const AnalyticalForceType& force,
                                    const AnalyticalDirichletDataType& dirichletData )
            : C_11_( C_11 ),
            C_12_( C_12 ),
            D_11_( D_11 ),
            D_12_( D_12 ),
            force_( force ),
            dirichletData_( dirichletData )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~DiscreteStokesModelDefault()
        {}

        /**
         *  \brief  returns true
         *
         *  since problem has a \f$\hat{u}_{\sigma}\f$ contribution
         *  \return true
         **/
        bool hasVelocitySigmaFlux() const
        {
            return true;
        }

        /**
         *  \brief  returns true
         *
         *  since problem has a \f$\hat{u}_{p}\f$ contribution
         *  \return true
         **/
        bool hasVelocityPressureFlux() const
        {
            return true;
        }

        /**
         *  \brief  returns true
         *
         *  since problem has a \f$\hat{p}\f$ contribution
         *  \return true
         **/
        bool hasPressureFlux() const
        {
            return true;
        }

        /**
         *  \brief  returns true
         *
         *  since problem has a \f$\hat{\sigma}\f$ contribution
         *  \return true
         **/
        bool hasSigmaFlux() const
        {
            return true;
        }

        /**
         *  \brief  returns true
         *
         *  since problem has a \f$f\f$ contribution
         *  \return true
         **/
        bool hasForce() const
        {
            return true;
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{\sigma}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{u}_{\sigma}^{U}:=\{\{u\}\}+\underline{\left[\left[u\right]\right]}\cdot C_{12}\f$
         *
         *  and
         *  - \f$\hat{u}_{\sigma}^{RHS}:=0\f$.
         *
         *  \attention  Assumption: \f$n_{-}=-1\cdot n_{+}\f$
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *              faceiterator
         *  \param[in]  time
         *              global time
         *  \param[in]  x
         *              point to evaluate at (on the face)
         *  \param[in]  uInner
         *              value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *              value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[out] uContribInner
         *              \f$\hat{u}_{\sigma}^{U}\f$ (seen from the inside)
         *  \param[out] uContribOuter
         *              \f$\hat{u}_{\sigma}^{U}\f$ (seen from the outside)
         *  \param[out] rhsContribInner
         *              \f$\hat{u}_{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param[out] rhsContribOuter
         *              \f$\hat{u}_{\sigma}^{RHS}\f$ (seen from the outside)
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
         *  \brief  implementation of \f$\hat{u}_{\sigma}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{u}_{\sigma}^{U}:=0\f$
         *
         *  and
         *  - \f$\hat{u}_{\sigma}^{RHS}:=g_{D}\f$.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{\sigma}^{U}\f$
         *  \param[out]  rhsContribInner
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
            // some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType global = it->intersectionSelfLocal().global( x );

            // contribution to u vector ( from inside entity )
            uContribInner = 0.0;

            // contribution to rhs ( from inside entity )
            dirichletData_.evaluate( global,  rhsContribInner );
        }

        /**
         *  \brief  implementation of \f$\hat{u}_{p}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{u}_{p}^{U}:=\{\{u\}\}+D_{12}\cdot\underline{\left[\left[u\right]\right]}\f$,
         *  - \f$\hat{u}_{p}^{P}:=D_{11}\left[\left[p\right]\right]\f$
         *
         *  and
         *  - \f$\hat{u}_{p}^{RHS}:=0\f$.
         *
         *  \attention  Assumption: \f$n_{-}=-1\cdot n_{+}\f$
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the inside)
         *  \param[out]  uContribOuter
         *          \f$\hat{u}_{p}^{U}\f$ (seen from the outside)
         *  \param[out]  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the inside)
         *  \param[out]  pContribOuter
         *          \f$\hat{u}_{p}^{P}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{u}_{p}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \brief  implementation of \f$\hat{u}_{p}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{u}_{p}^{U}:=0\f$,
         *  - \f$\hat{u}_{p}^{P}:=0\f$
         *
         *  and
         *  - \f$\hat{u}_{p}^{RHS}:=g_{D}\f$.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$
         *  \param[out]  uContribInner
         *          \f$\hat{u}_{p}^{U}\f$
         *  \param[out]  pContribInner
         *          \f$\hat{u}_{p}^{P}\f$
         *  \param[out]  rhsContribInner
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
         *  \brief  implementation of \f$\hat{p}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{p}^{P}:=\{\{p\}\}-D_{12}\cdot\left[\left[p\right]\right]\f$
         *
         *  and
         *  - \f$\hat{p}^{RHS}:=0\f$.
         *
         *  \attention  Assumption: \f$n_{-}=-1\cdot n_{+}\f$
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  pOuter
         *          value of \f$p\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  pContribInner
         *          \f$\hat{p}^{P}\f$ (seen from the inside)
         *  \param[out]  pContribOuter
         *          \f$\hat{p}^{P}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{p}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \brief  implementation of \f$\hat{p}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{p}^{P}:=p^{+}\f$
         *
         *  and
         *  - \f$\hat{p}^{RHS}:=0\f$.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  pInner
         *          value of \f$p\f$ in \f$x\f$ (seen from the inside)
         *  \param[out]  pContribInner
         *          \f$\hat{p}^{P}\f$ (seen from the inside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{p}^{RHS}\f$ (seen from the inside)
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
         *  \brief  implementation of \f$\hat{\sigma}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{\sigma}^{\sigma}:=\{\{\sigma\}\}+\left[\left[\sigma\right]\right]\otimes C_{12}\f$,
         *  - \f$\hat{\sigma}^{U}:=-C_{11}\underline{\left[\left[u\right]\right]}\f$
         *
         *  and
         *  - \f$\hat{\sigma}^{RHS}:=0\f$.
         *
         *  \attention  Assumption: \f$n_{-}=-1\cdot n_{+}\f$
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  uOuter
         *          value of \f$u\f$ in \f$x\f$ (seen from the outside)
         *  \param[in]  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  sigmaOuter
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the outside)
         *  \param[out]  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the inside)
         *  \param[out]  sigmaContribOuter
         *          \f$\hat{\sigma}^{\sigma}\f$ (seen from the outside)
         *  \param[out]  uContribInner
         *          \f$\hat{\sigma}^{U}\f$ (seen from the inside)
         *  \param[out]  uContribOuter
         *          \f$\hat{\sigma}^{U}\f$ (seen from the outside)
         *  \param[out]  rhsContribInner
         *          \f$\hat{\sigma}^{RHS}\f$ (seen from the inside)
         *  \param[out]  rhsContribOuter
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
         *  \brief  implementation of \f$\hat{\sigma}\f$
         *
         *  Under the assumption of linearity (see DiscreteStokesModelInterface)
         *  this flux returns
         *  - \f$\hat{\sigma}^{\sigma}:=\sigma^{+}\f$,
         *  - \f$\hat{\sigma}^{U}:=-C_{11}u^{+}\otimes n^{+}\f$
         *
         *  and
         *  - \f$\hat{\sigma}^{RHS}:=C_{11}g_{D}\otimes n^{+}\f$.
         *
         *  \tparam FaceDomainType
         *          domain type on given face
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face)
         *  \param[in]  uInner
         *          value of \f$u\f$ in \f$x\f$ (seen from the inside)
         *  \param[in]  sigmaInner
         *          value of \f$\sigma\f$ in \f$x\f$ (seen from the inside)
         *  \param[out]  sigmaContribInner
         *          \f$\hat{\sigma}^{\sigma}\f$
         *  \param[out]  uContribInner
         *          \f$\hat{\sigma}^{U}\f$
         *  \param[out]  rhsContribInner
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
         *  \brief  implementation of \f$f\f$.
         *
         *          Calls the implementation of the user.
         *
         *  \tparam DomainType
         *          domain type in entity
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at
         *  \param[out]  forceContrib
         *          value of \f$f\f$ in \f$x\f$
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
