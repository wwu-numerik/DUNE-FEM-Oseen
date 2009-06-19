/**
 *  \file   discretestokesmodelinterface.hh
 *  \brief  contains a class DiscreteStokesModelInterface
 *          and a class DiscreteStokesModelDefault with traits class
 *          DiscreteStokesModelDefaultTraits.
 **/
#ifndef DUNE_DISCRESTOKESTEMODELINTERFACE_HH
#define DUNE_DISCRESTOKESTEMODELINTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include <dune/stokes/discretefunctionspacepair.hh>

#ifndef NLOG
    #include <dune/stuff/printing.hh>
    #include <dune/stuff/logging.hh>
#endif

#include <dune/stuff/misc.hh>

namespace Dune
{

/**
 *  \brief  Interface class for stokes problem definition in the LDG context.
 *
 *          A discretization of the stokes problem using LDG methods consist of
 *          several dependent classes.
 *          The Dune::DiscreteStokesFunctionSpaceWrapper is a wrapper class to
 *          combine the discrete function spaces of the velocity and the
 *          pressure into one discrete function space. The
 *          Dune::DiscreteStokesFunctionWrapper accordingly is a wrapper class
 *          to combine the velocity and the pressure itself into one discrete
 *          function.
 *          The assembling of the system matrices and right hand sides takes
 *          place in the Dune::StokesPass. This class also solves for the
 *          velocity and the pressure (which have to be wrapped inside a
 *          Dune::DiscreteStokesFunctionWrapper).
 *          The problem-dependent data (force terms, boundary data, fluxes,
 *          viscosity) are implemented in a discrete model, which should be
 *          derived from the Dune::DiscreteStokesModelInterface interface class
 *          (see Dune::DiscreteStokesModelDefault for example).
 *
 *          A discrete model implementation of the user should be derived from
 *          this interface class to be compatible with the Dune::StokesPass.
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
 *          \f$V\f$, \f$Q\f$ are discrete function spaces.
 *          (For a detailed description see B. Cockburn, G. Kanschat,
 *          D. Schötzau, C. Schwab: <EM>Local Discontinuous Galerkin Methods
 *          for the Stokes System</EM> (2000)).\n
 *          The fluxes \f$\hat{u}_{\sigma}\f$, \f$\hat{\sigma}\f$,
 *          \f$\hat{p}\f$, \f$\hat{u}_{p}\f$ in the corresponding surface
 *          integrals are implemented in the methods velocitySigmaFlux(),
 *          sigmaFlux(), pressureFlux(), velocityPressureFlux().
 *          If the face in consideration is on the boundary of \f$\Omega\f$, the
 *          computation is done by velocitySigmaBoundaryFlux(),
 *          sigmaBoundaryFlux(), pressureBoundaryFlux() and
 *          velocityPressureBoundaryFlux().\n
 *          The analytical fluxes \f$\hat{\sigma}\f$ and \f$\hat{u}_{p}\f$
 *          each depend on several functions. This is realized by overloading
 *          sigmaFlux() and velocityPressureFlux(). All the boundary fluxes need
 *          to be overloaded as well.\n
 *          Since the fluxes are linear, they can be decomposed as follows:\n
 *          - \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$\n
 *            \f$
 *                  \hat{u}_{\sigma}(u) = \hat{u}_{\sigma}^{U^{+}}
 *                  + \hat{u}_{\sigma}^{U^{-}}
 *                  + \hat{u}_{\sigma}^{RHS}
 *             \f$
 *            - for inner faces
 *              - \f$\hat{u}_{\sigma}^{U^{+}}\f$ and
 *                \f$\hat{u}_{\sigma}^{U^{-}}\f$ are implemented in
 *                velocitySigmaFlux() (const IntersectionIteratorType& it, const
 *                double time, const FaceDomainType& x, const Side side, const
 *                VelocityRangeType& u, VelocityRangeType& uReturn), where
 *                <b>side</b> determines, whether \f$\hat{u}_{\sigma}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{\sigma}^{U^{-}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{u}_{\sigma}^{U^{+}}\f$ is implemented in
 *                velocitySigmaBoundaryFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x,
 *                const VelocityRangeType& <b>u</b>, VelocityRangeType&
 *                <b>uReturn</b> )
 *              - \f$\hat{u}_{\sigma}^{RHS}\f$ is implemented in
 *                velocitySigmaBoundaryFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x,
 *                VelocityRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{u}_{p}(u,p):\Omega\rightarrow R^{d}\f$\n
 *            \f$
 *                  \hat{u}_{p}(u,p) = \hat{u}_{p}^{U^{+}}
 *                  + \hat{u}_{p}^{U^{-}}
 *                  + \hat{u}_{p}^{P^{+}}
 *                  + \hat{u}_{p}^{P^{-}}
 *                  + \hat{u}_{p}^{RHS}
 *            \f$
 *            - for inner faces
 *              - \f$\hat{u}_{p}^{U^{+}}\f$ and \f$\hat{u}_{p}^{U^{-}}\f$ are
 *                implemented in velocityPressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const VelocityRangeType&
 *                <b>u</b>, VelocityRangeType& <b>uReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{u}_{p}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{p}^{U^{-}}\f$
 *                (side=outside) is returned
 *              - \f$\hat{u}_{p}^{P^{+}}\f$ and \f$\hat{u}_{p}^{P^{+}}\f$ are
 *                implemented by velocityPressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const PressureRangeType&
 *                <b>p</b>, VelocityRangeType& <b>pReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{u}_{p}^{P^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{p}^{P^{+}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{u}_{p}^{U^{+}}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const VelocityRangeType& <b>u</b>,
 *                VelocityRangeType& <b>uReturn</b> )
 *              - \f$\hat{u}_{p}^{P^{+}}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const PressureRangeType& <b>p</b>,
 *                VelocityRangeType& <b>pReturn</b> )
 *              - \f$\hat{u}_{p}^{RHS}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, VelocityRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{p}(p):\Omega\rightarrow R\f$\n
 *            \f$
 *                  \hat{p}(p) = \hat{p}^{P^{+}}
 *                  + \hat{p}^{P^{-}}
 *                  + \hat{p}^{RHS}
 *            \f$
 *            - for inner faces
 *              - \f$\hat{p}^{P^{+}}\f$ and \f$\hat{p}^{P^{-}}\f$ are
 *                implemented in pressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const PressureRangeType&
 *                p, PressureRangeType& pReturn ), where
 *                <b>side</b> determines, whether \f$\hat{p}^{P^{+}}\f$
 *                (side=inside) or \f$\hat{p}^{P^{-}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{p}^{P^{+}}\f$ is implemented in pressureBoundaryFlux() (
 *                const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const PressureRangeType& <b>p</b>,
 *                PressureRangeType& <b>pReturn</b> )
 *              - \f$\hat{p}^{RHS}\f$ is implemented in pressureBoundaryFlux() (
 *                const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, PressureRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$\n
 *            \f$
 *                  \hat{\sigma}(u,\sigma) = \hat{\sigma}^{U^{+}}
 *                  + \hat{\sigma}^{U^{-}}
 *                  + \hat{\sigma}^{\sigma^{+}}
 *                  + \hat{\sigma}^{\sigma^{+}}
 *                  + \hat{\sigma}^{RHS}
 *            \f$
 *            - for inner faces
 *              - \f$\hat{\sigma}^{U^{+}}\f$ and \f$\hat{\sigma}^{U^{-}}\f$ are
 *                implemented in sigmaFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x, const Side
 *                side, const VelocityRangeType& <b>u</b>, SigmaRangeType&
 *                <b>uReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{\sigma}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{\sigma}^{U^{-}}\f$
 *                (side=outside) is returned
 *              - \f$\hat{\sigma}^{\sigma^{+}}\f$ and
 *                \f$\hat{\sigma}^{\sigma^{-}}\f$ are implemented in sigmaFlux()
 *                ( const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const SigmaRangeType&
 *                <b>sigma</b>, SigmaRangeType& <b>sigmaReturn</b> ), where
 *                <b>side</b> determines, whether
 *                \f$\hat{\sigma}^{\sigma^{+}}\f$ (side=inside) or
 *                \f$\hat{\sigma}^{\sigma^{-}}\f$ (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{\sigma}^{U^{+}}\f$ is implemented in
 *                sigmaBoundaryFlux() ( const IntersectionIteratorType& it,
 *                const double time, const FaceDomainType& x, const
 *                VelocityRangeType& <b>u</b>, SigmaRangeType& <b>uReturn</b> )
 *              - \f$\hat{\sigma}^{\sigma^{+}}\f$ is implemented in
 *                sigmaBoundaryFlux() ( const IntersectionIteratorType& it,
 *                const double time, const FaceDomainType& x, const
 *                SigmaRangeType& <b>sigma</b>, SigmaRangeType&
 *                <b>sigmaReturn</b> )
 *              - \f$\hat{\sigma}^{RHS}\f$ is implemented in sigmaBoundaryFlux()
 *                ( const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, SigmaRangeType& <b>rhsReturn</b> )
 *
 *  \tparam DiscreteStokesModelTraits
 *          traits class defined by the user, should provide all types needed
 *          by this interface
 **/
template < class DiscreteStokesModelTraits >
class DiscreteStokesModelInterface
{
    public:

        //! traits class defined by the user
        typedef DiscreteStokesModelTraits
            Traits;

    private:

        //! implementation type for CRTP
        typedef typename Traits::DiscreteModelType
            DiscreteModelType;

    public:

        //! volume quadrature type to be used in pass
        typedef typename Traits::VolumeQuadratureType
            VolumeQuadratureType;

        //! face quadrature type to be used in pass
        typedef typename Traits::FaceQuadratureType
            FaceQuadratureType;

        //! discrete function space wrapper type
        typedef typename Traits::DiscreteStokesFunctionSpaceWrapperType
            DiscreteStokesFunctionSpaceWrapperType;

        //! discrete function wrapper type
        typedef typename Traits::DiscreteStokesFunctionWrapperType
            DiscreteStokesFunctionWrapperType;

    private:

        //! discrete function type for the velocity
        typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! discrete function space type for the velocity
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! function space type for the velocity
        typedef typename DiscreteVelocityFunctionSpaceType::FunctionSpaceType
            VelocityFunctionSpaceType;

    public:

        //! discrete function type for sigma
        typedef typename Traits::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

    private:

        //! discrete function space type for sigma
        typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
            DiscreteSigmaFunctionSpaceType;

        //! function space type for sigma
        typedef typename DiscreteSigmaFunctionSpaceType::FunctionSpaceType
            SigmaFunctionSpaceType;

        //! discrete function type for the pressure
        typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! discrete function space type for the pressure
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! function space type for the pressure
        typedef typename DiscretePressureFunctionSpaceType::FunctionSpaceType
            PressureFunctionSpaceType;

    protected:

        //! function type for analytical force
        typedef typename Traits::AnalyticalForceType
            AnalyticalForceType;

        //! function type for analytical dirichlet data
        typedef typename Traits::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

    private:

        //! coordinate type (world coordinates)
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            DomainType;

    protected:

        //! vector type of the velocity's discrete function space's range
        typedef typename DiscreteVelocityFunctionSpaceType::RangeType
            VelocityRangeType;

        //! matrix type of sigmas' discrete functions space's range
        typedef typename DiscreteSigmaFunctionSpaceType::RangeType
            SigmaRangeType;

        //! vector type of the pressure's discrete function space's range
        typedef typename DiscretePressureFunctionSpaceType::RangeType
            PressureRangeType;

    private:

        //! type of GridPart
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;

        //! type of the grid
        typedef typename GridPartType::GridType
            GridType;

    protected:

        //! intersection iterator of the grid
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;

    private:

        //! element type (codim 0 entity) of the grid
        typedef typename GridType::template Codim<0>::Entity
            EntityType;

    public:

        //! polynomial order for the discrete sigma function space
        static const int sigmaSpaceOrder = Traits::sigmaSpaceOrder;
        //! polynomial order for the discrete velocity function space
        static const int velocitySpaceOrder = Traits::velocitySpaceOrder;
        //! polynomial order for the discrete pressure function space
        static const int pressureSpaceOrder = Traits::pressureSpaceOrder;

        /**
         *  \brief  constructor
         *
         *  does nothing
         **/
        DiscreteStokesModelInterface()
        {}

        /**
         *  \brief  destructor
         *
         *  does nothing
         **/
        ~DiscreteStokesModelInterface()
        {}

        /**
         *  \brief  contains "inside" and "outside", needed for the fluxes
         **/
        enum Side{
            inside,
            outside
        };

        /**
         *  \brief  Returns true if problem has a flux contribution of type
         *          \f$\hat{u}_{\sigma}\f$.
         *          Calls the implementation of the derived class.
         *  \attention  If you let this method return true, make sure to
         *              implement <b>all</b> versions of velocitySigmaFlux() and
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
         *              implement <b>all</b> versions of velocityPressureFlux()
         *              and velocityPressureBoundaryFlux() as well.
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
         *              implement <b>all</b> versions of pressureFlux() and
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
         *              implement <b>all</b> versions of sigmaFlux() and
         *              sigmaBoundaryFlux() as well.
         **/
        bool hasSigmaFlux() const
        {
            return asImp().hasSigmaFlux();
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
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{U^{+}}\f$ and
         *          \f$\hat{u}_{\sigma}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{u}_{\sigma}^{U^{+}}\f$, if (side==inside) or
         *          \f$\hat{u}_{\sigma}^{U^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const Side side,
                                const VelocityRangeType& u,
                                VelocityRangeType& uReturn) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaFlux(  it,
                                            time,
                                            x,
                                            side,
                                            u,
                                            uReturn) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{U^{+}}\f$ for a face
         *          on the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{u}_{\sigma}^{U^{+}}\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const VelocityRangeType& u,
                                        VelocityRangeType& uReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaBoundaryFlux(  it,
                                                    time,
                                                    x,
                                                    u,
                                                    uReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{RHS}\f$ for a face
         *          on the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  rhsReturn
         *          \f$\hat{u}_{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        VelocityRangeType& rhsReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocitySigmaBoundaryFlux(  it,
                                                    time,
                                                    x,
                                                    rhsReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{U^{+}}\f$ and
         *          \f$\hat{u}_{p}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{u}_{p}^{U^{+}}\f$, if (side==inside) or
         *          \f$\hat{u}_{p}^{U^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const Side side,
                                    const VelocityRangeType& u,
                                    VelocityRangeType& uReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureFlux(   it,
                                                time,
                                                x,
                                                side,
                                                u,
                                                uReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{P^{+}}\f$ and
         *          \f$\hat{u}_{p}^{P^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>p</i> and
         *          <i>pReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  p
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  pReturn
         *          \f$\hat{u}_{p}^{P^{+}}\f$, if (side==inside) or
         *          \f$\hat{u}_{p}^{P^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const Side side,
                                    const PressureRangeType& p,
                                    VelocityRangeType& pReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureFlux(   it,
                                                time,
                                                x,
                                                side,
                                                p,
                                                pReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{U^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{u}_{p}^{U^{+}}\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const VelocityRangeType& u,
                                    VelocityRangeType& uReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureBoundaryFlux(   it,
                                                        time,
                                                        x,
                                                        u,
                                                        uReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{P^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  p
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out] pReturn
         *          \f$\hat{u}_{p}^{P^{+}}\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& p,
                                    VelocityRangeType& pReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureBoundaryFlux(   it,
                                                        time,
                                                        x,
                                                        p,
                                                        pReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{RHS}\f$ for a face on
         *          the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out] rhsReturn
         *          \f$\hat{u}_{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    VelocityRangeType& rhsReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().velocityPressureBoundaryFlux(   it,
                                                        time,
                                                        x,
                                                        rhsReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{P^{+}}\f$ and
         *          \f$\hat{p}^{P^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>p</i> and
         *          <i>pReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  p
         *          value of \f$p\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  pReturn
         *          \f$\hat{p}^{P^{+}}\f$, if (side==inside) or
         *          \f$\hat{p}^{P^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void pressureFlux(  const IntersectionIteratorType& it,
                            const double time,
                            const FaceDomainType& x,
                            const Side side,
                            const PressureRangeType& p,
                            PressureRangeType& pReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureFlux(   it,
                                        time,
                                        x,
                                        side,
                                        p,
                                        pReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{P^{+}}\f$ for a face on the
         *          boundary of \f$\Omega\f$. Calls the implementation of the
         *          derived class.
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  p
         *          value of \f$p\f$ in \f$x\f$
         *  \param[out]  pReturn
         *          \f$\hat{p}^{P^{+}}\f$
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& p,
                                    PressureRangeType& pReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureBoundaryFlux(   it,
                                                time,
                                                x,
                                                p,
                                                pReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{RHS}\f$ for a face on the
         *          boundary of \f$\Omega\f$. Calls the implementation of the
         *          derived class.
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  rhsReturn
         *          \f$\hat{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    PressureRangeType& rhsReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().pressureBoundaryFlux(   it,
                                                time,
                                                x,
                                                rhsReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{U^{+}}\f$ and
         *          \f$\hat{\sigma}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{U^{+}}\f$, if (side==inside) or
         *          \f$\hat{\sigma}^{U^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const Side side,
                        const VelocityRangeType& u,
                        SigmaRangeType& uReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaFlux(  it,
                                    time,
                                    x,
                                    side,
                                    u,
                                    uReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{\sigma^{+}}\f$ and
         *          \f$\hat{\sigma}^{\sigma^{-}}\f$ for a face inside
         *          \f$\Omega\f$. Calls the implementation of the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>sigma</i> and
         *          <i>sigmaReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  sigma
         *          value of \f$\sigma\f$ in \f$x\f$, interpreted once as seen
         *          from the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  sigmaReturn
         *          \f$\hat{\sigma}^{\sigma^{+}}\f$, if (side==inside) or
         *          \f$\hat{\sigma}^{\sigma^{-}}\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const Side side,
                        const SigmaRangeType& sigma,
                        SigmaRangeType& sigmaReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaFlux(  it,
                                    time,
                                    x,
                                    side,
                                    sigma,
                                    sigmaReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{U^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{U^{+}}\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& u,
                                SigmaRangeType& uReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaBoundaryFlux(  it,
                                            time,
                                            x,
                                            u,
                                            uReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{\sigma^{+}}\f$ for a face
         *          on the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  sigma
         *          value of \f$\sigma\f$ in \f$x\f$
         *  \param[out]  sigmaReturn
         *          \f$\hat{\sigma}^{\sigma^{+}}\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const SigmaRangeType& sigma,
                                SigmaRangeType& sigmaReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaBoundaryFlux(  it,
                                            time,
                                            x,
                                            sigma,
                                            sigmaReturn ) );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{RHS}\f$ for a face
         *          on the boundary of \f$\Omega\f$. Calls the implementation of
         *          the derived class.
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                SigmaRangeType& rhsReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().sigmaBoundaryFlux(  it,
                                            time,
                                            x,
                                            rhsReturn ) );
        }

        /**
         *  \brief  Implementation of \f$f\f$.
         *          Calls the implementation of the derived class.
         *
         *  \tparam DomainType
         *          domain type in entity (codim 0)
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (in world coordinates)
         *  \param[out]  forceReturn
         *          value of \f$f\f$ in \f$x\f$
         **/
        template < class DomainType >
        void force( const double time,
                    const DomainType& x,
                    VelocityRangeType& forceReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().force(  time,
                                x,
                                forceReturn ) );
        }

        /**
         *  \brief  Returns the viscosity \f$\mu\f$ of the fluid.
         *          Calls the implementation of the derived class.
         *  \return \f$\mu\f$
         **/
        double viscosity() const
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().viscosity() );
            return asImp().viscosity();
        }

    private:
        //! for CRTP trick
        DiscreteModelType& asImp()
        {
            return static_cast< DiscreteModelType& >(*this);
        }

        //! for CRTP trick
        const DiscreteModelType& asImp() const
        {
            return static_cast< const DiscreteModelType& >(*this);
        }

};

// forward declaration
template < class DiscreteStokesModelDefaultTraitsImp >
class DiscreteStokesModelDefault;

/**
 *  \brief  Traits class for DiscreteStokesModelDefault
 **/
template < class GridPartImp, class AnalyticalForceImp, class AnalyticalDirichletDataImp, int gridDim, int polOrder >
class DiscreteStokesModelDefaultTraits
{
    public:

        //! for CRTP trick
        typedef DiscreteStokesModelDefault < DiscreteStokesModelDefaultTraits >
            DiscreteModelType;

        //! we use caching quadratures for the entities
        typedef Dune::CachingQuadrature< GridPartImp, 0 >
            VolumeQuadratureType;

        //! we use caching quadratures for the faces
        typedef Dune::CachingQuadrature< GridPartImp, 1 >
            FaceQuadratureType;

        //! polynomial order for the discrete sigma function space
        static const int sigmaSpaceOrder = polOrder;
        //! polynomial order for the discrete velocity function space
        static const int velocitySpaceOrder = polOrder;
        //! polynomial order for the discrete pressure function space
        static const int pressureSpaceOrder = polOrder;

    private:

        //! function space type for the velocity
        typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
            VelocityFunctionSpaceType;

        //! discrete function space type for the velocity
        typedef Dune::DiscontinuousGalerkinSpace<   VelocityFunctionSpaceType,
                                                    GridPartImp,
                                                    velocitySpaceOrder >
            DiscreteVelocityFunctionSpaceType;

        //! function space type for the pressure
        typedef Dune::FunctionSpace< double, double, gridDim, 1 >
            PressureFunctionSpaceType;

        //! discrete function space type for the pressure
        typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
                                                    GridPartImp,
                                                    pressureSpaceOrder >
            DiscretePressureFunctionSpaceType;

    public:

        //! discrete function space wrapper type
        typedef Dune::DiscreteStokesFunctionSpaceWrapper< Dune::DiscreteStokesFunctionSpaceWrapperTraits<
                    DiscreteVelocityFunctionSpaceType,
                    DiscretePressureFunctionSpaceType > >
            DiscreteStokesFunctionSpaceWrapperType;

    private:

        //! discrete function type for the velocity
        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
            DiscreteVelocityFunctionType;

        //! discrete function type for the pressure
        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
            DiscretePressureFunctionType;

    public:

        //! discrete function wrapper type
        typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
                    DiscreteStokesFunctionSpaceWrapperType,
                    DiscreteVelocityFunctionType,
                    DiscretePressureFunctionType > >
            DiscreteStokesFunctionWrapperType;

    private:

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
                                                    sigmaSpaceOrder >
            DiscreteSigmaFunctionSpaceType;

    public:

        //! discrete function type for sigma
        typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
            DiscreteSigmaFunctionType;

        //! function type for the analytical force
        typedef Dune::Function< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType, AnalyticalForceImp >
            AnalyticalForceType;

        //! function type for the analytical dirichlet data
        typedef Dune::Function< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType, AnalyticalDirichletDataImp >
            AnalyticalDirichletDataType;

        /**
         *  \name   types needed for the pass
         *  \{
         **/
        //! return type of the pass
        typedef DiscreteStokesFunctionWrapperType
            DestinationType;
        /**
         *  \}
         **/

};


/**
 *  \brief  A default implementation of a discrete stokes model.
 *
 *          Implements the fluxes needed for the LDG method
 *          (see Dune::DiscreteStokesModelInterface for details).\n
 *          The fluxes \f$\hat{u}_{\sigma}\f$, \f$\hat{\sigma}\f$,
 *          \f$\hat{p}\f$ and \f$\hat{u}_{p}\f$ are implemented as proposed in
 *          B. Cockburn, G. Kanschat, D. Schötzau, C. Schwab: <EM>Local
 *          Discontinuous Galerkin Methodsfor the Stokes System</EM> (2000).\n\n
 *          To use this model, a user has to implement the analytical force
 *          \f$f\f$ and the dirichlet data \f$g_{D}\f$ as a Dune::Function
 *          (only the method evaluate( arg, ret ) is needed) and specify the
 *          types of this functions as template arguments for the traits class
 *          DiscreteStokesModelDefaultTraits.\n
 *
 *          <b>Notation:</b> Given simplices \f$T_{+}\f$ and
 *          \f$T_{-}\f$ and a face \f$\varepsilon\f$ between them, the values
 *          of a function \f$u\f$ on the face \f$\varepsilon\f$ are denoted by \f$u^{+}\f$,
 *          if seen from \f$T_{+}\f$ and \f$u^{-}\f$, if seen from \f$T_{-}\f$.
 *          The outer normals of \f$T_{+,-}\f$ in a given point on
 *          the face \f$\varepsilon\f$ are denoted by \f$n_{+,-}\f$,
 *          accordingly.\n
 *
 *          We define the <b>mean values</b>\n
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
 *
 *          <b>Attention:</b> For reasons of simplicity the assumtion \f$n^{-}=-1\cdot n^{+}\f$ is used.
 *          This may be not true for nonconforming grids.\n
 *
 *          With this notation at hand the fluxes can de described as
 *          - \f$\hat{u}_{\sigma}:\Omega\rightarrow R^{d}\f$ for an inner face
 *              \f[
 *                  \hat{u}_{\sigma}(u):=\{\{u\}\}-\underline{\left[\left[u\right]\right]}\cdot C_{12}\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{u}_{\sigma}:\Omega\rightarrow R^{d}\f$ for a boundary face
 *              \f[
 *                  \hat{u}_{\sigma}(u):=g_{D}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{\sigma}:\Omega\rightarrow R^{d\times d}\f$ for an inner face
 *              \f[
 *                  \hat{\sigma}(u,\sigma):=\{\{\sigma\}\}-C_{11}\underline{\left[\left[u\right]\right]}-\left[\left[\sigma\right]\right]\otimes C_{12}\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{\sigma}:\Omega\rightarrow R^{d\times d}\f$ for a boundary face
 *              \f[
 *                  \hat{\sigma}(u,\sigma):=\sigma^{+}-C_{11}\left(u^{+}-g_{D}\right)\otimes n^{+}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{p}:\Omega\rightarrow R\f$ for an inner face
 *              \f[
 *                  \hat{p}(p):=\{\{p\}\}-D_{12}\left[\left[p\right]\right]\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{p}:\Omega\rightarrow R\f$ for a boundary face
 *              \f[
 *                  \hat{p}(p):=p^{+}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *          - \f$\hat{u}_{p}:\Omega\rightarrow R^{d}\f$ for an inner face
 *              \f[
 *                  \hat{u}_{p}(u,p):=\{\{u\}\}+D_{11}\left[\left[p\right]\right]+D_{12}\left[\left[u\right]\right]\quad\quad\varepsilon\in\mathcal{E}_{I},
 *              \f]
 *          - \f$\hat{u}_{p}:\Omega\rightarrow R^{d}\f$ for a boundary face
 *              \f[
 *                  \hat{u}_{p}(u,p):=g_{D}\quad\quad\varepsilon\in\mathcal{E}_{D},
 *              \f]
 *
 *          where \f$C_{11},\;\;D_{11}\in R\f$ are the stability coefficients
 *          and \f$C_{12},\;\;D_{12}\in R^{d}\f$ are the coefficients
 *          concerning efficiency and accuracy.\n
 *
 *          These fluxes are then decomposed into several numerical fluxes (see
 *          Dune::DiscreteStokesModelInterface for details):\n
 *
 *          \f{tabular}{l||l}
 *              on $\mathcal{E}_{I}$ & on $\mathcal{E}_{I}$ \\
 *                  \hline\hline
 *              $\boldsymbol{\hat{u}_{\sigma}^{U^{+}}(u)} := \frac{1}{2} u + \left( u \otimes n^{+} \right) \cdot C_{12}$
 *                  & $\boldsymbol{\hat{u}_{\sigma}^{U^{+}}(u)} := 0$ \\
 *              $\boldsymbol{\hat{u}_{\sigma}^{U^{-}}(u)} := \frac{1}{2} u + \left( u \otimes n^{-} \right) \cdot C_{12}$
 *                  & $\boldsymbol{\hat{u}_{\sigma}^{RHS}} := g_{D}$ \\
 *                  \hline
 *              $\boldsymbol{\hat{u}_{p}^{U^{+}}(u)} := \frac{1}{2} u + D_{12} u \cdot n^{+}$
 *                  & $\boldsymbol{\hat{u}_{p}^{U^{+}}(u)} := 0$ \\
 *              $\boldsymbol{\hat{u}_{p}^{U^{-}}(u)} := \frac{1}{2} u + D_{12} u \cdot n^{-}$
 *                  & $\quad$ \\
 *              $\boldsymbol{\hat{u}_{p}^{P^{+}}(p)} := D_{11} p n^{+}$
 *                  & $\boldsymbol{\hat{u}_{p}^{P^{+}}(p)} := 0$ \\
 *              $\boldsymbol{\hat{u}_{p}^{P^{-}}(p)} := D_{11} p n^{-}$
 *                  & $\boldsymbol{\hat{u}_{p}^{RHS}} := g_{D}$\\
 *                  \hline
 *              $\boldsymbol{\hat{p}^{P^{+}}(p)} := \frac{1}{2} p - p D_{12} \cdot n^{+}$
 *                  & $\boldsymbol{\hat{p}^{P^{+}}(p)} := p$ \\
 *              $\boldsymbol{\hat{p}^{P^{-}}(p)} := \frac{1}{2} p - p D_{12} \cdot n^{-}$
 *                  & $\boldsymbol{\hat{p}^{RHS}} := 0$ \\
 *                  \hline
 *              $\boldsymbol{\hat{\sigma}^{U^{+}}(u)} := -C_{11} u \otimes n^{+}$
 *                  & $\boldsymbol{\hat{\sigma}^{U^{+}}(u)} := -C_{11} u \otimes n^{+}$ \\
 *              $\boldsymbol{\hat{\sigma}^{U^{-}}(u)} := -C_{11} u \otimes n^{-}$
 *                  & $\quad$ \\
 *              $\boldsymbol{\hat{\sigma}^{\sigma^{+}}(u)} := \frac{1}{2} \sigma - \left( \sigma \cdot n^{+} \right)$
 *                  & $\boldsymbol{\hat{\sigma}^{\sigma^{+}}(u)} := \sigma$ \\
 *              $\boldsymbol{\hat{\sigma}^{\sigma^{-}}(u)} := \frac{1}{2} \sigma - \left( \sigma \cdot n^{-} \right)$
 *                  & $\boldsymbol{\hat{\sigma}^{RHS}} := -C_{11} g_{D} \otimes n^{+}$
 *          \f}
 *
 *          The implementation is as follows:\n
 *
 *          - \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$\n
 *            - for inner faces
 *              - \f$\hat{u}_{\sigma}^{U^{+}}\f$ and
 *                \f$\hat{u}_{\sigma}^{U^{-}}\f$ are implemented in
 *                velocitySigmaFlux() (const IntersectionIteratorType& it, const
 *                double time, const FaceDomainType& x, const Side side, const
 *                VelocityRangeType& u, VelocityRangeType& uReturn), where
 *                <b>side</b> determines, whether \f$\hat{u}_{\sigma}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{\sigma}^{U^{-}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{u}_{\sigma}^{U^{+}}\f$ is implemented in
 *                velocitySigmaBoundaryFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x,
 *                const VelocityRangeType& <b>u</b>, VelocityRangeType&
 *                <b>uReturn</b> )
 *              - \f$\hat{u}_{\sigma}^{RHS}\f$ is implemented in
 *                velocitySigmaBoundaryFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x,
 *                VelocityRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{u}_{p}(u,p):\Omega\rightarrow R^{d}\f$\n
 *            - for inner faces
 *              - \f$\hat{u}_{p}^{U^{+}}\f$ and \f$\hat{u}_{p}^{U^{-}}\f$ are
 *                implemented in velocityPressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const VelocityRangeType&
 *                <b>u</b>, VelocityRangeType& <b>uReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{u}_{p}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{p}^{U^{-}}\f$
 *                (side=outside) is returned
 *              - \f$\hat{u}_{p}^{P^{+}}\f$ and \f$\hat{u}_{p}^{P^{+}}\f$ are
 *                implemented by velocityPressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const PressureRangeType&
 *                <b>p</b>, VelocityRangeType& <b>pReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{u}_{p}^{P^{+}}\f$
 *                (side=inside) or \f$\hat{u}_{p}^{P^{+}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{u}_{p}^{U^{+}}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const VelocityRangeType& <b>u</b>,
 *                VelocityRangeType& <b>uReturn</b> )
 *              - \f$\hat{u}_{p}^{P^{+}}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const PressureRangeType& <b>p</b>,
 *                VelocityRangeType& <b>pReturn</b> )
 *              - \f$\hat{u}_{p}^{RHS}\f$ is implemented in
 *                velocityPressureBoundaryFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, VelocityRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{p}(p):\Omega\rightarrow R\f$\n
 *            - for inner faces
 *              - \f$\hat{p}^{P^{+}}\f$ and \f$\hat{p}^{P^{-}}\f$ are
 *                implemented in pressureFlux() ( const
 *                IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const PressureRangeType&
 *                p, PressureRangeType& pReturn ), where
 *                <b>side</b> determines, whether \f$\hat{p}^{P^{+}}\f$
 *                (side=inside) or \f$\hat{p}^{P^{-}}\f$
 *                (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{p}^{P^{+}}\f$ is implemented in pressureBoundaryFlux() (
 *                const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const PressureRangeType& <b>p</b>,
 *                PressureRangeType& <b>pReturn</b> )
 *              - \f$\hat{p}^{RHS}\f$ is implemented in pressureBoundaryFlux() (
 *                const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, PressureRangeType& <b>rhsReturn</b> )
 *          - \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$\n
 *            - for inner faces
 *              - \f$\hat{\sigma}^{U^{+}}\f$ and \f$\hat{\sigma}^{U^{-}}\f$ are
 *                implemented in sigmaFlux() ( const IntersectionIteratorType&
 *                it, const double time, const FaceDomainType& x, const Side
 *                side, const VelocityRangeType& <b>u</b>, SigmaRangeType&
 *                <b>uReturn</b> ), where
 *                <b>side</b> determines, whether \f$\hat{\sigma}^{U^{+}}\f$
 *                (side=inside) or \f$\hat{\sigma}^{U^{-}}\f$
 *                (side=outside) is returned
 *              - \f$\hat{\sigma}^{\sigma^{+}}\f$ and
 *                \f$\hat{\sigma}^{\sigma^{-}}\f$ are implemented in sigmaFlux()
 *                ( const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, const Side side, const SigmaRangeType&
 *                <b>sigma</b>, SigmaRangeType& <b>sigmaReturn</b> ), where
 *                <b>side</b> determines, whether
 *                \f$\hat{\sigma}^{\sigma^{+}}\f$ (side=inside) or
 *                \f$\hat{\sigma}^{\sigma^{-}}\f$ (side=outside) is returned
 *            - for faces on the boundary of \f$\Omega\f$
 *              - \f$\hat{\sigma}^{U^{+}}\f$ is implemented in
 *                sigmaBoundaryFlux() ( const IntersectionIteratorType& it,
 *                const double time, const FaceDomainType& x, const
 *                VelocityRangeType& <b>u</b>, SigmaRangeType& <b>uReturn</b> )
 *              - \f$\hat{\sigma}^{\sigma^{+}}\f$ is implemented in
 *                sigmaBoundaryFlux() ( const IntersectionIteratorType& it,
 *                const double time, const FaceDomainType& x, const
 *                SigmaRangeType& <b>sigma</b>, SigmaRangeType&
 *                <b>sigmaReturn</b> )
 *              - \f$\hat{\sigma}^{RHS}\f$ is implemented in sigmaBoundaryFlux()
 *                ( const IntersectionIteratorType& it, const double time, const
 *                FaceDomainType& x, SigmaRangeType& <b>rhsReturn</b> )
 **/
template < class DiscreteStokesModelDefaultTraitsImp >
class DiscreteStokesModelDefault : public DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
{
    private:

        //! interface class
        typedef DiscreteStokesModelInterface< DiscreteStokesModelDefaultTraitsImp >
            BaseType;

        //! \copydoc Dune::DiscreteStokesModelInterface::IntersectionIteratorType
        typedef typename BaseType::IntersectionIteratorType
            IntersectionIteratorType;

    public:

        //! \copydoc Dune::DiscreteStokesModelInterface::VolumeQuadratureType
        typedef typename BaseType::VolumeQuadratureType
            VolumeQuadratureType;

        //! \copydoc Dune::DiscreteStokesModelInterface::FaceQuadratureType
        typedef typename BaseType::FaceQuadratureType
            FaceQuadratureType;

        //! \copydoc Dune::DiscreteStokesModelInterface::DiscreteStokesFunctionSpaceWrapperType
        typedef typename BaseType::DiscreteStokesFunctionSpaceWrapperType
            DiscreteStokesFunctionSpaceWrapperType;

        //! \copydoc Dune::DiscreteStokesModelInterface::DiscreteStokesFunctionSpaceWrapperType
        typedef typename BaseType::DiscreteStokesFunctionWrapperType
            DiscreteStokesFunctionWrapperType;

        //! \copydoc Dune::DiscreteStokesModelInterface::DiscreteSigmaFunctionType
        typedef typename BaseType::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

        //! \copydoc Dune::DiscreteStokesModelInterface::sigmaSpaceOrder
        static const int sigmaSpaceOrder
            = BaseType::sigmaSpaceOrder;
        //! \copydoc Dune::DiscreteStokesModelInterface::velocitySpaceOrder
        static const int velocitySpaceOrder
            = BaseType::velocitySpaceOrder;
        //! \copydoc Dune::DiscreteStokesModelInterface::pressureSpaceOrder
        static const int pressureSpaceOrder
            = BaseType::pressureSpaceOrder;

    private:

        //! codim 0 entity pointer type
        typedef typename IntersectionIteratorType::EntityPointer
            EntityPointer;

        //! codim 0 entity type
        typedef typename IntersectionIteratorType::Entity
            EntityType;

        //! geometry type of codim 0 entity
        typedef typename EntityType::Geometry
            EntityGeometryType;

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

    public:

        //! \copydoc Dune::DiscreteStokesModelInterface::Side
        typedef enum BaseType::Side
            Side;

        /**
         *  \brief  constructor
         *
         *  sets the coefficients and analytical data
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
         *  \param[in]  viscosity
         *          viscosity of the fluid
         **/
        DiscreteStokesModelDefault( const int C_11,
                                    const int C_12,
                                    const int D_11,
                                    const int D_12,
                                    const AnalyticalForceType& force,
                                    const AnalyticalDirichletDataType& dirichletData,
                                    const double viscosity )
            : viscosity_( viscosity ),
            C_11_( C_11 ),
            C_12_( C_12 ),
            D_11_( D_11 ),
            D_12_( D_12 ),
            force_( force ),
            dirichletData_( dirichletData )
        {}

        /**
         *  \brief  destructor
         *
         *  does nothing
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
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{U^{+}}\f$ and
         *          \f$\hat{u}_{\sigma}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{u}_{\sigma}^{U^{+}}(u) = \frac{1}{2} u + \left( u \otimes n^{+} \right) \cdot C_{12}\f$
         *          - \f$\hat{u}_{\sigma}^{U^{-}}(u) = \frac{1}{2} u + \left( u \otimes n^{-} \right) \cdot C_{12}\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{u}_{\sigma}^{U^{+}}(u)\f$, if (side==inside) or
         *          \f$\hat{u}_{\sigma}^{U^{-}}(u)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocitySigmaFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const Side side,
                                const VelocityRangeType& u,
                                VelocityRangeType& uReturn) const
        {
            // some preparations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType C_12( 1.0 );
            C_12 *= getStabScalar( x, it , C_12_ );

            // contribution to u vector ( from inside entity )
            if ( side == BaseType::inside ) {
                SigmaRangeType u_plus_tensor_n_plus = dyadicProduct( u, outerNormal );
                VelocityRangeType u_plus_tensor_n_plus_times_c_12( 0.0 );
                u_plus_tensor_n_plus.mv( C_12, u_plus_tensor_n_plus_times_c_12 );
                uReturn = u;
                uReturn *= 0.5;
                uReturn += u_plus_tensor_n_plus_times_c_12;
            }
            // contribution to u vector ( from outside entity )
            else if ( side == BaseType::outside ) {
                //some preparations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // calculation
                SigmaRangeType u_minus_tensor_n_minus = dyadicProduct( u, innerNormal );
                VelocityRangeType u_minus_tensor_n_minus_times_c_12( 0.0 );
                u_minus_tensor_n_minus.mv( C_12, u_minus_tensor_n_minus_times_c_12 );
                uReturn = u;
                uReturn *= 0.5;
                uReturn += u_minus_tensor_n_minus_times_c_12;
            }
            else {
                LOGIC_ERROR
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{U^{+}}\f$ for a face
         *          on the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{u}_{\sigma}^{U^{+}}(u) = 0\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{u}_{\sigma}^{U^{+}}(u)\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        const VelocityRangeType& u,
                                        VelocityRangeType& uReturn ) const
        {
            // contribution to u vector ( from inside entity )
            uReturn = 0.0;
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{\sigma}^{RHS}\f$ for a face
         *          on the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{u}_{\sigma}^{RHS} = g_{D}\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{\sigma}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  rhsReturn
         *          \f$\hat{u}_{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocitySigmaBoundaryFlux( const IntersectionIteratorType& it,
                                        const double time,
                                        const FaceDomainType& x,
                                        VelocityRangeType& rhsReturn ) const
        {
            // some preparations
            const EntityPointer entityPt = it.inside();
            const EntityType& entity = *entityPt;
            const EntityGeometryType& geometry = entity.geometry();
            const VelocityRangeType xIntersectionGlobal = it->intersectionSelfLocal().global( x );
            const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
            // contribution to rhs ( from inside entity )
            rhsReturn = 0.0;
            dirichletData_.evaluate( xWorld,  rhsReturn );
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{U^{+}}\f$ and
         *          \f$\hat{u}_{p}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{u}_{p}^{U^{+}}(u) = \frac{1}{2} u + D_{12} u \cdot n^{+}\f$
         *          - \f$\hat{u}_{p}^{U^{-}}(u) = \frac{1}{2} u + D_{12} u \cdot n^{-}\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{u}_{p}^{U^{+}}(u)\f$, if (side==inside) or
         *          \f$\hat{u}_{p}^{U^{-}}(u)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const Side side,
                                    const VelocityRangeType& u,
                                    VelocityRangeType& uReturn ) const
        {
            //some preparations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType D_12( 1.0 );
            D_12 *= getStabScalar( x, it , D_12_ );

            // contribution to u vector ( from inside entity )
            if ( side == BaseType::inside ) {
                const double u_plus_times_n_plus = u * outerNormal;
                VelocityRangeType d_12_times_u_plus_times_n_plus = D_12;
                d_12_times_u_plus_times_n_plus *= u_plus_times_n_plus;
                uReturn = u;
                uReturn *= 0.5;
                uReturn += d_12_times_u_plus_times_n_plus;
            }
            // contribution to u vector ( from outside entity )
            else if ( side == BaseType::outside ) {
                // some preparations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // calculations
                const double u_minus_times_n_minus = u * innerNormal;
                VelocityRangeType d_12_times_u_minus_times_n_minus = D_12;
                d_12_times_u_minus_times_n_minus *= u_minus_times_n_minus;
                uReturn = u;
                uReturn *= 0.5;
                uReturn += d_12_times_u_minus_times_n_minus;
            }
            else {
                LOGIC_ERROR
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{P^{+}}\f$ and
         *          \f$\hat{u}_{p}^{P^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{u}_{p}^{P^{+}}(p) = D_{11} p n^{+}\f$
         *          - \f$\hat{u}_{p}^{P^{-}}(p) = D_{11} p n^{-}\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>p</i> and
         *          <i>pReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  p
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  pReturn
         *          \f$\hat{u}_{p}^{P^{+}}(p)\f$, if (side==inside) or
         *          \f$\hat{u}_{p}^{P^{-}}(p)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void velocityPressureFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const Side side,
                                    const PressureRangeType& p,
                                    VelocityRangeType& pReturn ) const
        {
            //some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            double D_11( 1.0 );
            D_11 *= getStabScalar( x, it , D_11_ );

            // contribution to p vector ( from inside entity )
            if ( side == BaseType::inside ) {
                const double d_11_times_p_plus = D_11 * p;
                pReturn = outerNormal;
                pReturn *= d_11_times_p_plus;
            }
            // contribution to p vector ( from outside entity )
            else if ( side == BaseType::outside ) {
                //some preperations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // calculations
                const double d_11_times_p_minus = D_11 * p;
                pReturn = innerNormal;
                pReturn *= d_11_times_p_minus;
            }
            else {
                LOGIC_ERROR
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{U^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{u}_{p}^{U^{+}}(u) = 0\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{u}_{p}^{U^{+}}(u)\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const VelocityRangeType& u,
                                    VelocityRangeType& uReturn ) const
        {
            // contribution to u vector
            uReturn = 0.0;
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{P^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{u}_{p}^{P^{+}}(p) = 0\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  p
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out] pReturn
         *          \f$\hat{u}_{p}^{P^{+}}(p)\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& p,
                                    VelocityRangeType& pReturn ) const
        {
            // contribution to p vector
            pReturn = 0.0;
        }

        /**
         *  \brief  Implementation of \f$\hat{u}_{p}^{RHS}\f$ for a face on
         *          the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{u}_{p}^{RHS} = g_{D}\f$
         *
         *          For the docomposition of
         *          \f$\hat{u}_{p}(u):\Omega\rightarrow R^{d}\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out] rhsReturn
         *          \f$\hat{u}_{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void velocityPressureBoundaryFlux(
                                    const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    VelocityRangeType& rhsReturn ) const
        {
            //some preparations
            const EntityPointer entityPt = it.inside();
            const EntityType& entity = *entityPt;
            const EntityGeometryType& geometry = entity.geometry();
            const VelocityRangeType xIntersectionGlobal = it->intersectionSelfLocal().global( x );
            const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
            // contribution to rhs
            rhsReturn = 0.0;
            dirichletData_.evaluate( xWorld, rhsReturn );
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{P^{+}}\f$ and
         *          \f$\hat{p}^{P^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{p}^{P^{+}}(p) = \frac{1}{2} p - p D_{12} \cdot n^{+}\f$
         *          - \f$\hat{p}^{P^{-}}(p) = \frac{1}{2} p - p D_{12} \cdot n^{-}\f$
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>p</i> and
         *          <i>pReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  p
         *          value of \f$p\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  pReturn
         *          \f$\hat{p}^{P^{+}}(p)\f$, if (side==inside) or
         *          \f$\hat{p}^{P^{-}}(p)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void pressureFlux(  const IntersectionIteratorType& it,
                            const double time,
                            const FaceDomainType& x,
                            const Side side,
                            const PressureRangeType& p,
                            PressureRangeType& pReturn ) const
        {
            // some preperations
            VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType D_12( 1.0 );
            D_12 *= getStabScalar( x, it , D_12_ );

            // contribution to p vector ( from inside entity )
            if ( side == BaseType::inside ) {
                const double d_12_times_n_plus = D_12 * outerNormal;
                pReturn = 0.5 * p;
                pReturn -= p * d_12_times_n_plus;
            }
            // contribution to p vector ( from outside entity )
            else if ( side == BaseType::outside) {
                // some preperations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // claculations
                const double d_12_times_n_minus = D_12 * innerNormal;
                pReturn = 0.5 * p;
                pReturn -= p * d_12_times_n_minus;
            }
            else {
                LOGIC_ERROR
            }

        }

        template < class LocalPoint >
        double getStabScalar( const LocalPoint& x , const IntersectionIteratorType& it, const double param ) const
        {
            if ( param == -9 ) {
                return 0.0;
            }
            else if ( param == 0 ) {
                return 1.0;
            }
            else {
//                return std::pow( it.intersectionGlobal().integrationElement( x ), param );
                return std::pow( getLenghtOfIntersection( it ), param );
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{P^{+}}\f$ for a face on the
         *          boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{p}^{P^{+}}(p) = p\f$
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  p
         *          value of \f$p\f$ in \f$x\f$
         *  \param[out]  pReturn
         *          \f$\hat{p}^{P^{+}}(p)\f$
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    const PressureRangeType& p,
                                    PressureRangeType& pReturn ) const
        {
            // contribution to p vector
            pReturn = p;
        }

        /**
         *  \brief  Implementation of \f$\hat{p}^{RHS}\f$ for a face on the
         *          boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{p}^{RHS} = 0\f$
         *
         *          For the docomposition of
         *          \f$\hat{p}(p):\Omega\rightarrow R\f$, see the
         *          documentation of the Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  rhsReturn
         *          \f$\hat{p}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void pressureBoundaryFlux(  const IntersectionIteratorType& it,
                                    const double time,
                                    const FaceDomainType& x,
                                    PressureRangeType& rhsReturn ) const
        {
            // contribution to rhs
            rhsReturn = 0.0;
        }


        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{U^{+}}\f$ and
         *          \f$\hat{\sigma}^{U^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{\sigma}^{U^{+}}(u) = -C_{11} u \otimes n^{+}\f$
         *          - \f$\hat{\sigma}^{U^{-}}(u) = -C_{11} u \otimes n^{-}\f$
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>u</i> and
         *          <i>uReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$, interpreted once as seen from
         *          the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{U^{+}}(u)\f$, if (side==inside) or
         *          \f$\hat{\sigma}^{U^{-}}(u)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const Side side,
                        const VelocityRangeType& u,
                        SigmaRangeType& uReturn ) const
        {
            // some preparations
            const VelocityRangeType outerNormal = it.unitOuterNormal( x );
            double C_11( 1.0 );
            C_11 *= getStabScalar( x, it , C_11_ );

            // contribution to u vector ( from inside entity )
            if ( side == BaseType::inside ) {
                uReturn = dyadicProduct( u, outerNormal );
                uReturn *= ( -1.0 * C_11 );
            }
            // contribution to u vector ( from outside entity )
            else if ( side == BaseType::outside ) {
                // some preparations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // calculations
                uReturn = dyadicProduct( u, innerNormal );
                uReturn *= ( -1.0 * C_11 );
            }
            else {
                LOGIC_ERROR
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{\sigma^{+}}\f$ and
         *          \f$\hat{\sigma}^{\sigma^{-}}\f$ for a face inside
         *          \f$\Omega\f$.
         *
         *          Implements\n
         *          - \f$\hat{\sigma}^{\sigma^{+}}(\sigma) = \frac{1}{2} \sigma - \left( \sigma \cdot n^{+} \right) \otimes C_{12}\f$
         *          - \f$\hat{\sigma}^{\sigma^{+}}(\sigma) = \frac{1}{2} \sigma - \left( \sigma \cdot n^{-} \right) \otimes C_{12}\f$
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  side
         *          determines the interpretation of <i>sigma</i> and
         *          <i>sigmaReturn</i>\n
         *          legal values are DiscreteStokesModelInterface::inside and
         *          DiscreteStokesModelInterface::outside
         *  \param[in]  sigma
         *          value of \f$\sigma\f$ in \f$x\f$, interpreted once as seen
         *          from the inside entity (side==inside), once as seen from the
         *          neighbouring entity (side==outside)
         *  \param[out]  sigmaReturn
         *          \f$\hat{\sigma}^{\sigma^{+}}(\sigma)\f$, if (side==inside) or
         *          \f$\hat{\sigma}^{\sigma^{-}}(\sigma)\f$, if (side==outside)
         **/
        template < class FaceDomainType >
        void sigmaFlux( const IntersectionIteratorType& it,
                        const double time,
                        const FaceDomainType& x,
                        const Side side,
                        const SigmaRangeType& sigma,
                        SigmaRangeType& sigmaReturn ) const
        {
            // some preparations
            const VelocityRangeType outerNormal = it.unitOuterNormal( x );
            VelocityRangeType C_12( 1.0 );
            C_12 *= getStabScalar( x, it , C_12_ );

            // contribution to sigma vector ( from inside entity )
            if ( side == BaseType::inside ) {
                VelocityRangeType sigma_plus_times_n_plus( 0.0 );
                sigma.mv( outerNormal, sigma_plus_times_n_plus );
                const SigmaRangeType
                    sigma_plus_times_n_plus_times_c_12 =
                        dyadicProduct( sigma_plus_times_n_plus, C_12 );
                sigmaReturn = sigma;
                sigmaReturn *= 0.5;
                sigmaReturn -= sigma_plus_times_n_plus_times_c_12;
            }
            // contribution to sigma vector ( from outside entity )
            else if ( side == BaseType::outside ) {
                // some preparations
                VelocityRangeType innerNormal = outerNormal;
                innerNormal *= -1.0;
                // calculations
                VelocityRangeType sigma_minus_times_n_minus( 0.0 );
                sigma.mv( innerNormal, sigma_minus_times_n_minus );
                const SigmaRangeType
                    sigma_minus_times_n_minus_times_c_12 =
                        dyadicProduct( sigma_minus_times_n_minus, C_12 );
                sigmaReturn = sigma;
                sigmaReturn *= 0.5;
                sigmaReturn -= sigma_minus_times_n_minus_times_c_12;
            }
            else {
                LOGIC_ERROR
            }
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{U^{+}}\f$ for a face on
         *          the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{\sigma}^{U^{+}}(u) = -C_{11} u \otimes n^{+}\f$
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  u
         *          value of \f$u\f$ in \f$x\f$
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{U^{+}}(u)\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const VelocityRangeType& u,
                                SigmaRangeType& uReturn ) const
        {
            // some preparations
            const VelocityRangeType outerNormal = it.unitOuterNormal( x );
            double C_11( 1.0 );
            C_11 *= getStabScalar( x, it , C_11_ );

            // contribution to u vector
            uReturn = dyadicProduct( u, outerNormal );
            uReturn *= ( -1.0 * C_11 );
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{\sigma^{+}}\f$ for a face
         *          on the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{\sigma}^{\sigma^{+}}(\sigma) = \sigma\f$
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[in]  sigma
         *          value of \f$\sigma\f$ in \f$x\f$
         *  \param[out]  sigmaReturn
         *          \f$\hat{\sigma}^{\sigma^{+}}(\sigma)\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                const SigmaRangeType& sigma,
                                SigmaRangeType& sigmaReturn ) const
        {
            // contribution to sigma vector
            sigmaReturn = sigma;
        }

        /**
         *  \brief  Implementation of \f$\hat{\sigma}^{RHS}\f$ for a face
         *          on the boundary of \f$\Omega\f$.
         *
         *          Implements \f$\hat{\sigma}^{RHS} = C_{11} g_{D} \otimes n^{+}\f$
         *
         *          For the docomposition of
         *          \f$\hat{\sigma}(u,\sigma):\Omega\rightarrow R^{d\times d}\f$
         *          , see the documentation of the
         *          Dune::DiscreteStokesModelInterface.
         *
         *  \tparam FaceDomainType
         *          domain type on given face (codim 1 coordinates)
         *  \param[in]  it
         *          faceiterator
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (on the face) in face coordiantes
         *          (codim 1) (eg. as returned by Dune::CachingPointList::localPoint())
         *  \param[out]  uReturn
         *          \f$\hat{\sigma}^{RHS}\f$
         **/
        template < class FaceDomainType >
        void sigmaBoundaryFlux( const IntersectionIteratorType& it,
                                const double time,
                                const FaceDomainType& x,
                                SigmaRangeType& rhsReturn ) const
        {
            // some preparations
            const EntityPointer entityPt = it.inside();
            const EntityType& entity = *entityPt;
            const EntityGeometryType& geometry = entity.geometry();
            const VelocityRangeType xIntersectionGlobal = it->intersectionSelfLocal().global( x );
            const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
            const VelocityRangeType outerNormal = it.unitOuterNormal( x );
            double C_11( 1.0 );
            C_11 *= getStabScalar( x, it , C_11_ );

            // contribution to rhs
            VelocityRangeType gD( 0.0 );
            dirichletData_.evaluate( xWorld, gD );
            rhsReturn = dyadicProduct( gD, outerNormal );
            rhsReturn *= C_11;
        }

        /**
         *  \brief  Implementation of \f$f\f$.
         *
         *          Evaluates the analytical force given by the constructor
         *
         *  \tparam DomainType
         *          domain type in entity (codim 0)
         *  \param[in]  time
         *          global time
         *  \param[in]  x
         *          point to evaluate at (in world coordinates)
         *  \param[out]  forceReturn
         *          value of \f$f\f$ in \f$x\f$
         **/
        template < class DomainType >
        void force( const double time,
                    const DomainType& x,
                    VelocityRangeType& forceReturn ) const
        {
            force_.evaluate( x, forceReturn );
        }

        /**
         *  \brief  Returns the viscosity \f$\mu\f$ of the fluid.
         *  \return \f$\mu\f$
         **/
        double viscosity() const
        {
            return viscosity_;
        }

    private:

        const double viscosity_;
        const int C_11_, C_12_, D_11_, D_12_;
        const AnalyticalForceType& force_;
        const AnalyticalDirichletDataType& dirichletData_;

        /**
         *  \brief  dyadic product
         *
         *          Implements \f$\left(arg_{1} \otimes arg_{2}\right)_{i,j}:={arg_{1}}_{i} {arg_{2}}_{j}\f$
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
        }

        /**
         *  \brief  calculates length of given intersection in world coordinates
         *  \tparam IntersectionIteratorType
         *          IntersectionIteratorType
         *  \param[in]  intIt
         *          intersection
         *  \return length of intersection
         **/
        template < class IntersectionIteratorType >
        double getLenghtOfIntersection( const IntersectionIteratorType& intIt ) const
        {
            typedef typename IntersectionIteratorType::Geometry
                IntersectionGeometryType;
            const IntersectionGeometryType& intersectionGeoemtry = intIt.intersectionGlobal();
            assert( intersectionGeoemtry.corners() == 2 );
            typedef typename IntersectionIteratorType::ctype
                ctype;
            const int dimworld = IntersectionIteratorType::dimensionworld;
            typedef Dune::FieldVector< ctype, dimworld >
                DomainType;
            const DomainType cornerOne = intersectionGeoemtry[0];
            const DomainType cornerTwo = intersectionGeoemtry[1];
            const DomainType difference = cornerOne - cornerTwo;
            return difference.two_norm();
        }

};

}; // end of namespace Dune

#endif // end of discretestokesmodelinterface.hh
