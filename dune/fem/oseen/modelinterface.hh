/**
 *  \file   discretestokesmodelinterface.hh
 *  \brief  contains a class DiscreteOseenModelInterface
 *          and a class DiscreteOseenModelDefault with traits class
 *          DiscreteOseenModelDefaultTraits.
 **/
#ifndef DUNE_DISCRESTOKESTEMODELINTERFACE_HH
#define DUNE_DISCRESTOKESTEMODELINTERFACE_HH

#include <dune/common/fvector.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/fem/oseen/functionspacewrapper.hh>
#include <dune/fem/oseen/boundaryinfo.hh>
#include <dune/fem/oseen/stab_coeff.hh>

#ifndef NLOG
    #include <dune/stuff/common/print.hh>
    #include <dune/stuff/common/logging.hh>
#endif

#include <dune/stuff/common/misc.hh>

#include <algorithm>

// include this file after all other includes because some of them might undef
// the macros we want to use
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{
/**
 *  \brief  Interface class for stokes problem definition in the LDG context.
 *
 *          A discretization of the stokes problem using LDG methods consist of
 *          several dependent classes.
 *          The Dune::DiscreteOseenFunctionSpaceWrapper is a wrapper class to
 *          combine the discrete function spaces of the velocity and the
 *          pressure into one discrete function space. The
 *          Dune::DiscreteOseenFunctionWrapper accordingly is a wrapper class
 *          to combine the velocity and the pressure itself into one discrete
 *          function.
 *          The assembling of the system matrices and right hand sides takes
 *          place in the Dune::OseenLDGMethod. This class also solves for the
 *          velocity and the pressure (which have to be wrapped inside a
 *          Dune::DiscreteOseenFunctionWrapper).
 *          The problem-dependent data (force terms, boundary data, fluxes,
 *          viscosity) are implemented in a discrete model, which should be
 *          derived from the Dune::DiscreteOseenModelInterface interface class
 *          (see Dune::DiscreteOseenModelDefault for example).
 *
 *          A discrete model implementation of the user should be derived from
 *          this interface class to be compatible with the Dune::OseenLDGMethod.
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
 *  \tparam DiscreteOseenModelTraits
 *          traits class defined by the user, should provide all types needed
 *          by this interface
 **/
template < class DiscreteOseenModelTraits >
class DiscreteOseenModelInterface
{
    public:

        //! traits class defined by the user
        typedef DiscreteOseenModelTraits
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
        typedef typename Traits::DiscreteOseenFunctionSpaceWrapperType
            DiscreteOseenFunctionSpaceWrapperType;

        //! discrete function wrapper type
        typedef typename Traits::DiscreteOseenFunctionWrapperType
            DiscreteOseenFunctionWrapperType;

    private:

        //! discrete function type for the velocity
        typedef typename DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! discrete function space type for the velocity
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;   

    public:
		//! function space type for the velocity
		typedef typename DiscreteVelocityFunctionSpaceType::FunctionSpaceType
			VelocityFunctionSpaceType;

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
        typedef typename DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType
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
         *  \brief  contains "inside" and "outside", needed for the fluxes
         **/
        enum Side{
            inside,
            outside
        };

        const StabilizationCoefficients& getStabilizationCoefficients() const {
            return asImp().getStabilizationCoefficients();
        }

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

        bool isGeneralized() const
        {
            return asImp().isGeneralized();
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

        template < class IntersectionType, class DomainType >
        void dirichletData( const IntersectionType& intIt,
                            const double time,
                            const DomainType& x,
                            VelocityRangeType& dirichletDataReturn ) const
        {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().dirichletData(  intIt,
                                        time,
                                        x,
                                        dirichletDataReturn ) );
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

        /**
         *  \brief  constant for generalized stokes
         *
         *  \todo   doc
         **/
        double alpha() const
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().alpha() );
            return asImp().alpha();
        }

		//! ZAF
		double convection_scaling() const
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().convection_scaling() );
            return asImp().convection_scaling();
        }

		//! ZAF
		double pressure_gradient_scaling() const
        {
            CHECK_INTERFACE_IMPLEMENTATION( asImp().pressure_gradient_scaling() );
            return asImp().pressure_gradient_scaling();
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
template < class DiscreteOseenModelDefaultTraitsImp >
class DiscreteOseenModelDefault;

/**
 *  \brief  Traits class for DiscreteOseenModelDefault
 **/
template < class GridImp, template <class > class AnalyticalForceImp, class AnalyticalDirichletDataTraits,
            int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder,
           template< class, class, int, template<class> class BaseFunctionStorageImp = Dune::CachingStorage > class GalerkinSpaceImp = Dune::DiscontinuousGalerkinSpace>
class DiscreteOseenModelDefaultTraits
{
    public:
        //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
        typedef Dune::DGAdaptiveLeafGridPart< GridImp >
            GridPartType;
        //! for CRTP trick
        typedef DiscreteOseenModelDefault < DiscreteOseenModelDefaultTraits >
            DiscreteModelType;

        //! we use caching quadratures for the entities
        typedef Dune::CachingQuadrature< GridPartType, 0 >
            VolumeQuadratureType;

        //! we use caching quadratures for the faces
        typedef Dune::CachingQuadrature< GridPartType, 1 >
            FaceQuadratureType;

        //! polynomial order for the discrete sigma function space
        static const int sigmaSpaceOrder = sigmaOrder;
        //! polynomial order for the discrete velocity function space
        static const int velocitySpaceOrder = velocityOrder;
        //! polynomial order for the discrete pressure function space
        static const int pressureSpaceOrder = pressureOrder;

//    private:

        //! function space type for the velocity
        typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
            VelocityFunctionSpaceType;

        //! discrete function space type for the velocity
        typedef GalerkinSpaceImp<   VelocityFunctionSpaceType,
                                                    GridPartType,
                                                    velocitySpaceOrder >
            DiscreteVelocityFunctionSpaceType;

        //! function space type for the pressure
        typedef Dune::FunctionSpace< double, double, gridDim, 1 >
            PressureFunctionSpaceType;

        //! discrete function space type for the pressure
        typedef GalerkinSpaceImp<   PressureFunctionSpaceType,
                                                    GridPartType,
                                                    pressureSpaceOrder >
            DiscretePressureFunctionSpaceType;

    public:

        //! discrete function space wrapper type
        typedef Dune::DiscreteOseenFunctionSpaceWrapper< Dune::DiscreteOseenFunctionSpaceWrapperTraits<
                    DiscreteVelocityFunctionSpaceType,
                    DiscretePressureFunctionSpaceType > >
            DiscreteOseenFunctionSpaceWrapperType;

    private:
        //! function space type for sigma
        typedef Dune::MatrixFunctionSpace<  double,
                                            double,
                                            gridDim,
                                            gridDim,
                                            gridDim >
            SigmaFunctionSpaceType;

        //! discrete function space type for sigma
        typedef GalerkinSpaceImp<   SigmaFunctionSpaceType,
                                                    GridPartType,
                                                    sigmaSpaceOrder >
            DiscreteSigmaFunctionSpaceType;

        //! discrete function type for the velocity
        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
            DiscreteVelocityFunctionType;

        //! discrete function type for the pressure
        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
            DiscretePressureFunctionType;

    public:
        //! discrete function type for sigma
        typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
            DiscreteSigmaFunctionType;

        //! discrete function wrapper type
        typedef Dune::DiscreteOseenFunctionWrapper< Dune::DiscreteOseenFunctionWrapperTraits<
                    DiscreteOseenFunctionSpaceWrapperType,
                    DiscreteVelocityFunctionType,
                    DiscretePressureFunctionType > >
            DiscreteOseenFunctionWrapperType;

        //! function type for the analytical force
        typedef AnalyticalForceImp<VelocityFunctionSpaceType>
            AnalyticalForceType;

        //! function type for the analytical dirichlet data
        typedef typename AnalyticalDirichletDataTraits::template Implementation<VelocityFunctionSpaceType,GridPartType >
				AnalyticalDirichletDataTraitsImplementation;
		typedef typename AnalyticalDirichletDataTraitsImplementation::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

		typedef DiscreteVelocityFunctionType
			ExtraDataDiscreteFunctionType;
        /**
         *  \name   types needed for the pass
         *  \{
         **/
        //! return type of the pass
        typedef DiscreteOseenFunctionWrapperType
            DestinationType;
        /**
         *  \}
         **/

};

} // end of namespace Dune

#endif // end of discretestokesmodelinterface.hh

/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

