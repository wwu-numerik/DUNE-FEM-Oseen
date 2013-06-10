#ifndef MODELDEFAULT_HH
#define MODELDEFAULT_HH

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

#include "modelinterface.hh"

namespace Dune
{

/**
 *  \brief  A default implementation of a discrete stokes model.
 *
 *          Implements the fluxes needed for the LDG method
 *          (see Dune::DiscreteOseenModelInterface for details).\n
 *          The fluxes \f$\hat{u}_{\sigma}\f$, \f$\hat{\sigma}\f$,
 *          \f$\hat{p}\f$ and \f$\hat{u}_{p}\f$ are implemented as proposed in
 *          B. Cockburn, G. Kanschat, D. Schötzau, C. Schwab: <EM>Local
 *          Discontinuous Galerkin Methodsfor the Stokes System</EM> (2000).\n\n
 *          To use this model, a user has to implement the analytical force
 *          \f$f\f$ and the dirichlet data \f$g_{D}\f$ as a Dune::Fem::Function
 *          (only the method evaluate( arg, ret ) is needed) and specify the
 *          types of this functions as template arguments for the traits class
 *          DiscreteOseenModelDefaultTraits.\n
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
 *          Dune::DiscreteOseenModelInterface for details):\n
 *
 *          \f {tabular} {l||l}
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
template < class DiscreteOseenModelTraitsImp >
class DiscreteOseenModelDefault : public DiscreteOseenModelInterface< DiscreteOseenModelTraitsImp >
{
    private:

        //! interface class
        typedef DiscreteOseenModelInterface< DiscreteOseenModelTraitsImp >
            BaseType;

        //! \copydoc Dune::DiscreteOseenModelInterface::IntersectionIteratorType
        typedef typename BaseType::IntersectionIteratorType
            IntersectionIteratorType;

    public:

        //! \copydoc Dune::DiscreteOseenModelInterface::VolumeQuadratureType
        typedef typename BaseType::VolumeQuadratureType
            VolumeQuadratureType;

        //! \copydoc Dune::DiscreteOseenModelInterface::FaceQuadratureType
        typedef typename BaseType::FaceQuadratureType
            FaceQuadratureType;

        //! \copydoc Dune::DiscreteOseenModelInterface::DiscreteOseenFunctionSpaceWrapperType
        typedef typename BaseType::DiscreteOseenFunctionSpaceWrapperType
            DiscreteOseenFunctionSpaceWrapperType;

        //! \copydoc Dune::DiscreteOseenModelInterface::DiscreteOseenFunctionSpaceWrapperType
        typedef typename BaseType::DiscreteOseenFunctionWrapperType
            DiscreteOseenFunctionWrapperType;

        //! \copydoc Dune::DiscreteOseenModelInterface::DiscreteSigmaFunctionType
        typedef typename BaseType::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

        //! \copydoc Dune::DiscreteOseenModelInterface::sigmaSpaceOrder
        static const int sigmaSpaceOrder
            = BaseType::sigmaSpaceOrder;
        //! \copydoc Dune::DiscreteOseenModelInterface::velocitySpaceOrder
        static const int velocitySpaceOrder
            = BaseType::velocitySpaceOrder;
        //! \copydoc Dune::DiscreteOseenModelInterface::pressureSpaceOrder
        static const int pressureSpaceOrder
            = BaseType::pressureSpaceOrder;
        //! type of analytical force (usually Dune::Fem::Function)
        typedef typename BaseType::AnalyticalForceType
            AnalyticalForceType;
    private:


        //! Vector type of the velocity's discrete function space's range
        typedef typename BaseType::VelocityRangeType
            VelocityRangeType;

        //! Matrix type of the sigma's discrete function space's range
        typedef typename BaseType::SigmaRangeType
            SigmaRangeType;

        //! Vector type of the pressure's discrete function space's range
        typedef typename BaseType::PressureRangeType
            PressureRangeType;



        //! type of analytical dirichlet data (usually Dune::Fem::Function)
        typedef typename BaseType::AnalyticalDirichletDataType
            AnalyticalDirichletDataType;

    public:

        //! \copydoc Dune::DiscreteOseenModelInterface::Side
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
        DiscreteOseenModelDefault( const StabilizationCoefficients stab_coeff_in,
                                    const AnalyticalForceType force_in,
                                    const AnalyticalDirichletDataType dirichletData_in,
                                    const double viscosity_in = 1.0,
                                    const double alpha_in = 0.0,
                                   const double convection_scaling_in = 1.0,
                                   const double pressure_gradient_scaling_in = 1.0
                                   )
            : viscosity_( viscosity_in ),
            alpha_( alpha_in ),
            convection_scaling_( convection_scaling_in ),
            pressure_gradient_scaling_( pressure_gradient_scaling_in ),
            stabil_coeff_( stab_coeff_in ),
            force_( force_in ),
            dirichletData_( dirichletData_in )
        {
//            if ( !isGeneralized() ) {
//                if ( ( alpha_ < 0.0 ) || ( alpha_ > 0.0 ) ) {
//                    assert( !"isGeneralized() returns false, but alpha is not zero!" );
//                }
//            }
        }

        /**
         *  \brief  destructor
         *
         *  does nothing
         **/
    virtual ~DiscreteOseenModelDefault()
        {}

        const StabilizationCoefficients& getStabilizationCoefficients() const {
            return stabil_coeff_;
        }

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
        void force( const double /*time*/,
                    const DomainType& x,
                    VelocityRangeType& forceReturn ) const
        {
            force_.evaluate( x, forceReturn );
        }

        template < class IntersectionType, class DomainType >
        void dirichletData( const IntersectionType& intIt,
                            const double /*time*/,
                            const DomainType& x,
                            VelocityRangeType& dirichletDataReturn ) const
        {
            assert( ( !intIt.neighbor() && intIt.boundary() ) || !"this intersection does not lie on the boundary" );
            dirichletData_.evaluate( x, dirichletDataReturn, intIt );
        }

        /**
         *  \brief  Returns the viscosity \f$\mu\f$ of the fluid.
         *  \return \f$\mu\f$
         **/
        double viscosity() const
        {
            return viscosity_;
        }

        /**
         *  \brief  constant for generalized stokes
         *
         *  \todo   doc
         **/
        double alpha() const
        {
            return alpha_;
        }

        bool isGeneralized() const
        {
            return false;
        }

        //! ZAF
        double convection_scaling() const
        {
            return convection_scaling_;
        }

        //! ZAF
        double pressure_gradient_scaling() const
        {
            return pressure_gradient_scaling_;
        }

        const AnalyticalForceType& forceF() const
        {
            return force_;
        }

    private:

        const double viscosity_;
        const double alpha_;
        const double convection_scaling_;
        const double pressure_gradient_scaling_;
        StabilizationCoefficients stabil_coeff_;
        const AnalyticalForceType force_;
        const AnalyticalDirichletDataType dirichletData_;

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

        //! avoid code duplication by doing calculations for C_1X and D_1X here
        template < class LocalPoint >
        double getStabScalar( const LocalPoint& /*x*/ , const IntersectionIteratorType& it, const std::string coeffName ) const
        {
            const StabilizationCoefficients::PowerType  power   = stabil_coeff_.Power   ( coeffName );
            const StabilizationCoefficients::FactorType factor  = stabil_coeff_.Factor  ( coeffName );
            if ( power == StabilizationCoefficients::invalid_power ) {
                return 0.0;
            }
//                return std::pow( it.geometry().integrationElement( x ), param );
            return factor * std::pow( getLenghtOfIntersection( *it ), power );
        }
};

} // end of namespace Dune

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
#endif // MODELDEFAULT_HH
