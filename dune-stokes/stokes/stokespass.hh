/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

#ifndef NLOG // if we want logging, should be removed in the end
    #include "../src/stuff.hh"
    #include "../src/logging.hh"
#endif

#include <cmath>

namespace Dune
{
template <  class DiscreteModelImp,
            class PreviousPassImp,
            int PassID = 0 >
class StokesPass
    : public Pass < DiscreteModelImp, PreviousPassImp, PassID >
{
    public:

        //! base type
        typedef Pass < DiscreteModelImp, PreviousPassImp, PassID >
            BaseType;

        //! previous pass type
        typedef PreviousPassImp
            PreviousPassType;

        //! discrete model type
        typedef DiscreteModelImp
            DiscreteModelType;

        //! volume quadrature type
        typedef typename DiscreteModelType::VolumeQuadratureType
            VolumeQuadratureType;

        //! face quadrature type
        typedef typename DiscreteModelType::FaceQuadratureType
            FaceQuadratureType;

        //! type of discrete function space wrapper
        typedef typename DiscreteModelType::DiscreteStokesFunctionSpaceWrapperType
            DiscreteStokesFunctionSpaceWrapperType;

        //! discrete function wrapper type
        typedef typename DiscreteModelType::DiscreteStokesFunctionWrapperType
            DiscreteStokesFunctionWrapperType;

        //! discrete function type for the velocity
        typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! discrete function space type for the velocity
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! discrete function type for sigma
        typedef typename DiscreteModelType::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

        //! discrete function space type for sigma
        typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
            DiscreteSigmaFunctionSpaceType;

        //! discrete fucntion type for the pressure
        typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! discrete function space type for the pressure
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! Coordinate type (world coordinates)
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            WorldCoordinateType;

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

        //! Intersection iterator of the gridpart
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;

        //! entity iterator of the gridpart
        typedef typename GridPartType::template Codim< 0 >::IteratorType
            EntityIteratorType;

        //! type of the grid
        typedef typename GridPartType::GridType
            GridType;

        //! type of codim 0 entity
        typedef typename GridType::template Codim< 0 >::Entity
            EntityType;

        //! polynomial order for the discrete sigma function space
        static const int sigmaSpaceOrder = DiscreteModelType::sigmaSpaceOrder;
        //! polynomial order for the discrete velocity function space
        static const int velocitySpaceOrder = DiscreteModelType::velocitySpaceOrder;
        //! polynomial order for the discrete pressure function space
        static const int pressureSpaceOrder = DiscreteModelType::pressureSpaceOrder;


        /**
         *  \name typedefs for interface compliance
         *  \{
         **/
        typedef typename BaseType::DestinationType
            DestinationType;

        typedef typename BaseType::DomainType
            DomainType;

        typedef typename BaseType::RangeType
            RangeType;

        typedef typename BaseType::TotalArgumentType
            TotalArgumentType;
        /**
         *  \}
         **/


        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        StokesPass( PreviousPassType& prevPass,
                    DiscreteModelType& discreteModel,
                    GridPartType& gridPart,
                    DiscreteStokesFunctionSpaceWrapperType& spaceWrapper )
                : BaseType( prevPass ),
                discreteModel_( discreteModel ),
                gridPart_( gridPart ),
                spaceWrapper_( spaceWrapper ),
                velocitySpace_( spaceWrapper.discreteVelocitySpace() ),
                pressureSpace_( spaceWrapper.discretePressureSpace() ),
                sigmaSpace_( gridPart )
        {}

        /**
         *  \brief  empty constructor
         **/
        StokesPass()
        {}

        virtual void apply( const DomainType &arg, RangeType &dest) const
        {

            // viscosity
            double mu = discreteModel_.viscosity();

            // functions
            DiscreteVelocityFunctionType& velocity = dest.discreteVelocity();
            DiscretePressureFunctionType& pressure = dest.discretePressure();
            DiscreteSigmaFunctionType sigma( "sigma", sigmaSpace_ );

            // local functions
            typedef typename DiscreteVelocityFunctionType::LocalFunctionType
                LocalDiscreteVelocityFunctionType;
            typedef typename DiscretePressureFunctionType::LocalFunctionType
                LocalDiscretePressureFunctionType;
            typedef typename DiscreteSigmaFunctionType::LocalFunctionType
                LocalDiscreteSigmaFunctionType;

            // matrices
            // M\in R^{M\times M}
            typedef SparseRowMatrixObject< DiscreteSigmaFunctionSpaceType, DiscreteSigmaFunctionSpaceType >
                MmatrixType;
            MmatrixType Mmatrix( sigmaSpace_, sigmaSpace_ );
            Mmatrix.reserve();
            // W\in R^{M\times L}
            typedef SparseRowMatrixObject< DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
                WmatrixType;
            WmatrixType Wmatrix( sigmaSpace_, velocitySpace_ );
            Wmatrix.reserve();
            // X\in R^{L\times M}
            typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType >
                XmatrixType;
            XmatrixType Xmatrix( velocitySpace_, sigmaSpace_ );
            Xmatrix.reserve();
            // Y\in R^{L\times L}
            typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
                YmatrixType;
            YmatrixType Ymatrix( velocitySpace_, velocitySpace_ );
            Ymatrix.reserve();
            // Z\in R^{L\times K}
            typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscretePressureFunctionSpaceType >
                ZmatrixType;
            ZmatrixType Zmatrix( velocitySpace_, pressureSpace_ );
            Zmatrix.reserve();
            // E\in R^{K\times L}
            typedef SparseRowMatrixObject< DiscretePressureFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
                EmatrixType;
            EmatrixType Ematrix( pressureSpace_, velocitySpace_ );
            Ematrix.reserve();
            // R\in R^{K\times K}
            typedef SparseRowMatrixObject< DiscretePressureFunctionSpaceType, DiscretePressureFunctionSpaceType >
                RmatrixType;
            RmatrixType Rmatrix( pressureSpace_, pressureSpace_ );
            Rmatrix.reserve();

            // local matrices
            // M\in R^{M\times M}
            typedef typename MmatrixType::LocalMatrixType
                LocalMmatrixType;
            // W\in R^{M\times L}
            typedef typename WmatrixType::LocalMatrixType
                LocalWmatrixType;
            // X\in R^{L\times M}
            typedef typename XmatrixType::LocalMatrixType
                LocalXmatrixType;
            // Y\in R^{L\times L}
            typedef typename YmatrixType::LocalMatrixType
                LocalYmatrixType;
            // Z\in R^{L\times K}
            typedef typename ZmatrixType::LocalMatrixType
                LocalZmatrixType;
            // E\in R^{K\times L}
            typedef typename EmatrixType::LocalMatrixType
                LocalEmatrixType;
            // R\in R^{K\times K}
            typedef typename RmatrixType::LocalMatrixType
                LocalRmatrixType;

            // right hand sides
            typedef SparseRowMatrix< double > RHSType;
            // H_{1}\in R^{M}
            RHSType H1rhs( sigmaSpace_.size(), 1, 1 );
            // H_{2}\in R^{L}
            RHSType H2rhs( velocitySpace_.size(), 1, 1 );
            // H_{3}\in R^{K}
            RHSType H3rhs( pressureSpace_.size(), 1, 1 );

            // base functions
            // of type sigma
            typedef typename DiscreteSigmaFunctionSpaceType::BaseFunctionSetType
                SigmaBaseFunctionSetType;
            // of type u
            typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
                VelocityBaseFunctionSetType;
            // of type p
            typedef typename DiscretePressureFunctionSpaceType::BaseFunctionSetType
                PressureBaseFunctionSetType;

            // eps
            const double eps = 1.0e-14;

#ifndef NLOG
            // logging stuff
            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();
            int infoLogState = Logger().GetStreamFlags( Logging::LOG_INFO );
            int debugLogState = Logger().GetStreamFlags( Logging::LOG_DEBUG );
            bool output = false;
            const int outputEntity = 0;
            int entityNR = 0;
            infoStream << "\nthis is StokesPass::apply()" << std::endl;
            infoStream << "- starting gridwalk" << std::endl;
            Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif

            // walk the grid
            EntityIteratorType entityItEnd = velocitySpace_.end();
            for ( EntityIteratorType entityIt = velocitySpace_.begin(); entityIt != entityItEnd; ++entityIt ) {

                // entity and geometry
                EntityType& entity = *entityIt;
                typedef typename EntityType::Geometry
                    EntityGeometryType;

                const EntityGeometryType& geometry = entity.geometry();

                // local functions
                LocalDiscreteVelocityFunctionType localVelocity = velocity.localFunction( entity );
                LocalDiscretePressureFunctionType localPressure = pressure.localFunction( entity );
                LocalDiscreteSigmaFunctionType localSigma = sigma.localFunction( entity );

                // local matrices for the volume integral
                LocalMmatrixType localMmatrixElement = Mmatrix.localMatrix( entity, entity );
                LocalWmatrixType localWmatrixElement = Wmatrix.localMatrix( entity, entity );
                LocalXmatrixType localXmatrixElement = Xmatrix.localMatrix( entity, entity );
                LocalYmatrixType localYmatrixElement = Ymatrix.localMatrix( entity, entity );
                LocalZmatrixType localZmatrixElement = Zmatrix.localMatrix( entity, entity );
                LocalEmatrixType localEmatrixElement = Ematrix.localMatrix( entity, entity );
                LocalRmatrixType localRmatrixElement = Rmatrix.localMatrix( entity, entity );

                // get basefunctionsets
                SigmaBaseFunctionSetType sigmaBaseFunctionSetElement = sigmaSpace_.baseFunctionSet( entity );
                VelocityBaseFunctionSetType velocityBaseFunctionSetElement = velocitySpace_.baseFunctionSet( entity );
                PressureBaseFunctionSetType pressureBaseFunctionSetElement = pressureSpace_.baseFunctionSet( entity );
                const int numSigmaBaseFunctionsElement = sigmaBaseFunctionSetElement.numBaseFunctions();
                const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.numBaseFunctions();
                const int numPressureBaseFunctionsElement = pressureBaseFunctionSetElement.numBaseFunctions();

                // get quadrature
                VolumeQuadratureType volumeQuadratureElement( entity, ( 2 * sigmaSpaceOrder ) + 1 );
#ifndef NLOG
                if ( outputEntity == entityNR ) output = true;
                if ( output ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                debugStream << "  - entity " << outputEntity << std::endl;
                debugStream << "  - numSigmaBaseFunctionsElement: " << numSigmaBaseFunctionsElement << std::endl;
                debugStream << "  - numVelocityBaseFunctionsElement: " << numVelocityBaseFunctionsElement << std::endl;
                debugStream << "  - numPressureBaseFunctionsElement: " << numPressureBaseFunctionsElement << std::endl;
                debugStream << "  - calculating E" << std::endl;
                debugStream << "    =============" << std::endl;
                bool Moutput = false;
                bool Xoutput = false;
                bool Zoutput = false;
                bool Eoutput = false;
                bool Routput = false;
                // we want logging at the following base functions
                const int logBaseI = 0;
                const int logBaseJ = 0;
                Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                // calculate volume integrals

                // (M)_{i,j} = \int_{T}\tau_{i}:\tau_{j}dx
                // we build M^-1 in fact, because M should be diagonal, and
                // inversion is pretty easy
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {

                        double M_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Moutput = true;
//                        if ( output && Moutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            WorldCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // calculate \tau_{i}:\tau_{j}
                            SigmaRangeType tau_i( 0.0 );
                            SigmaRangeType tau_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                            sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                            double tau_i_times_tau_j = colonProduct( tau_i, tau_j );
                            // calculate M_i_j
                            M_i_j += elementVolume *
                                        integrationWeight *
                                        tau_i_times_tau_j;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
                            Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
                            debugStream << "\n      - colonProduct( tau_i, tau_j ): " << tau_i_times_tau_j << std::endl;
                            debugStream << "      - M_i_j: " << M_i_j << std::endl;
#endif
                        } // done sum over quadrature points

                        // if small, should be zero
                        if ( fabs( M_i_j ) < eps ) {
                            M_i_j = 0.0;
                        }
                        // else invert
                        else {
                            M_i_j = 1.0 / M_i_j;
                        }

                        // add to matrix
                        localMmatrixElement.add( i, j, M_i_j );
                    }
#ifndef NLOG
                        Moutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                } // done calculating M

                // (W)_{i,j}=\int_{\partial T}\hat{u}^{U}(v_{j})\cdot\tau_{i}\cdot n_{T}ds
                //           -\int_{T}v_{j}\cdot(\nabla \cdot \tau_{i})dx
                // we only compute the contribution of the 2nd integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double W_i_j = 0.0;
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            WorldCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // calculate \tau_{i}:\tau_{j}
                            SigmaRangeType tau_i( 0.0 );
                            VelocityRangeType v_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            VelocityRangeType divergence_of_tau_i = sigmaDivergenceOf( tau_i );
                            double v_j_times_divergence_of_tau_i = v_j * divergence_of_tau_i;
                            W_i_j += elementVolume *
                                        integrationWeight *
                                        v_j_times_divergence_of_tau_i * -1.0;
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( W_i_j ) < eps ) {
                            W_i_j = 0.0;
                        }
                        // add to matrix
                        localWmatrixElement.add( i, j, W_i_j );
                    }
                } // done calculationg W

                // \mu(X)_{i,j}=\mu\int_{T}\tau_{j}:\nabla v_{i} dx
                //              -\mu\int_{\partial}v_{i}\cdot\hat{\sigma}^{\sigma}(\tau_{j})\cdot n_{T}ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double X_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
//                        if ( output && Xoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            WorldCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // calculate \tau_{j}:\nabla v_{i}
                            SigmaRangeType gradient_of_v_i( 0.0 );
                            SigmaRangeType tau_j( 0.0 );
                            velocityBaseFunctionSetElement.jacobian( i, x, gradient_of_v_i );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_j );
                            double tau_j_times_gradient_v_i =
                                colonProduct( tau_j, gradient_of_v_i );
                            X_i_j += elementVolume *
                                        integrationWeight *
                                        tau_j_times_gradient_v_i *
                                        mu;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldVector( gradient_of_v_i, "gradient of v_i", debugStream, "      " );
                            Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
                            debugStream << "\n      - colonProduct( tau_j, grad v_i ): " << tau_j_times_gradient_v_i << std::endl;
                            debugStream << "      - X_i_j: " << X_i_j << std::endl;
#endif
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( X_i_j ) < eps ) {
                            X_i_j = 0.0;
                        }
                        // add to matrix
                        localXmatrixElement.add( i, j, X_i_j );
#ifndef NLOG
                        Xoutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculating X

                // (Z)_{i,j}=\int_{T} -q_{j}\cdot(\naba\cdot v_i) dx
                //           +\int_{\partial T}\hat{p}^{P}(q_j)\cdot v_{i}\cdot n_{T} ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                        double Z_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
//                        if ( output && Zoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            WorldCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // calculate q_{j}\cdot(\nabla\cdot v_i)
                            SigmaRangeType gradient_of_v_i( 0.0 );
                            PressureRangeType q_j( 0.0 );
                            velocityBaseFunctionSetElement.jacobian( i, x, gradient_of_v_i );
                            double divergence_of_v_i = velocityDivergenceOutOfGradient( gradient_of_v_i );
                            pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                            double q_j_times_divergence_of_v_i = q_j * divergence_of_v_i;
                            Z_i_j += -1.0 *
                                        elementVolume *
                                        integrationWeight *
                                        q_j_times_divergence_of_v_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( gradient_of_v_i, "gradient_of_v_i", debugStream, "      " );
                            Stuff::printFieldVector( q_j, "q_j", debugStream, "      " );
                            debugStream << "\n      - q_j_times_divergence_of_v_i: " << q_j_times_divergence_of_v_i << std::endl;
                            debugStream << "      - Z_i_j: " << Z_i_j << std::endl;
#endif
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( Z_i_j ) < eps ) {
                            Z_i_j = 0.0;
                        }
                        // add to matrix
                        localZmatrixElement.add( i, j, Z_i_j );
                    }
#ifndef NLOG
                        Zoutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                } // done calculating Z

                // (E)_{i,j}=\int_{T} -v_{j}\cdot(\naba q_i) dx
                //           +\int_{\partial T}\hat{u}^{U}(v_j)\cdot n_{T} q_{i}ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double E_i_j = 0.0;
#ifndef NLOG
                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
                        if ( output && Eoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            WorldCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // calculate v_{j}\cdot(\nabla q_i)
                            typename DiscretePressureFunctionSpaceType::JacobianRangeType jacobian_of_q_i( 0.0 );
                            VelocityRangeType v_j( 0.0 );
                            pressureBaseFunctionSetElement.jacobian( i, x, jacobian_of_q_i );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            VelocityRangeType gradient_of_q_i( jacobian_of_q_i[0] );
                            double v_j_times_gradient_of_q_i = v_j * gradient_of_q_i;
                            E_i_j += -1.0 *
                                        elementVolume *
                                        integrationWeight *
                                        v_j_times_gradient_of_q_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldVector( gradient_of_q_i, "gradient_of_q_i", debugStream, "      " );
                            Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                            debugStream << "\n      - v_j_times_gradient_of_q_i: " << v_j_times_gradient_of_q_i << std::endl;
                            debugStream << "      - E_i_j: " << E_i_j << std::endl;
#endif
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( E_i_j ) < eps ) {
                            E_i_j = 0.0;
                        }
                        // add to matrix
                        localEmatrixElement.add( i, j, E_i_j );
                    }
#ifndef NLOG
                        Zoutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                } // done calculating E

#ifndef NLOG
                if ( output ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                debugStream << "  - done calculating E" << std::endl;
                debugStream << "    ==================" << std::endl;
                Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
                output = false;
                ++entityNR;
#endif

            } // done walking the grid
#ifndef NLOG
            infoStream << "- gridwalk done" << std::endl;
            Ematrix.matrix().print( std::cout );
            Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // return to original state
#endif
        } // end of apply

        virtual void compute( const TotalArgumentType &arg, DestinationType &dest) const
        {}

        virtual void allocateLocalMemory()
        {}

    private:
        DiscreteModelType& discreteModel_;
        GridPartType& gridPart_;
        DiscreteStokesFunctionSpaceWrapperType& spaceWrapper_;
        DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscretePressureFunctionSpaceType& pressureSpace_;
        DiscreteSigmaFunctionSpaceType sigmaSpace_;

        /**
         *  \todo   doc
         **/
        double colonProduct(    const SigmaRangeType& arg1,
                                const SigmaRangeType& arg2 ) const
        {
            assert( arg1.rowdim() == arg2.coldim() );
            double ret = 0.0;
            // iterators
            typedef typename SigmaRangeType::ConstRowIterator
                ConstRowIteratorType;
            typedef typename SigmaRangeType::row_type::ConstIterator
                ConstIteratorType;
            ConstRowIteratorType arg1RowItEnd = arg1.end();
            ConstRowIteratorType arg2RowItEnd = arg2.end();
            ConstRowIteratorType arg2RowIt = arg2.begin();
            for (   ConstRowIteratorType arg1RowIt = arg1.begin();
                    arg1RowIt != arg1RowItEnd, arg2RowIt != arg2RowItEnd;
                    ++arg1RowIt, ++arg2RowIt ) {
                ConstIteratorType row1ItEnd = arg1RowIt->end();
                ConstIteratorType row2ItEnd = arg2RowIt->end();
                ConstIteratorType row2It = arg2RowIt->begin();
                for (   ConstIteratorType row1It = arg1RowIt->begin();
                        row1It != row1ItEnd, row2It != row2ItEnd;
                        ++row1It, ++row2It ) {
                    ret += *row1It * *row2It;
                }
            }
            return ret;
        }

        /**
         *  \todo   doc
         **/
        VelocityRangeType sigmaDivergenceOf( const SigmaRangeType& arg ) const
        {
            VelocityRangeType ret( 0.0 );
            typedef typename SigmaRangeType::ConstRowIterator
                ArgConstRowIteratorType;
            typedef typename SigmaRangeType::row_type::ConstIterator
                ArgConstIteratorType;
            typedef typename VelocityRangeType::Iterator
                RetIteratorType;
            ArgConstRowIteratorType argRowItEnd = arg.end();
            RetIteratorType retItEnd = ret.end();
            RetIteratorType retIt = ret.begin();
            for (   ArgConstRowIteratorType argRowIt = arg.begin();
                    argRowIt != argRowItEnd, retIt != retItEnd;
                    ++argRowIt, ++retIt ) {
                ArgConstIteratorType argItEnd = argRowIt->end();
                for (   ArgConstIteratorType argIt = argRowIt->begin();
                        argIt != argItEnd;
                        ++argIt ) {
                            *retIt += *argIt;
                }
            }
            return ret;
        }

        /**
         *  \todo   doc
         **/
        double velocityDivergenceOutOfGradient( const SigmaRangeType& arg ) const
        {
            double ret( 0.0 );
            typedef typename SigmaRangeType::ConstRowIterator
                ConstRowIteratorType;
            ConstRowIteratorType rowItEnd = arg.end();
            int i = 0;
            for (   ConstRowIteratorType rowIt = arg.begin();
                    rowIt != rowItEnd;
                    ++rowIt ) {
                ret += (*rowIt)[i];
                ++i;
            }
        }




};

}
#endif  // end of stokespass.hh




//            typedef SparseRowMatrix<double> MatrixType;
//            typedef MatrixType B_OperatorType;
//            typedef MatrixType B_Transposed_OperatorType;
//            typedef MatrixType C_OperatorType;
//

//        private:
//            enum { VelocityDimDomain = VelocityDiscreteFunctionSpaceType::DimDomain };
//            enum { VelocityDimRange = VelocityDiscreteFunctionSpaceType::DimRange };
//            enum { VelocityPolOrder = VelocityDiscreteFunctionSpaceType::polOrd };
//
//            //no need for diff velo/press typedefs, these will always be the same
//            typedef typename VelocityDiscreteFunctionSpaceType::GridType GridType;
//            typedef typename VelocityDiscreteFunctionSpaceType::GridPartType GridPartType;
//
//            typedef typename VelocityDiscreteFunctionSpaceType::DomainFieldType VelocityDomainFieldType;
//            typedef typename VelocityDiscreteFunctionSpaceType::RangeFieldType VelocityRangeFieldType;
//
//            typedef typename DiscreteModelType::VolumeQuadratureType VolumeQuadratureType;
//            typedef typename DiscreteModelType::FaceQuadratureType FaceQuadratureType;
//
//            //afaics these are ony ever used internally
//            typedef MatrixFunctionSpace<VelocityDomainFieldType, VelocityRangeFieldType, VelocityDimDomain, VelocityDimRange, VelocityDimRange > GradientSpaceType;
//            typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, VelocityPolOrder> DiscreteGradientSpaceType;
//
//            typedef SparseRowMatrixObject<DiscreteGradientSpaceType,VelocityDiscreteFunctionSpaceType> VelocityGradMatType;
//            typedef typename VelocityGradMatType::LocalMatrixType LocalVelocityGradMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,DiscreteGradientSpaceType> VelocityDivMatType;
//            typedef typename VelocityDivMatType::LocalMatrixType LocalVelocityDivMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> VelocityStabMatType;
//            typedef typename VelocityStabMatType::LocalMatrixType LocalVelocityStabMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,PressureDiscreteFunctionSpaceType> PressureGradMatType;
//            typedef typename PressureGradMatType::LocalMatrixType LocalPressureGradMatType;
//
//            typedef SparseRowMatrixObject<PressureDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> PressureDivMatType;
//            typedef typename PressureDivMatType::LocalMatrixType LocalPressureDivMatType;
//
//            typedef SparseRowMatrixObject<PressureDiscreteFunctionSpaceType,PressureDiscreteFunctionSpaceType> PressureStabMatType;
//            typedef typename PressureStabMatType::LocalMatrixType LocalPressureStabMatType;
//
//
//        public:
//
//            StokesPass( PreviousPassType& prev_pass,
//                        const DiscreteFunctionSpaceType& velo_space,
//                        const PressureDiscreteFunctionSpaceType& press_space,
//                        const DiscreteModelType& disc_model)
//                : BaseType( prev_pass, velo_space ),
//                velo_space_( velo_space ),
//                press_space_( press_space ),
//                grad_space_( press_space.gridPart() ),
//                disc_model_( disc_model )
//            {}
//
//            const B_OperatorType& Get_B_Operator() const { return b_op_; }
//            const B_Transposed_OperatorType& Get_B_Transposed_Operator() const { return b_transp_op_; }
//            const C_OperatorType& Get_C_Operator() const { return c_op_; }
//
//            MatrixType& systemMatrix()
//            {}
//
//            const VelocityDiscreteFunctionType& rhs1() const
//            {}
//
//
//            virtual void applyLocal(EntityType& entity) const
//            {
//                    //- typedefs
//                typedef typename VelocityDiscreteFunctionSpaceType::IndexSetType IndexSetType;
//                typedef typename VelocityDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
//
//                typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
//                typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
//
//                typedef typename PressureDiscreteFunctionSpaceType::IndexSetType PressureIndexSetType;
//                typedef typename PressureDiscreteFunctionSpaceType::BaseFunctionSetType PressureBaseFunctionSetType;
//                //- statements
//
//                typedef typename DestinationType::LocalFunctionType LocalFuncType;
//                VolumeQuadratureType volQuad(entity, 2); //! \todo get order from model too
//
//
//                double massVolElInv;
//                double vol = volumeElement(entity, volQuad,massVolElInv);
//
//                const BaseFunctionSetType& bsetEn = velo_space_.baseFunctionSet(entity);
//                const int numDofs = bsetEn.numBaseFunctions();
//
//                const GradientBaseFunctionSetType& gradbsetEn = grad_space_.baseFunctionSet(entity);
//                const int gradientNumDofs = gradbsetEn.numBaseFunctions();
//
//                const PressureBaseFunctionSetType& pressurebsetEn = press_space_.baseFunctionSet(entity);
//                const int pressureNumDofs = pressurebsetEn.numBaseFunctions();
//                LocalVelocityGradMatType en_grad = gradMatrix_.localMatrix(entity,entity);
//                LocalVelocityDivMatType en_div= divMatrix_.localMatrix(entity,entity);
//                LocalVelocityStabMatType en_stab=stabMatrix_.localMatrix(entity,entity);
//                LocalPressureGradMatType en_pressGrad=pressureGradMatrix_.localMatrix(entity,entity);
//                LocalPressureDivMatType en_pressDiv=pressureDivMatrix_.localMatrix(entity,entity);
//                LocalPressureStabMatType en_pressStab=pressureStabMatrix_.localMatrix(entity,entity);


                //compute volume integral contributions
//                for (int l = 0; l < quadNop ; ++l)
//                {
//                } // end volume integral contributions


                //compute surface integral contributions
//                IntersectionIterator endnit = entity.ileafend();
//                double dtLocal = 0.0;
//                double minvol = vol;
//                for (IntersectionIterator nit = entity.ileafbegin(); nit != endnit; ++nit)
//                {
//                    //	  int twistSelf = twistUtil_.twistInSelf(nit);
//                    FaceQuadratureType faceQuadInner(spc_.gridPart(),nit, faceQuadOrd_,
//                    FaceQuadratureType::INSIDE);
//
//                    if (nit.neighbor())
//                    {
//                    } //end if inner
//
//                    if (nit.boundary())
//                    {
//                    } // end if boundary
//                }
                //end compute surface integral contributions
            //}//end void applyLocal(EntityType& entity) const


        //private:
            //will prolly be generated on-the-fly, not stored as members
//            B_OperatorType b_op_;
//            B_Transposed_OperatorType b_transp_op_;
//            C_OperatorType c_op_;
//
//            //more dummies
//            // VelocityDiscreteFunctionType& rhs;
//
//
//            VelocityGradMatType gradMatrix_;
//            VelocityDivMatType divMatrix_;
//            VelocityStabMatType stabMatrix_;
//            PressureGradMatType pressureGradMatrix_;
//            PressureDivMatType pressureDivMatrix_;
//            PressureStabMatType pressureStabMatrix_;
//
