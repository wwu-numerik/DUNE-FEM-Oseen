/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/stokes/saddlepoint_inverse_operator.hh>


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
        //! own type
        typedef StokesPass< DiscreteModelImp, PreviousPassImp, PassID >
            ThisType;

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

        //! Coordinate type on the element
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            ElementCoordinateType;

        //! Coordinate type on an intersection
        typedef typename FaceQuadratureType::LocalCoordinateType
            IntersectionCoordinateType;

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

        //! type of the used solver
        typedef SaddlepointInverseOperator< ThisType >
            InvOpType;

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

        //! used in Postprocessing to get refs to gridparts, spaces
        const DiscreteStokesFunctionSpaceWrapperType& GetFunctionSpaceWrapper() const
        {
            return spaceWrapper_;
        }

        virtual void apply( const DomainType &arg, RangeType &dest) const
        {

            // viscosity
            const double mu = discreteModel_.viscosity();

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

            // twist utility, needed to evaluate basefunctions from both sides
            // of a face
            typedef typename Dune::TwistUtility< GridType >
                TwistUtilityType;
            TwistUtilityType twistUtility( gridPart_.grid() );

            // eps
            const double eps = 1.0e-14;

#ifndef NLOG
            // logging stuff
            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();
            int infoLogState = Logger().GetStreamFlags( Logging::LOG_INFO );
            int debugLogState = Logger().GetStreamFlags( Logging::LOG_DEBUG );
            bool entityOutput = false;
            bool intersectionOutput = false;
            const int outputEntity = 1;
            const int outputIntersection = 2;
            int entityNR = 0;
            int intersectionNR = 0;
            int numberOfBoundaryIntersections = 0;
            int numberOfInnerIntersections = 0;
            const bool Mprint = true;
            const bool Wprint = false;
            const bool Xprint = false;
            const bool Zprint = false;
            const bool Eprint = false;
            const bool Rprint = false;
            const bool H1print = false;
            const bool H2print = false;
            const bool H3print = false;
            infoStream << "\nthis is StokesPass::apply()" << std::endl;
            infoStream << "- starting gridwalk" << std::endl;
            Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif

            // walk the grid
            EntityIteratorType entityItEnd = velocitySpace_.end();
            for ( EntityIteratorType entityIt = velocitySpace_.begin(); entityIt != entityItEnd; ++entityIt ) {

                // entity and geometry
                const EntityType& entity = *entityIt;
                typedef typename EntityType::Geometry
                    EntityGeometryType;
                const EntityGeometryType& geometry = entity.geometry();

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
                if ( outputEntity == entityNR ) entityOutput = true;
                if ( entityOutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                debugStream << "  - entity " << outputEntity << std::endl;
                debugStream << "  - numSigmaBaseFunctionsElement: " << numSigmaBaseFunctionsElement << std::endl;
                debugStream << "  - numVelocityBaseFunctionsElement: " << numVelocityBaseFunctionsElement << std::endl;
                debugStream << "  - numPressureBaseFunctionsElement: " << numPressureBaseFunctionsElement << std::endl;
                debugStream << "  - start calculations on entity" << std::endl;
                debugStream << "    ============================" << std::endl;
                bool Moutput = false;
                bool Woutput = false;
                bool Xoutput = false;
                bool Zoutput = false;
                bool Eoutput = false;
                bool Routput = false;
                bool H1output = false;
                bool H2output = false;
                bool H3output = false;
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
                        if ( entityOutput && Moutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = M ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
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
                            debugStream << "\n      - tau_i_times_tau_j: " << tau_i_times_tau_j << std::endl;
                            debugStream << "      - M_" << i << "_" << j << "+=: " << M_i_j << std::endl;
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
#ifndef NLOG
                        Moutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculating M

                // (W)_{i,j}=\int_{\partial T}\hat{u}^{U}\cdot\tau_{i}\cdot n_{T}ds
                //           -\int_{T}v_{j}\cdot(\nabla \cdot \tau_{i})dx
                // we only compute the contribution of the 2nd integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double W_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
                        if ( entityOutput && Woutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = W =========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
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
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
                            Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                            debugStream << "\n      - v_j_times_divergence_of_tau_i: " << v_j_times_divergence_of_tau_i << std::endl;
                            debugStream << "      - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
#endif
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( W_i_j ) < eps ) {
                            W_i_j = 0.0;
                        }
                        // add to matrix
                        localWmatrixElement.add( i, j, W_i_j );
#ifndef NLOG
                        Woutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculationg W

                // \mu(X)_{i,j}=\mu\int_{T}\tau_{j}:\nabla v_{i} dx
                //              -\mu\int_{\partial}v_{i}\cdot\hat{\sigma}^{\sigma}\cdot n_{T}ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double X_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
                        if ( entityOutput && Xoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = mu X =====================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
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
                            debugStream << "      - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
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
                //           +\int_{\partial T}\hat{p}^{P}\cdot v_{i}\cdot n_{T} ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                        double Z_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
                        if ( entityOutput && Zoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = Z ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
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
                            debugStream << "      - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
#endif
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( Z_i_j ) < eps ) {
                            Z_i_j = 0.0;
                        }
                        // add to matrix
                        localZmatrixElement.add( i, j, Z_i_j );
#ifndef NLOG
                        Zoutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculating Z

                // (E)_{i,j}=\int_{T} -v_{j}\cdot(\naba q_i) dx
                //           +\int_{\partial T}\hat{u}^{U}\cdot n_{T} q_{i}ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double E_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
                        if ( entityOutput && Eoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = E ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
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
                            debugStream << "      - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
#endif
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( E_i_j ) < eps ) {
                            E_i_j = 0.0;
                        }
                        // add to matrix
                        localEmatrixElement.add( i, j, E_i_j );
#ifndef NLOG
                        Eoutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculating E

                // (H2)_{j}=\int_{T}f\cdot v_{j} dx
                //          +\mu(\int_{\partial T}v_{j}\cdot\hat{\sigma}^{RHS} \cdot n_{T})
                //          -\int_{\partial T}\hat{p}^{RHS}\cdot v_j \cdot n_{T} ds
                // we only compute the contribution of the 1st integral here,
                // the surface integral will be computed later on
                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                    double H2_j = 0.0;
#ifndef NLOG
//                    if ( ( j == logBaseJ ) ) H2output = true;
                    if ( entityOutput && H2output ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                    debugStream << "    = H2 =======================" << std::endl;
                    debugStream << "    basefunction " << " " << j << std::endl;
                    debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                    // sum over all quadratur points
                    for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                        // get x
                        ElementCoordinateType x = volumeQuadratureElement.point( quad );
                        // get the integration factor
                        double elementVolume = geometry.integrationElement( x );
                        // get the quadrature weight
                        double integrationWeight = volumeQuadratureElement.weight( quad );
                        // calculate f\cdot v_j
                        VelocityRangeType v_j( 0.0 );
                        VelocityRangeType f( 0.0 );
                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                        discreteModel_.force( 0.0, x, f );
                        double f_times_v_j = f * v_j;
                        H2_j += elementVolume *
                                    integrationWeight *
                                    f_times_v_j;
#ifndef NLOG
                        debugStream << "    - quadPoint " << quad;
                        Stuff::printFieldVector( x, "x", debugStream, "      " );
                        debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                        debugStream << "      - integrationWeight: " << integrationWeight;
                        Stuff::printFieldVector( f, "f", debugStream, "      " );
                        Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                        debugStream << "\n      - f_times_v_j: " << f_times_v_j << std::endl;
                        debugStream << "      - H2_" << j << "+=: " << H2_j << std::endl;
#endif
                    } // done sum over all quadrature points
                    // if small, should be zero
                    if ( fabs( H2_j ) < eps ) {
                        H2_j = 0.0;
                    }
                    // add to matrix
                    H2rhs.set( velocitySpace_.mapToGlobal( entity, j ), 0, H2_j );
#ifndef NLOG
                    H2output = false;
                    Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                } // done calculating H2

                // walk the neighbours
                IntersectionIteratorType intItEnd = gridPart_.iend( entity );
                for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
                        intIt != intItEnd;
                        ++intIt ) {
#ifndef NLOG
                    if ( ( outputIntersection == intersectionNR ) && entityOutput ) intersectionOutput = true;
                    if ( intersectionOutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                    debugStream << "    - intersection " << intersectionNR << std::endl;
                    debugStream << "    - start calculations on intersection" << std::endl;
                    debugStream << "      ==================================" << std::endl;
                    Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    // get intersection informations, seen from the inside
                    typedef typename IntersectionIteratorType::LocalGeometry
                        IntersectionGeometryType;
                    const IntersectionGeometryType& intersectionGeometryElement = intIt.intersectionSelfLocal();
                    // get intersection quadrature, seen from inside
                    FaceQuadratureType faceQuadratureElement(   gridPart_,
                                                                intIt,
                                                                ( 2 * sigmaSpaceOrder ) + 1,
                                                                FaceQuadratureType::INSIDE );
                    // if we are inside the grid
                    if ( intIt.neighbor() && !intIt.boundary() ) {
                        // get neighbour
                        const typename IntersectionIteratorType::EntityPointer neighbourPtr = intIt.outside();
                        const EntityType& neighbour = *neighbourPtr;

                        // local matrices for the surface integral
                        LocalMmatrixType localMmatrixNeighbour = Mmatrix.localMatrix( entity, neighbour );
                        LocalWmatrixType localWmatrixNeighbour = Wmatrix.localMatrix( entity, neighbour );
                        LocalXmatrixType localXmatrixNeighbour = Xmatrix.localMatrix( entity, neighbour );
                        LocalYmatrixType localYmatrixNeighbour = Ymatrix.localMatrix( entity, neighbour );
                        LocalZmatrixType localZmatrixNeighbour = Zmatrix.localMatrix( entity, neighbour );
                        LocalEmatrixType localEmatrixNeighbour = Ematrix.localMatrix( entity, neighbour );
                        LocalRmatrixType localRmatrixNeighbour = Rmatrix.localMatrix( entity, neighbour );

                        // get basefunctionsets
                        SigmaBaseFunctionSetType sigmaBaseFunctionSetNeighbour = sigmaSpace_.baseFunctionSet( neighbour );
                        VelocityBaseFunctionSetType velocityBaseFunctionSetNeighbour = velocitySpace_.baseFunctionSet( neighbour );
                        PressureBaseFunctionSetType pressureBaseFunctionSetNeighbour = pressureSpace_.baseFunctionSet( neighbour );
                        const int numSigmaBaseFunctionsNeighbour = sigmaBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numVelocityBaseFunctionsNeighbour = velocityBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numPressureBaseFunctionsNeighbour = pressureBaseFunctionSetNeighbour.numBaseFunctions();

                        // get intersection quadrature, seen from outside
                        FaceQuadratureType faceQuadratureNeighbour( gridPart_,
                                                                    intIt,
                                                                    ( 2 * sigmaSpaceOrder ) + 1,
                                                                    FaceQuadratureType::OUTSIDE );

                        // compute the surface integrals

                    } // done with those inside the grid

                    // if we are on the boundary of the grid
                    if ( !intIt.neighbor() && intIt.boundary() ) {
                        ++numberOfInnerIntersections;
                    } // done with those on the boundary
#ifndef NLOG
                    if ( intersectionOutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                    debugStream << "    - done calculations on intersection" << std::endl;
                    debugStream << "      =================================" << std::endl;
                    Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
                    intersectionOutput = false;
                    ++intersectionNR;
#endif
                } // done walking the neighbours

#ifndef NLOG
                if ( entityOutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                debugStream << "  - done calculations on entity" << std::endl;
                debugStream << "    ===========================" << std::endl;
                Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
                entityOutput = false;
                ++entityNR;
#endif
            } // done walking the grid
#ifndef NLOG
            infoStream << "- gridwalk done" << std::endl;
            Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
            debugStream << "  found " << entityNR << " entities," << std::endl;
            debugStream << "  found " << intersectionNR << " intersections," << std::endl;
            debugStream << "        " << numberOfBoundaryIntersections << " intersections inside and" << std::endl;
            debugStream << "        " << numberOfInnerIntersections << " intersections on the boundary." << std::endl;
            if ( Mprint || Wprint || Xprint || Zprint || Eprint || Rprint || H1print || H2print || H3print ) {
                debugStream << "- printing matrices" << std::endl;
                if ( Mprint ) {
                    debugStream << " - M ===============" << std::endl;
                    Mmatrix.matrix().print( std::cout );
                }
                if ( Wprint ) {
                    debugStream << " - W ===============" << std::endl;
                    Wmatrix.matrix().print( std::cout );
                }
                if ( Xprint ) {
                    debugStream << " - X ===============" << std::endl;
                    Xmatrix.matrix().print( std::cout );
                }
                if ( Zprint ) {
                    debugStream << " - Z ===============" << std::endl;
                    Zmatrix.matrix().print( std::cout );
                }
                if ( Eprint ) {
                    debugStream << " - E ===============" << std::endl;
                    Ematrix.matrix().print( std::cout );
                }
                if ( Rprint ) {
                    debugStream << " - R ===============" << std::endl;
                    Rmatrix.matrix().print( std::cout );
                }
                if ( H1print ) {
                    debugStream << " - H1 ==============" << std::endl;
                    H1rhs.print( std::cout );
                }
                if ( H2print ) {
                    debugStream << " - H2 ==============" << std::endl;
                    H2rhs.print( std::cout );
                }
                if ( H3print ) {
                    debugStream << " - H3 ==============" << std::endl;
                    H3rhs.print( std::cout );
                }
                debugStream << "- done printing matrices" << std::endl;
            }
            Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // return to original state
#endif
//            infoStream << "- build global matrices - " << std::endl;
//            typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
//                AmatrixType;
//            AmatrixType Amatrix( velocitySpace_, velocitySpace_ );
//            Amatrix.reserve();
//
//            XmatrixType neg_X_Minv_mat( velocitySpace_, sigmaSpace_ );
//            neg_X_Minv_mat.reserve();
//            Xmatrix.matrix().multiply( Mmatrix.matrix(), neg_X_Minv_mat.matrix() );
//            infoStream << "-    1te feritg- " << std::endl;
//            neg_X_Minv_mat.matrix().scale( -1 );
//
//            neg_X_Minv_mat.matrix().multiply( Wmatrix.matrix(), Amatrix.matrix() );
//            Ymatrix.matrix().add( Amatrix.matrix() );
////            Amatrix = Ymatrix;
//
//            Ematrix.matrix().scale( -1 );
//            Rmatrix.matrix().scale( -1 );
//
//            RHSType Fmat ( velocitySpace_.size(), 1, 1 );
//            neg_X_Minv_mat.matrix().scale ( mu );
//            neg_X_Minv_mat.matrix().multiply( H1rhs, Fmat );
////            H2rhs.add( Fmat );
//            //Fmat = H2rhs;
//
//            H3rhs.scale( -1 );
//            infoStream << "- build global matrices - done" << std::endl;
//
//            InvOpType op( *this, 1.0,1.0,1,1 );
//            op.solve( arg, dest, Ymatrix, Zmatrix, Ematrix, Rmatrix, H2rhs, H3rhs );


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
