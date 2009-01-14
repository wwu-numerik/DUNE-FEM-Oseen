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

        //! local coordinate type on an intersection
        typedef typename FaceQuadratureType::LocalCoordinateType
            LocalIntersectionCoordinateType;

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
                MInversMatrixType;
            MInversMatrixType MInversMatrix( sigmaSpace_, sigmaSpace_ );
            MInversMatrix.reserve();
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
            typedef typename MInversMatrixType::LocalMatrixType
                LocalMInversMatrixType;
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
            // H_{1}\in R^{M}
            DiscreteSigmaFunctionType H1rhs( "H1", sigmaSpace_ );
            H1rhs.clear();
            // H_{2}\in R^{L}
            DiscreteVelocityFunctionType H2rhs( "H2", velocitySpace_ );
            H2rhs.clear();
            // H_{3}\in R^{K}
            DiscretePressureFunctionType H3rhs( "H3", pressureSpace_ );
            H3rhs.clear();

            // local right hand sides
            // H_{1}\in R^{M}
            typedef typename DiscreteSigmaFunctionType::LocalFunctionType
                LocalH1rhsType;
            // H_{2}\in R^{L}
            typedef typename DiscreteVelocityFunctionType::LocalFunctionType
                LocalH2rhsType;
            // H_{3}\in R^{K}
            typedef typename DiscretePressureFunctionType::LocalFunctionType
                LocalH3rhsType;

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

            // twist utility, needed to twist inside and outside face quadratures
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
            const int outputEntity = 0;
            const int outputIntersection = 0;
            int entityNR = 0;
            int intersectionNR = 0;
            int numberOfBoundaryIntersections = 0;
            int numberOfInnerIntersections = 0;
            const bool Mprint = false;
            const bool Wprint = true;
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
                LocalMInversMatrixType localMInversMatrixElement = MInversMatrix.localMatrix( entity, entity );
                LocalWmatrixType localWmatrixElement = Wmatrix.localMatrix( entity, entity );
                LocalXmatrixType localXmatrixElement = Xmatrix.localMatrix( entity, entity );
                LocalYmatrixType localYmatrixElement = Ymatrix.localMatrix( entity, entity );
                LocalZmatrixType localZmatrixElement = Zmatrix.localMatrix( entity, entity );
                LocalEmatrixType localEmatrixElement = Ematrix.localMatrix( entity, entity );
                LocalRmatrixType localRmatrixElement = Rmatrix.localMatrix( entity, entity );

                // local right hand sides
                LocalH1rhsType LocalH1rhs = H1rhs.localFunction( entity );
                LocalH2rhsType LocalH2rhs = H2rhs.localFunction( entity );
                LocalH3rhsType LocalH3rhs = H3rhs.localFunction( entity );

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
                // calculate volume integrals on the entity

                // (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1}
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
                            double tau_j_times_tau_i = colonProduct( tau_j, tau_i );
                            // calculate M_i_j
                            M_i_j += elementVolume *
                                        integrationWeight *
                                        tau_j_times_tau_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                            debugStream << "      - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
                            Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
                            debugStream << "\n      - tau_j_times_tau_i: " << tau_j_times_tau_i << std::endl;
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
                        localMInversMatrixElement.add( i, j, M_i_j );
#ifndef NLOG
                        Moutput = false;
                        Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                    }
                } // done calculating M

                // (W)_{i,j} += \int_{T}v_{j}\cdot(\nabla\cdot\tau_{i})dx
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double W_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
                        if ( entityOutput && Woutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = W ========================" << std::endl;
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
                                        v_j_times_divergence_of_tau_i;
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

                // (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double X_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
                        if ( entityOutput && Xoutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                        debugStream << "    = X ========================" << std::endl;
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

                // (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx
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

                // (H2)_{j} += \int_{T}f\cdot v_{j}dx
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
                    LocalH2rhs[ j ] = H2_j;
#ifndef NLOG
                    H2output = false;
                    Logger().SetStreamFlags( Logging::LOG_DEBUG, Logging::LOG_NONE ); // disable logging
#endif
                } // done calculating H2

                // (E)_{i,j} += -\int_{T}v_{j}\cdot\nabla q_{i}dx
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

                        // get intersection geometry
                        typedef typename IntersectionIteratorType::LocalGeometry
                            IntersectionGeometryType;
                        const IntersectionGeometryType& intersectionGeoemtry = intIt.intersectionSelfLocal();

                        // local matrices for the surface integral
                        LocalMInversMatrixType localMInversMatrixNeighbour = MInversMatrix.localMatrix( entity, neighbour );
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

                        //                                                                                                               // we will call this one
                        // (W)_{i,j} += -\int_{\varepsilon\in \Epsilon_{I}^{T}}\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}ds // element integral
                        //           += -\int_{\varepsilon\in \Epsilon_{D}^{T}}\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}ds // neighbour integral
                        for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                            // compute the element integral
                            for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                double W_i_j = 0.0;
#ifndef NLOG
                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
                                if ( intersectionOutput && Woutput ) Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // enable logging
                                debugStream << "      = W element ======================" << std::endl;
                                debugStream << "      basefunctions " << i << " " << j << std::endl;
                                debugStream << "      volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                                // sum over all quadrature points
                                for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x
                                    ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                    // get the integration factor
                                    double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                    // get the quadrature weight
                                    double integrationWeight = faceQuadratureElement.weight( quad );
                                    // calculate \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
                                    SigmaRangeType tau_i( 0.0 );
                                    VelocityRangeType v_j( 0.0 );
                                    VelocityRangeType u_sigma_u_plus_flux( 0.0 );
                                    VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                    sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                                    velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                    discreteModel_.velocitySigmaFlux(   intIt,
                                                                        0.0,
                                                                        localX,
                                                                        DiscreteModelType::inside,
                                                                        v_j,
                                                                        u_sigma_u_plus_flux );
                                    VelocityRangeType tau_i_times_n_t( 0.0 );
                                    tau_i.mv( outerNormal, tau_i_times_n_t );
                                    const double flux_times_tau_i_times_n_t = u_sigma_u_plus_flux * tau_i_times_n_t;
                                    W_i_j += -1.0 *
                                                elementVolume *
                                                integrationWeight *
                                                flux_times_tau_i_times_n_t;
#ifndef NLOG
                                    debugStream << "      - quadPoint " << quad;
                                    Stuff::printFieldVector( x, "x", debugStream, "        " );
                                    Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                    Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                    debugStream << "\n      - elementVolume: " << elementVolume << std::endl;
                                    debugStream << "      - integrationWeight: " << integrationWeight;
                                    Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
                                    Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                    Stuff::printFieldVector( u_sigma_u_plus_flux, "u_sigma_u_plus_flux", debugStream, "        " );
                                    Stuff::printFieldVector( tau_i_times_n_t, "tau_i_times_n_t", debugStream, "        " );
                                    debugStream << "\n        - flux_times_tau_i_times_n_t: " << flux_times_tau_i_times_n_t << std::endl;
                                    debugStream << "        - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
#endif
                                } // done sum over all quadrature points
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
                            } // done with the element Integral
                            for ( int j = 0; j < numVelocityBaseFunctionsNeighbour; ++j ) {
                            } // done with the neighbour integral
                        } // done calculating W

                        ++numberOfInnerIntersections;

                    } // done with those inside the grid

                    // if we are on the boundary of the grid
                    if ( !intIt.neighbor() && intIt.boundary() ) {
                        ++numberOfBoundaryIntersections;
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
                intersectionNR = 0;
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
            debugStream << "        " << numberOfInnerIntersections << " intersections inside and" << std::endl;
            debugStream << "        " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
            if ( Mprint || Wprint || Xprint || Zprint || Eprint || Rprint || H1print || H2print || H3print ) {
                debugStream << "- printing matrices" << std::endl;
                if ( Mprint ) {
                    debugStream << " - M ==============" << std::endl;
                    debugStream.Log( &MInversMatrixType::MatrixType::print,  MInversMatrix.matrix() );
                }
                if ( Wprint ) {
                    debugStream << " - W ==============" << std::endl;
                    debugStream.Log( &WmatrixType::MatrixType::print,  Wmatrix.matrix() );
                }
                if ( Xprint ) {
                    debugStream << " - X ==============" << std::endl;
                    debugStream.Log( &XmatrixType::MatrixType::print,  Xmatrix.matrix() );
                }
                if ( Zprint ) {
                    debugStream << " - Z ==============" << std::endl;
                    debugStream.Log( &ZmatrixType::MatrixType::print,  Zmatrix.matrix() );
                }
                if ( Eprint ) {
                    debugStream << " - E ==============" << std::endl;
                    debugStream.Log( &EmatrixType::MatrixType::print,  Ematrix.matrix() );
                }
                if ( Rprint ) {
                    debugStream << " - R ==============" << std::endl;
                    debugStream.Log( &RmatrixType::MatrixType::print,  Rmatrix.matrix() );
                }
                if ( H1print ) {
                    debugStream << " - H1 =============" << std::endl;
                    debugStream.Log( &DiscreteSigmaFunctionType::print, H1rhs );
                }
                if ( H2print ) {
                    debugStream << " - H2 =============" << std::endl;
                    debugStream.Log( &DiscreteVelocityFunctionType::print, H2rhs );
                }
                if ( H3print ) {
                    debugStream << " - H3 =============" << std::endl;
                    debugStream.Log( &DiscretePressureFunctionType::print, H3rhs );
                }
                debugStream << "- done printing matrices" << std::endl;
            }
            Logger().SetStreamFlags( Logging::LOG_DEBUG, debugLogState ); // return to original state
#endif

            InvOpType op( *this, 1.0,1.0,1,1 );
            op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs );


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
