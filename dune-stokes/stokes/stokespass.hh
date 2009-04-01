/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/stokes/saddlepoint_inverse_operator.hh>


#ifndef NLOG // if we want logging, should be removed in the end
    #include "../src/stuff.hh"
    #include "../src/logging.hh"
#endif

#include "../src/profiler.hh"

#ifdef CHEAT
    #if (POLORDER > 0)
        #warning "gradient reconstruction disabled because polorder > 0"
        #undef CHEAT
    #else
        #include "discretegradientpass.hh"
    #endif
#endif

namespace Dune
{

/**
 *  \brief  StokesPass
 *
 *             more
 *  \todo   doc
 **/
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
        static const int sigmaSpaceOrder
            = DiscreteModelType::sigmaSpaceOrder;
        //! polynomial order for the discrete velocity function space
        static const int velocitySpaceOrder
            = DiscreteModelType::velocitySpaceOrder;
        //! polynomial order for the discrete pressure function space
        static const int pressureSpaceOrder
            = DiscreteModelType::pressureSpaceOrder;


        /**
         *  \name typedefs needed for interface compliance
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

        /**
         *  \todo doc
         *  \attention  think about quadrature orders
         **/
        virtual void apply( const DomainType &arg, RangeType &dest) const
        {
            profiler().StartTiming("Pass");
            profiler().StartTiming("Pass -- ASSEMBLE");
            // viscosity
            const double mu = discreteModel_.viscosity();

            // matrices
            // M\in R^{M\times M}
            typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
                                            DiscreteSigmaFunctionSpaceType >
                MInversMatrixType;
            MInversMatrixType MInversMatrix( sigmaSpace_, sigmaSpace_ );
            MInversMatrix.reserve();
MInversMatrixType MMatrix( sigmaSpace_, sigmaSpace_ );
MMatrix.reserve();
            // W\in R^{M\times L}
            typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
                                            DiscreteVelocityFunctionSpaceType >
                WmatrixType;
            WmatrixType Wmatrix( sigmaSpace_, velocitySpace_ );
            Wmatrix.reserve();
            // X\in R^{L\times M}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
                                            DiscreteSigmaFunctionSpaceType >
                XmatrixType;
            XmatrixType Xmatrix( velocitySpace_, sigmaSpace_ );
            Xmatrix.reserve();
            // Y\in R^{L\times L}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
                                            DiscreteVelocityFunctionSpaceType >
                YmatrixType;
            YmatrixType Ymatrix( velocitySpace_, velocitySpace_ );
            Ymatrix.reserve();
            // Z\in R^{L\times K}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
                                            DiscretePressureFunctionSpaceType >
                ZmatrixType;
            ZmatrixType Zmatrix( velocitySpace_, pressureSpace_ );
            Zmatrix.reserve();
            // E\in R^{K\times L}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
                                            DiscreteVelocityFunctionSpaceType >
                EmatrixType;
            EmatrixType Ematrix( pressureSpace_, velocitySpace_ );
            Ematrix.reserve();
            // R\in R^{K\times K}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
                                            DiscretePressureFunctionSpaceType >
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

#ifdef CHEAT
            DiscreteSigmaFunctionType exactH1rhs( "exact H1", sigmaSpace_ );
            DiscreteSigmaFunctionType tmpH1rhs( "tmp H1", sigmaSpace_ );
            exactH1rhs.clear();
            tmpH1rhs.clear();

            DiscreteVelocityFunctionType exactH2rhs( "exact H2", velocitySpace_ );
            DiscreteVelocityFunctionType tmpH2rhs( "tmp H2", velocitySpace_ );
            exactH2rhs.clear();
            tmpH2rhs.clear();

            DiscretePressureFunctionType exactH3rhs( "exact H3", pressureSpace_ );
            DiscretePressureFunctionType tmpH3rhs( "tmp H3", pressureSpace_ );
            exactH3rhs.clear();
            tmpH3rhs.clear();
#endif
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
            typedef typename SigmaBaseFunctionSetType::JacobianRangeType
                SigmaJacobianRangeType;
            // of type u
            typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
                VelocityBaseFunctionSetType;
            // of type p
            typedef typename DiscretePressureFunctionSpaceType::BaseFunctionSetType
                PressureBaseFunctionSetType;

            // eps
            const double eps = Parameters().getParam( "eps", 1.0e-14 );

#ifndef NLOG
            // logging stuff
            Logging::LogStream& infoStream = Logger().Info();
//            Logging::LogStream& debugStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();
            bool entityOutput = false;
            bool intersectionOutput = false;
            const int outputEntity = 0;
            const int outputIntersection = -1;
            int entityNR = 0;
            int intersectionNR = 0;
            int numberOfEntities = 0;
            int numberOfIntersections = 0;
            int numberOfBoundaryIntersections = 0;
            int numberOfInnerIntersections = 0;
            const bool Mprint = Parameters().getParam( "Mprint", false );
            const bool Wprint = Parameters().getParam( "Wprint", false );
            const bool Xprint = Parameters().getParam( "Xprint", false );
            const bool Yprint = Parameters().getParam( "Yprint", false );
            const bool Zprint = Parameters().getParam( "Zprint", false );
            const bool Eprint = Parameters().getParam( "Eprint", false );
            const bool Rprint = Parameters().getParam( "Rprint", false );
            const bool H1print = Parameters().getParam( "H1print", false );
            const bool H2print = Parameters().getParam( "H2print", false );
            const bool H3print = Parameters().getParam( "H3print", false );
            const bool allOutput = Parameters().getParam( "allOutput", false );
            int fivePercentOfEntities = 0;
            int fivePercents = 0;
            infoStream << "this is StokesPass::apply()" << std::endl;

#ifdef CHEAT
            // init the gradient pass
            typedef Dune::DiscreteGradientModelTraits<  GridPartType,
                                                        DiscreteSigmaFunctionType >
                DiscreteGradientModelTraitsImp;

            typedef Dune::DiscreteGradientModel< DiscreteGradientModelTraitsImp >
                DiscreteGradientModelType;


            typedef Dune::StartPass< DiscreteVelocityFunctionType, -1 >
                StartPassType;

            typedef Dune::DiscreteGradientPass< DiscreteGradientModelType,
                                                StartPassType,
                                                0 >
                DiscreteGradientPassType;
            DiscreteGradientModelType discreteGradientModel;
            StartPassType startPass;
            DiscreteGradientPassType discreteGradientPass(  startPass,
                                                            discreteGradientModel,
                                                            gridPart_ );
            // get the exact solution out of the argument
            DiscreteVelocityFunctionType discreteExactVelocity = arg.discreteVelocity();
            DiscreteSigmaFunctionType computedVelocityGradient( "computed velocity gradient", sigmaSpace_ );
            computedVelocityGradient.clear();
            discreteGradientPass.apply( discreteExactVelocity, computedVelocityGradient );
            DiscretePressureFunctionType discreteExactPressure = arg.discretePressure();
		double cg_min,cg_max;
		Stuff::getMinMaxOfDiscreteFunction( computedVelocityGradient, cg_min, cg_max );
		Stuff::oneLinePrint( infoStream, computedVelocityGradient );
		Stuff::oneLinePrint( infoStream, arg.discreteVelocity() );
		Stuff::oneLinePrint( infoStream, arg.discretePressure() );
#endif

            // do an empty grid walk to get informations
            double maxGridWidth( 0.0 );
            EntityIteratorType entityItEndLog = velocitySpace_.end();
            for (   EntityIteratorType entityItLog = velocitySpace_.begin();
                    entityItLog != entityItEndLog;
                    ++entityItLog ) {
                const EntityType& entity = *entityItLog;
                // count entities
                ++numberOfEntities;
                // walk the intersections
                IntersectionIteratorType intItEnd = gridPart_.iend( entity );
                for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
                        intIt != intItEnd;
                        ++intIt ) {
                    // count intersections
                    ++numberOfIntersections;
                    maxGridWidth = Stuff::getLenghtOfIntersection( intIt ) > maxGridWidth ? Stuff::getLenghtOfIntersection( intIt ) : maxGridWidth;
                    // if we are inside the grid
                    if ( intIt.neighbor() && !intIt.boundary() ) {
                        // count inner intersections
                        ++numberOfInnerIntersections;
                    }
                    // if we are on the boundary of the grid
                    if ( !intIt.neighbor() && intIt.boundary() ) {
                        // count boundary intersections
                        ++numberOfBoundaryIntersections;
                    }
                }
            }
            if ( numberOfEntities > 19 ) {
                infoStream << "found " << numberOfEntities << " entities," << std::endl;
                infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
                infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
                infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
                infoStream << "maxGridWidth is " << maxGridWidth << std::endl;
                infoStream << "- starting gridwalk" << std::endl;
                fivePercentOfEntities = int( std::floor(double(numberOfEntities) / double(20)));
                infoStream << "  [ assembling         ]" << std::endl;
                infoStream << "  [";
            } else {
                infoStream << "found " << numberOfEntities << " entities," << std::endl;
                infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
                infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
                infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
                infoStream << "maxGridWidth is " << maxGridWidth << std::endl;
                infoStream << "- starting gridwalk" << std::endl;
            }
#endif
            // walk the grid
            EntityIteratorType entityItEnd = velocitySpace_.end();
            for (   EntityIteratorType entityIt = velocitySpace_.begin();
                    entityIt != entityItEnd;
                    ++entityIt ) {

                // get entity and geometry
                const EntityType& entity = *entityIt;
                typedef typename EntityType::Geometry
                    EntityGeometryType;
                const EntityGeometryType& geometry = entity.geometry();

                // get local matrices for the volume integral
                LocalMInversMatrixType localMInversMatrixElement = MInversMatrix.localMatrix( entity, entity );
                LocalMInversMatrixType localMMatrixElement = MMatrix.localMatrix( entity, entity );
                LocalWmatrixType localWmatrixElement = Wmatrix.localMatrix( entity, entity );
                LocalXmatrixType localXmatrixElement = Xmatrix.localMatrix( entity, entity );
                LocalYmatrixType localYmatrixElement = Ymatrix.localMatrix( entity, entity );
                LocalZmatrixType localZmatrixElement = Zmatrix.localMatrix( entity, entity );
                LocalEmatrixType localEmatrixElement = Ematrix.localMatrix( entity, entity );
                LocalRmatrixType localRmatrixElement = Rmatrix.localMatrix( entity, entity );

                // get local right hand sides
                LocalH1rhsType LocalH1rhs = H1rhs.localFunction( entity );
                LocalH2rhsType LocalH2rhs = H2rhs.localFunction( entity );
                LocalH3rhsType LocalH3rhs = H3rhs.localFunction( entity );

                // get basefunctionsets
                const SigmaBaseFunctionSetType sigmaBaseFunctionSetElement = sigmaSpace_.baseFunctionSet( entity );
                const VelocityBaseFunctionSetType velocityBaseFunctionSetElement = velocitySpace_.baseFunctionSet( entity );
                const PressureBaseFunctionSetType pressureBaseFunctionSetElement = pressureSpace_.baseFunctionSet( entity );
                const int numSigmaBaseFunctionsElement = sigmaBaseFunctionSetElement.numBaseFunctions();
                const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.numBaseFunctions();
                const int numPressureBaseFunctionsElement = pressureBaseFunctionSetElement.numBaseFunctions();

                // get quadrature
                const VolumeQuadratureType volumeQuadratureElement( entity,
                                                                    ( 2 * pressureSpaceOrder ) + 1 );
//                const VolumeQuadratureType volumeQuadratureElement( entity,
//                                                                    3 );
#ifndef NLOG
                if ( numberOfEntities > 19 ) {
                    if ( ( entityNR % fivePercentOfEntities ) == 0 ) {
                        if ( fivePercents < 20 ) {
                            infoStream << "=";
                            ++fivePercents;
                            infoStream.Flush();
                        }
                    }
                }
                debugStream.Suspend(); // disable logging
//                if ( outputEntity == entityNR ) entityOutput = true;
                if ( allOutput ) entityOutput = true;
                if ( entityOutput ) debugStream.Resume(); // enable logging
                debugStream << "  - numSigmaBaseFunctionsElement: " << numSigmaBaseFunctionsElement << std::endl;
                debugStream << "  - numVelocityBaseFunctionsElement: " << numVelocityBaseFunctionsElement << std::endl;
                debugStream << "  - numPressureBaseFunctionsElement: " << numPressureBaseFunctionsElement << std::endl;
                debugStream << "  - == start calculations on entity " << entityNR << std::endl;
                debugStream << "    - entity " << entityNR << " has " << geometry.corners() << " corners:";
                for ( int i = 0; i < geometry.corners(); ++i ) {
                    const VelocityRangeType corner = geometry[i];
                    Stuff::printFieldVector( corner, "corner_"+Stuff::toString(i), debugStream, "      " );
                }
                debugStream << std::endl;
                bool Moutput = false;
                bool Woutput = false;
                bool Xoutput = false;
                bool Youtput = false;
                bool Zoutput = false;
                bool Eoutput = false;
                bool Routput = false;
                bool H1output = false;
                bool H2output = false;
                bool H3output = false;
                // we want logging at the following base functions
                const int logBaseI = 0;
                const int logBaseJ = 0;
                debugStream.Suspend(); // disable logging
#endif
                // compute volume integrals

                //                                                     // we will call this one
                // (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double M_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Moutput = true;
//                        if ( allOutput ) Moutput = true;
                        if ( entityOutput && Moutput ) debugStream.Resume(); // enable logging
                        debugStream << "    = M ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute \tau_{i}:\tau_{j}
                            SigmaRangeType tau_i( 0.0 );
                            SigmaRangeType tau_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                            sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                            const double tau_j_times_tau_i = colonProduct( tau_j, tau_i );
                            // compute M_i_j
                            M_i_j += elementVolume
                                * integrationWeight
                                * tau_j_times_tau_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
                            Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
                            debugStream << "\n        - tau_j_times_tau_i: " << tau_j_times_tau_i << std::endl;
                            debugStream << "        - M_" << i << "_" << j << "+=: " << M_i_j << std::endl;
#endif
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( M_i_j ) < eps ) {
                            M_i_j = 0.0;
localMMatrixElement.add( i, j, M_i_j );
                        } // else invert
                        else {
localMMatrixElement.add( i, j, M_i_j );
                            M_i_j = 1.0 / M_i_j;
                        }
                        // add to matrix
                        localMInversMatrixElement.add( i, j, M_i_j );
#ifndef NLOG
                        Moutput = false;
                        debugStream.Suspend(); // disable logging
#endif
                    }
                } // done computing Minvs' volume integral

                //                                                        // we will call this one
                // (W)_{i,j} += \int_{T}v_{j}\cdot(\nabla\cdot\tau_{i})dx // W's volume integral
                //                                                        // see also "W's entitity surface integral", "W's neighbour surface integral" and "W's boundary integral" below
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double W_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
//                        if ( allOutput ) Woutput = true;
                        if ( entityOutput && Woutput ) debugStream.Resume(); // enable logging
                        debugStream << "    = W ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute \tau_{i}:\tau_{j}
                            SigmaRangeType tau_i( 0.0 );
                            VelocityRangeType v_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                            SigmaJacobianRangeType gradient_of_tau_i( 0.0 );
                            sigmaBaseFunctionSetElement.jacobian( i, x, gradient_of_tau_i );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            const VelocityRangeType divergence_of_tau_i = sigmaDivergenceOutOfGradient( gradient_of_tau_i );
                            const double v_j_times_divergence_of_tau_i = v_j * divergence_of_tau_i;
                            W_i_j += elementVolume
                                * integrationWeight
                                * v_j_times_divergence_of_tau_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
                            Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                            debugStream << "\n        - v_j_times_divergence_of_tau_i: " << v_j_times_divergence_of_tau_i << std::endl;
                            debugStream << "        - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
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
                        debugStream.Suspend(); // disable logging
#endif
                    }
                } // done computing W's volume integral

                //                                                  // we will call this one
                // (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx // X's volume integral
                //                                                  // see also "X's entitity surface integral", "X's neighbour surface integral" and "X's boundary integral" below
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double X_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
//                        if ( allOutput ) Xoutput = true;
                        if ( entityOutput && Xoutput ) debugStream.Resume(); // enable logging
                        debugStream << "    = X ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadrature points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute \tau_{j}:\nabla v_{i}
                            SigmaRangeType gradient_of_v_i( 0.0 );
                            SigmaRangeType tau_j( 0.0 );
                            velocityBaseFunctionSetElement.jacobian( i, x, gradient_of_v_i );
                            sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                            const double tau_j_times_gradient_v_i =
                                colonProduct( tau_j, gradient_of_v_i );
                            X_i_j += elementVolume
                                * integrationWeight
                                * mu
                                * tau_j_times_gradient_v_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldVector( gradient_of_v_i, "gradient of v_i", debugStream, "      " );
                            Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
                            debugStream << "\n        - tau_j_times_gradient_v_i: " << tau_j_times_gradient_v_i << std::endl;
                            debugStream << "        - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
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
                        debugStream.Suspend(); // disable logging
#endif
                    }
                } // done computing X's volume integral

                //                                                  // we will call this one
                // (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx // Z's volume integral
                //                                                  // see also "Z's entitity surface integral", "Z's neighbour surface integral" and "Z's boundary integral" below
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                        double Z_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
//                        if ( allOutput ) Zoutput = true;
                        if ( entityOutput && Zoutput ) debugStream.Resume(); // enable logging
                        debugStream << "    = Z ========================" << std::endl;
                        debugStream << "    basefunctions " << i << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute q_{j}\cdot(\nabla\cdot v_i)
                            SigmaRangeType gradient_of_v_i( 0.0 );
                            PressureRangeType q_j( 0.0 );
                            velocityBaseFunctionSetElement.jacobian( i, x, gradient_of_v_i );
                            const double divergence_of_v_i = velocityDivergenceOutOfGradient( gradient_of_v_i );
                            pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                            const double q_j_times_divergence_of_v_i = q_j * divergence_of_v_i;
                            Z_i_j += -1.0
                                * elementVolume
                                * integrationWeight
                                * q_j_times_divergence_of_v_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldMatrix( gradient_of_v_i, "gradient_of_v_i", debugStream, "      " );
                            Stuff::printFieldVector( q_j, "q_j", debugStream, "      " );
                            debugStream << "\n        - q_j_times_divergence_of_v_i: " << q_j_times_divergence_of_v_i << std::endl;
                            debugStream << "        - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
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
                        debugStream.Suspend(); // disable logging
#endif
                    }
                } // done computing Z's volume integral

                //                                    // we will call this one
                // (H2)_{j} += \int_{T}f\cdot v_{j}dx // H2's volume integral
                //                                    // see also "H2's boundary integral" further down
                if ( discreteModel_.hasForce() ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double H2_j = 0.0;
#ifndef NLOG
    //                    if ( ( j == logBaseJ ) ) H2output = true;
//                        if ( allOutput ) H2output = true;
                        if ( entityOutput && H2output ) debugStream.Resume(); // enable logging
                        debugStream << "    = H2 =======================" << std::endl;
                        debugStream << "    basefunction " << " " << j << std::endl;
                        debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                        // sum over all quadratur points
                        for ( int quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            const VelocityRangeType xWorld = geometry.global( x );

                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute f\cdot v_j
                            VelocityRangeType v_j( 0.0 );
                            VelocityRangeType f( 0.0 );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            discreteModel_.force( 0.0, xWorld, f );
                            const double f_times_v_j = f * v_j;
                            H2_j += elementVolume
                                * integrationWeight
                                * f_times_v_j;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            Stuff::printFieldVector( xWorld, "xWorld", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldVector( f, "f", debugStream, "      " );
                            Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                            debugStream << "\n        - f_times_v_j: " << f_times_v_j << std::endl;
                            debugStream << "        - H2_" << j << "+=: " << H2_j << std::endl;
#endif
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( H2_j ) < eps ) {
                            H2_j = 0.0;
                        }
                        // add to rhs
                        LocalH2rhs[ j ] += H2_j;
#ifndef NLOG
                        H2output = false;
                        debugStream.Suspend(); // disable logging
#endif
                    } // done computing H2's volume integral
                }

                //                                                // we will call this one
                // (E)_{i,j} += -\int_{T}v_{j}\cdot\nabla q_{i}dx // E's volume integral
                //                                                // see also "E's entitity surface integral", "E's neighbour surface integral" and "E's boundary integral" below
                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double E_i_j = 0.0;
#ifndef NLOG
//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
//                        if ( allOutput ) Eoutput = true;
                        if ( entityOutput && Eoutput ) debugStream.Resume(); // enable logging
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
                            // compute v_{j}\cdot(\nabla q_i)
                            typename DiscretePressureFunctionSpaceType::JacobianRangeType jacobian_of_q_i( 0.0 );
                            VelocityRangeType v_j( 0.0 );
                            pressureBaseFunctionSetElement.jacobian( i, x, jacobian_of_q_i );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            VelocityRangeType gradient_of_q_i( jacobian_of_q_i[0] );
                            double v_j_times_gradient_of_q_i = v_j * gradient_of_q_i;
                            E_i_j += -1.0
                                * elementVolume
                                * integrationWeight
                                * v_j_times_gradient_of_q_i;
#ifndef NLOG
                            debugStream << "    - quadPoint " << quad;
                            Stuff::printFieldVector( x, "x", debugStream, "      " );
                            debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                            debugStream << "        - integrationWeight: " << integrationWeight;
                            Stuff::printFieldVector( gradient_of_q_i, "gradient_of_q_i", debugStream, "      " );
                            Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
                            debugStream << "\n        - v_j_times_gradient_of_q_i: " << v_j_times_gradient_of_q_i << std::endl;
                            debugStream << "        - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
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
                        debugStream.Suspend(); // disable logging
#endif
                    }
                } // done computing E's volume integral

                // walk the intersections
                IntersectionIteratorType intItEnd = gridPart_.iend( entity );
                for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
                        intIt != intItEnd;
                        ++intIt ) {
#ifndef NLOG
//                    if ( ( outputIntersection == intersectionNR ) && entityOutput ) intersectionOutput = true;
                    if ( entityOutput ) intersectionOutput = true;
                    if ( intersectionOutput ) debugStream.Resume(); // enable logging
                    debugStream << "    - ==== start calculations on intersection " << intersectionNR << std::endl;
#endif

                    // get intersection geometry
                    typedef typename IntersectionIteratorType::Geometry
                        IntersectionGeometryType;
//                    const IntersectionGeometryType& intersectionGeoemtry = intIt.intersectionSelfLocal();
                    const IntersectionGeometryType& intersectionGeoemtry = intIt.intersectionGlobal();
#ifndef NLOG
                    // get corners
                    debugStream << "      - intersection " << intersectionNR << " has " << intersectionGeoemtry.corners() << " corners:";
                    for ( int i = 0; i < intersectionGeoemtry.corners(); ++i ) {
                        const VelocityRangeType corner = intersectionGeoemtry[i];
                        Stuff::printFieldVector( corner, "corner_"+Stuff::toString(i), debugStream, "        " );
                    }
                    debugStream << "\n        length of intersection " << intersectionNR << " is " << Stuff::getLenghtOfIntersection( intIt ) << std::endl;
                    debugStream.Suspend(); // disable logging
#endif

                    // get intersection quadrature, seen from inside
                    const FaceQuadratureType faceQuadratureElement( gridPart_,
                                                                    intIt,
                                                                    ( 2 * pressureSpaceOrder ) + 1,
                                                                    FaceQuadratureType::INSIDE );
//                    const FaceQuadratureType faceQuadratureElement( gridPart_,
//                                                                    intIt,
//                                                                    3,
//                                                                    FaceQuadratureType::INSIDE );

                    // if we are inside the grid
                    if ( intIt.neighbor() && !intIt.boundary() ) {
                        // get neighbour
                        const typename IntersectionIteratorType::EntityPointer neighbourPtr = intIt.outside();
                        const EntityType& neighbour = *neighbourPtr;

                        // get local matrices for the surface integrals
                        LocalMInversMatrixType localMInversMatrixNeighbour = MInversMatrix.localMatrix( entity, neighbour );
                        LocalWmatrixType localWmatrixNeighbour = Wmatrix.localMatrix( entity, neighbour );
                        LocalXmatrixType localXmatrixNeighbour = Xmatrix.localMatrix( entity, neighbour );
                        LocalYmatrixType localYmatrixNeighbour = Ymatrix.localMatrix( entity, neighbour );
                        LocalZmatrixType localZmatrixNeighbour = Zmatrix.localMatrix( entity, neighbour );
                        LocalEmatrixType localEmatrixNeighbour = Ematrix.localMatrix( entity, neighbour );
                        LocalRmatrixType localRmatrixNeighbour = Rmatrix.localMatrix( entity, neighbour );

                        // get basefunctionsets
                        const SigmaBaseFunctionSetType sigmaBaseFunctionSetNeighbour = sigmaSpace_.baseFunctionSet( neighbour );
                        const VelocityBaseFunctionSetType velocityBaseFunctionSetNeighbour = velocitySpace_.baseFunctionSet( neighbour );
                        const PressureBaseFunctionSetType pressureBaseFunctionSetNeighbour = pressureSpace_.baseFunctionSet( neighbour );
                        const int numSigmaBaseFunctionsNeighbour = sigmaBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numVelocityBaseFunctionsNeighbour = velocityBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numPressureBaseFunctionsNeighbour = pressureBaseFunctionSetNeighbour.numBaseFunctions();

                        // get intersection quadrature, seen from outside
                        const FaceQuadratureType faceQuadratureNeighbour(   gridPart_,
                                                                            intIt,
                                                                            ( 2 * pressureSpaceOrder ) + 1,
                                                                            FaceQuadratureType::OUTSIDE );
//                        const FaceQuadratureType faceQuadratureNeighbour(   gridPart_,
//                                                                            intIt,
//                                                                            3,
//                                                                            FaceQuadratureType::OUTSIDE );

                        // compute the surface integrals

                        //                                                                                                               // we will call this one
                        // (W)_{i,j} += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's element surface integral
                        //           += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's neighbour surface integral
                        //                                                                                                               // see also "W's boundary integral" below
                        //                                                                                                               // and "W's volume integral" above
                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
                            for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                                // compute W's element surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double W_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
    //                                if ( allOutput ) Woutput = true;
                                    if ( intersectionOutput && Woutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = W element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x in codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
                                        SigmaRangeType tau_i( 0.0 );
                                        VelocityRangeType v_j( 0.0 );
                                        VelocityRangeType u_sigma_u_plus_flux( 0.0 );
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
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
                                        W_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * flux_times_tau_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_sigma_u_plus_flux, "u_sigma_u_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( tau_i_times_n_t, "tau_i_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - flux_times_tau_i_times_n_t: " << flux_times_tau_i_times_n_t << std::endl;
                                        debugStream << "          - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
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
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing W's element surface integral
                                // compute W's neighbour surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsNeighbour; ++j ) {
                                    double W_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
    //                                if ( allOutput ) Woutput = true;
                                    if ( intersectionOutput && Woutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = W neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x in codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{j}\cdot n_{T}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        SigmaRangeType tau_i( 0.0 );
                                        VelocityRangeType v_j( 0.0 );
                                        VelocityRangeType u_sigma_u_plus_flux( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                                        velocityBaseFunctionSetNeighbour.evaluate( j, x, v_j );
                                        discreteModel_.velocitySigmaFlux(   intIt,
                                                                            0.0,
                                                                            localX,
                                                                            DiscreteModelType::outside,
                                                                            v_j,
                                                                            u_sigma_u_plus_flux );
                                        VelocityRangeType tau_i_times_n_t( 0.0 );
                                        tau_i.mv( outerNormal, tau_i_times_n_t );
                                        const double flux_times_tau_i_times_n_t = u_sigma_u_plus_flux * tau_i_times_n_t;
                                        W_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * flux_times_tau_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_sigma_u_plus_flux, "u_sigma_u_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( tau_i_times_n_t, "tau_i_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - flux_times_tau_i_times_n_t: " << flux_times_tau_i_times_n_t << std::endl;
                                        debugStream << "          - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( W_i_j ) < eps ) {
                                        W_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localWmatrixNeighbour.add( i, j, W_i_j );
#ifndef NLOG
                                    Woutput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing W's neighbour surface integral
                            } // done computing W's surface integrals
                        }

                        //                                                                                                                   // we will call this one
                        // (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's element sourface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}ds // X's neighbour sourface integral
                        //                                                                                                                   // see also "X's boundary integral" below
                        //                                                                                                                   // and "X's volume integral" above
                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                // compute X's element sourface integral
                                for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                    double X_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
    //                                if ( allOutput ) Xoutput = true;
                                    if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = X element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      volumeQuadratureElement.nop() " << volumeQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x in codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                                        SigmaRangeType sigma_sigma_plus_flux( 0.0 );
                                        discreteModel_.sigmaFlux(   intIt,
                                                                    0.0,
                                                                    localX,
                                                                    DiscreteModelType::inside,
                                                                    tau_j,
                                                                    sigma_sigma_plus_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_sigma_plus_flux.mv( outerNormal, flux_times_n_t );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_flux_times_n_t = v_i * flux_times_n_t;
                                        X_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_i_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_sigma_plus_flux, "sigma_sigma_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_flux_times_n_t: " << v_i_times_flux_times_n_t << std::endl;
                                        debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( X_i_j ) < eps ) {
                                        X_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localXmatrixElement.add( i, j, X_i_j );
#ifndef NLOG
                                    Xoutput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing X's element sourface integral
                                // compute X's neighbour sourface integral
                                for ( int j = 0; j < numSigmaBaseFunctionsNeighbour; ++j ) {
                                    double X_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
    //                                if ( allOutput ) Xoutput = true;
                                    if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = X neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetNeighbour.evaluate( j, x, tau_j );
                                        SigmaRangeType sigma_sigma_minus_flux( 0.0 );
                                        discreteModel_.sigmaFlux(   intIt,
                                                                    0.0,
                                                                    localX,
                                                                    DiscreteModelType::outside,
                                                                    tau_j,
                                                                    sigma_sigma_minus_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_sigma_minus_flux.mv( outerNormal, flux_times_n_t );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_flux_times_n_t = v_i * flux_times_n_t;
                                        X_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_i_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_sigma_minus_flux, "sigma_sigma_minus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_flux_times_n_t: " << v_i_times_flux_times_n_t << std::endl;
                                        debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( X_i_j ) < eps ) {
                                        X_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localXmatrixNeighbour.add( i, j, X_i_j );
#ifndef NLOG
                                    Xoutput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing X's neighbour sourface integral
                            } // done computing X's sourface integrals
                        }

                        //                                                                                                         // we call this one
                        // (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}ds // Y's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}ds // Y's neighbour surface integral
                        //                                                                                                         // see also "Y's boundary integral" below
                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                // compute Y's element surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double Y_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Youtput = true;
    //                                if ( allOutput ) Youtput = true;
                                    if ( intersectionOutput && Youtput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = Y element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        SigmaRangeType sigma_u_plus_flux( 0.0 );
                                        discreteModel_.sigmaFlux(   intIt,
                                                                    0.0,
                                                                    localX,
                                                                    DiscreteModelType::inside,
                                                                    v_j,
                                                                    sigma_u_plus_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_u_plus_flux.mv( outerNormal, flux_times_n_t );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_flux_times_n_t = v_i * flux_times_n_t;
                                        Y_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_i_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_u_plus_flux, "sigma_u_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_flux_times_n_t: " << v_i_times_flux_times_n_t << std::endl;
                                        debugStream << "          - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Y_i_j ) < eps ) {
                                        Y_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localYmatrixElement.add( i, j, Y_i_j );
#ifndef NLOG
                                    Youtput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing Y's element surface integral
                                // compute Y's neighbour surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsNeighbour; ++j ) {
                                    double Y_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Youtput = true;
    //                                if ( allOutput ) Youtput = true;
                                    if ( intersectionOutput && Youtput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = Y neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetNeighbour.evaluate( j, x, v_j );
                                        SigmaRangeType sigma_u_minus_flux( 0.0 );
                                        discreteModel_.sigmaFlux(   intIt,
                                                                    0.0,
                                                                    localX,
                                                                    DiscreteModelType::outside,
                                                                    v_j,
                                                                    sigma_u_minus_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_u_minus_flux.mv( outerNormal, flux_times_n_t );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_flux_times_n_t = v_i * flux_times_n_t;
                                        Y_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_i_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_u_minus_flux, "sigma_u_minus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_flux_times_n_t: " << v_i_times_flux_times_n_t << std::endl;
                                        debugStream << "          - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Y_i_j ) < eps ) {
                                        Y_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localYmatrixNeighbour.add( i, j, Y_i_j );
#ifndef NLOG
                                    Youtput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing Y's neighbour surface integral
                            } // done computing Y's surface integrals
                        }

                        //                                                                                                  // we will call this one
                        // (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{-}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's neighbour surface integral
                        //                                                                                                  // see also "Z's boundary integral" below
                        //                                                                                                  // and "Z's volume integral" above
                        if ( discreteModel_.hasPressureFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                // compute Z's element surface integral
                                for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                    double Z_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
    //                                if ( allOutput ) Zoutput = true;
                                    if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = Z element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_n_t = v_i * outerNormal;
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        PressureRangeType p_p_plus_flux( 0.0 );
                                        discreteModel_.pressureFlux(    intIt,
                                                                        0.0,
                                                                        localX,
                                                                        DiscreteModelType::inside,
                                                                        q_j,
                                                                        p_p_plus_flux );
                                        Z_i_j += elementVolume
                                            * integrationWeight
                                            * p_p_plus_flux
                                            * v_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_n_t: " << v_i_times_n_t << std::endl;
                                        debugStream << "\n          - p_p_plus_flux: " << p_p_plus_flux << std::endl;
                                        debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
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
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing Z's element surface integral
                                // compute Z's neighbour surface integral
                                for ( int j = 0; j < numPressureBaseFunctionsNeighbour; ++j ) {
                                    double Z_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
    //                                if ( allOutput ) Zoutput = true;
                                    if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = Z neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_n_t = v_i * outerNormal;
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetNeighbour.evaluate( j, x, q_j );
                                        PressureRangeType p_p_minus_flux( 0.0 );
                                        discreteModel_.pressureFlux(    intIt,
                                                                        0.0,
                                                                        localX,
                                                                        DiscreteModelType::outside,
                                                                        q_j,
                                                                        p_p_minus_flux );
                                        Z_i_j += elementVolume
                                            * integrationWeight
                                            * p_p_minus_flux
                                            * v_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_n_t: " << v_i_times_n_t << std::endl;
                                        debugStream << "\n          - p_p_minus_flux: " << p_p_minus_flux << std::endl;
                                        debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Z_i_j ) < eps ) {
                                        Z_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localZmatrixNeighbour.add( i, j, Z_i_j );
#ifndef NLOG
                                    Zoutput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing Z's neighbour surface integral
                            } // done computing Z's surface integrals
                        }

                        //                                                                                                // we will call this one
                        // (E)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}ds // E's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}ds // E's neighbour surface integral
                        //                                                                                                // see also "E's boundary integral" below
                        //                                                                                                // and "E's volume integral" above
                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                // compute E's element surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double E_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
    //                                if ( allOutput ) Eoutput = true;
                                    if ( intersectionOutput && Eoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = E element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        VelocityRangeType u_p_u_plus_flux( 0.0 );
                                        discreteModel_.velocityPressureFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                DiscreteModelType::inside,
                                                                                v_j,
                                                                                u_p_u_plus_flux );
                                        const double flux_times_n_t = u_p_u_plus_flux * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        E_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_u_plus_flux, "u_p_u_plus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
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
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing E's element surface integral
                                // compute E's neighbour surface integral
                                for ( int j = 0; j < numVelocityBaseFunctionsNeighbour; ++j ) {
                                    double E_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
    //                                if ( allOutput ) Eoutput = true;
                                    if ( intersectionOutput && Eoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = E neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetNeighbour.evaluate( j, x, v_j );
                                        VelocityRangeType u_p_u_minus_flux( 0.0 );
                                        discreteModel_.velocityPressureFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                DiscreteModelType::outside,
                                                                                v_j,
                                                                                u_p_u_minus_flux );
                                        const double flux_times_n_t = u_p_u_minus_flux * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        E_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_u_minus_flux, "u_p_u_minus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( E_i_j ) < eps ) {
                                        E_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localEmatrixNeighbour.add( i, j, E_i_j );
#ifndef NLOG
                                    Eoutput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing E's neighbour surface integral
                            } // done computing E's surface integrals
                        }

                        //                                                                                                // we will call this one
                        // (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}ds // R's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}ds // R's neighbour surface integral
                        //                                                                                                // see also "R's boundary integral" below
                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                // compute R's element surface integral
                                for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                    double R_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Routput = true;
    //                                if ( allOutput ) Routput = true;
                                    if ( intersectionOutput && Routput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = R element ======================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        VelocityRangeType u_p_p_plus_flux( 0.0 );
                                        discreteModel_.velocityPressureFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                DiscreteModelType::inside,
                                                                                q_j,
                                                                                u_p_p_plus_flux );
                                        const double flux_times_n_t = u_p_p_plus_flux
                                            * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        R_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_p_plus_flux, "u_p_p_plus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - R_" << i << "_" << j << "+=: " << R_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( R_i_j ) < eps ) {
                                        R_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localRmatrixElement.add( i, j, R_i_j );
#ifndef NLOG
                                    Routput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing R's element surface integral
                                // compute R's neighbour surface integral
                                for ( int j = 0; j < numPressureBaseFunctionsNeighbour; ++j ) {
                                    double R_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Routput = true;
    //                                if ( allOutput ) Routput = true;
                                    if ( intersectionOutput && Routput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = R neighbour ====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetNeighbour.evaluate( j, x, q_j );
                                        VelocityRangeType u_p_p_minus_flux( 0.0 );
                                        discreteModel_.velocityPressureFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                DiscreteModelType::outside,
                                                                                q_j,
                                                                                u_p_p_minus_flux );
                                        const double flux_times_n_t = u_p_p_minus_flux
                                            * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        R_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_p_minus_flux, "u_p_p_minus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - R_" << i << "_" << j << "+=: " << R_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( R_i_j ) < eps ) {
                                        R_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localRmatrixNeighbour.add( i, j, R_i_j );
#ifndef NLOG
                                    Routput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                } // done computing R's neighbour surface integral
                            } // done computing R's surface integrals
                        }

                    } // done with those inside the grid

                    // if we are on the boundary of the grid
                    if ( !intIt.neighbor() && intIt.boundary() ) {
                        // compute the boundary integrals

                        //                                                                                                               // we wil call this one
                        // (W)_{i,j} += \int_{\varepsilon\in \Epsilon_{D}^{T}}-\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's boundary integral
                        //                                                                                                               // see also "W's volume integral", "W's element surface integral" and "W's neighbour surface integral" above
                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
                            for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double W_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
    //                                if ( allOutput ) Woutput = true;
                                    if ( intersectionOutput && Woutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = W boundary =====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        SigmaRangeType tau_i( 0.0 );
                                        VelocityRangeType v_j( 0.0 );
                                        VelocityRangeType u_sigma_u_plus_flux( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        discreteModel_.velocitySigmaBoundaryFlux(   intIt,
                                                                                    0.0,
                                                                                    localX,
                                                                                    v_j,
                                                                                    u_sigma_u_plus_flux );
                                        VelocityRangeType tau_i_times_n_t( 0.0 );
                                        tau_i.mv( outerNormal, tau_i_times_n_t );
                                        const double flux_times_tau_i_times_n_t = u_sigma_u_plus_flux * tau_i_times_n_t;
                                        W_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * flux_times_tau_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_sigma_u_plus_flux, "u_sigma_u_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( tau_i_times_n_t, "tau_i_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - flux_times_tau_i_times_n_t: " << flux_times_tau_i_times_n_t << std::endl;
                                        debugStream << "          - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
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
                                }
                            } // done computing W's boundary integral
                        }

                        //                                                                                                    // we will call this one
                        // (H1)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}ds // H1's boundary integral
                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
                            for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                double H1_j = 0.0;
#ifndef NLOG
    //                            if ( j == logBaseJ ) H1output = true;
//                                if ( allOutput ) H1output = true;
                                if ( intersectionOutput && H1output ) debugStream.Resume(); // enable logging
                                debugStream << "      = H1 boundary ====================" << std::endl;
                                debugStream << "      basefunction " << j << std::endl;
                                debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                // sum over all quadrature points
                                for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const VelocityRangeType globalX = geometry.global( x );
                                    const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // compute \hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}
                                    const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                    SigmaRangeType tau_j( 0.0 );
                                    sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                                    VelocityRangeType tau_j_times_n_t( 0.0 );
                                    tau_j.mv( outerNormal, tau_j_times_n_t );
                                    VelocityRangeType u_sigma_rhs_flux( 0.0 );
                                    discreteModel_. velocitySigmaBoundaryFlux(  intIt,
                                                                                0.0,
                                                                                localX,
                                                                                u_sigma_rhs_flux );
                                    const double flux_times_tau_j_times_n_t = u_sigma_rhs_flux * tau_j_times_n_t;
                                    H1_j += elementVolume
                                        * integrationWeight
                                        * flux_times_tau_j_times_n_t;
#ifndef NLOG
                                    debugStream << "      - quadPoint " << quad;
                                    Stuff::printFieldVector( x, "x", debugStream, "        " );
                                    Stuff::printFieldVector( localX, "localX", debugStream, "        " );
                                    Stuff::printFieldVector( globalX, "globalX", debugStream, "        " );
                                    Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
                                    debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                    debugStream << "        - integrationWeight: " << integrationWeight;
                                    Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
                                    Stuff::printFieldVector( u_sigma_rhs_flux, "u_sigma_rhs_flux", debugStream, "        " );
                                    Stuff::printFieldVector( tau_j_times_n_t, "tau_j_times_n_t", debugStream, "        " );
                                    debugStream << "\n        - flux_times_tau_j_times_n_t: " << flux_times_tau_j_times_n_t << std::endl;
                                    debugStream << "        - H1_" << j << "+=: " << H1_j << std::endl;
#endif
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H1_j ) < eps ) {
                                    H1_j = 0.0;
                                }
                                // add to rhs
                                LocalH1rhs[ j ] += H1_j;
#ifndef NLOG
                                H1output = false;
                                debugStream.Suspend(); // disable logging
#endif
                            } // done computing H1's boundary integral
                        }

                        //                                                                                                                   // we will call this one
                        // (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's boundary integral
                        //                                                                                                                   // see also "X's volume integral", "X's element surface integral" and "X's neighbour surface integral" above
                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                    double X_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
    //                                if ( allOutput ) Xoutput = true;
                                    if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = X boundary =====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                                        SigmaRangeType sigma_sigma_plus_flux( 0.0 );
                                        discreteModel_.sigmaBoundaryFlux(   intIt,
                                                                            0.0,
                                                                            localX,
                                                                            tau_j,
                                                                            sigma_sigma_plus_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_sigma_plus_flux.mv( outerNormal, flux_times_n_t );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_flux_times_n_t = v_i * flux_times_n_t;
                                        X_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_i_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_sigma_plus_flux, "sigma_sigma_plus_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_flux_times_n_t: " << v_i_times_flux_times_n_t << std::endl;
                                        debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
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
                            } // done computing X's boundary integral
                        }

                        //                                                                                                  // we will call this one
                        // (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's boundary integral
                        //                                                                                                  // see also "Z's volume integral", "Z's element surface integral" and "Z's neighbour surface integral" above
                        if ( discreteModel_.hasPressureFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                // compute the boundary integral
                                for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                    double Z_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
    //                                if ( allOutput ) Zoutput = true;
                                    if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = Z boundary =====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_n_t = v_i * outerNormal;
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        PressureRangeType p_p_plus_flux( 0.0 );
                                        discreteModel_.pressureBoundaryFlux(    intIt,
                                                                        0.0,
                                                                        localX,
                                                                        q_j,
                                                                        p_p_plus_flux );
                                        Z_i_j += elementVolume
                                            * integrationWeight
                                            * p_p_plus_flux
                                            * v_i_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        debugStream << "\n          - v_i_times_n_t: " << v_i_times_n_t << std::endl;
                                        debugStream << "\n          - p_p_plus_flux: " << p_p_plus_flux << std::endl;
                                        debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
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
                                    debugStream.Suspend(); // disable logging
#endif
                                }
                            } // done computing Z's boundary integral
                        }

                        //                                                                                                                 // we will call this one
                        // (H2)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left( \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}ds         // H2's 1st boundary integral
                        //                                                         -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}ds        \right) // H2's 2nd boundary integral
                        //                                                                                                                 // see also "H2's volume integral" above
                        if ( discreteModel_.hasSigmaFlux() && discreteModel_.hasPressureFlux() ) {
                            for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                double H2_j = 0.0;
#ifndef NLOG
    //                            if ( j == logBaseJ ) H2output = true;
//                                if ( allOutput ) H2output = true;
                                if ( intersectionOutput && H2output ) debugStream.Resume(); // enable logging
                                debugStream << "      = H2 boundary ====================" << std::endl;
                                debugStream << "      basefunction " << j << std::endl;
                                debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                // sum over all quadrature points
                                for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                    const VelocityRangeType globalX = geometry.global( x );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // prepare
                                    const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                    VelocityRangeType v_j( 0.0 );
                                    velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                    // compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
                                    if ( discreteModel_.hasSigmaFlux() ) {
                                        SigmaRangeType sigma_rhs_flux( 0.0 );
                                        discreteModel_.sigmaBoundaryFlux(   intIt,
                                                                            0.0,
                                                                            localX,
                                                                            sigma_rhs_flux );
                                        VelocityRangeType flux_times_n_t( 0.0 );
                                        sigma_rhs_flux.mv( outerNormal, flux_times_n_t );
                                        const double v_j_times_flux_times_n_t = v_j * flux_times_n_t;
                                        H2_j += elementVolume
                                            * integrationWeight
                                            * mu
                                            * v_j_times_flux_times_n_t;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "        " );
                                        Stuff::printFieldVector( globalX, "globalX", debugStream, "        " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldMatrix( sigma_rhs_flux, "sigma_rhs_flux", debugStream, "        " );
                                        Stuff::printFieldVector( flux_times_n_t, "flux_times_n_t", debugStream, "        " );
                                        debugStream << "\n        - v_j_times_flux_times_n_t: " << v_j_times_flux_times_n_t << std::endl;
                                        debugStream << "        - H2_" << j << "+=: " << H2_j;
#endif
                                    }
                                    // done computing H2's 1st boundary integral
                                    // compute -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}
                                    if ( discreteModel_.hasPressureFlux() ) {
                                        const double v_j_times_n_t = v_j * outerNormal;
                                        PressureRangeType p_rhs_flux( 0.0 );
                                        discreteModel_.pressureBoundaryFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                p_rhs_flux );
                                        const double flux_times_v_j_times_n_t = p_rhs_flux * v_j_times_n_t;
                                        H2_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * flux_times_v_j_times_n_t;
#ifndef NLOG
                                        Stuff::printFieldVector( p_rhs_flux, "p_rhs_flux", debugStream, "        " );
                                        debugStream << "\n        - flux_times_v_j_times_n_t: " << flux_times_v_j_times_n_t << std::endl;
                                        debugStream << "        - H2_" << j << "+=: " << H2_j << std::endl;
#endif
                                    }
                                    // done computing H2's 2nd boundary integral
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H2_j ) < eps ) {
                                    H2_j = 0.0;
                                }
                                // add to rhs
                                LocalH2rhs[ j ] += H2_j;
#ifndef NLOG
                                H2output = false;
                                debugStream.Suspend(); // disable logging
#endif
                            } // done computing H2's boundary integrals
                        }

                        //                                                                                               // we will call this one
                        // (E)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{p}^{U^{+}}(v_{j}\cdot n_{T}q_{i}ds // E's boundary integral
                        //                                                                                               // see also "E's volume integral", "E's element surface integral" and "E's neighbour surface integral" above
                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                // compute the boundary integral
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double E_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
    //                                if ( allOutput ) Eoutput = true;
                                    if ( intersectionOutput && Eoutput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = E boundary =====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        VelocityRangeType u_p_u_plus_flux( 0.0 );
                                        discreteModel_.velocityPressureBoundaryFlux(    intIt,
                                                                                0.0,
                                                                                localX,
                                                                                v_j,
                                                                                u_p_u_plus_flux );
                                        const double flux_times_n_t = u_p_u_plus_flux * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        E_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_u_plus_flux, "u_p_u_plus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
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
                                    debugStream.Suspend(); // disable logging
#endif
                                }
                            } // done computing E's boundary integral
                        }

                        //                                                                                                // we call this one
                        // (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{t}q_{i}ds // R's boundary integral
                        //                                                                                                // see also "R's element surface integral" and "R's neighbour surface integral" above
                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                    double R_i_j = 0.0;
#ifndef NLOG
    //                                if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Routput = true;
    //                                if ( allOutput ) Routput = true;
                                    if ( intersectionOutput && Routput ) debugStream.Resume(); // enable logging
                                    debugStream << "      = R boundary =====================" << std::endl;
                                    debugStream << "      basefunctions " << i << " " << j << std::endl;
                                    debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                    // sum over all quadrature points
                                    for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
                                        const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        VelocityRangeType u_p_p_plus_flux( 0.0 );
                                        discreteModel_.velocityPressureBoundaryFlux(    intIt,
                                                                                        0.0,
                                                                                        localX,
                                                                                        q_j,
                                                                                        u_p_p_plus_flux );
                                        const double flux_times_n_t = u_p_p_plus_flux
                                            * outerNormal;
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        const double flux_times_n_t_times_q_i = q_i * flux_times_n_t;
                                        R_i_j += elementVolume
                                            * integrationWeight
                                            * flux_times_n_t_times_q_i;
#ifndef NLOG
                                        debugStream << "      - quadPoint " << quad;
                                        Stuff::printFieldVector( x, "x", debugStream, "        " );
                                        Stuff::printFieldVector( localX, "localX", debugStream, "          " );
                                        Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
                                        debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                        debugStream << "        - integrationWeight: " << integrationWeight;
                                        Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
                                        Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                        Stuff::printFieldVector( u_p_p_plus_flux, "u_p_p_plus_flux", debugStream, "        " );
                                        debugStream << "\n          - flux_times_n_t: " << flux_times_n_t << std::endl;
                                        debugStream << "\n          - flux_times_n_t_times_q_i: " << flux_times_n_t_times_q_i << std::endl;
                                        debugStream << "          - R_" << i << "_" << j << "+=: " << R_i_j << std::endl;
#endif
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( R_i_j ) < eps ) {
                                        R_i_j = 0.0;
                                    }
                                    // add to matrix
                                    localRmatrixElement.add( i, j, R_i_j );
#ifndef NLOG
                                    Routput = false;
                                    debugStream.Suspend(); // disable logging
#endif
                                }
                            } // done computing R's boundary integral
                        }

                        //                                                                                        // we will call this one
                        // (H3)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}-\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}ds // H3's boundary integral
                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                double H3_j = 0.0;
#ifndef NLOG
    //                            if ( j == logBaseJ ) H3output = true;
//                                if ( allOutput ) H3output = true;
                                if ( intersectionOutput && H3output ) debugStream.Resume(); // enable logging
                                debugStream << "      = H3 boundary ====================" << std::endl;
                                debugStream << "      basefunction " << j << std::endl;
                                debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
#endif
                                // sum over all quadrature points
                                for ( int quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const LocalIntersectionCoordinateType localX = faceQuadratureElement.localPoint( quad );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeoemtry.integrationElement( localX );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // compute -\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}
                                    const VelocityRangeType outerNormal = intIt.unitOuterNormal( localX );
                                    VelocityRangeType u_p_rhs_flux( 0.0 );
                                    discreteModel_.velocityPressureBoundaryFlux( intIt,
                                                                                0.0,
                                                                                localX,
                                                                                u_p_rhs_flux );
                                    const double flux_times_n_t = u_p_rhs_flux
                                        * outerNormal;
                                    PressureRangeType q_j( 0.0 );
                                    pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                    const double flux_times_n_t_times_q_j = q_j
                                        * flux_times_n_t;
                                    H3_j += elementVolume
                                        * integrationWeight
                                        * flux_times_n_t_times_q_j;
#ifndef NLOG
                                    debugStream << "      - quadPoint " << quad;
                                    Stuff::printFieldVector( x, "x", debugStream, "        " );
                                    Stuff::printFieldVector( localX, "localX", debugStream, "        " );
                                    Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
                                    debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
                                    debugStream << "        - integrationWeight: " << integrationWeight;
                                    Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
                                    Stuff::printFieldVector( u_p_rhs_flux, "u_p_rhs_flux", debugStream, "        " );
                                    debugStream << "\n        - flux_times_n_t: " << flux_times_n_t << std::endl;
                                    debugStream << "\n        - flux_times_n_t_times_q_j: " << flux_times_n_t_times_q_j << std::endl;
                                    debugStream << "        - H3_" << j << "+=: " << H3_j << std::endl;
#endif
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H3_j ) < eps ) {
                                    H3_j = 0.0;
                                }
                                // add to rhs
                                LocalH3rhs[ j ] += H3_j;
#ifndef NLOG
                                H3output = false;
                                debugStream.Suspend(); // disable logging
#endif
                            } // done computing H3's boundary integral
                        }

                    } // done with those on the boundary
#ifndef NLOG
                    if ( intersectionOutput ) debugStream.Resume(); // enable logging
                    debugStream << "    - ==== done calculations on intersection " << intersectionNR << std::endl;
                    debugStream.Suspend(); // disable logging
                    intersectionOutput = false;
                    ++intersectionNR;
#endif
                } // done walking the neighbours

#ifndef NLOG
                intersectionNR = 0;
                if ( entityOutput ) debugStream.Resume(); // enable logging
                debugStream << "  - == done calculations on entity " << entityNR << std::endl;
                debugStream.Suspend(); // disable logging
                entityOutput = false;
                ++entityNR;
#endif
            } // done walking the grid
#ifndef NLOG

#ifdef CHEAT
            // compute the artificial right hand sides, should be the right ones
            // H1
//            debugStream = Logger().Info();

            debugStream.Resume();
            debugStream << "  - computing artificial right hand sides" << std::endl;
            debugStream << "    - H1" << std::endl;
            Wmatrix.matrix().apply( discreteExactVelocity, tmpH1rhs );
            exactH1rhs += tmpH1rhs;
            tmpH1rhs.clear();
            MMatrix.matrix().apply( computedVelocityGradient, tmpH1rhs );
            exactH1rhs += tmpH1rhs;
            tmpH1rhs.clear();
            // H2
            debugStream << "    - H2" << std::endl;
            Zmatrix.matrix().apply( discreteExactPressure, tmpH2rhs );
            exactH2rhs += tmpH2rhs;
            tmpH2rhs.clear();
            Ymatrix.matrix().apply( discreteExactVelocity, tmpH2rhs );
            exactH2rhs += tmpH2rhs;
            tmpH2rhs.clear();
            Xmatrix.matrix().apply( computedVelocityGradient, tmpH2rhs );
            exactH2rhs += tmpH2rhs;
            tmpH2rhs.clear();
            // H3
            debugStream << "    - H3" << std::endl;
            Rmatrix.matrix().apply( discreteExactPressure, tmpH3rhs );
            exactH3rhs += tmpH3rhs;
            tmpH3rhs.clear();
            Ematrix.matrix().apply( discreteExactVelocity, tmpH3rhs );
            exactH3rhs += tmpH3rhs;
            tmpH3rhs.clear();
            debugStream << "  - done computing artificial rihgt hand sides" << std::endl;
#endif


            if ( numberOfEntities > 19 ) {
                infoStream << "]";
            }
            infoStream << "\n- gridwalk done" << std::endl << std::endl;

            if ( Mprint || Wprint || Xprint || Yprint || Zprint || Eprint || Rprint || H1print || H2print || H3print ) {
                debugStream.Resume();
                debugStream << "- printing matrices" << std::endl;
                if ( Mprint ) {
                    debugStream << " - = Minvers ======" << std::endl;
                    debugStream.Log( &MInversMatrixType::MatrixType::print,  MInversMatrix.matrix() );
                }
if ( Mprint ) {
    debugStream << " - = M ============" << std::endl;
    debugStream.Log( &MInversMatrixType::MatrixType::print,  MMatrix.matrix() );
}
                if ( Wprint ) {
                    debugStream << " - = W ============" << std::endl;
                    debugStream.Log( &WmatrixType::MatrixType::print,  Wmatrix.matrix() );
                }
                if ( Xprint ) {
                    debugStream << " - = X ============" << std::endl;
                    debugStream.Log( &XmatrixType::MatrixType::print,  Xmatrix.matrix() );
                }
                if ( Yprint ) {
                    debugStream << " - = Y ============" << std::endl;
                    debugStream.Log( &YmatrixType::MatrixType::print,  Ymatrix.matrix() );
                }
                if ( Zprint ) {
                    debugStream << " - = Z ============" << std::endl;
                    debugStream.Log( &ZmatrixType::MatrixType::print,  Zmatrix.matrix() );
                }
                if ( Eprint ) {
                    debugStream << " - = E ============" << std::endl;
                    debugStream.Log( &EmatrixType::MatrixType::print,  Ematrix.matrix() );
                }
                if ( Rprint ) {
                    debugStream << " - = R ============" << std::endl;
                    debugStream.Log( &RmatrixType::MatrixType::print,  Rmatrix.matrix() );
                }
                if ( H1print ) {
//                    debugStream << " - = H1 ===========" << std::endl;
//                    debugStream.Log( &DiscreteSigmaFunctionType::print, H1rhs );
                    Stuff::oneLinePrint( debugStream, H1rhs );
                }
#ifdef CHEAT
                if ( H1print ) {
//                    debugStream << " - = exact H1 =====" << std::endl;
//                    debugStream.Log( &DiscreteSigmaFunctionType::print, exactH1rhs );
                    Stuff::oneLinePrint( debugStream, exactH1rhs );
                }
#endif
                if ( H2print ) {
//                    debugStream << " - = H2 ===========" << std::endl;
//                    debugStream.Log( &DiscreteVelocityFunctionType::print, H2rhs );
                    Stuff::oneLinePrint( debugStream, H2rhs );
                }
#ifdef CHEAT
                if ( H2print ) {
//                    debugStream << " - = exact H2 =====" << std::endl;
//                    debugStream.Log( &DiscreteVelocityFunctionType::print, exactH2rhs );
                    Stuff::oneLinePrint( debugStream, exactH2rhs );
                }
#endif
                if ( H3print ) {
//                    debugStream << " - = H3 ===========" << std::endl;
//                    debugStream.Log( &DiscretePressureFunctionType::print, H3rhs );
                    Stuff::oneLinePrint( debugStream, H3rhs );
                }
#ifdef CHEAT
                if ( H3print ) {
//                    debugStream << " - = exact H3 =====" << std::endl;
//                    debugStream.Log( &DiscretePressureFunctionType::print, exactH3rhs );
                    Stuff::oneLinePrint( debugStream, exactH3rhs );
                }
#endif
                debugStream << std::endl;
            }

            double H1Min = 0.0;
            double H1Max = 0.0;
            double H2Min = 0.0;
            double H2Max = 0.0;
            double H3Min = 0.0;
            double H3Max = 0.0;
            Stuff::getMinMaxOfDiscreteFunction( H1rhs, H1Min, H1Max );
            Stuff::getMinMaxOfDiscreteFunction( H2rhs, H2Min, H2Max );
            Stuff::getMinMaxOfDiscreteFunction( H3rhs, H3Min, H3Max );
            double H1avg = Stuff::getFuncAvg( H1rhs );
            double H2avg = Stuff::getFuncAvg( H2rhs );
            double H3avg = Stuff::getFuncAvg( H3rhs );

#ifdef CHEAT
            double exactH1Min = 0.0;
            double exactH1Max = 0.0;
            double exactH2Min = 0.0;
            double exactH2Max = 0.0;
            double exactH3Min = 0.0;
            double exactH3Max = 0.0;
            Stuff::getMinMaxOfDiscreteFunction( exactH1rhs, exactH1Min, exactH1Max );
            Stuff::getMinMaxOfDiscreteFunction( exactH2rhs, exactH2Min, exactH2Max );
            Stuff::getMinMaxOfDiscreteFunction( exactH3rhs, exactH3Min, exactH3Max );
            double eH1avg = Stuff::getFuncAvg( exactH1rhs );
            double eH2avg = Stuff::getFuncAvg( exactH2rhs );
            double eH3avg = Stuff::getFuncAvg( exactH3rhs );
            double tol = Parameters().getParam( "diff-tolerance", 0.01 );
            double h1_diffs = Stuff::getNumDiffDofs( H1rhs, exactH1rhs, tol );
            double h2_diffs = Stuff::getNumDiffDofs( H2rhs, exactH2rhs, tol );
            double h3_diffs = Stuff::getNumDiffDofs( H3rhs, exactH3rhs, tol );
            Dune::L2Norm< GridPartType > l2 ( exactH3rhs.space().gridPart() );
            tmpH1rhs.assign( exactH1rhs );
            tmpH1rhs -= H1rhs;
            double h1_err = l2.norm( tmpH1rhs );
            tmpH2rhs.assign( exactH2rhs );
            tmpH2rhs -= H2rhs;
            double h2_err = l2.norm( tmpH2rhs );
            tmpH3rhs.assign( exactH3rhs );
            tmpH3rhs -= H3rhs;
            double h3_err = l2.norm( tmpH3rhs );
#endif
            debugStream << "- printing infos" << std::endl
                        << "  - H1" << std::endl
                        << "    min: " << H1Min << std::endl
                        << "    max: " << H1Max << std::endl
                        << "    avg: " << H1avg << std::endl
#ifdef CHEAT
                        << "  - exact H1" << std::endl
                        << "    min: " << exactH1Min << std::endl
                        << "    max: " << exactH1Max << std::endl
                        << "    avg: " << eH1avg << std::endl
                        << "    err: " << h1_err << std::endl
                        << "   #dif: " << h1_diffs << " / " << H1rhs.size()
                                << " (" << 100 * h1_diffs / (double) H1rhs.size() << "%)" << std::endl
#endif
                        << "  - H2" << std::endl
                        << "    min: " << H2Min << std::endl
                        << "    max: " << H2Max << std::endl
                        << "    avg: " << H2avg << std::endl
#ifdef CHEAT
                        << "  - exact H2" << std::endl
                        << "    min: " << exactH2Min << std::endl
                        << "    max: " << exactH2Max << std::endl
                        << "    avg: " << eH2avg << std::endl
                        << "    err: " << h2_err << std::endl
                        << "   #dif: " << h2_diffs << " / " << H2rhs.size()
                                << " (" << 100 * h2_diffs / (double) H2rhs.size() << "%)" << std::endl
#endif
                        << "  - H3" << std::endl
                        << "    min: " << H3Min << std::endl
                        << "    max: " << H3Max << std::endl
                        << "    avg: " << H3avg << std::endl
#ifdef CHEAT
                        << "  - exact H3" << std::endl
                        << "    min: " << exactH3Min << std::endl
                        << "    max: " << exactH3Max << std::endl
                        << "    avg: " << eH3avg << std::endl
                        << "    err: " << h3_err << std::endl
                        << "   #dif: " << h3_diffs << " / " << H3rhs.size()
                                << " (" << 100 * h3_diffs / (double) H3rhs.size() << "%)" << std::endl
                        << " - " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl
#endif
                        << std::endl;
#endif

            // do the matlab logging stuff
            Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
            Stuff::printSparseRowMatrixMatlabStyle( MInversMatrix.matrix(), "M_invers", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Wmatrix.matrix(), "W", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Xmatrix.matrix(), "X", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Ymatrix.matrix(), "Y", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Zmatrix.matrix(), "Z", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Ematrix.matrix(), "E", matlabLogStream );
            Stuff::printSparseRowMatrixMatlabStyle( Rmatrix.matrix(), "R", matlabLogStream );
            Stuff::printDiscreteFunctionMatlabStyle( H1rhs, "H1", matlabLogStream );
            Stuff::printDiscreteFunctionMatlabStyle( H2rhs, "H2", matlabLogStream );
            Stuff::printDiscreteFunctionMatlabStyle( H3rhs, "H3", matlabLogStream );
//            matlabLogStream << "\nA = Y - X * M_invers * W;" << std::endl;
//            matlabLogStream << "B = Z;" << std::endl;
//            matlabLogStream << "B_T = - E;" << std::endl;
//            matlabLogStream << "C = R;" << std::endl;
//            matlabLogStream << "F = H2 - X * M_invers * H1;" << std::endl;
//            matlabLogStream << "G = - H3;" << std::endl;
//            matlabLogStream << "%A_invers = inv( A );" << std::endl;
//            matlabLogStream << "%schur_S = B_T * A_invers * B + C;" << std::endl;
//            matlabLogStream << "%schur_f = B_T * A_invers * F - G;" << std::endl;
//            matlabLogStream << "%p = schur_S \\ schur_f;" << std::endl;
//            matlabLogStream << "%u = A_invers * ( F - B * p );\n" << std::endl;
//            matlabLogStream << "%fprintf(1, 'Condition A: %d\\n', cond( A ) );\n" << std::endl;
//            matlabLogStream << "%fprintf(1, 'Condition S: %d\\n', cond( schur_S ) );\n" << std::endl;
//#endif
            profiler().StopTiming("Pass -- ASSEMBLE");
            profiler().StartTiming("Pass -- SOLVER");
            InvOpType op;
#ifndef CHEAT
            op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs );
#else
            if ( Parameters().getParam( "use-cheat", true ) ) {
                infoStream << "solving with cheated rhs" << std::endl;
                op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, exactH1rhs, exactH2rhs, exactH3rhs );
            }
            else {
                infoStream << "solving with \"normal\" rhs" << std::endl;
                op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs );
            }
#endif
            if ( Parameters().getParam( "solution-print", true ) ) {
                Stuff::oneLinePrint( infoStream, dest.discretePressure() );
                Stuff::oneLinePrint( infoStream, dest.discreteVelocity() );
            }
            profiler().StopTiming("Pass -- SOLVER");
            profiler().StopTiming("Pass");
#ifndef NLOG
//            debugStream.Resume();
//            Stuff::oneLinePrint( debugStream, dest.discretePressure() );
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
         *  \tparam C   SigmaJacobianRangeType
         **/
        template < class C >
        VelocityRangeType sigmaDivergenceOutOfGradient( const C& arg ) const
        {
            VelocityRangeType ret( 0.0 );
            const int dim = arg[0].dim();
            for ( int i = 0; i < dim; ++i ) {
                for ( int j = 0; j < dim; ++j ) {
                    const VelocityRangeType gradientRow = arg[( 3 * i ) + j];
                    ret[i] += gradientRow[ i ];
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
            return ret;
        }

};

}
#endif  // end of stokespass.hh
