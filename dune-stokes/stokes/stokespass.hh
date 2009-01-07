/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>

#include <dune/stokes/saddlepoint_inverse_operator.hh>

#ifndef BLUB
    #include "../src/stuff.hh" // should be removed in the end
    #include "../src/logging.hh" // should be removed in the end
#endif

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

        virtual void apply( const DomainType &arg, RangeType &dest) const
        {

#ifndef BLUB
            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();
#endif

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
            const double eps = 1.0e-15;

#ifndef BLUB
            bool output = false;
            int outputEntity = 0;
            infoStream << "\n== starting gridwalk" << std::endl;
#endif

            // walk the grid
            EntityIteratorType entityItEnd = velocitySpace_.end();
            for ( EntityIteratorType entityIt = velocitySpace_.begin(); entityIt != entityItEnd; ++entityIt ) {

#ifndef BLUB
                if ( outputEntity == 0 ) output = true;
#endif

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

#ifndef BLUB
                if ( output ) {
                    debugStream << "\nnumSigmaBaseFunctionsElement: " << numSigmaBaseFunctionsElement << std::endl;
                    debugStream << "\nnumVelocityBaseFunctionsElement: " << numVelocityBaseFunctionsElement << std::endl;
                    debugStream << "\nnumPressureBaseFunctionsElement: " << numPressureBaseFunctionsElement << std::endl;
                }
#endif

                // calculate volume integrals
                // (M)_{i,j} = \int_{T}\tau_{i}:\tau_{j}dx
                // we build M^-1 in fact, besauce M should be diagnal, and
                // inversion is pretty easy
#ifndef BLUB
                if ( output ) {
                    debugStream << "\ncalculating M" << std::endl;
                    debugStream << "=============" << std::endl;
                }
                bool Moutput = false;
#endif
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
#ifndef BLUB
                        if ( ( i == 0 ) && ( j == 8 ) ) Moutput = true;
#endif
                        double M_i_j = 0.0;
                        // get quadrature
                        VolumeQuadratureType volumeQuadrature( entity, ( 2 * sigmaSpaceOrder ) + 1 );
                        // sum over all quadrature points
#ifndef BLUB
                        if ( output && Moutput ) debugStream << "Basefunctions " << i << " " << j << std::endl;
                        if ( output && Moutput ) debugStream << "volumeQuadrature.nop() " << volumeQuadrature.nop() << std::endl;
#endif
                        for ( int quad = 0; quad < volumeQuadrature.nop(); ++quad ) {
#ifndef BLUB
                            if ( output && Moutput ) debugStream << "-quadPoint " << quad;
#endif
                            WorldCoordinateType x = volumeQuadrature.point( quad );
#ifndef BLUB
                            if ( output && Moutput ) Stuff::printFieldVector( x, "x", debugStream );
#endif
                            double elementVolume = geometry.integrationElement( x );
#ifndef BLUB
                            if ( output && Moutput ) debugStream << "\nelementVolume: " << elementVolume << std::endl;
#endif
                            double integrationWeight = volumeQuadrature.weight( quad );
#ifndef BLUB
                            if ( output && Moutput ) debugStream << "integrationWeight: " << integrationWeight;
#endif
                            // calculate \tau_{i}:\tau_{j}
                            SigmaRangeType tau_i( 0.0 );
                            SigmaRangeType tau_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
#ifndef BLUB
                            if ( output && Moutput ) Stuff::printFieldMatrix( tau_i, "tau_i", debugStream );
#endif
                            sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                            double tau_i_times_tau_j = colonProduct( tau_i, tau_j );

#ifndef BLUB
                            if ( output && Moutput ) Stuff::printFieldMatrix( tau_j, "tau_j", debugStream );
                            if ( output && Moutput ) debugStream << "\ncolonProduct( tau_i, tau_j ): " << tau_i_times_tau_j << std::endl;
#endif
                            M_i_j += elementVolume *
                                        integrationWeight *
                                        tau_i_times_tau_j;
#ifndef BLUB
                            if ( output && Moutput ) debugStream << "M_i_j: " << M_i_j << std::endl;
#endif
                        }
#ifndef BLUB
                        Moutput = false;
#endif
                        // if small, should be zero
                        if ( M_i_j <= eps ) {
                            M_i_j = 0.0;
                        }
                        // else invert
                        else {
                            M_i_j = 1.0 / M_i_j;
                        }

                        // add to matrix
                        localMmatrixElement.add( i, j, M_i_j );
                    }
                } // done calculating M
#ifndef BLUB
                if ( output ) {
                    debugStream << "\ndone calculating M" << std::endl;
                    debugStream << "==================" << std::endl;
                }
                output = false;
                ++outputEntity;
            }
#endif
#ifndef BLUB
            infoStream << "\n== gridwalk done." << std::endl;
#endif
#ifndef BLUB
            debugStream << "\nMmatrix" << std::endl;
            Mmatrix.matrix().print( std::cout );
#endif
            // build global matrices
            typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
                AmatrixType;
            AmatrixType Amatrix( velocitySpace_, velocitySpace_ );
            Amatrix.reserve();
            typedef SparseRowMatrixObject< DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
                Tmp_matrixType;
            Tmp_matrixType tmp( sigmaSpace_, velocitySpace_ );
            tmp.reserve();



        //    Mmatrix.matrix().multiply( Wmatrix.matrix(), Amatrix.matrix() );

            InvOpType op( *this, 1.0,1.0,1,1 );
            op.solve( arg, dest, Amatrix, Amatrix, tmp );



        }

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

        template < class FieldMatrixType >
        double colonProduct(    const FieldMatrixType& arg1,
                                const FieldMatrixType& arg2 ) const
        {
            assert( arg1.rowdim() == arg2.coldim() );
            double ret = 0.0;
            // iterators
            typedef typename FieldMatrixType::ConstRowIterator
                ConstRowIterator;
            typedef typename FieldMatrixType::row_type::ConstIterator
                ConstIterator;
            ConstRowIterator arg1RowItEnd = arg1.end();
            ConstRowIterator arg2RowItEnd = arg2.end();
            ConstRowIterator arg2RowIt = arg2.begin();
            for (   ConstRowIterator arg1RowIt = arg1.begin();
                    arg1RowIt != arg1RowItEnd, arg2RowIt != arg2RowItEnd;
                    ++arg1RowIt, ++arg2RowIt ) {
                ConstIterator row1ItEnd = arg1RowIt->end();
                ConstIterator row2ItEnd = arg2RowIt->end();
                ConstIterator row2It = arg2RowIt->begin();
                for (   ConstIterator row1It = arg1RowIt->begin();
                        row1It != row1ItEnd, row2It != row2ItEnd;
                        ++row1It, ++row2It ) {
                    ret += *row1It * *row2It;
                }
            }
        }


};

}
#endif  // end of stokespass.hh
