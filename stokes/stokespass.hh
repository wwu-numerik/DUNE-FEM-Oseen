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

#include <dune/stokes/solver/solvercaller.hh>
#include <dune/stokes/integrators.hh>

#include <dune/common/stdstreams.hh>
#include <dune/stuff/matrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stuff/progressbar.hh>

#include <omp.h>

template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
	struct StencilType {
		template <class T>
		static int nonZerosEstimate( T& rangeSpace_ ) {
			return 100;
		}
	};
};


#ifndef NLOG // if we want logging, should be removed in the end
    #include <dune/stuff/printing.hh>
    #include <dune/stuff/misc.hh>
    #include <dune/stuff/logging.hh>
#endif

#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

#include <dune/stuff/profiler.hh>

namespace Dune
{

template < class DiscreteModelImp >
struct StokesTraits
{
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

	typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
		VelocityJacobianRangeType;

	//! vector type of sigmas' discrete functions space's range
	typedef typename DiscreteSigmaFunctionSpaceType::RangeType
		SigmaRangeType;

	typedef typename DiscreteSigmaFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
		SigmaJacobianRangeType;

	//! Vector type of the pressure's discrete function space's range
	typedef typename DiscretePressureFunctionSpaceType::RangeType
		PressureRangeType;

	typedef typename DiscretePressureFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
		PressureJacobianRangeType;

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

	//! polynomial order for the discrete sigma function space
	static const int sigmaSpaceOrder
		= DiscreteModelType::sigmaSpaceOrder;
	//! polynomial order for the discrete velocity function space
	static const int velocitySpaceOrder
		= DiscreteModelType::velocitySpaceOrder;
	//! polynomial order for the discrete pressure function space
	static const int pressureSpaceOrder
		= DiscreteModelType::pressureSpaceOrder;

	//! the stab coeff. for sigma is a vector field, paramterized by the element's normal
	typedef StabilizationCoefficients::C12< VelocityRangeType >
		C12;
};

/**
 *  \brief  StokesPass
 *
 *  \todo   doc
 **/
template <  class DiscreteModelImp,
            class PreviousPassImp,
			int PassID = 0,
		  template <class> class TraitsImp = StokesTraits >
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

		typedef TraitsImp< DiscreteModelType >
			Traits;

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

		//! when requested we store \f$ \vardelta u, \nabla p (u \cdot \nabla ) u\f$ in this struct after the solver
		struct RhsDatacontainer {
			typename Traits::DiscreteVelocityFunctionType velocity_laplace;
			typename Traits::DiscreteVelocityFunctionType pressure_gradient;
			typename Traits::DiscreteSigmaFunctionType velocity_gradient;
			typename Traits::DiscreteVelocityFunctionType convection;

			RhsDatacontainer( const typename Traits::DiscreteVelocityFunctionSpaceType& space,
							  const typename Traits::DiscreteSigmaFunctionSpaceType& sigma_space)
				: velocity_laplace( "velocity_laplace", space ),
				pressure_gradient( "pressure_gradient", space ),
				velocity_gradient( "velocity_gradient", sigma_space ),
				convection( "convection", space )
			{}
			void scale( double factor ) {
				velocity_laplace	*= factor;
				pressure_gradient	*= factor;
				velocity_gradient	*= factor;
				convection			*= factor;
			}

		};

        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        StokesPass( PreviousPassType& prevPass,
                    DiscreteModelType& discreteModel,
					typename Traits::GridPartType& gridPart,
					const typename Traits::DiscreteStokesFunctionSpaceWrapperType& spaceWrapper,
					const typename Traits::DiscreteVelocityFunctionType& beta,
					const bool do_oseen_discretization )//! \todo move to model
            : BaseType( prevPass ),
            discreteModel_( discreteModel ),
            gridPart_( gridPart ),
            spaceWrapper_( spaceWrapper ),
            velocitySpace_( spaceWrapper.discreteVelocitySpace() ),
            pressureSpace_( spaceWrapper.discretePressureSpace() ),
			sigmaSpace_( gridPart ),
			beta_( beta ),
			do_oseen_discretization_( do_oseen_discretization )
        {}

        /**
         *  \brief  empty constructor
         **/
        StokesPass()
        {}

        //! used in Postprocessing to get refs to gridparts, spaces
		const typename Traits::DiscreteStokesFunctionSpaceWrapperType& GetFunctionSpaceWrapper() const
        {
            return spaceWrapper_;
        }


		void apply( const DomainType &arg, RangeType &dest ) const
		{
			apply<RhsDatacontainer,typename Traits::DiscreteSigmaFunctionType>( arg, dest, 0,0 );
		}

		template < class RhsDatacontainerType >
		void apply( const DomainType &arg, RangeType &dest,RhsDatacontainerType* rhs_datacontainer ) const
		{
			apply<RhsDatacontainerType,typename Traits::DiscreteSigmaFunctionType>( arg, dest, rhs_datacontainer,0 );
		}
        /**
         *  \todo doc
         *  \attention  think about quadrature orders
         **/
		template < class RhsDatacontainerType, class ExactSigmaType >
		void apply( const DomainType &arg, RangeType &dest, RhsDatacontainerType* rhs_datacontainer, const ExactSigmaType* sigma_exact ) const
        {

            // profiler information
            profiler().StartTiming("Pass -- ASSEMBLE");

            // entity and geometry types
			typedef typename Traits::EntityType::Geometry
                EntityGeometryType;
            typedef typename Dune::FieldMatrix< typename EntityGeometryType::ctype,
                                                EntityGeometryType::coorddimension,
                                                EntityGeometryType::mydimension >
                JacobianInverseTransposedType;

            // viscosity
			const double viscosity = discreteModel_.viscosity();

            // generalized stokes alpha
            const double alpha = discreteModel_.alpha();

            // matrices
            // M\in R^{M\times M}
			typedef typename Traits::DiscreteSigmaFunctionSpaceType
				DiscreteSigmaFunctionSpaceType;
			typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
											DiscreteSigmaFunctionSpaceType,
											MatrixTraits<DiscreteSigmaFunctionSpaceType,DiscreteSigmaFunctionSpaceType> >
                MInversMatrixType;
            MInversMatrixType MInversMatrix( sigmaSpace_, sigmaSpace_ );
            MInversMatrix.reserve();
            assert( MInversMatrix.matrix().rows() == MInversMatrix.matrix().cols() );
            // W\in R^{M\times L}
			typedef typename Traits::DiscreteVelocityFunctionSpaceType
				DiscreteVelocityFunctionSpaceType;
            typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType> >
                WmatrixType;
            WmatrixType Wmatrix( sigmaSpace_, velocitySpace_ );
            Wmatrix.reserve();
            // X\in R^{L\times M}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteSigmaFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType> >
                XmatrixType;
            XmatrixType Xmatrix( velocitySpace_, sigmaSpace_ );
            Xmatrix.reserve();
            // Y\in R^{L\times L}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                YmatrixType;
            YmatrixType Ymatrix( velocitySpace_, velocitySpace_ );
            Ymatrix.reserve();
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                OmatrixType;
			OmatrixType Omatrix( velocitySpace_, velocitySpace_ );
			Omatrix.reserve();
            // Z\in R^{L\times K}
			typedef typename Traits::DiscretePressureFunctionSpaceType
				DiscretePressureFunctionSpaceType;
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscretePressureFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscretePressureFunctionSpaceType> >
                ZmatrixType;
            ZmatrixType Zmatrix( velocitySpace_, pressureSpace_ );
            Zmatrix.reserve();
            // E\in R^{K\times L}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscretePressureFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                EmatrixType;
            EmatrixType Ematrix( pressureSpace_, velocitySpace_ );
            Ematrix.reserve();
            // R\in R^{K\times K}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
											DiscretePressureFunctionSpaceType,
											MatrixTraits<DiscretePressureFunctionSpaceType,DiscretePressureFunctionSpaceType> >
                RmatrixType;
            RmatrixType Rmatrix( pressureSpace_, pressureSpace_ );
            Rmatrix.reserve();

            // local matrices
            // M\in R^{M\times M}
            typedef typename MInversMatrixType::LocalMatrixType
                LocalMInversMatrixType;
            // X\in R^{L\times M}
            typedef typename XmatrixType::LocalMatrixType
                LocalXmatrixType;
            // Y\in R^{L\times L}
            typedef typename YmatrixType::LocalMatrixType
                LocalYmatrixType;
			typedef typename OmatrixType::LocalMatrixType
				LocalOmatrixType;
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
			typename Traits::DiscreteSigmaFunctionType H1rhs( "H1", sigmaSpace_ );
            H1rhs.clear();
            // H_{2}\in R^{L}
			typename Traits::DiscreteVelocityFunctionType H2rhs( "H2", velocitySpace_ );
            H2rhs.clear();
			typename Traits::DiscreteVelocityFunctionType H2_O_rhs( "H2_O", velocitySpace_ );
			H2_O_rhs.clear();
            // H_{3}\in R^{K}
			typename Traits::DiscretePressureFunctionType H3rhs( "H3", pressureSpace_ );
            H3rhs.clear();

            // local right hand sides
            // H_{1}\in R^{M}
			typedef typename Traits::DiscreteSigmaFunctionType::LocalFunctionType
                LocalH1rhsType;
            // H_{2}\in R^{L}
			typedef typename Traits::DiscreteVelocityFunctionType::LocalFunctionType
                LocalH2rhsType;
            // H_{3}\in R^{K}
			typedef typename Traits::DiscretePressureFunctionType::LocalFunctionType
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

            // eps
            const double eps = Parameters().getParam( "eps", 1.0e-14 );
			const double convection_scaling = discreteModel_.convection_scaling();
			const double pressure_gradient_scaling = discreteModel_.pressure_gradient_scaling();

			Logger().Info() << boost::format( "pressure_gradient/convection scaling: %e | %e\npass viscosity: %e\n")
								% convection_scaling
								% pressure_gradient_scaling
								% viscosity;

#ifndef NLOG
            // logging stuff
            Logging::LogStream& infoStream = Logger().Info();
            int entityNR = 0;
            int numberOfEntities = 0;
            int numberOfIntersections = 0;
            int numberOfBoundaryIntersections = 0;
            int numberOfInnerIntersections = 0;
            infoStream << "this is StokesPass::apply()" << std::endl;

            // do an empty grid walk to get informations
            double maxGridWidth( 0.0 );
			typename Traits::EntityIteratorType entityItEndLog = velocitySpace_.end();
			for (   typename Traits::EntityIteratorType entityItLog = velocitySpace_.begin();
                    entityItLog != entityItEndLog;
                    ++entityItLog ) {
				const typename Traits::EntityType& entity = *entityItLog;
                // count entities
                ++numberOfEntities;
                // walk the intersections
				typename Traits::IntersectionIteratorType intItEnd = gridPart_.iend( entity );
				for (   typename Traits::IntersectionIteratorType intIt = gridPart_.ibegin( entity );
                        intIt != intItEnd;
                        ++intIt ) {
                    // count intersections
                    ++numberOfIntersections;
					maxGridWidth = std::max( Stuff::getLenghtOfIntersection( *intIt ), maxGridWidth );
                    // if we are inside the grid
					if ( intIt->neighbor() && !intIt->boundary() ) {
                        // count inner intersections
                        ++numberOfInnerIntersections;
                    }
                    // if we are on the boundary of the grid
					if ( !intIt->neighbor() && intIt->boundary() ) {
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
                infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
                infoStream << "- starting gridwalk" << std::endl;
            } else {
                infoStream << "found " << numberOfEntities << " entities," << std::endl;
                infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
                infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
                infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
                infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
                infoStream << "- starting gridwalk" << std::endl;
            }
			infoStream.Suspend();
#endif

            // walk the grid

			Stokes::Integrators::W< Traits >
					w_integrator( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_ );
			w_integrator.apply( Wmatrix );

			typename Traits::EntityIteratorType entityItEnd = velocitySpace_.end();
			typename Traits::EntityIteratorType entityIt = velocitySpace_.begin();
			infoStream.Resume();
            for (   Stuff::SimpleProgressBar<Logging::LogStream> progress( (numberOfEntities-1), infoStream, 40 );
                    entityIt != entityItEnd;
					++entityIt,++entityNR,++progress ) {

                // get entity and geometry
				const typename Traits::EntityType& entity = *entityIt;
				const EntityGeometryType& geometry = entity.geometry();

                // get local matrices for the volume integral
                LocalMInversMatrixType localMInversMatrixElement = MInversMatrix.localMatrix( entity, entity );
                LocalXmatrixType localXmatrixElement = Xmatrix.localMatrix( entity, entity );
                LocalYmatrixType localYmatrixElement = Ymatrix.localMatrix( entity, entity );
				LocalOmatrixType localOmatrixElement = Omatrix.localMatrix( entity, entity );
                LocalZmatrixType localZmatrixElement = Zmatrix.localMatrix( entity, entity );
                LocalEmatrixType localEmatrixElement = Ematrix.localMatrix( entity, entity );
                LocalRmatrixType localRmatrixElement = Rmatrix.localMatrix( entity, entity );

                // get local right hand sides
                LocalH1rhsType localH1rhs = H1rhs.localFunction( entity );
                LocalH2rhsType localH2rhs = H2rhs.localFunction( entity );
				LocalH2rhsType localH2_O_rhs = H2_O_rhs.localFunction( entity );
                LocalH3rhsType localH3rhs = H3rhs.localFunction( entity );

                // get basefunctionsets
				const SigmaBaseFunctionSetType sigmaBaseFunctionSetElement = sigmaSpace_.baseFunctionSet( entity );
				const VelocityBaseFunctionSetType velocityBaseFunctionSetElement = velocitySpace_.baseFunctionSet( entity );
				const PressureBaseFunctionSetType pressureBaseFunctionSetElement = pressureSpace_.baseFunctionSet( entity );
                const int numSigmaBaseFunctionsElement = sigmaBaseFunctionSetElement.numBaseFunctions();
                const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.numBaseFunctions();
                const int numPressureBaseFunctionsElement = pressureBaseFunctionSetElement.numBaseFunctions();

                // get quadrature
				const typename Traits::VolumeQuadratureType volumeQuadratureElement( entity,
																	( 4 * Traits::pressureSpaceOrder ) + 1 );

				typedef typename Traits::ElementCoordinateType
					ElementCoordinateType;
				typedef typename Traits::SigmaRangeType
					SigmaRangeType;
				typedef typename Traits::VelocityRangeType
					VelocityRangeType;
				typedef typename Traits::PressureRangeType
					PressureRangeType;
				typedef typename Traits::VelocityJacobianRangeType
					VelocityJacobianRangeType;
				typedef typename Traits::PressureJacobianRangeType
					PressureJacobianRangeType;

				// compute volume integrals
                //                                                     // we will call this one
                // (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
                for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double M_i_j = 0.0;
                        // sum over all quadrature points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
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
							const double tau_i_times_tau_j = Stuff::colonProduct( tau_i, tau_j );
                            // compute M_i_j
                            M_i_j += elementVolume
                                * integrationWeight
                                * tau_i_times_tau_j;
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( M_i_j ) < eps ) {
                            M_i_j = 0.0;
                        } // else invert
                        else {
                            M_i_j = 1.0 / M_i_j;
                            // add to matrix
                            localMInversMatrixElement.add( i, j, M_i_j );
                        }
                    }
                } // done computing Minvs' volume integral

                //                                                  // we will call this one
                // (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx // X's volume integral
                //                                                  // see also "X's entitity surface integral", "X's neighbour surface integral" and "X's boundary integral" below
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                        double X_i_j = 0.0;
                        // sum over all quadrature points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute \tau_{j}:\nabla v_{i}
                            SigmaRangeType tau_j( 0.0 );
                            sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                            const double gradient_of_v_i_times_tau_j = velocityBaseFunctionSetElement.evaluateGradientSingle( i, entity, x, tau_j );
                            X_i_j += elementVolume
                                * integrationWeight
                                * gradient_of_v_i_times_tau_j;
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( X_i_j ) < eps ) {
                            X_i_j = 0.0;
                        }
                        else
                            // add to matrix
                            localXmatrixElement.add( i, j, X_i_j );
                    }
                } // done computing X's volume integral

                //
                // (Y)
                //
//                if ( discreteModel_.isGeneralized() ) 
                {
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double Y_i_j = 0.0;
                        // sum over all quadrature points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute \tau_{j}:\nabla v_{i}
                            VelocityRangeType v_i( 0.0 );
                            velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                            VelocityRangeType v_j( 0.0 );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            const double v_i_times_v_j = v_i * v_j;
                            Y_i_j += elementVolume
                                * integrationWeight
                                * alpha
                                * v_i_times_v_j;
                        } // done sum over quadrature points
                        // if small, should be zero
                        if ( fabs( Y_i_j ) < eps ) {
                            Y_i_j = 0.0;
                        }
                        else
                            // add to matrix
                            localYmatrixElement.add( i, j, Y_i_j );
                    }
                } // done computing Y's volume integral
                }

				for ( int i = 0; (i < numVelocityBaseFunctionsElement ) && do_oseen_discretization_; ++i ) {
					for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
						double O_i_j = 0.0;
						double O_i_j_d = 0.0;
						// sum over all quadrature points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
							// get x
							const ElementCoordinateType x = volumeQuadratureElement.point( quad );
							// get the integration factor
							const double elementVolume = geometry.integrationElement( x );
							// get the quadrature weight
							const double integrationWeight = volumeQuadratureElement.weight( quad );
							//calc u_h * \nabla * (v \tensor \beta )
							VelocityRangeType v_i( 0.0 );
							velocityBaseFunctionSetElement.evaluate( i, x, v_i );
							VelocityRangeType v_j( 0.0 );
							velocityBaseFunctionSetElement.evaluate( j, x, v_j );
							VelocityRangeType beta_eval;
							beta_.localFunction( entity ).evaluate( x, beta_eval );

							VelocityJacobianRangeType v_i_jacobian;
							velocityBaseFunctionSetElement.jacobian( i, x, v_i_jacobian );
							VelocityJacobianRangeType v_j_jacobian;
							velocityBaseFunctionSetElement.jacobian( j, x, v_j_jacobian );
							VelocityJacobianRangeType beta_jacobian;
							const typename Traits::DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
									beta_.localFunction( entity );
							beta_lf.jacobian( x, beta_jacobian );


							VelocityRangeType divergence_of_v_i_tensor_beta;
							for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
								double row_result = 0;
								for ( size_t m = 0; m < beta_eval.dim(); ++m ) {
									row_result += beta_jacobian[l][m] * v_i[l] + v_i_jacobian[l][m] * beta_eval[l];
								}
								divergence_of_v_i_tensor_beta[l] = row_result;
							}
							for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
								assert( !isnan(divergence_of_v_i_tensor_beta[l]) );
							}

						   divergence_of_v_i_tensor_beta[0] = beta_eval[0] * v_i_jacobian[0][0]
								   + v_i[0] * beta_jacobian[0][0]
								   + beta_eval[0] * v_i_jacobian[1][1]
								   + v_i[1] * beta_jacobian[0][1];
						   divergence_of_v_i_tensor_beta[1] = beta_eval[1] * v_i_jacobian[0][0]
								   + v_i[0] * beta_jacobian[1][0]
								   + beta_eval[1] * v_i_jacobian[1][1]
								   + v_i[1] * beta_jacobian[1][1];


							const double u_h_times_divergence_of_beta_v_j_tensor_beta =
									v_j * divergence_of_v_i_tensor_beta;
							VelocityJacobianRangeType v_i_tensor_beta = dyadicProduct( v_i, beta_eval );
							const double ret = Stuff::colonProduct( v_i_tensor_beta, v_j_jacobian );

							O_i_j -= elementVolume
								* integrationWeight
								* convection_scaling
//								* u_h_times_divergence_of_beta_v_j_tensor_beta;
									* ret;
							O_i_j_d-= elementVolume
							        * integrationWeight
									* convection_scaling
									* u_h_times_divergence_of_beta_v_j_tensor_beta;
//										* ret;

						}
						if ( fabs( O_i_j ) < eps ) {
							O_i_j = 0.0;
						}
						else {
							// add to matrix
							localOmatrixElement.add( i, j, O_i_j );
//							double diff = O_i_j- O_i_j_d;
//							std::cout << boost::format( "DIFF %e\n") % diff;
						}
					}
				}

                //                                                  // we will call this one
                // (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx // Z's volume integral
                //                                                  // see also "Z's entitity surface integral", "Z's neighbour surface integral" and "Z's boundary integral" below
                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                        double Z_i_j = 0.0;
                        // sum over all quadratur points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute q_{j}\cdot(\nabla\cdot v_i)
                            PressureRangeType q_j( 0.0 );
                            pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                            const double divergence_of_v_i_times_q_j = velocityBaseFunctionSetElement.evaluateGradientSingle( i, entity, x, preparePressureRangeTypeForVelocityDivergence( q_j ) );
                            Z_i_j += -1.0
                                * elementVolume
                                * integrationWeight
                                * pressure_gradient_scaling
                                * divergence_of_v_i_times_q_j;
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( Z_i_j ) < eps ) {
                            Z_i_j = 0.0;
                        }
                        else
                            // add to matrix
                            localZmatrixElement.add( i, j, Z_i_j );
                    }
                } // done computing Z's volume integral

                //                                    // we will call this one
                // (H2)_{j} += \int_{T}f\cdot v_{j}dx // H2's volume integral
                //                                    // see also "H2's boundary integral" further down
                if ( discreteModel_.hasForce() ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double H2_j = 0.0;
                        // sum over all quadratur points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            const ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            const VelocityRangeType xWorld = geometry.global( x );
                            // get the integration factor
                            const double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            const double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute f\cdot v_j
                            VelocityRangeType v_j( 0.0 );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            VelocityRangeType f( 0.0 );
#if MODEL_PROVIDES_LOCALFUNCTION
										discreteModel_.forceF().localFunction(entity).evaluate( x, f );
#else
										discreteModel_.force( 0.0, xWorld, f );
#endif
                            const double f_times_v_j = f * v_j;
                            H2_j += elementVolume
                                * integrationWeight
                                * f_times_v_j;
                        } // done sum over all quadrature points
                        // if small, should be zero
                        if ( fabs( H2_j ) < eps ) {
                            H2_j = 0.0;
                        }
						else
							// add to rhs
							localH2rhs[ j ] += H2_j;
                    } // done computing H2's volume integral
                }

                //                                                // we will call this one
                // (E)_{i,j} += -\int_{T}v_{j}\cdot\nabla q_{i}dx // E's volume integral
                //                                                // see also "E's entitity surface integral", "E's neighbour surface integral" and "E's boundary integral" below
                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                    for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                        double E_i_j = 0.0;
                        // sum over all quadratur points
						for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
                            // get x
                            ElementCoordinateType x = volumeQuadratureElement.point( quad );
                            // get the integration factor
                            double elementVolume = geometry.integrationElement( x );
                            // get the quadrature weight
                            double integrationWeight = volumeQuadratureElement.weight( quad );
                            // compute v_{j}\cdot(\nabla q_i)
                            VelocityRangeType v_j( 0.0 );
                            velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                            PressureJacobianRangeType jacobian_of_q_i( 0.0 );
                            pressureBaseFunctionSetElement.jacobian( i, x, jacobian_of_q_i );
                            const VelocityRangeType gradient_of_q_i_untransposed( jacobian_of_q_i[0] );
                            const JacobianInverseTransposedType jacobianInverseTransposed = geometry.jacobianInverseTransposed( x );
                            VelocityRangeType gradient_of_q_i( 0.0 );
                            jacobianInverseTransposed.mv( gradient_of_q_i_untransposed, gradient_of_q_i );
                            const double gradient_of_q_i_times_v_j = gradient_of_q_i * v_j;
                            E_i_j += -1.0
                                * elementVolume
                                * integrationWeight
                                * gradient_of_q_i_times_v_j;
                        } // done sum over all quadrature points
                        // if small, should be zero
						if ( fabs( E_i_j ) < eps ) {
                            E_i_j = 0.0;
                        }
                        else
                            // add to matrix
                            localEmatrixElement.add( i, j, E_i_j );
                    }
                } // done computing E's volume integral

                // walk the intersections
				typename Traits::IntersectionIteratorType intItEnd = gridPart_.iend( entity );
				for (   typename Traits::IntersectionIteratorType intIt = gridPart_.ibegin( entity );
						intIt != intItEnd;
						++intIt ) {
					const typename Traits::IntersectionIteratorType::Intersection& intersection = *intIt;

                    // get intersection geometry
					typedef typename Traits::IntersectionIteratorType::Geometry
                        IntersectionGeometryType;
					const IntersectionGeometryType& intersectionGeometry = intersection.intersectionGlobal();
                    // get intersection quadrature, seen from inside
					const typename Traits::FaceQuadratureType faceQuadratureElement( gridPart_,
																	intersection,
																	( 4 * Traits::pressureSpaceOrder ) + 1,
																	Traits::FaceQuadratureType::INSIDE );

                    // get flux coefficients
					const double lengthOfIntersection = Stuff::getLenghtOfIntersection( intersection );
					const StabilizationCoefficients& stabil_coeff ( discreteModel_.getStabilizationCoefficients() );
                    const double C_11 = stabil_coeff.Factor("C11") * std::pow( lengthOfIntersection, stabil_coeff.Power("C11") );
                    const double D_11 = stabil_coeff.Factor("D11") * std::pow( lengthOfIntersection, stabil_coeff.Power("D11") );
					//we'll leave this on 0 for the time being so it does not generate any additional  penalty terms
//					const VelocityRangeType D_12(stabil_coeff.Factor("D12") );//TODO FIXME
					VelocityRangeType D_12( 1 );//TODO FIXME
					D_12 /= D_12.two_norm();
					D_12 *= stabil_coeff.Factor("D12");
//					VelocityRangeType E_11(stabil_coeff.Factor("E12"));

                    // if we are inside the grid
					if ( intersection.neighbor() && !intersection.boundary() ) {
                        // get neighbour
						//! DO NOT TRY TO DEREF outside() DIRECTLY
						const typename Traits::IntersectionIteratorType::EntityPointer neighbourPtr = intersection.outside();
						//there's some (copy ctor? / implicit type conversion) at play here which will make shit fail if you do
						const typename Traits::EntityType& neighbour = *neighbourPtr;
						// esp. these two line ARE NOT equiv. to const EntityType& neighbour = * (static_cast<const typename IntersectionIteratorType::EntityPointer>(intersection.outside()) );

                        // get local matrices for the surface integrals
                        LocalMInversMatrixType localMInversMatrixNeighbour = MInversMatrix.localMatrix( entity, neighbour );

                        LocalXmatrixType localXmatrixNeighbour = Xmatrix.localMatrix( entity, neighbour );
                        LocalYmatrixType localYmatrixNeighbour = Ymatrix.localMatrix( neighbour, entity );
						LocalOmatrixType localOmatrixNeighbour = Omatrix.localMatrix( neighbour, entity );
                        LocalZmatrixType localZmatrixNeighbour = Zmatrix.localMatrix( entity, neighbour );
                        LocalEmatrixType localEmatrixNeighbour = Ematrix.localMatrix( neighbour, entity );
                        LocalRmatrixType localRmatrixNeighbour = Rmatrix.localMatrix( entity, neighbour );

                        // get basefunctionsets
                        const SigmaBaseFunctionSetType sigmaBaseFunctionSetNeighbour = sigmaSpace_.baseFunctionSet( neighbour );
                        const VelocityBaseFunctionSetType velocityBaseFunctionSetNeighbour = velocitySpace_.baseFunctionSet( neighbour );
                        const PressureBaseFunctionSetType pressureBaseFunctionSetNeighbour = pressureSpace_.baseFunctionSet( neighbour );
                        const int numSigmaBaseFunctionsNeighbour = sigmaBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numVelocityBaseFunctionsNeighbour = velocityBaseFunctionSetNeighbour.numBaseFunctions();
                        const int numPressureBaseFunctionsNeighbour = pressureBaseFunctionSetNeighbour.numBaseFunctions();

                        // get intersection quadrature, seen from outside
						const typename Traits::FaceQuadratureType faceQuadratureNeighbour(   gridPart_,
																			intersection,
																			( 4 * Traits::pressureSpaceOrder ) + 1,
																			Traits::FaceQuadratureType::OUTSIDE );
						typedef typename Traits::LocalIntersectionCoordinateType
							LocalIntersectionCoordinateType;

                        // compute surface integrals



                        //                                                                                                                   // we will call this one
                        // (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's element sourface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}ds // X's neighbour sourface integral
                        //                                                                                                                   // see also "X's boundary integral" below
                        //                                                                                                                   // and "X's volume integral" above
//                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                // compute X's element sourface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                    double X_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x in codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );

										VelocityRangeType tau_j_times_normal( 0.0 );
										tau_j.mv( outerNormal, tau_j_times_normal );
										typename Traits::C12 c_12( outerNormal, discreteModel_.getStabilizationCoefficients() );
										SigmaRangeType tau_j_times_normal_dyadic_C12
												= dyadicProduct( tau_j_times_normal, c_12 );

										SigmaRangeType flux_value = tau_j;
										flux_value *= 0.5;
										flux_value -= tau_j_times_normal_dyadic_C12;

										VelocityRangeType flux_times_normal( 0.0 );
										flux_value.mv( outerNormal, flux_times_normal );

										const double v_i_times_flux_times_normal
												= velocityBaseFunctionSetElement.evaluateSingle( i, x, flux_times_normal );
										X_i_j -= elementVolume
											* integrationWeight
											* v_i_times_flux_times_normal;

                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( X_i_j ) < eps ) {
                                        X_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localXmatrixElement.add( i, j, X_i_j );
                                } // done computing X's element sourface integral
                                // compute X's neighbour sourface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
                                    double X_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
                                        const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetNeighbour.evaluate( j, xOutside, tau_j );

										VelocityRangeType tau_j_times_normal( 0.0 );
										tau_j.mv( outerNormal, tau_j_times_normal );
										typename Traits::C12 c_12( outerNormal, discreteModel_.getStabilizationCoefficients() );
										SigmaRangeType tau_j_times_normal_dyadic_C12
												= dyadicProduct( tau_j_times_normal, c_12 );

										SigmaRangeType flux_value = tau_j;
										flux_value *= 0.5;
										flux_value += tau_j_times_normal_dyadic_C12;

										VelocityRangeType flux_times_normal( 0.0 );
										flux_value.mv( outerNormal, flux_times_normal );

										const double v_i_times_flux_times_normal
												= velocityBaseFunctionSetNeighbour.evaluateSingle( i, xInside, flux_times_normal );
										X_i_j -= elementVolume
                                            * integrationWeight
											* v_i_times_flux_times_normal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( X_i_j ) < eps ) {
                                        X_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localXmatrixNeighbour.add( i, j, X_i_j );
                                } // done computing X's neighbour sourface integral
                            } // done computing X's sourface integrals
//                        }

                        //                                                                                                         // we call this one
                        // (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}ds // Y's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}ds // Y's neighbour surface integral
                        //                                                                                                         // see also "Y's boundary integral" below
//                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                // compute Y's element surface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                    double Y_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_v_j = v_i * v_j;
                                        Y_i_j += C_11
                                            * elementVolume
                                            * integrationWeight
                                            * v_i_times_v_j;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Y_i_j ) < eps ) {
                                        Y_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localYmatrixElement.add( i, j, Y_i_j );
                                } // done computing Y's element surface integral
                                // compute Y's neighbour surface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
                                    double Y_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
                                        const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetNeighbour.evaluate( i, xOutside, v_i );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
                                        const double v_i_times_v_j = v_i * v_j;
                                        Y_i_j += -1.0
                                            * C_11
                                            * elementVolume
                                            * integrationWeight
                                            * v_i_times_v_j;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Y_i_j ) < eps ) {
                                        Y_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localYmatrixNeighbour.add( i, j, Y_i_j );
                                } // done computing Y's neighbour surface integral
                            } // done computing Y's surface integrals
//                        }
							//                                                                                                         // we call this one
							// (O)_{i,j} += \int_{ // O's element surface integral
							//           += \int_{ // O's neighbour surface integral
							//                                                                                                         // see also "O's boundary integral" below

							for ( int j = 0; (j < numVelocityBaseFunctionsElement ) && do_oseen_discretization_; ++j ) {
								// compute O's element surface integral
								for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
									double O_i_j = 0.0;
									// sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
										// get x codim<0> and codim<1> coordinates
										const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
										const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
										const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
										const VelocityRangeType xWorld = geometry.global( xInside );
										const VelocityRangeType xWorld_Outside = geometry.global( xOutside );
										// get the integration factor
										const double elementVolume = intersectionGeometry.integrationElement( xLocal );
										// get the quadrature weight
										const double integrationWeight = faceQuadratureElement.weight( quad );
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
										VelocityRangeType v_i( 0.0 );
										velocityBaseFunctionSetElement.evaluate( i, xInside, v_i );
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );

										VelocityRangeType beta_eval;
										beta_.localFunction(entity).evaluate( xInside, beta_eval );
										const double beta_times_normal = beta_eval * outerNormal;
										//calc u^c_h \tensor beta * v \tensor n (self part), the flux value
										double c_s = (beta_times_normal) * 0.5;
										VelocityRangeType u_h = v_i;
										VelocityJacobianRangeType mean_value = dyadicProduct( u_h, beta_eval );
										mean_value *= 0.5;
										VelocityJacobianRangeType u_jump = dyadicProduct( v_i, outerNormal );
										u_jump *= c_s;
										VelocityJacobianRangeType flux_value = mean_value;
										flux_value += u_jump;

										// \int_{dK} flux_value : ( v_j \ctimes n ) ds
										VelocityJacobianRangeType v_i_tensor_n = dyadicProduct( v_j, outerNormal );
										double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );

										O_i_j += elementVolume
												* integrationWeight
												* convection_scaling
												* ret;
									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( O_i_j ) < eps ) {
										O_i_j = 0.0;
									}
									else {
										// add to matrix
//										std::cerr<< boost::format( "O face value (el) on entity %d: %e\n") % entityNR % O_i_j;
										localOmatrixElement.add( i, j, O_i_j );
									}
								} // done computing Y's element surface integral
								// compute O's neighbour surface integral
								for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
									double O_i_j = 0.0;
									// sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
										// get x codim<0> and codim<1> coordinates
										const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
										const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
										const VelocityRangeType xWorld = geometry.global( xInside );
										const VelocityRangeType xWorld_Outside = geometry.global( xOutside );
										const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
										// get the integration factor
										const double elementVolume = intersectionGeometry.integrationElement( xLocal );
										// get the quadrature weight
										const double integrationWeight = faceQuadratureNeighbour.weight( quad );
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );

										VelocityRangeType v_i( 0.0 );
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetNeighbour.evaluate( i, xOutside, v_i );
										velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
										VelocityRangeType beta_eval;
										beta_.localFunction(entity).evaluate( xInside, beta_eval );
										const double beta_times_normal =  ( beta_eval * outerNormal );

										//calc u^c_h \tensor beta * v \tensor n (self part), the flux value
										double c_s = (beta_times_normal) * 0.5;
										VelocityRangeType u_h = v_i;
										VelocityJacobianRangeType mean_value = dyadicProduct( u_h, beta_eval );
										mean_value *= 0.5;
										VelocityJacobianRangeType u_jump = dyadicProduct( v_i, outerNormal );
										u_jump *= c_s;
										VelocityJacobianRangeType flux_value = mean_value;
										flux_value += u_jump;

										// \int_{dK} flux_value : ( v_j \ctimes n ) ds
										VelocityJacobianRangeType v_i_tensor_n = dyadicProduct( v_j, outerNormal );
										double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );

										O_i_j += elementVolume
												* integrationWeight
												* convection_scaling
												* ret;
									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( O_i_j ) < eps ) {
										O_i_j = 0.0;
									}
									else {
										// add to matrix
//										std::cerr<< boost::format( "O face value (ne) on entity %d: %e\n") % entityNR % O_i_j;
										localOmatrixNeighbour.add( i, j, O_i_j );
									}
								} // done computing Y's neighbour surface integral
							} // done computing Y's surface integrals

                        //                                                                                                  // we will call this one
                        // (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{-}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's neighbour surface integral
                        //                                                                                                  // see also "Z's boundary integral" below
                        //                                                                                                  // and "Z's volume integral" above
//                        if ( discreteModel_.hasPressureFlux() ) {
                            for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                // compute Z's element surface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                    double Z_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        const double v_i_times_normal = v_i * outerNormal;

										const double p_factor = ( 0.5 - ( D_12 * outerNormal ) );// (0.5 p - p D_12 ) n ) <- p+
//										const double p_factor = ( 0.5 - ( 1 ) );// (0.5 p - p D_12 ) n ) <- p+
                                        const double q_j_times_v_i_times_normal =  q_j * v_i_times_normal;
										Z_i_j += p_factor
                                            * elementVolume
                                            * integrationWeight
											* pressure_gradient_scaling
                                            * q_j_times_v_i_times_normal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Z_i_j ) < eps ) {
                                        Z_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localZmatrixElement.add( i, j, Z_i_j );
                                } // done computing Z's element surface integral
                                // compute Z's neighbour surface integral
                                for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
                                    double Z_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
                                        const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, xInside, v_i );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetNeighbour.evaluate( j, xOutside, q_j );
                                        const double v_i_times_normal = v_i * outerNormal;

										const double p_factor = ( 0.5 + ( D_12 * outerNormal ) );// (0.5 p + p D_12 ) n ) <- p-
//										const double p_factor = ( 0.5 + ( 1 ) );// (0.5 p + p D_12 ) n ) <- p-
                                        const double q_j_times_v_i_times_normal = q_j * v_i_times_normal;
										Z_i_j += p_factor
                                            * elementVolume
                                            * integrationWeight
                                            * pressure_gradient_scaling
                                            * q_j_times_v_i_times_normal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Z_i_j ) < eps ) {
                                        Z_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localZmatrixNeighbour.add( i, j, Z_i_j );
                                } // done computing Z's neighbour surface integral
                            } // done computing Z's surface integrals
//                        }

                        //                                                                                                // we will call this one
                        // (E)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}ds // E's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}ds // E's neighbour surface integral
                        //                                                                                                // see also "E's boundary integral" below
                        //                                                                                                // and "E's volume integral" above
//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                // compute E's element surface integral
                                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                    double E_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
										VelocityRangeType flux_value = v_j;
										flux_value *= 0.5;
										const double v_j_times_outerNormal = v_j * outerNormal;
										VelocityRangeType jump = D_12;
										jump *= v_j_times_outerNormal;
										flux_value += jump;
										VelocityRangeType q_i_times_flux = flux_value;
										q_i_times_flux *= q_i;
										const double q_i_times_flux_times_outerNormal = q_i_times_flux * outerNormal;

										E_i_j += elementVolume
											* integrationWeight
											* q_i_times_flux_times_outerNormal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( E_i_j ) < eps ) {
                                        E_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localEmatrixElement.add( i, j, E_i_j );
                                } // done computing E's element surface integral
                                // compute E's neighbour surface integral
                                for ( int i = 0; i < numPressureBaseFunctionsNeighbour; ++i ) {
                                    double E_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
                                        const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetNeighbour.evaluate( i, xOutside, q_i );

										VelocityRangeType flux_value = v_j;
										flux_value *= 0.5;
										const double v_j_times_outerNormal = v_j * outerNormal;
										VelocityRangeType jump = D_12;
										jump *= v_j_times_outerNormal;
										flux_value += jump;
										VelocityRangeType q_i_times_flux = flux_value;
										q_i_times_flux *= q_i;
										const double q_i_times_flux_times_outerNormal = q_i_times_flux * outerNormal;

										E_i_j -= elementVolume
                                            * integrationWeight
											* q_i_times_flux_times_outerNormal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( E_i_j ) < eps ) {
                                        E_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localEmatrixNeighbour.add( i, j, E_i_j );
                                } // done computing E's neighbour surface integral
                            } // done computing E's surface integrals
//                        }

                        //                                                                                                // we will call this one
                        // (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}ds // R's element surface integral
                        //           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}ds // R's neighbour surface integral
                        //                                                                                                // see also "R's boundary integral" below
//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                // compute R's element surface integral
                                for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
                                    double R_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, x, q_i );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        const double q_i_times_q_j = q_i * q_j;
                                        R_i_j += D_11
                                            * elementVolume
                                            * integrationWeight
                                            * q_i_times_q_j;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( R_i_j ) < eps ) {
                                        R_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localRmatrixElement.add( i, j, R_i_j );
                                } // done computing R's element surface integral
                                // compute R's neighbour surface integral
                                for ( int i = 0; i < numPressureBaseFunctionsNeighbour; ++i ) {
                                    double R_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
                                        const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureNeighbour.weight( quad );
                                        // compute \hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetNeighbour.evaluate( j, xOutside, q_j );
                                        PressureRangeType q_i( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( i, xInside, q_i );
                                        const double q_i_times_q_j = q_i * q_j;
										R_i_j += -1.0
                                            * D_11
                                            * elementVolume
                                            * integrationWeight
                                            * q_i_times_q_j;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( R_i_j ) < eps ) {
                                        R_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localRmatrixNeighbour.add( i, j, R_i_j );
                                } // done computing R's neighbour surface integral
                            } // done computing R's surface integrals
//                        }
                    } // done with those inside the grid

					typedef typename Traits::LocalIntersectionCoordinateType
						LocalIntersectionCoordinateType;
                    // if we are on the boundary of the grid
					if ( !intersection.neighbor() && intersection.boundary() ) {
                        //                                                                                                    // we will call this one
                        // (H1)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}ds // H1's boundary integral
                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
                            for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                double H1_j = 0.0;
                                // sum over all quadrature points
								for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const VelocityRangeType xWorld = geometry.global( x );
                                    const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // compute \hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}
									const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                    SigmaRangeType tau_j( 0.0 );
                                    sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                                    VelocityRangeType tau_j_times_normal( 0.0 );
                                    tau_j.mv( outerNormal, tau_j_times_normal );
                                    VelocityRangeType gD( 0.0 );
									discreteModel_.dirichletData( intersection, 0.0, xWorld,  gD );
                                    const double gD_times_tau_j_times_normal = gD * tau_j_times_normal;
                                    H1_j += elementVolume
                                        * integrationWeight
										* viscosity
                                        * gD_times_tau_j_times_normal;
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H1_j ) < eps ) {
                                    H1_j = 0.0;
                                }
								else
									// add to rhs
									localH1rhs[ j ] += H1_j;
                            } // done computing H1's boundary integral
                        }

                        //                                                                                                                   // we will call this one
                        // (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's boundary integral
                        //                                                                                                                   // see also "X's volume integral", "X's element surface integral" and "X's neighbour surface integral" above
//                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
                                    double X_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        SigmaRangeType tau_j( 0.0 );
                                        sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        VelocityRangeType tau_times_normal( 0.0 );
                                        tau_j.mv( outerNormal, tau_times_normal );
                                        const double v_i_times_tau_times_normal = v_i * tau_times_normal;
                                        X_i_j += -1.0
                                            * elementVolume
                                            * integrationWeight
                                            * v_i_times_tau_times_normal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( X_i_j ) < eps ) {
                                        X_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localXmatrixElement.add( i, j, X_i_j );
                                }
                            } // done computing X's boundary integral
//                        }

                        //                                                                                                           // we will call this one
                        // (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}ds // Y's boundary integral
                        //                                                                                                           // see also "Y's element surface integral" and "Y's neighbour surface integral" above
//                        if ( discreteModel_.hasSigmaFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                    double Y_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute -\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_j( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        const double v_i_times_v_j = v_i * v_j;
                                        Y_i_j += C_11
                                            * elementVolume
                                            * integrationWeight
                                            * v_i_times_v_j;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Y_i_j ) < eps ) {
                                        Y_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localYmatrixElement.add( i, j, Y_i_j );
                                }
                            } // done computing Y's boundary integral
//                        }

							//                                                                                                           // we will call this one
							// (O)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}} STUFF n_{t}ds											// O's boundary integral
							//                                                                                                           // see also "O's element surface integral" and "Y's neighbour surface integral" above
							for ( int i = 0; (i < numVelocityBaseFunctionsElement ) && do_oseen_discretization_; ++i ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									double O_i_j = 0.0;
									// sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
										// get x codim<0> and codim<1> coordinates
										const ElementCoordinateType x = faceQuadratureElement.point( quad );
										const VelocityRangeType xWorld = geometry.global( x );
										const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
										// get the integration factor
										const double elementVolume = intersectionGeometry.integrationElement( xLocal );
										// get the quadrature weight
										const double integrationWeight = faceQuadratureElement.weight( quad );
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
										//calc u^c_h \tensor beta * v \tensor n
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, x, v_j );
										VelocityRangeType v_i( 0.0 );
										velocityBaseFunctionSetElement.evaluate( i, x, v_i );
										VelocityRangeType beta_eval;
										beta_.localFunction(entity).evaluate( x, beta_eval );
										const double beta_times_normal = beta_eval * outerNormal;

										double c_s;
										if ( beta_times_normal < 0 ) {
											c_s = beta_times_normal * 0.5;
										}
										else {
											c_s = - beta_times_normal * 0.5;
										}

										VelocityJacobianRangeType mean_value = dyadicProduct( v_i, beta_eval );
										mean_value *= 0.5;

										VelocityJacobianRangeType u_jump = dyadicProduct( v_i, outerNormal );
										u_jump *= c_s;

										VelocityJacobianRangeType flux_value = mean_value;
										flux_value += u_jump;

										VelocityJacobianRangeType v_i_tensor_n = dyadicProduct( v_j, outerNormal );
										double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );
										//inner edge (self)
										O_i_j += elementVolume
											* integrationWeight
											* convection_scaling
											* ret;

									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( O_i_j ) < eps ) {
										O_i_j = 0.0;
									}
									else {
										// add to matrix
//										std::cerr<< boost::format( "O face value (bnd) on entity %d: %e\n") % entityNR % O_i_j;
										localOmatrixElement.add( i, j, O_i_j );
									}
								}
							} // done computing O's boundary integral
                        //                                                                                                  // we will call this one
                        // (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's boundary integral
                        //                                                                                                  // see also "Z's volume integral", "Z's element surface integral" and "Z's neighbour surface integral" above
//                        if ( discreteModel_.hasPressureFlux() ) {
                            for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
                                // compute the boundary integral
                                for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                    double Z_i_j = 0.0;
                                    // sum over all quadrature points
									for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                        // get x codim<0> and codim<1> coordinates
                                        const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                        const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                        // get the integration factor
                                        const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                        // get the quadrature weight
                                        const double integrationWeight = faceQuadratureElement.weight( quad );
                                        // compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
										const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                        VelocityRangeType v_i( 0.0 );
                                        velocityBaseFunctionSetElement.evaluate( i, x, v_i );
                                        PressureRangeType q_j( 0.0 );
                                        pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                        const double v_i_times_normal = v_i * outerNormal;
                                        const double q_j_times_v_times_normal = q_j * v_i_times_normal;
                                        Z_i_j += elementVolume
                                            * integrationWeight
                                            * pressure_gradient_scaling
                                            * q_j_times_v_times_normal;
                                    } // done sum over all quadrature points
                                    // if small, should be zero
                                    if ( fabs( Z_i_j ) < eps ) {
                                        Z_i_j = 0.0;
                                    }
                                    else
                                        // add to matrix
                                        localZmatrixElement.add( i, j, Z_i_j );
                                }
                            } // done computing Z's boundary integral
//                        }

                        //                                                                                                                 // we will call this one
                        // (H2)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left( \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}ds         // H2's 1st boundary integral
                        //                                                         -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}ds        \right) // H2's 2nd boundary integral
                        //                                                                                                                 // see also "H2's volume integral" above
                        if ( discreteModel_.hasSigmaFlux() && discreteModel_.hasPressureFlux() ) {
                            for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
                                double H2_j = 0.0;
                                // sum over all quadrature points
								for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                    const VelocityRangeType globalX = geometry.global( x );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // prepare
									const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                    VelocityRangeType v_j( 0.0 );
                                    velocityBaseFunctionSetElement.evaluate( j, x, v_j );
                                    // compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
//                                    if ( discreteModel_.hasSigmaFlux() ) {
										const VelocityRangeType xIntersectionGlobal = intersection.intersectionSelfLocal().global( xLocal );
                                        const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
                                        VelocityRangeType gD( 0.0 );
										discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
                                        SigmaRangeType gD_times_normal( 0.0 );
                                        gD_times_normal = dyadicProduct( gD, outerNormal );
                                        VelocityRangeType gD_times_normal_times_normal( 0.0 );
                                        gD_times_normal.mv( outerNormal, gD_times_normal_times_normal );
                                        const double v_j_times_gD_times_normal_times_normal= v_j * gD_times_normal_times_normal;
                                        H2_j += C_11
                                            * elementVolume
                                            * integrationWeight
                                            * v_j_times_gD_times_normal_times_normal;
//                                    }
                                    // done computing H2's 1st boundary integral
                                    // compute -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}
//                                    if ( discreteModel_.hasPressureFlux() ) {
                                        const double v_j_times_normal = v_j * outerNormal;
                                        const double flux_times_v_j_times_n_t = 0.0 * v_j_times_normal;
                                        H2_j += -1.0
                                            * elementVolume
                                            * integrationWeight
											* pressure_gradient_scaling
                                            * flux_times_v_j_times_n_t;
//                                    }
                                    // done computing H2's 2nd boundary integral
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H2_j ) < eps ) {
                                    H2_j = 0.0;
                                }
								else
									// add to rhs
									localH2rhs[ j ] += H2_j;
                            } // done computing H2's boundary integrals
                        }

						//                                                                                                                 // we will call this one
						// (H2_O)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left(  \beta n_{T} g_D v_j ds        \right) // H2_O's boundary integral
						for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
							double H2_O_j = 0.0;
							// sum over all quadrature points
							for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
								// get x codim<0> and codim<1> coordinates
								const ElementCoordinateType x = faceQuadratureElement.point( quad );
								const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
												// get the integration factor
								const double elementVolume = intersectionGeometry.integrationElement( xLocal );
								// get the quadrature weight
								const double integrationWeight = faceQuadratureElement.weight( quad );
								// prepare
								const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
								VelocityRangeType v_j( 0.0 );
								velocityBaseFunctionSetElement.evaluate( j, x, v_j );
								// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
								const VelocityRangeType xIntersectionGlobal = intersection.intersectionSelfLocal().global( xLocal );
								const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
								VelocityRangeType gD( 0.0 );
								discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );

								VelocityRangeType beta_eval;
								beta_.localFunction(entity).evaluate( x, beta_eval );
								const double beta_times_normal = beta_eval * outerNormal;

								// u^c = 0.5 gD \otimes beta + Cs -gD \otimes n
								VelocityJacobianRangeType gD_tensor_beta = dyadicProduct( gD, beta_eval );
								gD_tensor_beta *= 0.5;
								double c_s;
								if ( beta_times_normal < 0 ) {
									c_s = beta_times_normal * 0.5;
								}
								else {
									c_s = - beta_times_normal * 0.5;
								}
								VelocityJacobianRangeType jump = dyadicProduct( gD, outerNormal );
								jump *= c_s;

								VelocityJacobianRangeType flux_value = gD_tensor_beta;
								flux_value += jump;

								VelocityJacobianRangeType v_j_tensor_n = dyadicProduct( v_j,  outerNormal );
								const double ret = Stuff::colonProduct( flux_value, v_j_tensor_n );
								H2_O_j += elementVolume
										* convection_scaling
										* integrationWeight
										* ret;
							}
							if ( fabs( H2_O_j ) < eps ) {
									 H2_O_j = 0.0;
							}
							else {
								localH2_O_rhs[ j ] += H2_O_j;
							}
						}

                        //                                                                                        // we will call this one
                        // (H3)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}-\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}ds // H3's boundary integral
//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
                            for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
                                double H3_j = 0.0;
                                // sum over all quadrature points
								for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
                                    // get x codim<0> and codim<1> coordinates
                                    const ElementCoordinateType x = faceQuadratureElement.point( quad );
                                    const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
                                    const VelocityRangeType xWorld = geometry.global( x );
                                    // get the integration factor
                                    const double elementVolume = intersectionGeometry.integrationElement( xLocal );
                                    // get the quadrature weight
                                    const double integrationWeight = faceQuadratureElement.weight( quad );
                                    // compute -\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}
									const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
                                    VelocityRangeType gD( 0.0 );
									discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
                                    const double gD_times_normal = gD * outerNormal;
                                    PressureRangeType q_j( 0.0 );
                                    pressureBaseFunctionSetElement.evaluate( j, x, q_j );
                                    const double q_j_times_gD_times_normal = q_j * gD_times_normal;
                                    H3_j += elementVolume
                                        * integrationWeight
                                        * q_j_times_gD_times_normal;
                                } // done sum over all quadrature points
                                // if small, should be zero
                                if ( fabs( H3_j ) < eps ) {
                                    H3_j = 0.0;
                                }
								else
									// add to rhs
									localH3rhs[ j ] += H3_j;
                            } // done computing H3's boundary integral
//                        }

                    } // done with those on the boundary
                } // done walking the neighbours
//#endif //no_surface_ints

            } // done walking the grid


//            // do the matlab logging stuff
			if ( Parameters().getParam( "save_matrices", false ) ) {
				Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
				Stuff::printSparseRowMatrixMatlabStyle( MInversMatrix.matrix(), "M_invers", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Wmatrix.matrix(), "W", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Omatrix.matrix(), "O", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Xmatrix.matrix(), "X", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Ymatrix.matrix(), "Y", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Zmatrix.matrix(), "Z", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Ematrix.matrix(), "E", matlabLogStream );
				Stuff::printSparseRowMatrixMatlabStyle( Rmatrix.matrix(), "R", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( H1rhs, "H1", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( H2rhs, "H2", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( H3rhs, "H3", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( H2_O_rhs, "H_O", matlabLogStream );
				matlabLogStream.Flush();
			}
//
//            // log local matrices
//            Stuff::GridWalk<GridPartType> gw( gridPart_ );
//            typedef Logging::MatlabLogStream
//                FunctorStream;
//            FunctorStream& functorStream = matlabLogStream;
//            Stuff::GridWalk<DiscreteVelocityFunctionSpaceType> gw( velocitySpace_ );
//            Stuff::LocalMatrixPrintFunctor< EmatrixType,FunctorStream> f_E ( Ematrix, functorStream, "E" );
//            Stuff::LocalMatrixPrintFunctor< WmatrixType,FunctorStream> f_W ( Wmatrix, functorStream, "W" );
//            Stuff::LocalMatrixPrintFunctor< XmatrixType,FunctorStream> f_X ( Xmatrix, functorStream, "X" );
//            Stuff::LocalMatrixPrintFunctor< YmatrixType,FunctorStream> f_Y ( Ymatrix, functorStream, "Y" );
//            Stuff::LocalMatrixPrintFunctor< ZmatrixType,FunctorStream> f_Z ( Zmatrix, functorStream, "Z" );
//            Stuff::LocalMatrixPrintFunctor< RmatrixType,FunctorStream> f_R ( Rmatrix, functorStream, "R" );
//                gw( f_W );
//                gw( f_X );
//                gw( f_Y );
//                gw( f_Z );
//                gw( f_E );
//                gw( f_R );
            // do profiling
            profiler().StopTiming("Pass -- ASSEMBLE");

			if ( Parameters().getParam( "outputMatrixPlots", false ) ) {
                Stuff::matrixToGnuplotFile( Ematrix.matrix(),       std::string( "mat_E.gnuplot")       );
                Stuff::matrixToGnuplotFile( Wmatrix.matrix(),       std::string( "mat_W.gnuplot")       );
                Stuff::matrixToGnuplotFile( Xmatrix.matrix(),       std::string( "mat_X.gnuplot")       );
                Stuff::matrixToGnuplotFile( Ymatrix.matrix(),       std::string( "mat_Y.gnuplot")       );
                Stuff::matrixToGnuplotFile( Zmatrix.matrix(),       std::string( "mat_Z.gnuplot")       );
                Stuff::matrixToGnuplotFile( Rmatrix.matrix(),       std::string( "mat_R.gnuplot")       );
                Stuff::matrixToGnuplotFile( MInversMatrix.matrix(), std::string( "mat_M.gnuplot")   );
            }

            profiler().StartTiming("Pass -- SOLVER");

			if ( Parameters().getParam( "paranoid_checks", false ) )
			{//paranoid checks
				assert( !Stuff::MatrixContainsNanOrInf( Omatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Ematrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Wmatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Xmatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Ymatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Zmatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( Rmatrix.matrix() ) );
				assert( !Stuff::MatrixContainsNanOrInf( MInversMatrix.matrix() ) );
				assert( !Stuff::FunctionContainsNanOrInf( H1rhs ) );
				assert( !Stuff::FunctionContainsNanOrInf( H2rhs ) );
				assert( !Stuff::FunctionContainsNanOrInf( H3rhs ) );
				assert( !Stuff::FunctionContainsNanOrInf( H2_O_rhs ) );
//				Zmatrix.matrix().scale( -1 );
//				assert( areTransposed( Zmatrix.matrix(), Ematrix.matrix() ));
//				Zmatrix.matrix().scale( -1 );
			}

			Logger().Info().Resume();
			Logger().Info() << "Solving system with " << dest.discreteVelocity().size() << " + " << dest.discretePressure().size() << " unknowns" << std::endl;

			// do solving
			if ( do_oseen_discretization_  ) {
				H2rhs += H2_O_rhs;
			}
			//this lets us switch between standalone oseen and reduced oseen in  thete scheme easily
			const bool use_reduced_solver = do_oseen_discretization_ && Parameters().getParam( "reduced_oseen_solver", false );
			typedef SolverCaller< ThisType >
				SolverCallerType;

			//Select which solver we want to use
			typename SolverCallerType::SolverID solver_ID = SolverCallerType::SaddlePoint_Solver_ID;
			if( !use_reduced_solver ) {
				if ( Parameters().getParam( "use_nested_cg_solver", false ) )
					solver_ID = SolverCallerType::NestedCG_Solver_ID;
				else if ( Parameters().getParam( "use_full_solver", false ) )
					solver_ID = SolverCallerType::FullSystem_Solver_ID;
			}
			else
				solver_ID = SolverCallerType::Reduced_Solver_ID;

			info_ = SolverCallerType::solve( dest, rhs_datacontainer, solver_ID, arg, Xmatrix, MInversMatrix,
											Ymatrix, Omatrix, Ematrix, Rmatrix,
											Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs, beta_ );

            // do profiling
            profiler().StopTiming("Pass -- SOLVER");

			if ( Parameters().getParam( "save_matrices", false ) ) {
				Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
				Stuff::printDiscreteFunctionMatlabStyle( dest.discreteVelocity(), "u_computed", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( dest.discretePressure(), "p_computed", matlabLogStream );
			}
        } // end of apply

        virtual void compute( const TotalArgumentType& /*arg*/, DestinationType& /*dest*/ ) const
        {}

        virtual void allocateLocalMemory()
        {}

#ifdef HAS_RUN_INFO
        void getRuninfo( RunInfo& info )
        {
			info.iterations_inner_avg = int( info_.iterations_inner_avg );
            info.iterations_inner_min = info_.iterations_inner_min;
            info.iterations_inner_max = info_.iterations_inner_max;
            info.iterations_outer_total = info_.iterations_outer_total;
            info.max_inner_accuracy = info_.max_inner_accuracy;
        }
#endif

    private:
        DiscreteModelType& discreteModel_;
		const typename Traits::GridPartType& gridPart_;
		const typename Traits::DiscreteStokesFunctionSpaceWrapperType& spaceWrapper_;
		const typename Traits::DiscreteVelocityFunctionSpaceType& velocitySpace_;
		const typename Traits::DiscretePressureFunctionSpaceType& pressureSpace_;
		typename Traits::DiscreteSigmaFunctionSpaceType sigmaSpace_;
		const typename Traits::DiscreteVelocityFunctionType& beta_;
		const bool do_oseen_discretization_;
        mutable SaddlepointInverseOperatorInfo info_;

        /**
         *  \brief  dyadic product
         *
         *          Implements \f$\left(arg_{1} \otimes arg_{2}\right)_{i,j}:={arg_{1}}_{i} {arg_{2}}_{j}\f$
         **/
		static typename Traits::SigmaRangeType dyadicProduct(
										const typename Traits::VelocityRangeType& arg1,
										const typename Traits::VelocityRangeType& arg2 )
        {
			typename Traits::SigmaRangeType ret( 0.0 );
			typedef typename Traits::SigmaRangeType::RowIterator
                MatrixRowIteratorType;
			typedef typename Traits::VelocityRangeType::ConstIterator
                ConstVectorIteratorType;
			typedef typename Traits::VelocityRangeType::Iterator
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

        // VelocityRangeType is expected to be a FieldVector,
        // SigmaJacobianRangeType to be a Matrixmapping and
        // SigmaJacobianRangeType[i] to be a FieldVector
        //! \todo   doc me
		static typename Traits::SigmaJacobianRangeType prepareVelocityRangeTypeForSigmaDivergence(
													const typename Traits::VelocityRangeType& arg )
        {
			typename Traits::SigmaJacobianRangeType ret( 0.0 );
            assert( arg.dim() == ret[0].dim() );
            for ( int i = 0; i < int(arg.dim()) ; ++i ) {
                for ( int j = 0; j < int(arg.dim()); ++j ) {
					typename Traits::VelocityRangeType row( 0.0 );
                    row[ j ] = arg[ i ];
                    ret[ i * arg.dim() + j ] = row;
                }
            }
            return ret;
        }

        //! \todo   doc me
		static typename Traits::VelocityJacobianRangeType preparePressureRangeTypeForVelocityDivergence(
													const typename Traits::PressureRangeType& arg )
        {
			typename Traits::VelocityJacobianRangeType ret( 0.0 );
            for ( unsigned int i = 0; i < ret[0].dim(); ++i ) {
				typename Traits::VelocityRangeType row( 0.0 );
                row[ i ] = arg;
                ret[ i ] = row;
            }
            return ret;
        }
};

}
#endif  // end of stokespass.hh
