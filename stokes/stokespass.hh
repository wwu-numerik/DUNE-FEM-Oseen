/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/stokes/defaulttraits.hh>
#include <dune/stokes/solver/solvercaller.hh>
#include <dune/stokes/integrators/all.hh>

#include <dune/common/stdstreams.hh>
#include <dune/stuff/matrix.hh>
#include <dune/stuff/tuple.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stuff/progressbar.hh>

template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
	struct StencilType {
		static int nonZerosEstimate( const ColSpaceImp& rangeSpace ) {
			return rangeSpace.maxNumLocalDofs() * 1.5f;
		}
	};
};

template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct DiagonalMatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
	struct StencilType {
		static int nonZerosEstimate( const ColSpaceImp& ) {
			return 1;
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
		/**
		 *  \brief  empty constructor
		 **/
		StokesPass()
		{}

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
			void clear() {
				velocity_laplace.clear();
				pressure_gradient.clear();
				velocity_gradient.clear();
				convection.clear();
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
			profiler().StartTiming("Pass_init");

            // matrices
            // M\in R^{M\times M}
			typedef typename Traits::DiscreteSigmaFunctionSpaceType
				DiscreteSigmaFunctionSpaceType;
			typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
											DiscreteSigmaFunctionSpaceType,
											DiagonalMatrixTraits<DiscreteSigmaFunctionSpaceType,DiscreteSigmaFunctionSpaceType> >
                MInversMatrixType;
			typedef Stokes::Integrators::M< MInversMatrixType, Traits >
				MInversMatrixIntegratorType;
            MInversMatrixType MInversMatrix( sigmaSpace_, sigmaSpace_ );
            MInversMatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( MInversMatrix, Logger().Dbg(), "M" );
            assert( MInversMatrix.matrix().rows() == MInversMatrix.matrix().cols() );
            // W\in R^{M\times L}
			typedef typename Traits::DiscreteVelocityFunctionSpaceType
				DiscreteVelocityFunctionSpaceType;
            typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType> >
                WmatrixType;
			typedef Stokes::Integrators::W< WmatrixType, Traits >
				WmatrixTypeIntegratorType;
            WmatrixType Wmatrix( sigmaSpace_, velocitySpace_ );
            Wmatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Wmatrix, Logger().Dbg(), "W" );
            // X\in R^{L\times M}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteSigmaFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType> >
                XmatrixType;
			typedef Stokes::Integrators::X< XmatrixType, Traits >
				XmatrixTypeIntegratorType;
            XmatrixType Xmatrix( velocitySpace_, sigmaSpace_ );
            Xmatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Xmatrix, Logger().Dbg(), "X" );
            // Y\in R^{L\times L}
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                YmatrixType;
			typedef Stokes::Integrators::Y< YmatrixType, Traits >
				YmatrixTypeIntegratorType;
            YmatrixType Ymatrix( velocitySpace_, velocitySpace_ );
            Ymatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Ymatrix, Logger().Dbg(), "Y" );
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                OmatrixType;
			typedef Stokes::Integrators::O< OmatrixType, Traits, typename Traits::DiscreteVelocityFunctionType >
				OmatrixTypeIntegratorType;
			OmatrixType Omatrix( velocitySpace_, velocitySpace_ );
			Omatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Omatrix, Logger().Dbg(), "O" );
            // Z\in R^{L\times K}
			typedef typename Traits::DiscretePressureFunctionSpaceType
				DiscretePressureFunctionSpaceType;
            typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
											DiscretePressureFunctionSpaceType,
											MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscretePressureFunctionSpaceType> >
                ZmatrixType;
			typedef Stokes::Integrators::Z< ZmatrixType, Traits >
				ZmatrixTypeIntegratorType;
            ZmatrixType Zmatrix( velocitySpace_, pressureSpace_ );
            Zmatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Zmatrix, Logger().Dbg(), "Z" );
            // E\in R^{K\times L}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
											DiscreteVelocityFunctionSpaceType,
											MatrixTraits<DiscretePressureFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
                EmatrixType;
			typedef Stokes::Integrators::E< EmatrixType, Traits >
				EmatrixTypeIntegratorType;
            EmatrixType Ematrix( pressureSpace_, velocitySpace_ );
            Ematrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Ematrix, Logger().Dbg(), "E" );
            // R\in R^{K\times K}
            typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
											DiscretePressureFunctionSpaceType,
											MatrixTraits<DiscretePressureFunctionSpaceType,DiscretePressureFunctionSpaceType> >
                RmatrixType;
			typedef Stokes::Integrators::R< RmatrixType, Traits >
				RmatrixTypeIntegratorType;
            RmatrixType Rmatrix( pressureSpace_, pressureSpace_ );
            Rmatrix.reserve();
//			Stuff::Matrix::printMemUsageObject( Rmatrix, Logger().Dbg(), "R" );

            // right hand sides
            // H_{1}\in R^{M}
			typename Traits::DiscreteSigmaFunctionType H1rhs( "H1", sigmaSpace_ );
			typedef Stokes::Integrators::H1< typename Traits::DiscreteSigmaFunctionType, Traits >
				H1_IntegratorType;
            H1rhs.clear();
            // H_{2}\in R^{L}
			typename Traits::DiscreteVelocityFunctionType H2rhs( "H2", velocitySpace_ );
			typedef Stokes::Integrators::H2< typename Traits::DiscreteVelocityFunctionType , Traits >
				H2_IntegratorType;
            H2rhs.clear();
			typename Traits::DiscreteVelocityFunctionType H2_O_rhs( "H2_O", velocitySpace_ );
			typedef Stokes::Integrators::H2_O< typename Traits::DiscreteVelocityFunctionType , Traits, typename Traits::DiscreteVelocityFunctionType >
				H2_O_IntegratorType;
			H2_O_rhs.clear();
            // H_{3}\in R^{K}
			typename Traits::DiscretePressureFunctionType H3rhs( "H3", pressureSpace_ );
			typedef Stokes::Integrators::H3< typename Traits::DiscretePressureFunctionType, Traits >
				H3_IntegratorType;
            H3rhs.clear();

			profiler().StopTiming("Pass_init");

			//because of the 9-element limit in dune tuples i have to split the assembly in two...
			typedef Tuple<	MInversMatrixIntegratorType,
							WmatrixTypeIntegratorType,
							XmatrixTypeIntegratorType,
							YmatrixTypeIntegratorType,
							OmatrixTypeIntegratorType,
							ZmatrixTypeIntegratorType,
							EmatrixTypeIntegratorType,
							RmatrixTypeIntegratorType,
							H1_IntegratorType,
							H2_IntegratorType,
							H2_O_IntegratorType,
							H3_IntegratorType >
				OseenIntegratorTuple;
			typedef Tuple<	MInversMatrixIntegratorType,
							WmatrixTypeIntegratorType,
							XmatrixTypeIntegratorType,
							YmatrixTypeIntegratorType,
							ZmatrixTypeIntegratorType,
							EmatrixTypeIntegratorType,
							RmatrixTypeIntegratorType,
							H1_IntegratorType,
							H2_IntegratorType,
							H3_IntegratorType >
				StokesIntegratorTuple;

			MInversMatrixIntegratorType m_integrator( MInversMatrix );
			WmatrixTypeIntegratorType	w_integrator( Wmatrix );
			XmatrixTypeIntegratorType	x_integrator( Xmatrix );
			YmatrixTypeIntegratorType	y_integrator( Ymatrix );
			OmatrixTypeIntegratorType	o_integrator( Omatrix, beta_ );
			ZmatrixTypeIntegratorType	z_integrator( Zmatrix );
			EmatrixTypeIntegratorType	e_integrator( Ematrix );
			RmatrixTypeIntegratorType	r_integrator( Rmatrix );
			H1_IntegratorType			h1_integrator( H1rhs );
			H2_IntegratorType			h2_integrator( H2rhs );
			H2_O_IntegratorType			h2_o_integrator( H2rhs, beta_ );
			H3_IntegratorType			h3_integrator( H3rhs );
			if ( do_oseen_discretization_ )
			{
				Stokes::Integrators::Coordinator< Traits, OseenIntegratorTuple >
						coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

				OseenIntegratorTuple tuple(	m_integrator, w_integrator, x_integrator, y_integrator,
										o_integrator, z_integrator, e_integrator, r_integrator,
										h1_integrator, h2_integrator,h2_o_integrator, h3_integrator );
				coordinator.apply( tuple );
			}
			else
			{
				Stokes::Integrators::Coordinator< Traits, StokesIntegratorTuple >
						coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

				StokesIntegratorTuple tuple(	m_integrator, w_integrator, x_integrator, y_integrator,
										z_integrator, e_integrator, r_integrator,
										h1_integrator, h2_integrator,h3_integrator );
				coordinator.apply( tuple );
			}


		#ifndef NDEBUG
			{
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


			if ( Parameters().getParam( "outputMatrixPlots", false ) ) {
                Stuff::matrixToGnuplotFile( Ematrix.matrix(),       std::string( "mat_E.gnuplot")       );
                Stuff::matrixToGnuplotFile( Wmatrix.matrix(),       std::string( "mat_W.gnuplot")       );
                Stuff::matrixToGnuplotFile( Xmatrix.matrix(),       std::string( "mat_X.gnuplot")       );
                Stuff::matrixToGnuplotFile( Ymatrix.matrix(),       std::string( "mat_Y.gnuplot")       );
                Stuff::matrixToGnuplotFile( Zmatrix.matrix(),       std::string( "mat_Z.gnuplot")       );
                Stuff::matrixToGnuplotFile( Rmatrix.matrix(),       std::string( "mat_R.gnuplot")       );
                Stuff::matrixToGnuplotFile( MInversMatrix.matrix(), std::string( "mat_M.gnuplot")   );
            }

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
			}
		#endif

			Logger().Info().Resume();
			Logger().Info() << "Solving system with " << dest.discreteVelocity().size() << " + " << dest.discretePressure().size() << " unknowns" << std::endl;

			// do solving

			//this lets us switch between standalone oseen and reduced oseen in  thete scheme easily
			const bool use_reduced_solver = (do_oseen_discretization_ && Parameters().getParam( "reduced_oseen_solver", false ))
					|| Parameters().getParam( "parabolic", false );
			typedef Stokes::SolverCaller< ThisType>
				SolverCallerType;
			typedef Stokes::SolverCaller< ThisType, SmartReconstruction >
				SmartSolverCallerType;

			//Select which solver we want to use
//			typename Stokes::Solver::SolverID solver_ID = do_oseen_discretization_
//					? Stokes::Solver::BiCg_Saddlepoint_Solver_ID
//					: Stokes::Solver::SaddlePoint_Solver_ID;
			typename Stokes::Solver::SolverID solver_ID = Stokes::Solver::BiCg_Saddlepoint_Solver_ID;

			if( !use_reduced_solver ) {
				if ( Parameters().getParam( "use_nested_cg_solver", false ) )
					solver_ID = Stokes::Solver::NestedCG_Solver_ID;
				else if ( Parameters().getParam( "use_full_solver", false ) )
					solver_ID = Stokes::Solver::FullSystem_Solver_ID;
			}
			else
				solver_ID = Stokes::Solver::Reduced_Solver_ID;

			if ( Parameters().getParam( "smart_reconstruction", false ) )
				info_ = SmartSolverCallerType::solve(dest, rhs_datacontainer, solver_ID, do_oseen_discretization_,
												arg, Xmatrix, MInversMatrix, Ymatrix, Omatrix, Ematrix,
												Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs, beta_ );
			else
				info_ = SolverCallerType::solve(dest, rhs_datacontainer, solver_ID, do_oseen_discretization_,
												arg, Xmatrix, MInversMatrix, Ymatrix, Omatrix, Ematrix,
												Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs, beta_ );

		#ifndef NDEBUG
			if ( Parameters().getParam( "save_matrices", false ) ) {
				Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
				Stuff::printDiscreteFunctionMatlabStyle( dest.discreteVelocity(), "u_computed", matlabLogStream );
				Stuff::printDiscreteFunctionMatlabStyle( dest.discretePressure(), "p_computed", matlabLogStream );
			}
			if ( Parameters().getParam( "paranoid_checks", false ) )
			{//paranoid checks
				assert( !Stuff::FunctionContainsNanOrInf( dest.discretePressure() ) );
				assert( !Stuff::FunctionContainsNanOrInf( dest.discreteVelocity() ) );
			}
		#endif
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

	public:
		void printInfo() const
		{
#ifndef NLOG
			Logging::LogStream& infoStream = Logger().Info();
			infoStream << boost::format( "pressure_gradient/convection scaling: %e | %e\npass viscosity: %e\n")
								% discreteModel_.convection_scaling()
								% discreteModel_.pressure_gradient_scaling()
								% discreteModel_.viscosity();
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
		}
};

}
#endif  // end of stokespass.hh
