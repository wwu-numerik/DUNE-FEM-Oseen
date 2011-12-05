/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <cmake_config.h>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/stokes/defaulttraits.hh>
#include <dune/stokes/solver/solvercaller.hh>
#include <dune/stokes/integrators/all.hh>

#include <dune/common/stdstreams.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/stuff/matrix.hh>
#include <dune/stuff/tuple.hh>

#ifdef STOKES_USE_ISTL
#   include <dune/fem/operator/matrix/istlmatrix.hh>
#else
#   include <dune/fem/operator/matrix/spmatrix.hh>
#endif

#include <dune/stuff/progressbar.hh>

template < class FunctionSpaceImp, class TimeProviderImp >
class VelocityConvection : public Dune::TimeFunction < FunctionSpaceImp , VelocityConvection< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
	public:
		typedef VelocityConvection< FunctionSpaceImp, TimeProviderImp >
			ThisType;
		typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;

		VelocityConvection(	const TimeProviderImp& timeprovider,
					const FunctionSpaceImp& space,
					const double parameter_a = M_PI /2.0 ,
					const double parameter_d = M_PI /4.0)
			: BaseType( timeprovider, space ),
			parameter_a_( parameter_a ),
			parameter_d_( parameter_d )
		{}

		~VelocityConvection()
		{}

		void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
		{
			const double x			= arg[0];
			const double y			= arg[1];
			ret[0] = std::pow(time,5.0)*2*x*y;
			ret[1] = std::pow(time,5.0)*y*y;;
		}

	private:
		static const int dim_ = FunctionSpaceImp::dimDomain ;
		const double parameter_a_;
		const double parameter_d_;
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

		//! when requested we store \f$ \varDelta u, \nabla p (u \cdot \nabla ) u\f$ in this struct after the solver
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
			do_oseen_discretization_( do_oseen_discretization ),
	      info_( SaddlepointInverseOperatorInfo() )
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
		void apply( const DomainType &arg, RangeType &dest, RhsDatacontainerType* rhs_datacontainer, const ExactSigmaType* /*sigma_exact*/ ) const
        {
            // profiler information
            profiler().StartTiming("Pass_init");

            const bool verbose_reserve = true;
            // matrices
			typedef typename Traits::DiscreteSigmaFunctionSpaceType
				DiscreteSigmaFunctionSpaceType;
            typedef typename Traits::DiscreteVelocityFunctionSpaceType
                DiscreteVelocityFunctionSpaceType;
            typedef typename Traits::DiscretePressureFunctionSpaceType
                DiscretePressureFunctionSpaceType;
            typedef Stokes::Integrators::Factory< Traits >
                Factory;

            // M\in R^{M\times M}
            typedef typename Factory::MInversMatrixIntegratorType
                MInversMatrixIntegratorType;
            auto MInversMatrix = Factory::matrix( sigmaSpace_, sigmaSpace_ );
            assert( MInversMatrix->matrix().rows() == MInversMatrix->matrix().cols() );
            // W\in R^{M\times L}
            auto Wmatrix = Factory::matrix( velocitySpace_, sigmaSpace_ );

            // X\in R^{L\times M}

            typedef typename Factory::XmatrixTypeIntegratorType
				XmatrixTypeIntegratorType;
            auto Xmatrix = Factory::matrix( sigmaSpace_ , velocitySpace_ );

            // Y\in R^{L\times L}
            typedef typename Factory::YmatrixTypeIntegratorType
				YmatrixTypeIntegratorType;
            auto Ymatrix = Factory::matrix( velocitySpace_, velocitySpace_ );

            typedef typename Factory::YmatrixTypeIntegratorType
				OmatrixTypeIntegratorType;
            auto Omatrix = Factory::matrix( velocitySpace_, velocitySpace_ );

            // Z\in R^{L\times K}
            typedef typename Factory::ZmatrixTypeIntegratorType
				ZmatrixTypeIntegratorType;
            auto Zmatrix = Factory::matrix( pressureSpace_, velocitySpace_ );

            // E\in R^{K\times L}
            typedef typename Factory::EmatrixTypeIntegratorType
				EmatrixTypeIntegratorType;
            auto Ematrix = Factory::matrix( velocitySpace_, pressureSpace_ );

            // R\in R^{K\times K}
            typedef typename Factory::RmatrixTypeIntegratorType
				RmatrixTypeIntegratorType;
            auto Rmatrix = Factory::matrix( pressureSpace_, pressureSpace_ );


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
            typedef tuple<	typename Factory::MInversMatrixIntegratorType,
                            typename Factory::WmatrixTypeIntegratorType,
                            typename Factory::XmatrixTypeIntegratorType,
                            typename Factory::YmatrixTypeIntegratorType,
                            typename Factory::OmatrixTypeIntegratorType,
                            typename Factory::ZmatrixTypeIntegratorType,
                            typename Factory::EmatrixTypeIntegratorType,
                            typename Factory::RmatrixTypeIntegratorType,
                            H1_IntegratorType,
							H2_IntegratorType,
							H2_O_IntegratorType,
							H3_IntegratorType >
				OseenIntegratorTuple;
            typedef tuple<	typename Factory::MInversMatrixIntegratorType,
                            typename Factory::WmatrixTypeIntegratorType,
                            typename Factory::XmatrixTypeIntegratorType,
                            typename Factory::YmatrixTypeIntegratorType,
                            typename Factory::ZmatrixTypeIntegratorType,
                            typename Factory::EmatrixTypeIntegratorType,
                            typename Factory::RmatrixTypeIntegratorType,
							H1_IntegratorType,
							H2_IntegratorType,
							H3_IntegratorType >
				StokesIntegratorTuple;

            auto m_integrator = Factory::integrator( MInversMatrix );
            auto w_integrator = Factory::integrator( Wmatrix );
            auto x_integrator = Factory::integrator( Xmatrix );
            auto y_integrator = Factory::integrator( Ymatrix );
            auto o_integrator = Factory::integratorO( Omatrix, beta_ );
            auto z_integrator = Factory::integrator( Zmatrix );
            auto e_integrator = Factory::integrator( Ematrix );
            auto r_integrator = Factory::integrator( Rmatrix );
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
				Logger().Dbg() << "Oseen disc\n" ;
			}
			else
			{
				Stokes::Integrators::Coordinator< Traits, StokesIntegratorTuple >
						coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

				StokesIntegratorTuple tuple(	m_integrator, w_integrator, x_integrator, y_integrator,
										z_integrator, e_integrator, r_integrator,
										h1_integrator, h2_integrator,h3_integrator );
				coordinator.apply( tuple );
				Logger().Dbg() << "Stokes disc\n" ;
			}

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
			typename Stokes::Solver::SolverID solver_ID = do_oseen_discretization_
					? Stokes::Solver::BiCg_Saddlepoint_Solver_ID
					: Stokes::Solver::SaddlePoint_Solver_ID;

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
				Stuff::Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
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
		void getRuninfo( Stuff::RunInfo& info )
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
			Stuff::Logging::LogStream& infoStream = Logger().Info();
			infoStream << boost::format( "pressure_gradient/convection scaling: %e | %e\npass viscosity: %e\n")
								% discreteModel_.pressure_gradient_scaling()
								% discreteModel_.convection_scaling()
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
