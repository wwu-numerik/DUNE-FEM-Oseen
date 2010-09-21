#ifndef DUNE_STOKES_PROBLEMS_HH
#define DUNE_STOKES_PROBLEMS_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/misc.hh>

//! collection of data functions
namespace StokesProblems {
	//! a simple problem, d'oh
	namespace Simple {
		#include "problems/simple.hh"
	}
	//! cockburn, see ref
	namespace Cockburn {
		#include "problems/cockburn.hh"
	}
	//! docme
	namespace Constant {
		#include "problems/constant.hh"
	}
	//! docme
	namespace Rotate {
		#include "problems/rotate.hh"
	}
	//! docme
	namespace Aorta {
		#include "problems/aorta.hh"
	}

	/**
	 *  \brief  a collection of some analytical functions solving a stokes problem
	 *
	 *  namely velocity, pressure, force term and dirichlet boundary data
	 *
	 *  \tparam gridDim
	 *          dimension of the grid
	 *
	 *  \todo   extensive docu with latex
		\note super useless
	 **/
	template < int griddim, class DiscreteFunctionWrapperImp  >
	class Container
	{
		public:
			static const unsigned int gridDim = griddim;
			static const bool hasMeaningfulAnalyticalSolution = PROBLEM_NAMESPACE::hasExactSolution;

			typedef DiscreteFunctionWrapperImp
				DiscreteFunctionWrapperType;

			typedef typename DiscreteFunctionWrapperType::DiscreteFunctionSpaceType
				DiscreteFunctionSpaceWrapperType;

			typedef typename DiscreteFunctionSpaceWrapperType
					::DiscreteVelocityFunctionSpaceType
					::FunctionSpaceType
			  VelocityFunctionSpaceType;

			typedef PROBLEM_NAMESPACE::Velocity< VelocityFunctionSpaceType >
				VelocityType;

			typedef typename DiscreteFunctionSpaceWrapperType
					::DiscretePressureFunctionSpaceType
					::FunctionSpaceType
				PressureFunctionSpaceType;

			typedef PROBLEM_NAMESPACE::Pressure< PressureFunctionSpaceType >
				PressureType;
			typedef PROBLEM_NAMESPACE::Force< VelocityFunctionSpaceType >
				ForceType;
			typedef PROBLEM_NAMESPACE::DirichletData< VelocityFunctionSpaceType >
				DirichletDataType;

		/**
		 *  \brief  constructor
		 *
		 *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
		 **/
		Container( const double viscosity, const DiscreteFunctionWrapperType& funcWrapper, const DirichletDataType& dirichlet )
			: velocity_( funcWrapper.discreteVelocity().space() ),
			  pressure_ ( funcWrapper.discretePressure().space() ),
			  force_( viscosity, funcWrapper.discreteVelocity().space() ),
			  dirichletData_( dirichlet )
		{
		}

		/**
		 *  \brief  destructor
		 *
		 *  doing nothing
		 **/
		~Container()
		{
		}

		/**
		 *  \brief  to get the velocity
		 *
		 *  \return velocity
		 **/
		const VelocityType& velocity() const
		{
			return velocity_;
		}

		/**
		 *  \brief  to get the pressure
		 *
		 *  \return pressure
		 **/
		const PressureType& pressure() const
		{
			return pressure_;
		}

		/**
		 *  \brief  to get the force term
		 *
		 *  \return force
		 **/
		const ForceType& force() const
		{
			return force_;
		}

		/**
		 *  \brief  to get the dirichlet boundary data
		 *
		 *  \return dirichlet boundary data
		 **/
		const DirichletDataType& dirichletData() const
		{
			return dirichletData_;
		}

		/**
		 *  \brief  a simple test of all class' functionalities
		 **/
		void testMe()
		{
			// some logstreams
			Logging::LogStream& infoStream = Logger().Info();
			infoStream << "testing class Problem..." << std::endl;
			//tests
			velocity_.testMe();
			pressure_.testMe();
	//        force_.testMe();
	//        dirichletData_.testMe();
			// happy
			infoStream << "...test passed!" << std::endl;
		}

		private:
			VelocityType velocity_;
			PressureType pressure_;
			ForceType force_;
			const DirichletDataType& dirichletData_;
	};
} // namespace StokesProblems

#endif // DUNE_STOKES_PROBLEMS_HH
