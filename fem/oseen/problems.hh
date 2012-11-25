#ifndef DUNE_OSEEN_PROBLEMS_HH
#define DUNE_OSEEN_PROBLEMS_HH

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
	//! cockburn, see ref
	namespace Generalized {
		#include "problems/generalized.hh"
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
	namespace TimeDisc {
		#include "problems/timedisc.hh"
	}

#ifndef PROBLEM_NAMESPACE
	#define PROBLEM_NAMESPACE StokesProblems::Cockburn
#endif

	/**
	 *  \brief  a collection of some analytical functions describing a stokes problem
	 *
	 *  namely velocity, pressure, force term and dirichlet boundary data
	 *
	 *  \tparam gridDim
	 *          dimension of the grid
	 *
		\note only ever used in \file postprocessor.hh
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

		private:
			VelocityType velocity_;
			PressureType pressure_;
			ForceType force_;
			const DirichletDataType& dirichletData_;
	};
} // namespace OseenProblems

#endif // DUNE_OSEEN_PROBLEMS_HH

/** Copyright (c) 2012, Rene Milk 
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

