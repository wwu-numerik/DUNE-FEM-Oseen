#ifndef DUNE_OSEEN_DEFAULTTRAITS_HH
#define DUNE_OSEEN_DEFAULTTRAITS_HH

#include <dune/oseen/stab_coeff.hh>

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
		typedef typename DiscreteModelType::DiscreteOseenFunctionSpaceWrapperType
			DiscreteOseenFunctionSpaceWrapperType;

		//! discrete function wrapper type
		typedef typename DiscreteModelType::DiscreteOseenFunctionWrapperType
			DiscreteOseenFunctionWrapperType;

		//! discrete function type for the velocity
		typedef typename DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType
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
		typedef typename DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType
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
}//end namespace Dune
#endif // DEFAULTTRAITS_HH

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

