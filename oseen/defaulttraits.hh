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
}//end namespace Dune
#endif // DEFAULTTRAITS_HH
