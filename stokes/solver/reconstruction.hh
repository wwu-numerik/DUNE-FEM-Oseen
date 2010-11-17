#ifndef DUNE_STOKES_RECONSTRUCTION_HH
#define DUNE_STOKES_RECONSTRUCTION_HH

#include <dune/stuff/functionadapter.hh>

namespace Dune {

template < class DiscreteVelocityFunctionRangeType >
struct ConvectiveTerm : public DiscreteVelocityFunctionRangeType {
	template < class DiscreteSigmaFunctionRangeType >
	ConvectiveTerm( const DiscreteVelocityFunctionRangeType beta,
					DiscreteSigmaFunctionRangeType du)
	{
		for ( size_t d = 0; d< beta.dim(); ++d )  {
			DiscreteVelocityFunctionRangeType j;
			for ( size_t i = 0; i< beta.dim(); ++i )  {
				j[i] = du(d,i);
			}
			(*this)[d] = beta * j;
		}
	}
};

template < class MatrixType, class OperandFunctionType, class ResultFunctionType >
static inline void naiveMatrixMult( const MatrixType& matrix, const OperandFunctionType& operand, ResultFunctionType& result )
{
	const int rows = matrix.rows();
    const int cols = matrix.columns();
    for ( int i = 0; i < rows; ++i ) {
		result[i] = 0;
        for ( int j = 0; j < cols; ++j ) {
			result[i] += matrix.get(i,j) * operand[j];
		}
	}
}

template < class MatrixType, class OperandFunctionType, class ResultFunctionType >
static inline void naiveMatrixMultAdd( const MatrixType& matrix, const OperandFunctionType& operand, ResultFunctionType& result )
{
	const int rows = matrix.rows();
    const int cols = matrix.columns();
    for ( int i = 0; i < rows; ++i ) {
		double row_val = 0;
        for ( int j = 0; j < cols; ++j ) {
			row_val += matrix.get(i,j) * operand[j];
		}
		result[i] += row_val;
	}
}

template <class MatrixObjectType, class DiscretePressureFunctionType, class PressureGradientDiscreteFunctionType>
void getPressureGradient( MatrixObjectType& matrix_object, const DiscretePressureFunctionType& pressure, PressureGradientDiscreteFunctionType& pressure_gradient )
{
	typedef typename DiscretePressureFunctionType::FunctionSpaceType
		SpaceType;
	typedef typename SpaceType::GridPartType
        GridPart;
	typedef typename GridPart::GridType::template Codim< 0 >::Entity
		EntityType;
    typedef typename GridPart::template Codim< 0 >::IteratorType
        EntityIteratorType;
    typedef typename GridPart::IntersectionIteratorType
        IntersectionIteratorType;
    typedef typename IntersectionIteratorType::EntityPointer
        EntityPointer;
	typedef typename MatrixObjectType::LocalMatrixType
		LocalMatrixType;
	typedef typename PressureGradientDiscreteFunctionType::LocalFunctionType
		PressureGradientLocalFunction;
	typedef typename DiscretePressureFunctionType::LocalFunctionType
		PressureLocalFunction;
	const SpaceType& space_ = pressure.space();
	const GridPart& gridPart_ = space_.gridPart();
	Logger().Err().Resume( 9001 );
	Stuff::LocalMatrixPrintFunctor< MatrixObjectType, Logging::LogStream > local_print( matrix_object, Logger().Err(), std::string("LOCAL Z" ) );
	Stuff::LocalFunctionVerbatimPrintFunctor< PressureGradientDiscreteFunctionType, Logging::LogStream > local_print_pressure_grad( pressure_gradient, Logger().Err() );
	Stuff::LocalFunctionVerbatimPrintFunctor< DiscretePressureFunctionType, Logging::LogStream > local_print_pressure( pressure, Logger().Err() );
	EntityIteratorType entityItEndLog = space_.end();
	for (   EntityIteratorType it = space_.begin();
			it != entityItEndLog;
			++it )
	{
		const EntityType& entity = *it;
		LocalMatrixType local_matrix = matrix_object.localMatrix( *it, *it );
		PressureGradientLocalFunction local_pressure_gradient = pressure_gradient.localFunction( * it );
		PressureLocalFunction local_pressure = pressure.localFunction( * it );
//				local_print(*it,*it,0,0);
//				local_print_pressure_grad(*it,*it,0,0);
//				local_print_pressure(*it,*it,0,0);
//				Logger().Err() << std::endl;
//				local_matrix.multiplyAdd( local_pressure, local_pressure_gradient );
		naiveMatrixMult( local_matrix, local_pressure, local_pressure_gradient );


		IntersectionIteratorType intItEnd = space_.gridPart().iend( *it );
		for (   IntersectionIteratorType intIt = space_.gridPart().ibegin( *it );
				intIt != intItEnd;
				++intIt ) {
			const typename IntersectionIteratorType::Intersection& intersection = *intIt;
			if ( intersection.neighbor() && !intersection.boundary() ) {
				const typename IntersectionIteratorType::EntityPointer neighbourPtr = intersection.outside();
                const EntityType& neighbour = *neighbourPtr;
				LocalMatrixType local_matrix_neighbour = matrix_object.localMatrix( entity, neighbour );
				PressureLocalFunction local_pressure_neighbour = pressure.localFunction( neighbour );
				naiveMatrixMultAdd( local_matrix_neighbour, local_pressure_neighbour, local_pressure_gradient );
			}
		}
	}
}


template < class DataContainerType, class DiscreteModelType >
struct BruteForceReconstruction {
	template < class DiscreteVelocityFunctionType, class GradientFunctionType >
	static void getConvection( const DiscreteVelocityFunctionType& beta, const GradientFunctionType& sigma, DiscreteVelocityFunctionType& convection)
	{
		typedef typename DiscreteVelocityFunctionType::FunctionSpaceType
			DiscreteVelocityFunctionSpaceType;
		typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPart;
		typedef typename GridPart::GridType::template Codim< 0 >::Entity
			EntityType;
        typedef typename GridPart::template Codim< 0 >::IteratorType
            EntityIteratorType;
		const DiscreteVelocityFunctionSpaceType& discrete_velocity_space = beta.space();
		convection.clear();
		EntityIteratorType entityItEnd = discrete_velocity_space.end();
		for (   EntityIteratorType entityIt = discrete_velocity_space.begin();
				entityIt != entityItEnd;
				++entityIt)
		{
			const EntityType& entity = *entityIt;
			typedef typename EntityType::Geometry
				EntityGeometryType;
			const EntityGeometryType& geo = entity.geometry();
			typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
				VelocityBaseFunctionSetType;

			typename DiscreteVelocityFunctionType::LocalFunctionType beta_local = beta.localFunction( *entityIt );
			typename DiscreteVelocityFunctionType::LocalFunctionType convection_local = convection.localFunction( *entityIt );
			typedef typename DiscreteModelType::VolumeQuadratureType
	            VolumeQuadratureType;
			const VolumeQuadratureType quad( entity, ( 2 * discrete_velocity_space.order() ) + 1 );
			const VelocityBaseFunctionSetType velocityBaseFunctionSetElement = discrete_velocity_space.baseFunctionSet( entity );
			const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.numBaseFunctions();
			const int quadNop = quad.nop();

			for(int qP = 0; qP < quadNop ; ++qP)
			{
				const typename DiscreteVelocityFunctionSpaceType::DomainType xLocal = quad.point(qP);

				const double intel = quad.weight(qP);
//					const double intel =
//						quad.weight(qP)* geo.integrationElement( xLocal ); // wrong!

				typename DiscreteVelocityFunctionSpaceType::DomainType
					xWorld = geo.global( xLocal );

				// evaluate function
				typename GradientFunctionType::RangeType
					sigma_eval;
				sigma.evaluate( xWorld, sigma_eval );

				typename DiscreteVelocityFunctionType::RangeType
					beta_eval;
				beta_local.evaluate( quad[qP], beta_eval );

				ConvectiveTerm< typename DiscreteVelocityFunctionType::RangeType >
				        c(beta_eval,sigma_eval);

				// do projection
				for(int i=0; i<numVelocityBaseFunctionsElement; ++i)
				{
					typename DiscreteVelocityFunctionType::RangeType phi (0.0);
					convection_local.baseFunctionSet().evaluate(i, quad[qP], phi);
					convection_local[i] +=  intel * ( c * phi );
				}
			}
		}
	}
	template <  class FunctionWrapperType,
				class XmatrixObjectType,
				class MInversMatrixObjectType,
				class YmatrixObjectType,
				class OmatrixObjectType,
				class EmatrixObjectType,
				class RmatrixObjectType,
				class ZmatrixObjectType,
				class WmatrixObjectType,
				class DiscreteSigmaFunctionType,
				class DiscreteVelocityFunctionType,
				class DiscretePressureFunctionType  >
	static void reconstruct( DataContainerType& rhs_datacontainer,
				const FunctionWrapperType& solution,
				const DiscreteVelocityFunctionType& beta,
				const XmatrixObjectType& /*Xmatrix*/,
				const MInversMatrixObjectType& /*MInversMatrix*/,
				const YmatrixObjectType& /*Ymatrix*/,
				const OmatrixObjectType& /*Omatrix*/,
				const EmatrixObjectType& /*Ematrix*/,
				const RmatrixObjectType& /*Rmatrix*/,
				const ZmatrixObjectType& /*Zmatrix*/,
				const WmatrixObjectType& /*Wmatrix*/,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& /*H2rhs*/,
				const DiscretePressureFunctionType& /*H3rhs*/ )
	{
		DiscreteVelocityFunctionType velocity_tmp1( "velocity_tmp1", solution.discreteVelocity().space() );
		DiscreteSigmaFunctionType sigma_tmp( "sigma_dummy", H1rhs.space() );

		Stuff::GradientAdapterFunction< DiscretePressureFunctionType, DiscreteVelocityFunctionType,Stuff::ProductFunctorMatrixVector >
				pressure_grad ( solution.discretePressure(), velocity_tmp1 );
		rhs_datacontainer.pressure_gradient.assign( pressure_grad );

		Stuff::GradientAdapterFunction< DiscreteVelocityFunctionType, DiscreteSigmaFunctionType,Stuff::ProductFunctorMatrices >
				grad ( solution.discreteVelocity(), sigma_tmp );
		rhs_datacontainer.velocity_gradient .assign( grad );

		Stuff::LaplaceAdapterFunction< DiscreteVelocityFunctionType, DiscreteSigmaFunctionType,Stuff::ProductFunctorMatrices >
		        laplace( solution.discreteVelocity(), sigma_tmp );
		rhs_datacontainer.velocity_laplace.assign( laplace );

		getConvection( beta, rhs_datacontainer.velocity_gradient, rhs_datacontainer.convection );
	}

};

template < class DataContainerType, class DiscreteModelType >
struct SmartReconstruction {
	template <  class FunctionWrapperType,
				class XmatrixObjectType,
				class MInversMatrixObjectType,
				class YmatrixObjectType,
				class OmatrixObjectType,
				class EmatrixObjectType,
				class RmatrixObjectType,
				class ZmatrixObjectType,
				class WmatrixObjectType,
				class DiscreteSigmaFunctionType,
				class DiscreteVelocityFunctionType,
				class DiscretePressureFunctionType  >
	static void reconstruct( DataContainerType& rhs_datacontainer,
				const FunctionWrapperType& solution,
				const DiscreteVelocityFunctionType& beta,
				const XmatrixObjectType& Xmatrix,
				const MInversMatrixObjectType& MInversMatrix,
				const YmatrixObjectType& Ymatrix,
				const OmatrixObjectType& Omatrix,
				const EmatrixObjectType& Ematrix,
				const RmatrixObjectType& Rmatrix,
				const ZmatrixObjectType& Zmatrix,
				const WmatrixObjectType& Wmatrix,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& H2rhs,
				const DiscretePressureFunctionType& H3rhs )
	{
		CompileTimeChecker<false> dysfunctionalCode;
		//				Zmatrix.apply( dest.discretePressure(), rhs_datacontainer->pressure_gradient );
		//				rhs_datacontainer->pressure_gradient *= Parameters().getParam("pressure_gradient_scale", 1);
		//				getPressureGradient( Zmatrix,  dest.discretePressure(),  rhs_datacontainer->pressure_gradient);

						// \sigma = M^{-1} ( H_1 - Wu )
		//				const double m_inv_scale = MInversMatrix.matrix()(0,0);
		//				rhs_datacontainer->velocity_gradient.assign( H1rhs );

		//				Wmatrix.apply( dest.discreteVelocity(), sigma_tmp );
		//				rhs_datacontainer->velocity_gradient -= sigma_tmp;
		////				if ( viscosity != 0.0f )
		////					rhs_datacontainer->velocity_gradient /= viscosity;//since mu is assmenled into both W and H1
		//				rhs_datacontainer->velocity_gradient *= m_inv_scale;

		//				Stuff::printFunctionMinMax( std::cout, H1rhs );


		//				Xmatrix.apply( rhs_datacontainer->velocity_gradient, velocity_tmp1 );
		//				Ymatrix.apply( dest.discreteVelocity(), rhs_datacontainer->velocity_laplace );
		//				rhs_datacontainer->velocity_laplace += velocity_tmp1;
		//				velocity_tmp1.assign( dest.discreteVelocity() );
		//				velocity_tmp1 *= alpha;
		//				rhs_datacontainer->velocity_laplace -= velocity_tmp1;
		////				Stuff::printFunctionMinMax( std::cout, rhs_datacontainer->velocity_laplace );
		//				const double laplace_scale = Parameters().getParam("laplace_scale", -1/viscosity);
		//				rhs_datacontainer->velocity_laplace *= laplace_scale;
		////				Stuff::printFunctionMinMax( std::cout, rhs_datacontainer->velocity_laplace );
		//				Logger().Dbg().Resume();
		//				Logger().Dbg() << boost::format( "laplace_scale: %f\n") % laplace_scale;

		//				rhs_datacontainer->convection.clear();
		//				Omatrix.apply( dest.discreteVelocity(), rhs_datacontainer->convection );
		//				rhs_datacontainer->convection += H2_O_rhs;
		//				if ( sigma_exact )
		//					getConvection( beta_, *sigma_exact, rhs_datacontainer->convection );
		//				else

		//				dest.discreteVelocity() += rhs_datacontainer->convection;
		//				Stuff::printFunctionMinMax( std::cout, rhs_datacontainer->convection );

		//				rhs_datacontainer->scale( 1 / std::sqrt(2) );

	}
};
}

#endif // DUNE_STOKES_RECONSTRUCTION_HH
