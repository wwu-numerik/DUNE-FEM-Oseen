#ifndef DUNE_OSEEN_RECONSTRUCTION_HH
#define DUNE_OSEEN_RECONSTRUCTION_HH

#include <dune/stuff/fem/functionadapter.hh>
#include <dune/fem/oseen/datacontainer.hh>

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
    const auto& space_ = pressure.space();
    DSC_LOG_ERROR.resume( 9001 );

    for (const auto& entity : space_ )
	{
        auto local_matrix = matrix_object.localMatrix(entity, entity);
        auto local_pressure_gradient = pressure_gradient.localFunction(entity);
        auto local_pressure = pressure.localFunction(entity);
		naiveMatrixMult( local_matrix, local_pressure, local_pressure_gradient );


        auto  intItEnd = space_.gridPart().iend(entity);
        for (   auto  intIt = space_.gridPart().ibegin(entity);
				intIt != intItEnd;
				++intIt ) {
            const auto & intersection = *intIt;
			if ( intersection.neighbor() && !intersection.boundary() ) {
                const auto neighbourPtr = intersection.outside();
                const auto & neighbour = *neighbourPtr;
                auto  local_matrix_neighbour = matrix_object.localMatrix( entity, neighbour );
                auto local_pressure_neighbour = pressure.localFunction( neighbour );
				naiveMatrixMultAdd( local_matrix_neighbour, local_pressure_neighbour, local_pressure_gradient );
			}
		}
	}
}


template < class DiscreteModelType >
struct BruteForceReconstruction {
	template < class DiscreteVelocityFunctionType, class GradientFunctionType >
	static void getConvection( const DiscreteVelocityFunctionType& beta, const GradientFunctionType& sigma, DiscreteVelocityFunctionType& convection)
	{
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
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
            const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.size();
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
                class DiscretePressureFunctionType,
                class AnyTraits  >
    static void reconstruct( Dune::Oseen::RhsDatacontainer<AnyTraits>& rhs_datacontainer,
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

        DSFe::GradientAdapterFunction< DiscretePressureFunctionType, DiscreteVelocityFunctionType,DSFe::ProductFunctorMatrixVector >
				pressure_grad ( solution.discretePressure(), velocity_tmp1 );
		rhs_datacontainer.pressure_gradient.assign( pressure_grad );

        DSFe::GradientAdapterFunction< DiscreteVelocityFunctionType, DiscreteSigmaFunctionType,DSFe::ProductFunctorMatrices >
				grad ( solution.discreteVelocity(), sigma_tmp );
		rhs_datacontainer.velocity_gradient .assign( grad );

        DSFe::LaplaceAdapterFunction< DiscreteVelocityFunctionType, DiscreteSigmaFunctionType,DSFe::ProductFunctorMatrices >
		        laplace( solution.discreteVelocity(), sigma_tmp );
		rhs_datacontainer.velocity_laplace.assign( laplace );

		getConvection( beta, rhs_datacontainer.velocity_gradient, rhs_datacontainer.convection );
	}

};

template < class DiscreteModelType >
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
                class DiscretePressureFunctionType,
                class AnyTraits>
    static void reconstruct( Dune::Oseen::RhsDatacontainer<AnyTraits>& rhs_datacontainer,
				const FunctionWrapperType& solution,
				const DiscreteVelocityFunctionType& /*beta*/,
				const XmatrixObjectType& Xmatrix,
				const MInversMatrixObjectType& MInversMatrix,
				const YmatrixObjectType& Ymatrix,
				const OmatrixObjectType& Omatrix,
				const EmatrixObjectType& /*Ematrix*/,
				const RmatrixObjectType& /*Rmatrix*/,
				const ZmatrixObjectType& Zmatrix,
				const WmatrixObjectType& Wmatrix,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& H2rhs,
				const DiscretePressureFunctionType& /*H3rhs*/ )
	{
//		CompileTimeChecker<false> dysfunctionalCode;
						Zmatrix.apply( solution.discretePressure(), rhs_datacontainer.pressure_gradient );
//						rhs_datacontainer.pressure_gradient *= DSC_CONFIG_GET("pressure_gradient_scale", 1);
//						getPressureGradient( Zmatrix,  solution.discretePressure(),  rhs_datacontainer.pressure_gradient);

						 //\sigma = M^{-1} ( H_1 - Wu )
						const double m_inv_scale = MInversMatrix(0,0);
						rhs_datacontainer.velocity_gradient.assign( H1rhs );

						DiscreteSigmaFunctionType sigma_tmp ("s_tmp", H1rhs.space() );
						DiscreteVelocityFunctionType velocity_tmp1 ("s_tmp", H2rhs.space() );

						Wmatrix.apply( solution.discreteVelocity(), sigma_tmp );
						rhs_datacontainer.velocity_gradient -= sigma_tmp;
		//				if ( viscosity != 0.0f )
		//					rhs_datacontainer.velocity_gradient /= viscosity;//since mu is assmenled into both W and H1
						rhs_datacontainer.velocity_gradient *= m_inv_scale;

//						DSC::printFunctionMinMax( std::cout, H1rhs );


						Xmatrix.apply( rhs_datacontainer.velocity_gradient, velocity_tmp1 );
						Ymatrix.apply( solution.discreteVelocity(), rhs_datacontainer.velocity_laplace );
						rhs_datacontainer.velocity_laplace += velocity_tmp1;
						velocity_tmp1.assign( solution.discreteVelocity() );
//						velocity_tmp1 *= alpha;
						rhs_datacontainer.velocity_laplace -= velocity_tmp1;
		//				DSC::printFunctionMinMax( std::cout, rhs_datacontainer.velocity_laplace );
//						const double laplace_scale = DSC_CONFIG_GET("laplace_scale", -1/viscosity);
//						rhs_datacontainer.velocity_laplace *= laplace_scale;
		//				DSC::printFunctionMinMax( std::cout, rhs_datacontainer.velocity_laplace );
//						DSC_LOG_DEBUG.resume();
//						DSC_LOG_DEBUG << boost::format( "laplace_scale: %f\n") % laplace_scale;

						rhs_datacontainer.convection.clear();
						Omatrix.apply( solution.discreteVelocity(), rhs_datacontainer.convection );
//						rhs_datacontainer.convection += H2_O_rhs;
//						if ( sigma_exact )
//							getConvection( beta_, *sigma_exact, rhs_datacontainer.convection );
//						else

//						dest.discreteVelocity() += rhs_datacontainer.convection;
//						DSC::printFunctionMinMax( std::cout, rhs_datacontainer.convection );

//						rhs_datacontainer.scale( 1 / std::sqrt(2) );

	}
};
}

#endif // DUNE_OSEEN_RECONSTRUCTION_HH

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

