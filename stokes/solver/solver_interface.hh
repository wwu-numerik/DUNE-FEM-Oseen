#ifndef DUNE_STOKES_SOLVER_INTERFACE_HH
#define DUNE_STOKES_SOLVER_INTERFACE_HH

// OEMBICGSQOp will NOT compile
#ifndef INNER_CG_SOLVERTYPE
	#define INNER_CG_SOLVERTYPE OEMCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
	#define OUTER_CG_SOLVERTYPE OEMCGOp
#endif


#ifdef USE_BFG_CG_SCHEME
	#include <utility>
	//< iteration no , < absLimit, residuum > >
	typedef std::pair<int,std::pair<double,double> >
		IterationInfo;
	#include <dune/stokes/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE DuneStokes
#else
	#include <dune/fem/solver/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE Dune
#endif

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stokes/solver/cghelper.hh>

#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>

#include <cmath>
#include <boost/utility.hpp>



namespace Dune {
//! utility struct used to expose runtime statistics
struct SaddlepointInverseOperatorInfo {
	double iterations_inner_avg;
	int iterations_inner_min;
	int iterations_inner_max;
	int iterations_outer_total;
	double max_inner_accuracy;
};

template < class MatrixObjectType >
class MatrixWrapper : boost::noncopyable {
	private:
		typedef typename MatrixObjectType::MatrixType
			MatrixType;
	public:
		MatrixWrapper( const MatrixObjectType& matrix_object )
			:matrix_object_( matrix_object ),
			cumulative_scale_factor_( 1.0 )
		{}

		~MatrixWrapper()
		{
			assert( cumulative_scale_factor_ != 0 );
			matrix_object_.matrix().scale( 1.0/cumulative_scale_factor_ );
		}

		template <class DiscFType, class DiscFuncType>
		void apply(const DiscFType &f, DiscFuncType &ret) const
		{
			matrix_object_.apply( f, ret );
		}
		//! return diagonal of (this * A * B)
		template <class DiscrecteFunctionType>
		void getDiag(const MatrixType& A, const MatrixType& B, DiscrecteFunctionType& rhs) const
		{
			matrix_object_.matrix().getDiag( A, B, rhs );
		}

		//! return diagonal of (this * A)
		template <class DiscrecteFunctionType>
		void getDiag(const MatrixType& A, DiscrecteFunctionType& rhs) const
		{
			matrix_object_.matrix().getDiag( A, rhs );
		}

		//! same as apply A * x = ret, used by OEM-Solvers
		template <class VECtype>
		void multOEM(const VECtype *x, VECtype * ret) const
		{
			matrix_object_.matrix().multOEM( x, ret );
		}

		//! calculates ret += A * x
		template <class VECtype>
		void multOEMAdd(const VECtype *x, VECtype * ret) const
		{
			matrix_object_.matrix().multOEMAdd( x, ret );
		}

		void scale( const double factor )
		{
			cumulative_scale_factor_ *= factor;
			matrix_object_.matrix().scale( factor );
		}

	private:
		const MatrixObjectType& matrix_object_;
		double cumulative_scale_factor_;
};

template<class SolverType>
struct SolverCaller {
	template <  class DomainType,
				class RangeType,
				class XmatrixObjectType,
				class MInversMatrixObjectType,
				class YmatrixObjectType,
				class EmatrixObjectType,
				class RmatrixObjectType,
				class ZmatrixObjectType,
				class WmatrixObjectType,
				class DiscreteSigmaFunctionType,
				class DiscreteVelocityFunctionType,
				class DiscretePressureFunctionType  >
	static SaddlepointInverseOperatorInfo solve( const DomainType& arg,
				RangeType& dest,
				const XmatrixObjectType& Xmatrix,
				const MInversMatrixObjectType& MInversMatrix,
				const YmatrixObjectType& Ymatrix,
				const YmatrixObjectType& Omatrix,
				const EmatrixObjectType& Ematrix,
				const RmatrixObjectType& Rmatrix,
				const ZmatrixObjectType& Zmatrix,
				const WmatrixObjectType& Wmatrix,
				const DiscreteSigmaFunctionType& H1rhs,
				const DiscreteVelocityFunctionType& H2rhs,
				const DiscretePressureFunctionType& H3rhs )
	{
		MatrixWrapper<XmatrixObjectType> X(Xmatrix);
		SolverType().solve(	arg, dest,
							X, MInversMatrix, Ymatrix,
							Omatrix, Ematrix, Rmatrix,
							Zmatrix, Wmatrix,
							H1rhs, H2rhs, H3rhs );
	}

};

} //end namespace Dune

#endif // DUNE_STOKES_SOLVER_INTERFACE_HH
