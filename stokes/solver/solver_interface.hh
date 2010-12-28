#ifndef DUNE_STOKES_SOLVER_INTERFACE_HH
#define DUNE_STOKES_SOLVER_INTERFACE_HH

// OEMBICGSQOp will NOT compile
#ifndef INNER_CG_SOLVERTYPE
	#define INNER_CG_SOLVERTYPE OEMCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
	#define OUTER_CG_SOLVERTYPE OEMCGOp
#endif


#if defined(USE_BFG_CG_SCHEME) || defined(FORCE_CUSTOM_SOLVER)
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

	SaddlepointInverseOperatorInfo()
		:iterations_inner_avg(-1.0f),iterations_inner_min(-1),
		iterations_inner_max(-1),iterations_outer_total(-1),
		max_inner_accuracy(-1.0f)
	{}
};

template < class MatrixObjectType >
class MatrixWrapper : boost::noncopyable {
	public:
		typedef typename MatrixObjectType::MatrixType
			MatrixType;
		typedef typename MatrixObjectType::MatrixType
			RealMatrixType;

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
		template <class DiscrecteFunctionType, class OtherMatrixType_A, class OtherMatrixType_B>
		void getDiag(const OtherMatrixType_A& A, const OtherMatrixType_B& B, DiscrecteFunctionType& rhs) const
		{
			matrix_object_.matrix().getDiag( A.matrix(), B.matrix(), rhs );
		}

		//! return diagonal of (this * A)
		template <class DiscrecteFunctionType, class OtherMatrixType_A>
		void getDiag(const OtherMatrixType_A& A, DiscrecteFunctionType& rhs) const
		{
			matrix_object_.matrix().getDiag( A, rhs );
		}

		template <class DiscrecteFunctionType>
		void addDiag(DiscrecteFunctionType& rhs) const
		{
			matrix_object_.matrix().addDiag( rhs );
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

		double operator ()(const size_t i, const size_t j ) const
		{
			return matrix_object_.matrix()(i,j);
		}

		size_t rows() const
		{
			return matrix_object_.matrix().rows();
		}

		size_t cols() const
		{
			return matrix_object_.matrix().cols();
		}

		void scale( const double factor )
		{
			cumulative_scale_factor_ *= factor;
			matrix_object_.matrix().scale( factor );
		}

		const MatrixType& matrix() const
		{
			return matrix_object_.matrix();
		}

	private:
		const MatrixObjectType& matrix_object_;
		double cumulative_scale_factor_;
};


} //end namespace Dune

#endif // DUNE_STOKES_SOLVER_INTERFACE_HH
