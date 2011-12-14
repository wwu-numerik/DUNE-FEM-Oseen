#ifndef DUNE_STOKES_SOLVER_INTERFACE_HH
#define DUNE_STOKES_SOLVER_INTERFACE_HH

#include "solver_defines.hh"

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stokes/solver/cghelper.hh>

#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/matrix.hh>
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

template < class MatrixPointerType >
class MatrixWrapper : boost::noncopyable {
    public:
		typedef typename MatrixPointerType::element_type::MatrixType
			MatrixType;
		typedef typename MatrixPointerType::element_type::MatrixType
			RealMatrixType;
        typedef typename MatrixPointerType::element_type
            WrappedMatrixObjectType;

        MatrixWrapper( const MatrixPointerType& matrix_object, std::string name )
			:matrix_pointer_( matrix_object ),
			cumulative_scale_factor_( 1.0 )
		{
            #if ! STOKES_USE_ISTL
			Stuff::Matrix::printMemUsage( matrix_pointer_->matrix(), Logger().Dbg(), name );
		    #endif
		}

		~MatrixWrapper()
		{
			assert( cumulative_scale_factor_ != 0 );
			matrix_pointer_->matrix().scale( 1.0/cumulative_scale_factor_ );
		}

//		template <class DiscFType, class DiscFuncType>
        void apply(const typename WrappedMatrixObjectType::RowDiscreteFunctionType &f,
                   typename WrappedMatrixObjectType::ColumnDiscreteFunctionType &ret) const
		{
            matrix_pointer_->apply( f, ret );
		}
		template <class ArgBlockVectorType, class DestBlockVectorType, class ArgDofImp, class DestDofImp>
		void apply(const Dune::StraightenBlockVector<ArgBlockVectorType,ArgDofImp> &f,
				 Dune::StraightenBlockVector<DestBlockVectorType,DestDofImp> &ret) const
		{
		    matrix_pointer_->multOEM( f, ret );
		}

        template <class ArgBlockType, class DestBlockType, class ArgAllocatorType, class DestAllocatorType>
        void apply(const Dune::BlockVector<ArgBlockType,ArgAllocatorType> &f,
                 Dune::BlockVector<DestBlockType,DestAllocatorType> &ret) const
        {
            matrix_pointer_->multOEM( f, ret );
        }
		//! return diagonal of (this * A * B)
		template <class DiscrecteFunctionType, class OtherMatrixType_A, class OtherMatrixType_B>
		void getDiag(const OtherMatrixType_A& A, const OtherMatrixType_B& B, DiscrecteFunctionType& rhs) const
		{
            #if STOKES_USE_ISTL
			assert( false );
		    #else
			matrix_pointer_->matrix().getDiag( A.matrix(), B.matrix(), rhs );
		    #endif
		}

		//! return diagonal of (this * A)
		template <class DiscrecteFunctionType, class OtherMatrixType_A>
		void getDiag(const OtherMatrixType_A& A, DiscrecteFunctionType& rhs) const
		{
            #if STOKES_USE_ISTL
			assert( false );
		    #else
			matrix_pointer_->matrix().getDiag( A, rhs );
		    #endif
		}

		template <class DiscrecteFunctionType>
		void addDiag(DiscrecteFunctionType& rhs) const
		{
            #if STOKES_USE_ISTL
			assert( false );
		    #else
			matrix_pointer_->matrix().addDiag( rhs );
		    #endif
		}

//		template < class T >
//		void applyAdd( const T& x, T& ret ) const
        template <class ArgBlockType, class DestBlockType, class ArgDType, class DestDType>
        void applyAdd(const Dune::BlockVector<ArgBlockType, ArgDType>& x,
                 Dune::BlockVector<DestBlockType, DestDType>& ret ) const
		{
            matrix_pointer_->applyAdd( x, ret );
		}

		//! same as apply A * x = ret, used by OEM-Solvers
		template <class VECtype, class VECtypeR >
		void multOEM(const VECtype* x, VECtypeR* ret) const
		{
            #if STOKES_USE_ISTL
			//matrix wrapped by ISTLMatrixObject isn't interface compatible with Sp matrix
			matrix_pointer_->multOEM( x, ret );
		    #else
			matrix_pointer_->matrix().multOEM( x, ret );
		    #endif
		}
		template <class VECtype, class VECtypeR>
		void multOEM_t(const VECtype* x, VECtypeR* ret) const
		{
            #if STOKES_USE_ISTL
			matrix_pointer_->multOEM_t( x, ret );
		    #else
			matrix_pointer_->matrix().multOEM_t( x, ret );
		    #endif
		}

		//! calculates ret += A * x
		template <class VECtype, class VECtypeR>
		void multOEMAdd(const VECtype* x, VECtypeR* ret) const
		{
            #if STOKES_USE_ISTL
			//VecType == StraightenBlockVector<>*
			matrix_pointer_->multOEMAdd( x, ret );
		    #else
			matrix_pointer_->matrix().multOEMAdd( x, ret );
		    #endif
		}

        template <class ArgDofStorageType, class DestDofStorageType, class ArgRangeFieldType, class DestRangeFieldType>
        void multOEMAdd(const Dune::StraightenBlockVector<ArgDofStorageType,ArgRangeFieldType> &x,
                 Dune::StraightenBlockVector<DestDofStorageType,DestRangeFieldType> &ret) const
        {
            matrix_pointer_->multOEMAdd( x, ret );
        }

        template <class ArgDofStorageType, class DestDofStorageType>
        void multOEMAdd(const Dune::BlockVector<ArgDofStorageType> &x,
                 Dune::BlockVector<DestDofStorageType> &ret) const
        {
            matrix_pointer_->multOEMAdd( x, ret );
        }

        template <class ArgDofStorageType, class DestDofStorageType>
        void multOEM(const Dune::BlockVector<ArgDofStorageType> &x,
                 Dune::BlockVector<DestDofStorageType> &ret) const
        {
            matrix_pointer_->multOEM( x, ret );
        }

		double operator ()(const size_t i, const size_t j ) const
		{
            #if STOKES_USE_ISTL
			return matrix_pointer_->matrix()(i,j);
		    #else
			return matrix_pointer_->matrix()(i,j);
		    #endif
		}

		size_t rows() const
		{
		    return matrix_pointer_->matrix().rows();
		}

		size_t cols() const
		{
		    return matrix_pointer_->matrix().cols();
		}

		void scale( const double factor )
		{
			cumulative_scale_factor_ *= factor;
			matrix_pointer_->matrix().scale( factor );
		}

		const MatrixType& matrix() const
		{
			return matrix_pointer_->matrix();
		}

	private:
		const MatrixPointerType& matrix_pointer_;
		double cumulative_scale_factor_;
};


} //end namespace Dune

#endif // DUNE_STOKES_SOLVER_INTERFACE_HH
