#ifndef DUNE_OSEEN_SOLVER_INTERFACE_HH
#define DUNE_OSEEN_SOLVER_INTERFACE_HH

#include "solver_defines.hh"

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/oseen/solver/cghelper.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/common/logging.hh>

#include <cmath>
#include <boost/utility.hpp>

#if STOKES_USE_ISTL
#   define STOKES_MATRIX_ACCESS assert(matrix_pointer_);*(matrix_pointer_)
#else
#   define STOKES_MATRIX_ACCESS assert(matrix_pointer_);matrix_pointer_->matrix()
#endif

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
		}

		~MatrixWrapper()
		{
			assert( cumulative_scale_factor_ != 0 );
			matrix_pointer_->matrix().scale( 1.0/cumulative_scale_factor_ );
		}

		template <class DiscFType, class DiscFuncType>
		void apply(const DiscFType &f, DiscFuncType &ret) const
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
            STOKES_MATRIX_ACCESS.multOEM( x, ret );
		}
		template <class VECtype, class VECtypeR>
		void multOEM_t(const VECtype* x, VECtypeR* ret) const
		{
            STOKES_MATRIX_ACCESS.multOEM_t( x, ret );
		}

		//! calculates ret += A * x
		template <class VECtype, class VECtypeR>
		void multOEMAdd(const VECtype* x, VECtypeR* ret) const
		{
            STOKES_MATRIX_ACCESS.multOEMAdd( x, ret );
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
                return matrix_pointer_->operator()(i,j);
            #else
                return  matrix_pointer_->matrix()(i,j);
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

#endif // DUNE_OSEEN_SOLVER_INTERFACE_HH

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

