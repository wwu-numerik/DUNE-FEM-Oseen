#ifndef SCHURKOMPLEMENT_ADAPTER_HH
#define SCHURKOMPLEMENT_ADAPTER_HH

#include <dune/fem/operator/2order/dgmatrixsetup.hh>

namespace Dune {

template <class MatrixImp>
class SchurkomplementOperatorAdapter
//        : public AssembledLinearOperator< MatrixImp,
//    typename MatrixImp :: RowDiscreteFunctionType :: DofStorageType,
//    typename MatrixImp :: ColDiscreteFunctionType :: DofStorageType>
{
    typedef SchurkomplementOperatorAdapter< MatrixImp >
        ThisType ;
//    typedef AssembledLinearOperator< MatrixImp,
//            typename MatrixImp :: RowDiscreteFunctionType :: DofStorageType,
//            typename MatrixImp :: ColDiscreteFunctionType :: DofStorageType>
//        BaseType;

  public:
//    enum { category=SolverCategory::sequential };
    typedef MatrixImp MatrixType;
    typedef DSC::PrecondionWrapperDummy<MatrixType> PreconditionAdapterType;

    typedef typename MatrixType :: RowDiscreteFunctionType RowDiscreteFunctionType;
    typedef typename MatrixType :: ColDiscreteFunctionType ColumnDiscreteFunctionType;

    typedef typename RowDiscreteFunctionType :: DiscreteFunctionSpaceType RowSpaceType;

    typedef typename ColumnDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;
    typedef ParallelScalarProduct<ColumnDiscreteFunctionType> ParallelScalarProductType;

    typedef typename RowDiscreteFunctionType:: DofStorageType X;
    typedef typename ColumnDiscreteFunctionType:: DofStorageType Y;

    //! export types
    typedef MatrixType  matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

  protected:
    MatrixType& matrix_;
    const RowSpaceType& rowSpace_;
    const ColSpaceType& colSpace_;

    ParallelScalarProductType scp_;

    PreconditionAdapterType preconditioner_;
    mutable double averageCommTime_;

  public:
    //! constructor: just store a reference to a matrix
    SchurkomplementOperatorAdapter (const SchurkomplementOperatorAdapter& org)
      : matrix_(org.matrix_)
      , rowSpace_(org.rowSpace_)
      , colSpace_(org.colSpace_)
      , scp_(colSpace_)
      , preconditioner_(org.preconditioner_)
      , averageCommTime_( org.averageCommTime_ )
    {}

    //! constructor: just store a reference to a matrix
    SchurkomplementOperatorAdapter (MatrixType& A,
                 const RowSpaceType& rowSpace,
                 const ColSpaceType& colSpace)
      : matrix_(A)
      , rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , scp_(colSpace_)
      , preconditioner_( PreconditionAdapterType() )
      , averageCommTime_( 0.0 )
    {}

    //! return communication time
    double averageCommTime() const
    {
      return averageCommTime_ ;
    }

    //! return reference to preconditioner
    PreconditionAdapterType& preconditionAdapter() { return preconditioner_; }
    const PreconditionAdapterType& preconditionAdapter() const { return preconditioner_; }

    //! return reference to preconditioner
    ParallelScalarProductType& scp() { return scp_; }
    const ParallelScalarProductType& scp() const { return scp_; }

    //! apply operator to x:  \f$ y = A(x) \f$

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
      // exchange data first
      communicate( x );

      // apply vector to matrix
      matrix_.multOEM(x,y);

      // delete non-interior
      scp_.deleteNonInterior( y );
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
      // exchange data first
      communicate( x );

      // apply matrix
      matrix_.usmv(alpha,x,y);

      // delete non-interior
      scp_.deleteNonInterior( y );
    }

    template < class T, class O >
    double residuum(const T& rhs, O& x) const
    {
      // exchange data
      communicate( x );
      T tmp( rhs );
      apply( x, tmp );
      tmp -= rhs;
      double res = tmp.two_norm();

      res = rowSpace_.grid().comm().sum( res );
      // return global sum of residuum
      return std::sqrt( res );
    }

    //! get matrix via *
    virtual const MatrixType& getmat () const
    {
      return matrix_;
    }

  protected:
    void communicate(const X& x) const
    {

      if( rowSpace_.grid().comm().size() <= 1 ) return ;

      Timer commTime;

      // create temporary discretet function object
      RowDiscreteFunctionType tmp ("DGParallelMatrixAdapter::communicate",
                   rowSpace_, x );

      // exchange data by copying
      rowSpace_.communicate( tmp );

      // accumulate communication time
      averageCommTime_ += commTime.elapsed();
    }
};

} // namespace Dune {

#endif // SCHURKOMPLEMENT_ADAPTER_HH

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

