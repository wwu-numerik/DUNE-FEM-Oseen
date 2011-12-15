/**
 *  \file   stuff.hh
 *  \brief  contains some stuff
 **/

#ifndef DUNE_STOKES_SOLVER_CGHELPER_HH
#define DUNE_STOKES_SOLVER_CGHELPER_HH

#include "solver_defines.hh"

#include <dune/stokes/solver/new_bicgstab.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/matrix.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>

namespace Dune {


/** \brief Operator to wrap Matrix-vector multiplication in inner CG algorithms
	\see InnerCGSolverWrapper
	\see SchurkomplementOperator
    multOEM method evaluates matrix vector multiplication\n
	\f$ A := -1 \cdot X  M^{-1} W x  + Yx \f$\n
**/
template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteSigmaFunctionType,
			class DiscreteVelocityFunctionType>
class MatrixA_Operator : public SOLVER_INTERFACE_NAMESPACE::PreconditionInterface
{
	public:
	typedef MatrixA_Operator<   WMatType,
					MMatType,
					XMatType,
					YMatType,
					DiscreteSigmaFunctionType,
					DiscreteVelocityFunctionType>
				ThisType;

	// if shit goes south wrt precond working check if this doesn't need to be OEmSolver instead of SOLVER_INTERFACE_NAMESPACE
	friend class Conversion<ThisType,SOLVER_INTERFACE_NAMESPACE::PreconditionInterface>;
    typedef IdentityMatrixObject<typename YMatType::WrappedMatrixObjectType>
		PreconditionMatrixBaseType;

#if STOKES_USE_ISTL
	typedef SchurkomplementOperatorAdapter< ThisType/*,
						typename DiscretePressureFunctionType::DiscreteFunctionSpaceType,
						typename DiscretePressureFunctionType::DiscreteFunctionSpaceType*/ >
	    MatrixAdapterType;
#endif // STOKES_USE_ISTL

	typedef DiscreteVelocityFunctionType RowDiscreteFunctionType;
	typedef DiscreteVelocityFunctionType ColDiscreteFunctionType;

	class PreconditionMatrix : public PreconditionMatrixBaseType {
		const ThisType& a_operator_;

		public:
			PreconditionMatrix( const ThisType& a_operator)
				: PreconditionMatrixBaseType( a_operator.space_, a_operator.space_ ),
				a_operator_( a_operator )
			{
				DiscreteVelocityFunctionType precondition_diagonal( "diag1", a_operator_.space_ );
#if ! STOKES_USE_ISTL
                //!TODO
				a_operator_.getDiag( precondition_diagonal );
#endif
				Stuff::invertFunctionDofs( precondition_diagonal );
				setMatrixDiag( PreconditionMatrixBaseType::matrix(), precondition_diagonal );
			}

	#ifdef USE_BFG_CG_SCHEME
			template <class VECtype>
            void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& /*info*/ ) const
			{
				multOEM(x,ret);
			}
	#endif
			template <class VecType>
			void multOEM( const VecType* tmp, VecType* dest ) const
			{
				precondition(tmp,dest);
			}

			template <class VecType>
			void precondition( const VecType* tmp, VecType* dest ) const
			{
				PreconditionMatrixBaseType::matrix().multOEM( tmp, dest );
			}

			bool rightPrecondition() const { return false; }
	};

    public:
		/** The operator needs the


        **/
        MatrixA_Operator ( const WMatType& w_mat,
                const MMatType& m_mat,
                const XMatType& x_mat,
                const YMatType& y_mat,
				const YMatType& o_mat,
                const typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType& sig_space,
                const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& space)
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            o_mat_(o_mat),
            sig_tmp1( "sig_tmp1", sig_space ),
            sig_tmp2( "sig_tmp2", sig_space ),
            space_(space),
            precondition_matrix_( *this )
        #if STOKES_USE_ISTL
            , adapter_( *this, space, space )
        #endif // STOKES_USE_ISTL
		{}

        ~MatrixA_Operator()
        {}

	  //! ret = ( ( X * ( -1* ( M_inv * ( W * x ) ) ) ) + ( Y + O ) * x ) )
	template <class VECtype>//FEM
	void multOEM(const VECtype* x, VECtype*  ret) const
        {
            w_mat_.multOEM( x, sig_tmp1.leakPointer() );
            m_mat_.apply( sig_tmp1, sig_tmp2 );//Stuff:DiagmUlt

            sig_tmp2 *= ( -1 );
            x_mat_.multOEM( sig_tmp2.leakPointer(), ret );
            y_mat_.multOEMAdd( x, ret );
	    o_mat_.multOEMAdd( x, ret );
        }

	//! ret = ( ( X * ( -1* ( M_inv * ( W * x ) ) ) ) + ( Y + O ) * x ) )
//    template <class BlockType,class A>//ISTL
    void multOEM(const typename RowDiscreteFunctionType:: DofStorageType& x, typename RowDiscreteFunctionType:: DofStorageType&  ret) const
	{
	    w_mat_.apply( x, sig_tmp1.blockVector() );
	    m_mat_.apply( sig_tmp1, sig_tmp2 );//Stuff:DiagmUlt

	    sig_tmp2 *= ( -1 );
	    x_mat_.apply( sig_tmp2.blockVector(), ret );
	    y_mat_.applyAdd( x, ret );
	    o_mat_.applyAdd( x, ret );
	}

#ifdef USE_BFG_CG_SCHEME
        template <class VECtype>
	void multOEM(const VECtype* x, VECtype*  ret, const IterationInfo& /*info*/ ) const
        {
            multOEM(x,ret);
        }
#endif

#if STOKES_USE_ISTL
    //! called by matrix adater in ISTL case
	template <class NonPointerLeakPointerType>
	void mv( const NonPointerLeakPointerType& x, NonPointerLeakPointerType& ret )
	{
	    multOEM( x, ret );
	}
    template <class BlockType>
    void usmv( const typename MatrixAdapterType::field_type alpha,
               const Dune::BlockVector<BlockType>& x,
               Dune::BlockVector<BlockType>& ret ) const
	{
        Dune::BlockVector<BlockType> tmp( ret );
        multOEM(x,tmp);
        tmp *= alpha;
        ret += tmp;
	}

    void applyAdd ( const typename MatrixAdapterType::field_type alpha,
                    const DiscreteVelocityFunctionType& x,
                    DiscreteVelocityFunctionType& ret ) const
    {
        applyAdd( alpha, x.blockVector(), ret.blockVector() );
    }
#endif // STOKES_USE_ISTL

    double ddotOEM(const double*v, const double* w) const
	{
	    DiscreteVelocityFunctionType V( "ddot V", space_, v );
	    DiscreteVelocityFunctionType W( "ddot W", space_, w );
	    return V.scalarProductDofs( W );
	}

	void apply( const DiscreteVelocityFunctionType& rhs, DiscreteVelocityFunctionType& dest ) const
	{
	    multOEM( rhs.leakPointer(), dest.leakPointer() )	;
	}

    ThisType& systemMatrix () { return *this; }
    const ThisType& systemMatrix () const { return *this; }
#if STOKES_USE_ISTL
    const MatrixAdapterType& matrixAdapter() const { return adapter_; }
#endif // STOKES_USE_ISTL
    const PreconditionMatrix& preconditionMatrix() const { return precondition_matrix_; }

    bool hasPreconditionMatrix () const
    {
        return Parameters().getParam( "innerPrecond", false );
    }

	//! return diagonal entries of this matrix
	template <class DiscFuncType>
	void getDiag(DiscFuncType &precondition_diagonal_) const
	{
		x_mat_.getDiag( m_mat_, w_mat_, precondition_diagonal_);
		precondition_diagonal_ *= -1;
		y_mat_.addDiag( precondition_diagonal_ );
		o_mat_.addDiag( precondition_diagonal_ );
	}

	protected:
        const WMatType& w_mat_;
        const MMatType& m_mat_;
        const XMatType& x_mat_;
        const YMatType& y_mat_;
	const YMatType& o_mat_;
        mutable DiscreteSigmaFunctionType sig_tmp1;
        mutable DiscreteSigmaFunctionType sig_tmp2;
	const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& space_;
	PreconditionMatrix precondition_matrix_;
#if STOKES_USE_ISTL
	MatrixAdapterType adapter_;
#endif // STOKES_USE_ISTL

};



/** \brief wraps the solution of inner CG iteration and exposes only the minimally needed interface to
		\ref SaddlepointInverseOperator or \ref NestedCgSaddlepointInverseOperator respectively

  **/
template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteSigmaFunctionType,
            class DiscreteVelocityFunctionType >
class InnerCGSolverWrapper {
    public:
        typedef MatrixA_Operator<   WMatType,
                                    MMatType,
                                    XMatType,
                                    YMatType,
                                    DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType>
                A_OperatorType;

		typedef SOLVER_NAMESPACE::INNER_CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
            CG_SolverType;
		#ifdef USE_BFG_CG_SCHEME
			typedef typename CG_SolverType::ReturnValueType
				ReturnValueType;
		#endif

		InnerCGSolverWrapper( const WMatType& w_mat,
                const MMatType& m_mat,
                const XMatType& x_mat,
                const YMatType& y_mat,
				const YMatType& o_mat,
                const typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType& sig_space,
				const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& space,
                const double relLimit,
                const double absLimit,
                const bool verbose )
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            sig_tmp1( "sig_tmp1", sig_space ),
            sig_tmp2( "sig_tmp2", sig_space ),
            a_op_( w_mat, m_mat, x_mat, y_mat, o_mat, sig_space, space ),
            cg_solver( a_op_,   relLimit,
                                absLimit,
                                2000, //inconsequential anyways
                                verbose )
        {}

		/** \brief this signature is called if the CG solver uses non-standard third arg to expose runtime info
			\see SaddlepointInverseOperator (when compiled with BFG scheme support)
			**/
	#ifdef USE_BFG_CG_SCHEME
		void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest, ReturnValueType& ret )
		{
			cg_solver.apply(arg,dest, ret);
		}
	#endif

		//! the standard function call
        void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest )
        {
            cg_solver.apply(arg,dest);
        }
        //! the standard function call


#ifdef USE_BFG_CG_SCHEME
		/** - needed when using the BFG scheme
			- only works with modified Dune Solvers exposing a setAbsoluteLimit of their own
			**/
        void setAbsoluteLimit( const double abs )
        {
            cg_solver.setAbsoluteLimit( abs );
        }
#endif
		const A_OperatorType& getOperator() const { return a_op_;}

    private:
//        const MMatType precond_;
        const WMatType& w_mat_;
        const MMatType& m_mat_;
        const XMatType& x_mat_;
        const YMatType& y_mat_;
        mutable DiscreteSigmaFunctionType sig_tmp1;
        mutable DiscreteSigmaFunctionType sig_tmp2;
        A_OperatorType a_op_;
        CG_SolverType cg_solver;
};

}

#endif // DUNE_STOKES_SOLVER_CGHELPER_HH
