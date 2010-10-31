/**
 *  \file   stuff.hh
 *  \brief  contains some stuff
 **/

#ifndef DUNE_STOKES_SOLVER_CGHELPER_HH
#define DUNE_STOKES_SOLVER_CGHELPER_HH

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
class MatrixA_Operator //: public OEMSolver::PreconditionInterface
	{

    public:

        typedef MatrixA_Operator<   WMatType,
                        MMatType,
                        XMatType,
                        YMatType,
                        DiscreteSigmaFunctionType,
						DiscreteVelocityFunctionType>
                    ThisType;
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
            precondition_matrix_( y_mat_.rows(), y_mat_.cols(), 10 ),
			precondition_diagonal_( "diag1", space ),
			precondition_matrix_invers( y_mat_.cols(), y_mat_.rows(), 10 )
        {
			x_mat_.getDiag( m_mat_, w_mat_, precondition_diagonal_);
			precondition_diagonal_ *= -1;
			y_mat_.addDiag( precondition_diagonal_ );
			o_mat_.addDiag( precondition_diagonal_ );
			setMatrixDiag( precondition_matrix_, precondition_diagonal_ );
			DiscreteVelocityFunctionType precondition_diagonal_inv( "diag_inv", space );
			precondition_diagonal_inv.assign( precondition_diagonal_ );
			Stuff::invertFunctionDofs( precondition_diagonal_inv );
			setMatrixDiag( precondition_matrix_invers, precondition_diagonal_inv );
		}

        ~MatrixA_Operator()
        {}

		//! ret = ( ( X * ( -1* ( M_inv * ( W * x ) ) ) ) + ( Y + O ) * x ) )
        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            sig_tmp1.clear();
            sig_tmp2.clear();

            w_mat_.multOEM( x, sig_tmp1.leakPointer() );
            m_mat_.apply( sig_tmp1, sig_tmp2 );//Stuff:DiagmUlt

            sig_tmp2 *= ( -1 );
            x_mat_.multOEM( sig_tmp2.leakPointer(), ret );
            y_mat_.multOEMAdd( x, ret );
			o_mat_.multOEMAdd( x, ret );
        }

#ifdef USE_BFG_CG_SCHEME
        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& info ) const
        {
            multOEM(x,ret);
        }
#endif
        double ddotOEM(const double*v, const double* w) const
		{
	        DiscreteVelocityFunctionType V( "ddot V", space_, v );
	        DiscreteVelocityFunctionType W( "ddot W", space_, w );
	        return V.scalarProductDofs( W );
		}


        ThisType& systemMatrix ()
        {
            return *this;
        }

        YMatType& preconditionMatrix()
        {
            return precondition_matrix_;
        }

        bool hasPreconditionMatrix () const
        {
            return true;
        }

        bool rightPrecondition() const
        {
            return false;
        }

        template <class VecType>
        void precondition( const VecType* tmp, VecType* dest ) const
        {
			assert( false );
			precondition_matrix_invers.multOEM( tmp, dest );
        }


    private:
        const WMatType& w_mat_;
        const MMatType& m_mat_;
        const XMatType& x_mat_;
        const YMatType& y_mat_;
		const YMatType& o_mat_;
        mutable DiscreteSigmaFunctionType sig_tmp1;
        mutable DiscreteSigmaFunctionType sig_tmp2;
		const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& space_;
		YMatType precondition_matrix_;
		DiscreteVelocityFunctionType precondition_diagonal_;
		YMatType precondition_matrix_invers;
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

#ifdef USE_BFG_CG_SCHEME
		/** - needed when using the BFG scheme
			- only works with modified Dune Solvers exposing a setAbsoluteLimit of their own
			**/
        void setAbsoluteLimit( const double abs )
        {
            cg_solver.setAbsoluteLimit( abs );
        }
#endif

    private:
        const MMatType precond_;
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
