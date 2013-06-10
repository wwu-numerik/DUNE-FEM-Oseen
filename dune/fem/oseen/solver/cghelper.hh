/**
 *  \file   stuff.hh
 *  \brief  contains some stuff
 **/

#ifndef DUNE_OSEEN_SOLVER_CGHELPER_HH
#define DUNE_OSEEN_SOLVER_CGHELPER_HH

#include <cmake_config.h>
#include <dune/fem/oseen/oemsolver/oemsolver.hh>
#include <dune/fem/oseen/solver/new_bicgstab.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/matrix_object.hh>
#include <dune/stuff/fem/functions/transform.hh>

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
    typedef DSC::IdentityMatrixObject<YMatType>
		PreconditionMatrixBaseType;

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
                //!TODO
				a_operator_.getDiag( precondition_diagonal );
                DSFe::invertFunctionDofs( precondition_diagonal );
                DSC::setMatrixDiag( PreconditionMatrixBaseType::matrix(), precondition_diagonal );
			}

			template <class VECtype>
            void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& /*info*/ ) const
			{
				multOEM(x,ret);
			}

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
        y_mat_.matrix().multOEMAdd( x, ret );
        o_mat_.matrix().multOEMAdd( x, ret );
    }

    template <class VECtype>
	void multOEM(const VECtype* x, VECtype*  ret, const IterationInfo& /*info*/ ) const
    {
        multOEM(x,ret);
    }

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

    const PreconditionMatrix& preconditionMatrix() const { return precondition_matrix_; }

    bool hasPreconditionMatrix () const
    {
        return DSC_CONFIG_GET( "innerPrecond", false );
    }

	//! return diagonal entries of this matrix
	template <class DiscFuncType>
	void getDiag(DiscFuncType &precondition_diagonal_) const
	{
        x_mat_.matrix().getDiag( m_mat_.matrix(), w_mat_.matrix(), precondition_diagonal_);
		precondition_diagonal_ *= -1;
        y_mat_.matrix().addDiag( precondition_diagonal_ );
        o_mat_.matrix().addDiag( precondition_diagonal_ );
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
        typedef typename CG_SolverType::ReturnValueType
            ReturnValueType;

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
		void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest, ReturnValueType& ret )
		{
			cg_solver.apply(arg,dest, ret);
		}

		//! the standard function call
        void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest )
        {
            cg_solver.apply(arg,dest);
        }
        //! the standard function call

		/** - needed when using the BFG scheme
			- only works with modified Dune Solvers exposing a setAbsoluteLimit of their own
			**/
        void setAbsoluteLimit( const double abs )
        {
            cg_solver.setAbsoluteLimit( abs );
        }

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

#endif // DUNE_OSEEN_SOLVER_CGHELPER_HH

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

