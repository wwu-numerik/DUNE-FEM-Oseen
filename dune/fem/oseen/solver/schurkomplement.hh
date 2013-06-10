#ifndef DUNE_OSEEN_SOLVER_SCHURKOMPLEMENT_HH
#define DUNE_OSEEN_SOLVER_SCHURKOMPLEMENT_HH

#include <cmake_config.h>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/fem/preconditioning.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/fem/functions/analytical.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include <dune/fem/oseen/oemsolver/preconditioning.hh>
#include <dune/fem/oseen/oemsolver/oemsolver.hh>

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

template < class SchurkomplementOperatorType >
class PreconditionOperatorDefault {
    typedef PreconditionOperatorDefault< SchurkomplementOperatorType >
        ThisType;
    public:
        typedef typename SchurkomplementOperatorType::DiscretePressureFunctionType
            RowDiscreteFunctionType;
        typedef typename RowDiscreteFunctionType::DiscreteFunctionSpaceType
            DomainSpaceType;
        typedef typename SchurkomplementOperatorType::DiscretePressureFunctionType
            ColDiscreteFunctionType;
        typedef typename ColDiscreteFunctionType::DiscreteFunctionSpaceType
            RangeSpaceType;
    private:
        const SchurkomplementOperatorType& sk_op_;
        const typename SchurkomplementOperatorType::A_PreconditionMatrix& a_precond_;
        mutable typename SchurkomplementOperatorType::DiscreteVelocityFunctionType velo_tmp;
        mutable typename SchurkomplementOperatorType::DiscreteVelocityFunctionType velo_tmp2;

        const typename SchurkomplementOperatorType::Z_MatrixType::DomainSpaceType& pressure_space_;
        const typename SchurkomplementOperatorType::E_MatrixType::DomainSpaceType& velocity_space_;

    public:

        PreconditionOperatorDefault( const typename SchurkomplementOperatorType::A_SolverType& a_solver,
                           const SchurkomplementOperatorType& sk_op,
                           const typename SchurkomplementOperatorType::E_MatrixType::DomainSpaceType& velocity_space,
                           const typename SchurkomplementOperatorType::Z_MatrixType::DomainSpaceType& pressure_space)
            : sk_op_( sk_op),
            a_precond_( a_solver.getOperator().preconditionMatrix() ),
            velo_tmp( "sdeio", pressure_space ),
            velo_tmp2( "2sdeio", pressure_space ),
            pressure_space_(pressure_space),
            velocity_space_(velocity_space)
        {}

        template <class VECtype>
        void multOEM(const VECtype* x, VECtype* ret, const IterationInfo& ) const
        {
            multOEM(x,ret);
        }

        template <class VECtype>
        void multOEM(const VECtype* x, VECtype* ret) const
        {
            sk_op_.z_mat_.matrix().multOEM( x,velo_tmp.leakPointer());
            a_precond_.apply( velo_tmp, velo_tmp2);
            sk_op_.e_mat_.matrix().multOEM( velo_tmp2.leakPointer(),ret);
            sk_op_.r_mat_.matrix().multOEMAdd( x, ret);
        }

        ThisType& systemMatrix () { return *this; }
        const ThisType& systemMatrix () const { return *this; }

        double ddotOEM(const double*v, const double* w) const
        {
            typename SchurkomplementOperatorType::DiscretePressureFunctionType V( "ddot_V2", velocity_space_, v );
            typename SchurkomplementOperatorType::DiscretePressureFunctionType W( "ddot_W1", velocity_space_, w );
            return V.scalarProductDofs( W );
        }
};

/** \brief Operator wrapping Matrix vector multiplication for
			matrix \f$ S :=  B_t * A^-1 * B + rhs3 \f$
			**/
template <  class A_SolverImp,
            class E_MatrixImp,
			class R_MatrixType,
            class Z_MatrixImp,
			class M_invers_MatrixType,
            class DiscreteVelocityFunctionImp ,
            class DiscretePressureFunctionImp>
class SchurkomplementOperator //: public StokesOEMSolver::PreconditionInterface
{

        typedef DiscreteVelocityFunctionImp DiscreteVelocityFunctionType;
        typedef DiscretePressureFunctionImp DiscretePressureFunctionType;
        typedef A_SolverImp A_SolverType;
        typedef Z_MatrixImp Z_MatrixType;
        typedef E_MatrixImp E_MatrixType;
		typedef SchurkomplementOperator <   A_SolverType,
                                            E_MatrixType,
											R_MatrixType,
                                            Z_MatrixType,
											M_invers_MatrixType,
											DiscreteVelocityFunctionType,
											DiscretePressureFunctionType>
				ThisType;
		typedef typename A_SolverType::A_OperatorType::PreconditionMatrix
		    A_PreconditionMatrix;
        typedef DSC::IdentityMatrixObject<R_MatrixType>
		    PreconditionMatrixBaseType;
        typedef PreconditionOperatorDefault< ThisType >
            PreconditionOperator;
        //if shit goes south wrt precond working check if this doesn't need to be OEmSolver instead of StokesOEMSolver
        friend class Conversion<ThisType,StokesOEMSolver::PreconditionInterface>;
        friend class PreconditionOperatorDefault< ThisType >;


        typedef DSFe::OperatorBasedPreconditioner< PreconditionOperator,
                                                    DuneStokes::OUTER_CG_SOLVERTYPE,
													DiscretePressureFunctionType>
			PreconditionMatrix;
		typedef DiscretePressureFunctionType RowDiscreteFunctionType;
		typedef DiscretePressureFunctionType ColDiscreteFunctionType;

        typedef typename A_SolverType::ReturnValueType
                ReturnValueType;
public:
        SchurkomplementOperator( A_SolverType& a_solver,
								const E_MatrixType& e_mat,
								const R_MatrixType& r_mat,
								const Z_MatrixType& z_mat,
								const M_invers_MatrixType& m_inv_mat,
                                const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& velocity_space,
                                const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space )
            : a_solver_(a_solver),
			e_mat_(e_mat),
			r_mat_(r_mat),
			z_mat_(z_mat),
			m_inv_mat_(m_inv_mat),
			tmp1 ( "schurkomplementoperator_tmp1", velocity_space ),
			tmp2 ( "schurkomplementoperator_tmp2", velocity_space ),
            do_bfg( DSC_CONFIG_GET( "do-bfg", true ) ),
            total_inner_iterations( 0 ),
			pressure_space_(pressure_space),
            precond_operator_( a_solver, *this, pressure_space, velocity_space),
            precond_( precond_operator_, pressure_space )
	{}

    double ddotOEM(const double*v, const double* w) const
	{
	    DiscretePressureFunctionType V( "ddot_V2", pressure_space_, v );
	    DiscretePressureFunctionType W( "ddot_W1", pressure_space_, w );
	    assert( !DSFe::FunctionContainsNanOrInf( V ) );
	    assert( !DSFe::FunctionContainsNanOrInf( W ) );
	    const double ret = V.scalarProductDofs( W );
	    assert( std::isfinite( ret ) );
	    return ret;
	}

	void multOEM(const double *x, double * ret) const
        {
			// ret = ( ( -E * ( A^-1 * ( Z * x ) ) ) + ( R * x ) )
			z_mat_.multOEM( x, tmp1.leakPointer() );

			tmp2.clear();//don't remove w/o result testing
			assert( !DSFe::FunctionContainsNanOrInf( tmp1 ) );

            auto& info = DSC_LOG_INFO;
			const int solverVerbosity = DSC_CONFIG_GET( "solverVerbosity", 0 );
            ReturnValueType cg_info;
            a_solver_.apply( tmp1, tmp2, cg_info );
            if ( solverVerbosity > 1 )
                info << "\t\t\t\t\t inner iterations: " << cg_info.first << std::endl;
            total_inner_iterations += cg_info.first;

			assert( !DSFe::FunctionContainsNanOrInf( tmp2 ) );
			tmp2 *= -1;
			e_mat_.multOEM( tmp2.leakPointer(), ret );
			assert( !DSFe::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
            r_mat_.matrix().multOEMAdd( x, ret );
			assert( !DSFe::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
        }

	void apply( const DiscretePressureFunctionType& arg, DiscretePressureFunctionType& ret ) const
	{
	    assert( !DSFe::FunctionContainsNanOrInf( arg ) );
	    typedef typename DiscretePressureFunctionType::DofStorageType D;
	    typedef typename DiscretePressureFunctionType::RangeFieldType R;
        multOEM( arg.leakPointer(), ret.leakPointer() );
	    assert( !DSFe::FunctionContainsNanOrInf( ret ) );
	}

    template <class VECtype>
    void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& info ) const
    {
        if ( do_bfg ) {
            static const double tau = DSC_CONFIG_GET( "bfg-tau", 0.1 );
            double limit = info.second.first;
            const double residuum = fabs(info.second.second);
            const int n = info.first;

            if ( n == 0 )
                limit = DSC_CONFIG_GET( "absLimit", 10e-12 );
            limit = tau * std::min( 1. , limit / std::min ( std::pow( residuum, n ) , 1.0 ) );

            a_solver_.setAbsoluteLimit( limit );
            DSC_LOG_INFO << "\t\t\t Set inner error limit to: "<< limit << std::endl;
        }
        multOEM( x, ret );
    }

    ThisType& systemMatrix () { return *this; }
    const ThisType& systemMatrix () const { return *this; }
    const PreconditionMatrix& preconditionMatrix() const { return precond_; }

    bool hasPreconditionMatrix () const
    {
        return DSC_CONFIG_GET( "outerPrecond", false );
    }

    bool rightPrecondition() const { return false; }
    long getTotalInnerIterations() const { return total_inner_iterations; }

    private:
		A_SolverType& a_solver_;
		const E_MatrixType& e_mat_;
		const R_MatrixType& r_mat_;
		const Z_MatrixType& z_mat_;
		const M_invers_MatrixType& m_inv_mat_;
        mutable DiscreteVelocityFunctionType tmp1;
        mutable DiscreteVelocityFunctionType tmp2;
        bool do_bfg;
        mutable long total_inner_iterations;
		const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space_;
		PreconditionOperator precond_operator_;
		PreconditionMatrix precond_;
};

} //end namespace Dune

#endif // DUNE_OSEEN_SOLVER_SCHURKOMPLEMENT_HH

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

