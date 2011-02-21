#ifndef DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH
#define DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH

#include <dune/stuff/matrix.hh>
#include <dune/stuff/preconditioning.hh>

namespace Dune {

/** \brief Operator wrapping Matrix vector multiplication for
			matrix \f$ S :=  B_t * A^-1 * B + rhs3 \f$
			**/
template <  class A_SolverType,
			class E_MatrixType,
			class R_MatrixType,
			class Z_MatrixType,
			class M_invers_MatrixType,
			class DiscreteVelocityFunctionType ,
            class DiscretePressureFunctionType>
class SchurkomplementOperator //: public SOLVER_INTERFACE_NAMESPACE::PreconditionInterface
{
    public:
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
		typedef IdentityMatrixObject<typename R_MatrixType::WrappedMatrixObjectType>
			PreconditionMatrixBaseType;
		friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;

		class PreconditionOperator {
			const ThisType& sk_op_;
			const A_PreconditionMatrix& a_precond_;
			mutable DiscreteVelocityFunctionType velo_tmp;
			mutable DiscreteVelocityFunctionType velo_tmp2;
			const typename Z_MatrixType::WrappedMatrixObjectType::RangeSpaceType& pressure_space_;

			public:
				PreconditionOperator( const A_SolverType& a_solver,
								   const ThisType& sk_op,
								   const typename E_MatrixType::WrappedMatrixObjectType::RangeSpaceType& velocity_space,
								   const typename Z_MatrixType::WrappedMatrixObjectType::RangeSpaceType& pressure_space)
					: sk_op_( sk_op),
					a_precond_( a_solver.getOperator().preconditionMatrix() ),
					velo_tmp( "sdeio", velocity_space ),
					velo_tmp2( "2sdeio", velocity_space ),
					pressure_space_(pressure_space)
				{}

#ifdef USE_BFG_CG_SCHEME
				template <class VECtype>
				void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& ) const
				{
					multOEM(x,ret);
				}
#endif
				template <class VECtype>
				void multOEM(const VECtype *x, VECtype * ret) const
				{
					sk_op_.z_mat_.matrix().multOEM( x,velo_tmp.leakPointer());
					a_precond_.apply( velo_tmp, velo_tmp2);
					sk_op_.e_mat_.matrix().multOEM( velo_tmp2.leakPointer(),ret);
					sk_op_.r_mat_.matrix().multOEMAdd( x, ret);
				}

				PreconditionOperator& systemMatrix () { return *this; }

				double ddotOEM(const double*v, const double* w) const
				{
					DiscretePressureFunctionType V( "ddot_V2", pressure_space_, v );
					DiscretePressureFunctionType W( "ddot_W1", pressure_space_, w );
					return V.scalarProductDofs( W );
				}
		};

		typedef Stuff::OperatorBasedPreconditioner< PreconditionOperator,
													SOLVER_NAMESPACE::OUTER_CG_SOLVERTYPE,
													DiscretePressureFunctionType>
			PreconditionMatrix;

#ifdef USE_BFG_CG_SCHEME
        typedef typename A_SolverType::ReturnValueType
                ReturnValueType;
#endif

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
            tmp1 ( "tmp1", velocity_space ),
            tmp2 ( "tmp2", velocity_space ),
            do_bfg( Parameters().getParam( "do-bfg", true ) ),
            total_inner_iterations( 0 ),
			pressure_space_(pressure_space),
			precond_operator_( a_solver, *this, velocity_space, pressure_space),
			precond_( precond_operator_, pressure_space_ )
		{}

        double ddotOEM(const double*v, const double* w) const
		{
			DiscretePressureFunctionType V( "ddot_V2", pressure_space_, v );
			DiscretePressureFunctionType W( "ddot_W1", pressure_space_, w );
	        return V.scalarProductDofs( W );
		}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
			// ret = ( ( -E * ( A^-1 * ( Z * x ) ) ) + ( R * x ) )
			z_mat_.multOEM( x, tmp1.leakPointer() );

#ifdef USE_BFG_CG_SCHEME
			Logging::LogStream& info = Logger().Info();
			const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
            ReturnValueType cg_info;
            a_solver_.apply( tmp1, tmp2, cg_info );
            if ( solverVerbosity > 1 )
                info << "\t\t\t\t\t inner iterations: " << cg_info.first << std::endl;

            total_inner_iterations += cg_info.first;
#else
			a_solver_.apply( tmp1, tmp2 );
#endif
			tmp2 *= -1;
			e_mat_.multOEM( tmp2.leakPointer(), ret );
			r_mat_.multOEMAdd( x, ret );
        }

		void apply( const DiscretePressureFunctionType& arg, DiscretePressureFunctionType& ret ) const
		{
			multOEM( arg.leakPointer(), ret.leakPointer() );
		}

#ifdef USE_BFG_CG_SCHEME
        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& info ) const
        {
            if ( do_bfg ) {
                static const double tau = Parameters().getParam( "bfg-tau", 0.1 );
                double limit = info.second.first;
                const double residuum = fabs(info.second.second);
                const int n = info.first;

                if ( n == 0 )
                    limit = Parameters().getParam( "absLimit", 10e-12 );
                limit = tau * std::min( 1. , limit / std::min ( std::pow( residuum, n ) , 1.0 ) );

                a_solver_.setAbsoluteLimit( limit );
                Logger().Info() << "\t\t\t Set inner error limit to: "<< limit << std::endl;
            }
            multOEM( x, ret );
        }
#endif

        ThisType& systemMatrix ()
        {
            return *this;
        }

		const PreconditionMatrix& preconditionMatrix() const
        {
            return precond_;
        }

        bool hasPreconditionMatrix () const
        {
			return Parameters().getParam( "outerPrecond", false );
        }

        bool rightPrecondition() const
        {
            return false;
        }

		long getTotalInnerIterations() const
        {
            return total_inner_iterations;
        }

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

#endif // DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH
