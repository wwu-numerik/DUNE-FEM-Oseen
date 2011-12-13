#ifndef DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH
#define DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH

#include "solver_defines.hh"

#include <dune/stuff/matrix.hh>
#include <dune/stuff/preconditioning.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/parametercontainer.hh>

#if HAVE_DUNE_ISTL
#   include <dune/istl/operators.hh>
#   include <dune/fem/operator/matrix/istlmatrix.hh>
#   include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#endif

namespace Dune {

template < class SchurkomplementOperatorType >
class PreconditionOperatorDefault {
    typedef PreconditionOperatorDefault< SchurkomplementOperatorType >
        ThisType;
    public:
    #if STOKES_USE_ISTL
        typedef SchurkomplementOperatorAdapter<ThisType>
            MatrixAdapterType;
    #endif // STOKES_USE_ISTL
        typedef typename SchurkomplementOperatorType::DiscretePressureFunctionType RowDiscreteFunctionType;
        typedef typename SchurkomplementOperatorType::DiscretePressureFunctionType ColDiscreteFunctionType;
    private:
        const SchurkomplementOperatorType& sk_op_;
        const typename SchurkomplementOperatorType::A_PreconditionMatrix& a_precond_;
        mutable typename SchurkomplementOperatorType::DiscreteVelocityFunctionType velo_tmp;
        mutable typename SchurkomplementOperatorType::DiscreteVelocityFunctionType velo_tmp2;
        const typename SchurkomplementOperatorType::Z_MatrixType::WrappedMatrixObjectType::DomainSpaceType& pressure_space_;
        const typename SchurkomplementOperatorType::E_MatrixType::WrappedMatrixObjectType::DomainSpaceType& velocity_space_;
    #if STOKES_USE_ISTL
        MatrixAdapterType matrix_adapter_;
    #endif // STOKES_USE_ISTL

    public:

        PreconditionOperatorDefault( const typename SchurkomplementOperatorType::A_SolverType& a_solver,
                           const SchurkomplementOperatorType& sk_op,
                           const typename SchurkomplementOperatorType::E_MatrixType::WrappedMatrixObjectType::DomainSpaceType& velocity_space,
                           const typename SchurkomplementOperatorType::Z_MatrixType::WrappedMatrixObjectType::DomainSpaceType& pressure_space)
            : sk_op_( sk_op),
            a_precond_( a_solver.getOperator().preconditionMatrix() ),
            velo_tmp( "sdeio", pressure_space ),
            velo_tmp2( "2sdeio", pressure_space ),
            pressure_space_(pressure_space),
            velocity_space_(velocity_space)
        #if STOKES_USE_ISTL
            ,matrix_adapter_( *this, pressure_space, pressure_space )
        #endif // STOKES_USE_ISTL
        {}

#ifdef USE_BFG_CG_SCHEME
        template <class VECtype>
        void multOEM(const VECtype* x, VECtype* ret, const IterationInfo& ) const
        {
            multOEM(x,ret);
        }
#endif
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

    #if STOKES_USE_ISTL
        MatrixAdapterType& matrixAdapter() { return matrix_adapter_; }
        const MatrixAdapterType& matrixAdapter() const { return matrix_adapter_; }

        template <class ArgDofStorageType, class DestDofStorageType, class ArgRangeFieldType, class DestRangeFieldType>
        void multOEM(const Dune::StraightenBlockVector<ArgDofStorageType,ArgRangeFieldType> &x,
                 Dune::StraightenBlockVector<DestDofStorageType,DestRangeFieldType> &ret) const
        {
            multOEM( x.blockVector(), ret.blockVector() );
        }
        template <class ArgDofStorageType, class DestDofStorageType>
        void multOEM(const Dune::BlockVector<ArgDofStorageType> &x,
                 Dune::BlockVector<DestDofStorageType> &ret) const
        {
            sk_op_.z_mat_.multOEM( x,velo_tmp.blockVector());
            a_precond_.apply( velo_tmp, velo_tmp2);
            sk_op_.e_mat_.multOEM( velo_tmp2.blockVector(),ret);
            sk_op_.r_mat_.multOEMAdd( x, ret);
        }

        //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
        template < class field_type, class Arg, class Dest >
        void usmv( const field_type alpha, const Arg& arg, Dest& dest ) const
        {
            Dest tmp( dest );
            multOEM( arg, tmp );
            tmp *= alpha;
            dest += tmp;
        }
    #endif // STOKES_USE_ISTL

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
class SchurkomplementOperator //: public SOLVER_INTERFACE_NAMESPACE::PreconditionInterface
{
    public:
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
		typedef SchurkomplementOperatorAdapter< ThisType >
		    MatrixAdapterType;
		typedef typename A_SolverType::A_OperatorType::PreconditionMatrix
		    A_PreconditionMatrix;
		typedef IdentityMatrixObject<typename R_MatrixType::WrappedMatrixObjectType>
		    PreconditionMatrixBaseType;
        typedef PreconditionOperatorDefault< ThisType >
            PreconditionOperator;
		//if shit goes south wrt precond working check if this doesn't need to be OEmSolver instead of SOLVER_INTERFACE_NAMESPACE
		friend class Conversion<ThisType,SOLVER_INTERFACE_NAMESPACE::PreconditionInterface>;
        friend class PreconditionOperatorDefault< ThisType >;


		typedef Stuff::OperatorBasedPreconditioner< PreconditionOperator,
													SOLVER_NAMESPACE::OUTER_CG_SOLVERTYPE,
													DiscretePressureFunctionType>
			PreconditionMatrix;
		typedef DiscretePressureFunctionType RowDiscreteFunctionType;
		typedef DiscretePressureFunctionType ColDiscreteFunctionType;

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
			tmp1 ( "schurkomplementoperator_tmp1", velocity_space ),
			tmp2 ( "schurkomplementoperator_tmp2", velocity_space ),
            do_bfg( Parameters().getParam( "do-bfg", true ) ),
            total_inner_iterations( 0 ),
			pressure_space_(pressure_space),
            precond_operator_( a_solver, *this, pressure_space, velocity_space),
            precond_( precond_operator_, pressure_space )
        #if STOKES_USE_ISTL
            , adapter_( *this, pressure_space, pressure_space )
        #endif // STOKES_USE_ISTL
	{}

    double ddotOEM(const double*v, const double* w) const
	{
	    DiscretePressureFunctionType V( "ddot_V2", pressure_space_, v );
	    DiscretePressureFunctionType W( "ddot_W1", pressure_space_, w );
	    assert( !Stuff::FunctionContainsNanOrInf( V ) );
	    assert( !Stuff::FunctionContainsNanOrInf( W ) );
	    const double ret = V.scalarProductDofs( W );
	    assert( std::isfinite( ret ) );
	    return ret;
	}

	void multOEM(const double *x, double * ret) const
        {
			// ret = ( ( -E * ( A^-1 * ( Z * x ) ) ) + ( R * x ) )
			z_mat_.multOEM( x, tmp1.leakPointer() );

			tmp2.clear();//don't remove w/o result testing
			assert( !Stuff::FunctionContainsNanOrInf( tmp1 ) );
#ifdef USE_BFG_CG_SCHEME
			Stuff::Logging::LogStream& info = Logger().Info();
			const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
            ReturnValueType cg_info;
            a_solver_.apply( tmp1, tmp2, cg_info );
            if ( solverVerbosity > 1 )
                info << "\t\t\t\t\t inner iterations: " << cg_info.first << std::endl;

            total_inner_iterations += cg_info.first;
#else
			a_solver_.apply( tmp1, tmp2 );
#endif
			assert( !Stuff::FunctionContainsNanOrInf( tmp2 ) );
			tmp2 *= -1;
			e_mat_.multOEM( tmp2.leakPointer(), ret );
			assert( !Stuff::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
			r_mat_.multOEMAdd( x, ret );
			assert( !Stuff::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
//			ret[0] = x[0];
        }
#if STOKES_USE_ISTL
    template <class ArgDofStorageType, class DestDofStorageType>
    void multOEM(const Dune::BlockVector<ArgDofStorageType> &x,
             Dune::BlockVector<DestDofStorageType> &ret) const
    {
        // ret = ( ( -E * ( A^-1 * ( Z * x ) ) ) + ( R * x ) )
        z_mat_.apply( x, tmp1.blockVector() );

        tmp2.clear();//don't remove w/o result testing
        assert( !Stuff::FunctionContainsNanOrInf( tmp1 ) );

        a_solver_.apply( tmp1, tmp2 );

        assert( !Stuff::FunctionContainsNanOrInf( tmp2 ) );
        tmp2 *= -1;
        e_mat_.apply( tmp2.blockVector(), ret );
        assert( !Stuff::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
        r_mat_.multOEMAdd( x, ret );
        assert( !Stuff::FunctionContainsNanOrInf( ret, pressure_space_.size() ) );
//			ret[0] = x[0];
    }
	template <class ArgDofStorageType, class DestDofStorageType, class ArgRangeFieldType, class DestRangeFieldType>
	void multOEM(const Dune::StraightenBlockVector<ArgDofStorageType,ArgRangeFieldType> &x,
             Dune::StraightenBlockVector<DestDofStorageType,DestRangeFieldType> &ret) const
	{
        multOEM( x.blockVector(), ret.blockVector() );
	}

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    template < class field_type, class Arg, class Dest >
    void usmv( const field_type alpha, const Arg& arg, Dest& dest ) const
    {
        Dest tmp( dest );
        multOEM( arg, tmp );
        tmp *= alpha;
        dest += tmp;
    }
#endif // STOKES_USE_ISTL

	void apply( const DiscretePressureFunctionType& arg, DiscretePressureFunctionType& ret ) const
	{
	    assert( !Stuff::FunctionContainsNanOrInf( arg ) );
	    typedef typename DiscretePressureFunctionType::DofStorageType D;
	    typedef typename DiscretePressureFunctionType::RangeFieldType R;
        multOEM( arg.leakPointer(), ret.leakPointer() );
	    assert( !Stuff::FunctionContainsNanOrInf( ret ) );
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

#if STOKES_USE_ISTL
    const MatrixAdapterType& matrixAdapter() const { return adapter_; }
#endif // STOKES_USE_ISTL
    ThisType& systemMatrix () { return *this; }
    const ThisType& systemMatrix () const { return *this; }
    const PreconditionMatrix& preconditionMatrix() const { return precond_; }

    bool hasPreconditionMatrix () const
    {
        return Parameters().getParam( "outerPrecond", false );
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
    #if STOKES_USE_ISTL
		MatrixAdapterType adapter_;
    #endif // STOKES_USE_ISTL
};

} //end namespace Dune

#endif // DUNE_STOKES_SOLVER_SCHURKOMPLEMENT_HH
