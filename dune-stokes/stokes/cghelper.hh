#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED

#include <dune/fem/solver/oemsolver.hh>
#include "../src/logging.hh"
#include "../src/stuff.hh"

namespace Dune {



template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteSigmaFunctionType >
class MatrixA_Operator {

    public:

        typedef MatrixA_Operator<   WMatType,
                        MMatType,
                        XMatType,
                        YMatType,
                        DiscreteSigmaFunctionType >
                    ThisType;

        MatrixA_Operator ( const WMatType& w_mat,
                const MMatType& m_mat,
                const XMatType& x_mat,
                const YMatType& y_mat,
                const typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType& sig_space )
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            sig_tmp1( "sig_tmp1", sig_space ),
            sig_tmp2( "sig_tmp2", sig_space )
        {}

        ~MatrixA_Operator()
        {}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            sig_tmp1.clear();
            sig_tmp2.clear();
            Logging::LogStream& dbg = Logger().Dbg();

			// ret = ( ( X * ( -1* ( M_inv * ( W * x ) ) ) ) + ( Y * x ) )
            w_mat_.multOEM( x, sig_tmp1.leakPointer() );
            m_mat_.apply( sig_tmp1, sig_tmp2 );//Stuff:DiagmUlt
            sig_tmp2 *= ( -1 );
            x_mat_.multOEM( sig_tmp2.leakPointer(), ret );
            y_mat_.multOEMAdd( x, ret );
			//
        }

        ThisType& systemMatrix ()
        {
            return *this;
        }

    private:
        const WMatType& w_mat_;
        const MMatType& m_mat_;
        const XMatType& x_mat_;
        const YMatType& y_mat_;
        mutable DiscreteSigmaFunctionType sig_tmp1;
        mutable DiscreteSigmaFunctionType sig_tmp2;
};


template <  class A_SolverType,
            class B_t_matrixType,
            class CmatrixType,
            class BmatrixType,
            class MmatrixType,
            class DiscreteVelocityFunctionType ,
            class DiscretePressureFunctionType>
class SchurkomplementOperator
{
    public:

        typedef SchurkomplementOperator <   A_SolverType,
                                            B_t_matrixType,
                                            CmatrixType,
                                            BmatrixType,
                                            MmatrixType,
                                            DiscreteVelocityFunctionType,
                                            DiscretePressureFunctionType>
                ThisType;

        SchurkomplementOperator( A_SolverType& a_solver,
                                const B_t_matrixType& b_t_mat,
                                const CmatrixType& c_mat,
                                const BmatrixType& b_mat,
                                const MmatrixType& m_mat,
                                const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& velocity_space,
                                const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space )
            : a_solver_(a_solver),
            b_t_mat_(b_t_mat),
            c_mat_(c_mat),
            b_mat_(b_mat),
            m_mat_(m_mat),
            tmp1 ( "tmp1", velocity_space ),
            tmp2 ( "tmp2", velocity_space )
        {}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            Logging::LogStream& dbg = Logger().Dbg();

            const bool solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );

            tmp1.clear();
//            Stuff::addScalarToFunc( tmp2, 0.1 );

			// ret = ( ( B_t * ( A^-1 * ( B * x ) ) ) + ( C * x ) )
//			Stuff::oneLinePrint( std::cout, tmp1 ) ;
            b_t_mat_.multOEM_t( x, tmp1.leakPointer() );
//            Stuff::oneLinePrint( std::cout, tmp1 ) ;
//            Stuff::printDoubleVec( std::cout, x, b_mat_.cols() );
            a_solver_.apply( tmp1, tmp2 );
            b_t_mat_.multOEM( tmp2.leakPointer(), ret );
            c_mat_.multOEMAdd( x, ret );
			//
        }


        ThisType& systemMatrix ()
        {
            return *this;
        }

    private:
        A_SolverType& a_solver_;
        const B_t_matrixType& b_t_mat_;
        const CmatrixType& c_mat_;
        const BmatrixType& b_mat_;
        const MmatrixType& m_mat_;
        mutable DiscreteVelocityFunctionType tmp1;
        mutable DiscreteVelocityFunctionType tmp2;
};

template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteSigmaFunctionType,
            class DiscreteVelocityFunctionType >
class A_SolverCaller {
    public:
        typedef MatrixA_Operator<   WMatType,
                                    MMatType,
                                    XMatType,
                                    YMatType,
                                    DiscreteSigmaFunctionType >
                A_OperatorType;

        typedef OEMBICGSTABOp< DiscreteVelocityFunctionType, A_OperatorType >
            CG_SolverType;
        typedef typename CG_SolverType::ReturnValueType
                SolverReturnType;

        A_SolverCaller( const WMatType& w_mat,
                const MMatType& m_mat,
                const XMatType& x_mat,
                const YMatType& y_mat,
                const typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType& sig_space,
                const double relLimit,
                const double absLimit,
                const bool verbose )
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            sig_tmp1( "sig_tmp1", sig_space ),
            sig_tmp2( "sig_tmp2", sig_space ),
            a_op_( w_mat, m_mat, x_mat, y_mat, sig_space ),
            cg_solver( a_op_,   relLimit,
                                absLimit,
                                2000, //inconsequential anyways
                                verbose )
        {}

        void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest, SolverReturnType& ret )
        {
            cg_solver.apply(arg,dest, ret);
        }

        void apply ( const DiscreteVelocityFunctionType& arg, DiscreteVelocityFunctionType& dest )
        {
            cg_solver.apply(arg,dest);
        }

        const A_OperatorType& getA_operator( )
        {
            return a_op_;
        }

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

#endif // INNERCG_HH_INCLUDED
