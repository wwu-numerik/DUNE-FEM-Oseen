#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED
//! using the memprovider from FEM currently results in assertion failed
#undef USE_MEMPROVIDER
#include <dune/fem/solver/oemsolver.hh>

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
                const DiscreteSigmaFunctionType& sigma_dummy )
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            sigma_dummy_(sigma_dummy)
        {}

        ~MatrixA_Operator()
        {}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            DiscreteSigmaFunctionType tmp ( "tmp", sigma_dummy_.space() );
            DiscreteSigmaFunctionType tmp2 ( "tmp", sigma_dummy_.space() );
            w_mat_.multOEM( x, tmp.leakPointer() );
            m_mat_.multOEM( tmp.leakPointer(), tmp2.leakPointer() );
            x_mat_.multOEM( tmp2.leakPointer(), ret );

            y_mat_.multOEMAdd( x, ret );
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
        const DiscreteSigmaFunctionType& sigma_dummy_;
};


template <  class XmatrixType,
            class M_inv_matrixType,
            class YmatrixType,
            class B_t_matrixType,
            class CmatrixType,
            class BmatrixType,
            class WmatrixType,
            class FFunctype,
            class GFunctype,
            class DiscreteSigmaFunctionType >
class SchurkomplementOperator
{
    private:
        typedef MatrixA_Operator< WmatrixType, M_inv_matrixType, XmatrixType, YmatrixType, DiscreteSigmaFunctionType >
            A_OperatorType;

    public:

        typedef SchurkomplementOperator <  XmatrixType,
                                         M_inv_matrixType,
                                         YmatrixType,
                                         B_t_matrixType,
                                         CmatrixType,
                                         BmatrixType,
                                         WmatrixType,
                                         FFunctype,
                                         GFunctype,
                                         DiscreteSigmaFunctionType >
                ThisType;

        SchurkomplementOperator(  const XmatrixType& x_mat,
                                const M_inv_matrixType& m_inv_mat,
                                const YmatrixType& y_mat,
                                const B_t_matrixType& b_t_mat,
                                const CmatrixType& c_mat,
                                const BmatrixType& b_mat,
                                const WmatrixType& w_mat,
                                const FFunctype& f_func,
                                const GFunctype& g_func,
                                const FFunctype& arg,
                                FFunctype& dest,
                                const DiscreteSigmaFunctionType& sigma_dummy )
            : x_mat_(x_mat),
            m_inv_mat_(m_inv_mat),
            y_mat_(y_mat),
            b_t_mat_(b_t_mat),
            c_mat_(c_mat),
            b_mat_(b_mat),
            w_mat_(w_mat),
            f_func_(f_func),
            g_func_(g_func),
            arg_(arg),
            dest_(dest),
            sigma_dummy_(sigma_dummy)
        {}

//        template <  class ArgDescreteFunctionType,
//                    class DestDescreteFunctionType >
//        void apply (    const ArgDescreteFunctionType& arg,
//                        DestDescreteFunctionType& dest )
//        {
//
//
//
//
//        }

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            A_OperatorType a_op( w_mat_, m_inv_mat_, x_mat_, y_mat_, sigma_dummy_ );

            typedef OEMCGOp< FFunctype, A_OperatorType >
                AufSolver;
            AufSolver auf_solver( a_op, 0.001, 0.01, 2000, 1 );

            FFunctype tmp ( "tmp", arg_.space() );
            FFunctype tmp2 ( "tmp2", arg_.space() );
            b_mat_.multOEM( x, tmp.leakPointer() );
            auf_solver( tmp, tmp2 );
            b_t_mat_.multOEM( tmp2.leakPointer(), ret );
            c_mat_.multOEMAdd( x, ret );

        }


        ThisType& systemMatrix ()
        {
            return *this;
        }

    private:
        const XmatrixType& x_mat_;
        const M_inv_matrixType& m_inv_mat_;
        const YmatrixType& y_mat_;
        const B_t_matrixType& b_t_mat_;
        const CmatrixType& c_mat_;
        const BmatrixType& b_mat_;
        const WmatrixType& w_mat_;
        const FFunctype& f_func_;
        const GFunctype& g_func_;
        const FFunctype& arg_;
        FFunctype& dest_;
        const DiscreteSigmaFunctionType& sigma_dummy_;
};

}

#endif // INNERCG_HH_INCLUDED
