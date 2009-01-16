#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED

#include <dune/fem/solver/oemsolver.hh>

namespace Dune {



template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteSigmaFunctionType >
class MultA {

    public:
        MultA ( const WMatType& w_mat,
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

        ~MultA()
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
            class DiscreteSigmaFunctionType  >
class SchurkomplementSolver
{
    private:
        typedef MultA< WmatrixType, M_inv_matrixType, XmatrixType, YmatrixType, DiscreteSigmaFunctionType >
            MultAType;

    public:
        SchurkomplementSolver(  const XmatrixType& x_mat,
                                const M_inv_matrixType& m_inv_mat,
                                const YmatrixType& y_mat,
                                const B_t_matrixType& b_t_mat,
                                const CmatrixType& c_mat,
                                const BmatrixType& b_mat,
                                const WmatrixType& w_mat,
                                const FFunctype& f_func,
                                const GFunctype& g_func,
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
            sigma_dummy_(sigma_dummy)
        {}

        template <  class ArgDescreteFunctionType,
                    class DestDescreteFunctionType >
        void apply (    const ArgDescreteFunctionType& arg,
                        DestDescreteFunctionType& dest )
        {
            MultAType a_op( w_mat_, m_inv_mat_, x_mat_, y_mat_, sigma_dummy_ );
            a_op.multOEM( arg.discreteVelocity().leakPointer(), dest.discreteVelocity().leakPointer() );

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
        const DiscreteSigmaFunctionType& sigma_dummy_;
};

}

#endif // INNERCG_HH_INCLUDED
