#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED

namespace Dune {

template <  class XmatrixType,
            class M_inv_matrixType,
            class YmatrixType,
            class B_t_matrixType,
            class CmatrixType,
            class BmatrixType,
            class WmatrixType,
            class FFunctype,
            class GFunctype >
class SchurkomplementSolver
{
    public:
        SchurkomplementSolver(  const XmatrixType& x_mat,
                                const M_inv_matrixType& m_inv_mat,
                                const YmatrixType& y_mat,
                                const B_t_matrixType& b_t_mat,
                                const CmatrixType& c_mat,
                                const BmatrixType& b_mat,
                                const WmatrixType& w_mat,
                                const FFunctype& f_func,
                                const GFunctype& g_func )
            : x_mat_(x_mat),
            m_inv_mat_(m_inv_mat),
            y_mat_(y_mat),
            b_t_mat_(b_t_mat),
            c_mat_(c_mat),
            b_mat_(b_mat),
            w_mat_(w_mat),
            f_func_(f_func),
            g_func_(g_func)
        {}

        template <  class ArgDescreteFunctionType,
                    class DestDescreteFunctionType >
        void apply (    const ArgDescreteFunctionType& arg,
                        DestDescreteFunctionType& dest )
        {}

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
};


template < class BMatType, class AMatType, class B_t_MatType, class CMatType >
class InnerCG {

    public:
        InnerCG (   const BMatType& b_mat,
                    const AMatType& a_mat,
                    const B_t_MatType& b_t_mat,
                    const CMatType& c_mat )
        {

        }

        ~InnerCG() {}

        template < class DiscreteFunctionType >
        void apply( const DiscreteFunctionType& arg,
                    const DiscreteFunctionType& dest )
        {

        }

    private:
        const BMatType& b_mat_;
        const AMatType& a_mat_;
        const B_t_MatType& b_t_mat_;
        const CMatType& c_mat_;

};

}

#endif // INNERCG_HH_INCLUDED
