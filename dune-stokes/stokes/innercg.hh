#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED

namespace Dune {

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
