/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include<dune/fem/operator/matrix/spmatrix.hh>


namespace Dune
{
    template <class VelocityDiscreteFunctionImp,
            class PressureDiscreteFunctionImp, class DiscreteModelImp, class PreviousPassImp, int PassID = 0 >
    class StokesPass : public LocalPass < DiscreteModelImp, PreviousPassImp, PassID >
    {

        public:

            //! tempalte repetions etc
            typedef LocalPass < DiscreteModelImp, PreviousPassImp, PassID > BaseType;
            typedef PreviousPassImp PreviousPassType;
            typedef VelocityDiscreteFunctionImp VelocityDiscreteFunctionType;
            typedef PressureDiscreteFunctionImp PressureDiscreteFunctionType;
            typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFuntionSpaceType;

            typedef SparseRowMatrix<double> MatrixType;
            typedef MatrixType B_OperatorType;
            typedef MatrixType B_Transposed_OperatorType;
            typedef MatrixType C_OperatorType;

            //typedefs for interaface complaince, not definite
            typedef typename BaseType::ArgumentType ArgumentType;
            typedef typename BaseType::DestinationType DestinationType;
            typedef typename BaseType::Entity EntityType;

            StokesPass(PreviousPassType& prev_pass, const DiscreteFuntionSpaceType& disc_space)
                : BaseType( prev_pass, disc_space )
                //rhs( "pass_rhs"
            {}

            const B_OperatorType& Get_B_Operator() const { return b_op_; }
            const B_Transposed_OperatorType& Get_B_Transposed_Operator() const { return b_transp_op_; }
            const C_OperatorType& Get_C_Operator() const { return c_op_; }

            MatrixType& systemMatrix(){}

            const VelocityDiscreteFunctionType& rhs1() const {}

            virtual void prepare( const ArgumentType& arg,
                                    DestinationType& dest) const
            {}

            virtual void finalize(const ArgumentType& arg,
                                    DestinationType& dest) const
            {}

            virtual void applyLocal(EntityType& en) const
            {}





        private:
            //will prolly be generated on-the-fly, not stored as members
            B_OperatorType b_op_;
            B_Transposed_OperatorType b_transp_op_;
            C_OperatorType c_op_;

            //more dummies
           // VelocityDiscreteFunctionType& rhs;

    };
}
#endif  // end of stokespass.hh
