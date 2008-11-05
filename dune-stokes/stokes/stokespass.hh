/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include<dune/fem/operator/matrix/spmatrix.hh>


namespace Dune
{
    class StokesPass
    {

        public:
            typedef SparseRowMatrix<double> MatrixType;
            typedef MatrixType B_OperatorType;
            typedef MatrixType B_Transposed_OperatorType;
            typedef MatrixType C_OperatorType;

            StokesPass(){}

            B_OperatorType& Get_B_Operator() { return b_op_; }
            B_Transposed_OperatorType& Get_B_Transposed_Operator() { return b_transp_op_; }
            C_OperatorType& Get_C_Operator() { return c_op_; }

        private:
            //will prolly be generated on-the-fly, not stored as members
            B_OperatorType b_op_;
            B_Transposed_OperatorType b_transp_op_;
            C_OperatorType c_op_;

    };
}
#endif  // end of stokespass.hh
