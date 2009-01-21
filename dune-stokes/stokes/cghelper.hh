#ifndef INNERCG_HH_INCLUDED
#define INNERCG_HH_INCLUDED
//! using the memprovider from FEM currently results in assertion failed
#define USE_MEMPROVIDER
#include <dune/fem/solver/oemsolver.hh>
#include "../src/logging.hh"
#include "../src/stuff.hh"

namespace Dune {



template <  class WMatType,
            class MMatType,
            class XMatType,
            class YMatType,
            class DiscreteVelocityFunctionType >
class MatrixA_Operator {

    public:

        typedef MatrixA_Operator<   WMatType,
                        MMatType,
                        XMatType,
                        YMatType,
                        DiscreteVelocityFunctionType >
                    ThisType;

        MatrixA_Operator ( const WMatType& w_mat,
                const MMatType& m_mat,
                const XMatType& x_mat,
                const YMatType& y_mat,
                const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& space )
            :  w_mat_(w_mat),
            m_mat_(m_mat),
            x_mat_(x_mat),
            y_mat_(y_mat),
            tmp1( "tmp1", space ),
            tmp2( "tmp2", space )
        {}

        ~MatrixA_Operator()
        {}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            tmp1.clear();
            tmp2.clear();
            Logging::LogStream& dbg = Logger().Dbg();

//            dbg << "\n Auf x: " ;
//            Stuff::printDoubleVec( dbg, x, w_mat_.cols() );
//            dbg << std::endl ;
            w_mat_.multOEM( x, tmp1.leakPointer() );
            m_mat_.multOEM( tmp1.leakPointer(), tmp2.leakPointer() );//Stuff:DiagmUlt
            tmp2 *= ( -1 );
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
        mutable DiscreteVelocityFunctionType tmp1;
        mutable DiscreteVelocityFunctionType tmp2;
};


template <  class A_OperatorType,
            class B_t_matrixType,
            class CmatrixType,
            class BmatrixType,
            class DiscreteVelocityFunctionType ,
            class DiscretePressureFunctionType>
class SchurkomplementOperator
{
    public:

        typedef SchurkomplementOperator <   A_OperatorType,
                                            B_t_matrixType,
                                            CmatrixType,
                                            BmatrixType,
                                            DiscreteVelocityFunctionType,
                                            DiscretePressureFunctionType>
                ThisType;

        SchurkomplementOperator( A_OperatorType& a_op,
                                const B_t_matrixType& b_t_mat,
                                const CmatrixType& c_mat,
                                const BmatrixType& b_mat,
                                const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& velocity_space,
                                const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space )
            : a_op_(a_op),
            b_t_mat_(b_t_mat),
            c_mat_(c_mat),
            b_mat_(b_mat),
            tmp1 ( "tmp1", velocity_space ),
            tmp2 ( "tmp2", velocity_space )
        {}

        template <class VECtype>
        void multOEM(const VECtype *x, VECtype * ret) const
        {
            Logging::LogStream& dbg = Logger().Dbg();
            typedef OEMBICGSTABOp< DiscreteVelocityFunctionType, A_OperatorType >
                AufSolver;
            AufSolver auf_solver( a_op_, redEps, absLimit, 2000, 0 );

            tmp1.clear();
            tmp2.clear();
//            dbg << "\n SK x: " ;
//            Stuff::printDoubleVec( dbg, x, b_mat_.cols() );
//            dbg << std::endl ;
            b_mat_.multOEM( x, tmp1.leakPointer() );
//            dbg.Log (&DiscreteVelocityFunctionType::print , tmp1 ) ;
//            Stuff::oneLinePrint( dbg, tmp1 );
//            dbg << "begin: inner Ax=f" << std::endl;
            auf_solver( tmp1, tmp2 );

//            dbg << "end: inner Ax=f" << std::endl;
            b_t_mat_.multOEM( tmp2.leakPointer(), ret );
            c_mat_.multOEMAdd( x, ret );

        }


        ThisType& systemMatrix ()
        {
            return *this;
        }

    private:
         A_OperatorType& a_op_;
        const B_t_matrixType& b_t_mat_;
        const CmatrixType& c_mat_;
        const BmatrixType& b_mat_;
        mutable DiscreteVelocityFunctionType tmp1;
        mutable DiscreteVelocityFunctionType tmp2;
};

}

#endif // INNERCG_HH_INCLUDED
