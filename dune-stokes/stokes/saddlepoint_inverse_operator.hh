#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */


#define CG_SOLVERTYPE OEMCGOp
#ifndef CG_SOLVERTYPE
    #define CG_SOLVERTYPE OEMBICGSTABOp
#endif

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stokes/cghelper.hh>

#include "../src/logging.hh"
#include "../src/parametercontainer.hh"
#include "../src/stuff.hh" //DiagonalMult


namespace Dune {
  //!CG Verfahren fuer Sattelpunkt Problem
  //
  /**   \brief
        \tparam StokesPassImp discrete function types etc. get extracted from this
   */

    template < class StokesPassImp >
    class SaddlepointInverseOperator
{
  private:

    typedef StokesPassImp StokesPassType;

    typedef typename StokesPassType::DiscreteStokesFunctionWrapperType
        DiscreteStokesFunctionWrapperType;

    typedef typename StokesPassType::DomainType
        DomainType;

    typedef typename StokesPassType::RangeType
        RangeType;

    typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
        PressureDiscreteFunctionType;
    typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
        VelocityDiscreteFunctionType;


  public:
    /** \todo Please doc me!
     * \brief Constructor:
	*
	  **/
    SaddlepointInverseOperator()
	{}

    /** takes raw matrices from pass
    */
    template <  class XmatrixObjectType,
                class MmatrixObjectType,
                class YmatrixObjectType,
                class EmatrixObjectType,
                class RmatrixObjectType,
                class ZmatrixObjectType,
                class WmatrixObjectType,
                class DiscreteSigmaFunctionType,
                class DiscreteVelocityFunctionType,
                class DiscretePressureFunctionType  >
    void solve( const DomainType& arg,
                RangeType& dest,
                XmatrixObjectType& Xmatrix,
                MmatrixObjectType& Mmatrix,
                YmatrixObjectType& Ymatrix,
                EmatrixObjectType& Ematrix,
                RmatrixObjectType& Rmatrix,
                ZmatrixObjectType& Zmatrix,
                WmatrixObjectType& Wmatrix,
                const DiscreteSigmaFunctionType& rhs1,
                DiscreteVelocityFunctionType& rhs2,
                DiscretePressureFunctionType& rhs3 ) const
    {
        Logging::LogStream& logDebug = Logger().Dbg();
        Logging::LogStream& logError = Logger().Err();
        Logging::LogStream& logInfo = Logger().Info();

		// relative min. error at which cg-solvers will abort
        const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
		// aboslute min. error at which cg-solvers will abort
        const double absLimit = Parameters().getParam( "absLimit", 1e-3 );
        const bool solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );

        logInfo << "Begin SaddlePointInverseOperator " << std::endl;

        logDebug.Resume();
        //get some refs for more readability
        PressureDiscreteFunctionType& pressure = dest.discretePressure();
        VelocityDiscreteFunctionType& velocity = dest.discreteVelocity();

        typedef typename  XmatrixObjectType::MatrixType
            XmatrixType;
        typedef typename  MmatrixObjectType::MatrixType
            MmatrixType;
        typedef typename  YmatrixObjectType::MatrixType
            YmatrixType;
        typedef typename  EmatrixObjectType::MatrixType
            B_t_matrixType;                             //! renamed
        typedef typename  RmatrixObjectType::MatrixType
            CmatrixType;                                //! renamed
        typedef typename  ZmatrixObjectType::MatrixType
            BmatrixType;                                //! renamed
        typedef typename  WmatrixObjectType::MatrixType
            WmatrixType;

        XmatrixType& x_mat      = Xmatrix.matrix();
        MmatrixType& m_inv_mat  = Mmatrix.matrix();
        YmatrixType& y_mat      = Ymatrix.matrix();
        B_t_matrixType& b_t_mat = Ematrix.matrix(); //! renamed
        CmatrixType& c_mat      = Rmatrix.matrix(); //! renamed
        BmatrixType& b_mat      = Zmatrix.matrix(); //! renamed
        WmatrixType& w_mat      = Wmatrix.matrix();

        b_t_mat.scale( -1 ); //since B_t = -E

        DiscretePressureFunctionType& g_func = rhs3;
        g_func *= ( -1 ); //since G = -H_3

        logInfo << " \n\tbegin calc new_f,f_func " << std::endl;
        //Stuff::DiagonalMult( m_inv_mat, rhs1 ); //calc m_inv * H_1 "in-place"
        DiscreteSigmaFunctionType m_tmp ( "m_tom", rhs1.space() );
        DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
		f_func.clear();
		m_tmp.clear();

		// f_func = ( ( -1 * ( X * ( M_inv * rhs1 ) ) ) + rhs2 )
        m_inv_mat.apply( rhs1, m_tmp );
        x_mat.apply( m_tmp, f_func );
        f_func *= -1;
        f_func += rhs2;
		//

        typedef A_SolverCaller< WmatrixType,
                                MmatrixType,
                                XmatrixType,
                                YmatrixType,
                                DiscreteSigmaFunctionType,
                                DiscreteVelocityFunctionType >
            A_Solver;
        A_Solver a_solver( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space() );

//		typedef CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
//		F_Solver;
//        F_Solver f_solver( a_op, relLimit, absLimit, 2000, solverVerbosity );
        DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
		DiscretePressureFunctionType new_f ( "new_f", g_func.space() );
        tmp_f.clear();
		new_f.clear();

		// new_f := ( B * A^-1 * f_func ) + g_func
        a_solver.apply( f_func, tmp_f );
		b_t_mat.apply( tmp_f, new_f );
        new_f -= g_func;
		//
        logInfo << " \n\tend calc new_f,f_func " << std::endl;

        typedef SchurkomplementOperator<    A_Solver,
                                            B_t_matrixType,
                                            CmatrixType,
                                            BmatrixType,
                                            MmatrixType,
                                            DiscreteVelocityFunctionType,
                                            DiscretePressureFunctionType >
                Sk_Operator;

        typedef CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
                Sk_Solver;

        logInfo << " \n\tbegin S*p=new_f " << std::endl;
        Sk_Operator sk_op(  a_solver, b_t_mat, c_mat, b_mat, m_inv_mat,
                            velocity.space(), pressure.space() );
        Sk_Solver sk_solver( sk_op, relLimit, absLimit, 2000, solverVerbosity );
        pressure.clear();
//        Stuff::addScalarToFunc( pressure, 0.1 );
		// p = S^-1 * new_f = ( B_t * A^-1 * B + rhs3 )^-1 * new_f
		sk_solver( new_f, pressure );
		//
		logInfo << "\n\tend  S*p=new_f" << std::endl;

        DiscreteVelocityFunctionType Bp_temp ( "Bp_temp", velocity.space() );
        Bp_temp.clear();
		// velocity = A^-1 * ( ( -1 * ( B * pressure ) ) + f_func )
		b_mat.apply( pressure, Bp_temp );
        Bp_temp *= ( -1 );
        Bp_temp += f_func;
        a_solver.apply ( Bp_temp, velocity );

        logInfo << "\nEnd SaddlePointInverseOperator " << std::endl;
    }

  };

} // end namespace Dune

#endif

