#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */


#define CG_SOLVERTYPE OEMBICGSTABOp
#ifndef CG_SOLVERTYPE
    #define CG_SOLVERTYPE OEMCGOp
#endif

const double redEps = 1e-6;
const double absLimit = 1e-6;

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
//! using the memprovider from FEM currently results in assertion failed
//#undef USE_MEMPROVIDER
#include <dune/stokes/cghelper.hh>

#include "../src/logging.hh"
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

//    typedef typename StokesPassType::MatrixType
//        MatrixType;
    typedef typename StokesPassType::DiscreteStokesFunctionWrapperType
        DiscreteStokesFunctionWrapperType;

    typedef typename StokesPassType::DomainType
        DomainType;

    typedef typename StokesPassType::RangeType
        RangeType;

//    typedef typename StokesPassType::PressureDiscreteFunctionType::DiscreteFunctionSpaceType
//        PressureDiscreteSpaceType;
//    typedef typename StokesPassType::VelocityDiscreteFunctionType::DiscreteFunctionSpaceType
//        VelocityDiscreteSpaceType;
    typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
        PressureDiscreteFunctionType;
    typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
        VelocityDiscreteFunctionType;


  public:
    /** \todo Please doc me!
     * \brief Constructor:
    * aufSolver is the InverseOperator for Solving the elliptic Problem A^-1
    * velocity,rhs1 is  stored as member,no better idea
	  **/
    SaddlepointInverseOperator(
            const StokesPassType& pass,
            double redEps,
			double absLimit,
			int maxIter,
			int verbose
			)
      : pass_( pass ),
        error_reduction_per_step_ ( redEps ),
        epsilon_ ( absLimit ),
        max_iterations_ (maxIter ),
        verbosity_ ( verbose )
        //aufSolver_(aufSolver)
    {

    }

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
                DiscreteSigmaFunctionType& rhs1,
                DiscreteVelocityFunctionType& rhs2,
                DiscretePressureFunctionType& rhs3 ) const
    {
        Logging::LogStream& logDebug = Logger().Dbg();
        Logging::LogStream& logError = Logger().Err();
        Logging::LogStream& logInfo = Logger().Info();
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
        c_mat.scale( -1 ); //since C = -R

        DiscretePressureFunctionType& g_func = rhs3;
        g_func *= ( -1 ); //since G = -H_3

        //Stuff::DiagonalMult( m_inv_mat, rhs1 ); //calc m_inv * H_1 "in-place"
        DiscreteSigmaFunctionType m_tmp ( "m_tom", rhs1.space() );
        m_inv_mat.apply( rhs1, m_tmp );


        DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
        f_func.clear();
        x_mat.apply( m_tmp, f_func );
        f_func *= -1;
        f_func += rhs2;
        Stuff::oneLinePrint( logDebug, f_func );
        typedef MatrixA_Operator< WmatrixType, MmatrixType, XmatrixType, YmatrixType, DiscreteVelocityFunctionType >
            A_OperatorType;
        A_OperatorType a_op( w_mat, m_inv_mat, x_mat, y_mat, velocity.space() );


        typedef CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
                F_Solver;
        logInfo << " \nstart f solver " << std::endl;
        F_Solver f_solver( a_op, redEps, absLimit, 2000, 0 );

        // new_f := B * A^-1 * f_func + g_func
        DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
        tmp_f.clear();
        f_solver( f_func, tmp_f );
        logInfo << " \nend f solver " << std::endl;

        DiscretePressureFunctionType new_f ( "new_f", g_func.space() );
        b_t_mat.apply( tmp_f, new_f );
        new_f -= g_func;


        typedef SchurkomplementOperator<A_OperatorType,
                                        B_t_matrixType,
                                        CmatrixType,
                                        BmatrixType,
                                        DiscreteVelocityFunctionType,
                                        DiscretePressureFunctionType >
                Sk_Operator;

        typedef OEMBICGSTABOp< DiscretePressureFunctionType, Sk_Operator >
                Sk_Solver;
        logInfo << " \nbegin SK solver " << std::endl;
        Sk_Operator sk_op(  a_op, b_t_mat, c_mat, b_mat,
                            velocity.space(), pressure.space() );
        Sk_Solver sk_solver( sk_op, redEps, absLimit, 2000, 0 );
        pressure.assign( new_f );
        Stuff::addScalarToFunc( pressure, 0.0001 );
        sk_solver( new_f, pressure );
        logInfo << " \nend SK solver " << std::endl;

        //schurkomplement Operator ist: B2 * A_inv * B1 + C

        //schritt 1
        //löse S * u_0 = F - B1 * p_0
      //  typedef InnerCG<

        //schritt 2
        //setzte:
        //  r_0     = -B2 * u_0 + C * p_0 + G
        //  d_0     = r_0
        //  delta_0 = ( r_0, r_0 )

        //schritt 3:
        //wenn d_0 == 0 dann GOTO end

        //schritt 4:
        //löse S * X_m = B1 * d_m

        //schritt 5:
        //setze:
        //  h_m         = B2 * X_m + C * p_m
        //  rho_m       = delta_m /  (h_m,d_m)
        //  p_m+1       = p_m - rho_m + d_m
        //  u_m+1       = u_m + rho_m * X_m
        //  r_m+1       = r_m - rho_m * h_m
        //  delta_m+1   = (r_m+1,r_m+1)

        //schritt 6:
        //wenn delta_m+1 == 0 dann GOTO end

        //schritt 7:
        //setze:
        //  gamma_m     = delta_m+1 / delta_m
        //  d_m+1       = r_m - gamma_m * d_m

        //schritt 8:
        //GOTO schritt 4

        logInfo << "End SaddlePointInverseOperator " << std::endl;
    }

//  where do i need this?
//    VelocityDiscreteFunctionType& velocity() { return velocity_; }

  private:

    const StokesPassType& pass_;
    // reduce error each step by
    double error_reduction_per_step_;
//
//     minial error to reach
    double epsilon_;
//
//     number of maximal iterations
    int max_iterations_;
//
//     level of output
    int verbosity_ ;
//
//    the CGSolver for A^-1
   // const EllipticInverseOperatorType& aufSolver_;
//


  };

} // end namespace Dune

#endif

