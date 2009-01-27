#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */


//#define CG_SOLVERTYPE OEMGMRESOp
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
            double redEps,
			double absLimit,
			int maxIter,
			int verbose
			)
      : error_reduction_per_step_ ( redEps ),
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
                const DiscreteSigmaFunctionType& rhs1,
                DiscreteVelocityFunctionType& rhs2,
                DiscretePressureFunctionType& rhs3 ) const
    {
        Logging::LogStream& logDebug = Logger().Dbg();
        Logging::LogStream& logError = Logger().Err();
        Logging::LogStream& logInfo = Logger().Info();

        const double redEps = Parameters().getParam( "redEps", 1e-4 );
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
//        c_mat.scale(  ); //since C = -R

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

        typedef MatrixA_Operator<   WmatrixType,
                                    MmatrixType,
                                    XmatrixType,
                                    YmatrixType,
                                    DiscreteVelocityFunctionType,
                                    DiscreteSigmaFunctionType >
            A_OperatorType;
        A_OperatorType a_op( w_mat, m_inv_mat, x_mat, y_mat, velocity.space(), rhs1.space() );

        typedef CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
                F_Solver;
        logInfo << " \nstart f solver " << std::endl;
        F_Solver f_solver( a_op, redEps, absLimit, 2000, solverVerbosity );

        // new_f := B * A^-1 * f_func + g_func
        DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
        tmp_f.clear();

        f_solver( f_func, tmp_f );
        logInfo << " \nend f solver " << std::endl;

        DiscretePressureFunctionType new_f ( "new_f", g_func.space() );
        b_t_mat.apply( tmp_f, new_f );
        new_f -= g_func;

        double delta = std::numeric_limits<double>::max();
        double tau = Parameters().getParam( "tau", 0.1 );
        int i = 0;
        int max_iter = Parameters().getParam( "maxSKIterations", 10 );

        typedef SchurkomplementOperator<A_OperatorType,
                                        B_t_matrixType,
                                        CmatrixType,
                                        BmatrixType,
                                        DiscreteVelocityFunctionType,
                                        DiscretePressureFunctionType >
                Sk_Operator;

        typedef CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
                Sk_Solver;

        Sk_Operator sk_op(  a_op, b_t_mat, c_mat, b_mat,
                            velocity.space(), pressure.space() );
        Sk_Solver sk_solver( sk_op, redEps, absLimit, 2000, solverVerbosity+41 );
        pressure.clear();
//        Stuff::addScalarToFunc( pressure, 1 );
        if ( max_iter < 0 ){
            logInfo << "cg solver" << std::endl;
            sk_solver( new_f, pressure );
            logInfo << "\nend cg solver" << std::endl;
        }

        typedef CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
                U_Solver;

        typedef CG_SOLVERTYPE< DiscretePressureFunctionType, CmatrixType >
                C_Solver;

        DiscreteVelocityFunctionType Bp_temp ( "Bp_temp", velocity.space() );

        DiscretePressureFunctionType pressure_tmp ( "press_tmp", pressure.space() );
        DiscretePressureFunctionType pressure_tmp2 ( "press_tmp2", pressure.space() );
        DiscretePressureFunctionType pressure_tmp3 ( "press_tmp3", pressure.space() );


        DiscretePressureFunctionType pressure_last( "p_k-1", pressure.space() );
        DiscreteVelocityFunctionType velocity_last( "u_k-1", velocity.space() );
        DiscreteVelocityFunctionType velocity_tmp( "u_k-1", velocity.space() );

        pressure_last.assign( pressure );


        while ( i < max_iter ) {
            logInfo << " \n SK solver Iteration: " << i << std::endl;
            pressure_last.assign( pressure );
            pressure.clear();
            pressure_tmp.clear();
            velocity_last.assign( velocity );
            velocity_last.clear();

            //p_k+1 = S * p_k - C*p + C^-1 ( f_new  + G ) + G
//            Stuff::printDoubleVec( logInfo, pressure_last.leakPointer(), b_mat.cols() );
            sk_op.multOEM( pressure_last.leakPointer(), pressure_tmp.leakPointer() ); //S*p_k
            pressure_tmp *= tau; // tau * S *p_k
            pressure_tmp += pressure_last; // tau * S *p_k + p_k


            pressure_tmp2.clear();
            pressure_tmp2.assign( new_f );
            pressure_tmp2 *= tau; // tau ( B_t * A^-1  * B - G)
            pressure_tmp2 *= -1; // tau ( - B_t * A^-1  * B + G)

            pressure_tmp += pressure_tmp2;

            pressure.assign( pressure_tmp );


            // u_k+1 = A^-1 * ( f - B_t * p_k )


            i++;

        }

        Bp_temp.clear();
        b_mat.apply( pressure_last, Bp_temp );
        Bp_temp *= ( -1 );
        Bp_temp += f_func;
        U_Solver u_solver( a_op, redEps, absLimit, 2000, solverVerbosity );
        u_solver ( Bp_temp, velocity );
//
//        logInfo << "End SaddlePointInverseOperator " << std::endl;
    }

//  where do i need this?
//    VelocityDiscreteFunctionType& velocity() { return velocity_; }

  private:

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

/**
        DiscretePressureFunctionType residuum_tmp ( "resi_tmp", pressure.space() );
        DiscretePressureFunctionType residuum ( "resi", pressure.space() );
        DiscretePressureFunctionType d_m ( "d_m", pressure.space() );
        DiscretePressureFunctionType X_m ( "X_m", pressure.space() );
       DiscreteVelocityFunctionType h_m( "h_m", velocity.space() );

//  //schritt 2
//        //setzte:
//        //  r_0     = -B2 * u_0 + C * p_0 + G
//        //  d_0     = r_0
//        //  delta_0 = ( r_0, r_0 ) residuum.clear();
//        residuum_tmp.clear();
//        c_mat.apply( pressure, residuum );
//        residuum += g_func;
//        b_t_mat.apply( velocity, residuum_tmp );
//        residuum -= residuum_tmp;
//        d_m.assign( residuum );
//        delta = residuum.scalarProductDofs( residuum );
//        logInfo << "SK:: residuum-norm: " << delta << std::endl ;
//
//        //schritt 3:
//        //wenn d_0 == 0 dann GOTO end
//
//
//        while ( i < max_iter ) {
//            //schritt 4:
//            //löse S * X_m = B1 * d_m
//            b_mat.apply( d_m, velocity_tmp );
//            sk_solver( velocity_tmp, X_m );
//
//            //schritt 5:
//            //setze:
//            //  h_m         = B2 * X_m + C * p_m
//            //  rho_m       = delta_m /  (h_m,d_m)
//            //  p_m+1       = p_m - rho_m + d_m
//            //  u_m+1       = u_m + rho_m * X_m
//            //  r_m+1       = r_m - rho_m * h_m
//            //  delta_m+1   = (r_m+1,r_m+1)
//            c_mat.apply( pressure, velocity_tmp );
//            b_t_mat.apply( X_m, h_m );
//            h_m += velocity_tmp;
//            double rho_m = delta * 100 ;/// h_m.scalarProductDofs( d_m );
//            pressure_tmp.assign( d_m );
//            pressure_tmp *= rho_m;
//            pressure -= pressure_tmp;
//
//            velocity_tmp.assign( h_m );
//            velocity_tmp *= rho_m;
//            velocity += velocity_tmp;
//
//            pressure_tmp.assign( X_m );
//            pressure_tmp *= rho_m;
//            residuum -= pressure_tmp;
//            double delta_next = residuum.scalarProductDofs( residuum );
//            //schritt 6:
//            //wenn delta_m+1 == 0 dann GOTO end
//            if ( delta_next < absLimit )
//                break;
//
//            //schritt 7:
//            //setze:
//            //  gamma_m     = delta_m+1 / delta_m
//            //  d_m+1       = r_m - gamma_m * d_m
//            double gamma = delta_next / delta;
//            delta = delta_next;
//            d_m *= -1 * gamma;
//            d_m += residuum;
//            //schritt 8:
//            //GOTO schritt 4
//            i++;
//        }
*/
