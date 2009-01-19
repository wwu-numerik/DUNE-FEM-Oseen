#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
//! using the memprovider from FEM currently results in assertion failed
#undef USE_MEMPROVIDER
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
			//const EllipticInverseOperatorType& aufSolver
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

        Stuff::DiagonalMult( m_inv_mat, rhs1 ); //calc m_inv * H_1 "in-place"
        DiscreteVelocityFunctionType f_func( "f_func", velocity.space() );
        f_func.clear();
        x_mat.apply( rhs1, f_func );
        f_func += rhs2;

        typedef MatrixA_Operator< WmatrixType, MmatrixType, XmatrixType, YmatrixType, DiscreteSigmaFunctionType >
            A_OperatorType;
        A_OperatorType a_op( w_mat, m_inv_mat, x_mat, y_mat, rhs1 );


        typedef SchurkomplementOperator<A_OperatorType,
                                        XmatrixType,
                                        MmatrixType,
                                        YmatrixType,
                                        B_t_matrixType,
                                        CmatrixType,
                                        BmatrixType,
                                        WmatrixType,
                                        DiscreteVelocityFunctionType,
                                        DiscretePressureFunctionType,
                                        DiscreteSigmaFunctionType >
                Sk_Operator;

        typedef OEMCGOp< DiscreteVelocityFunctionType, Sk_Operator >
                Sk_Solver;

        Sk_Operator sk_op(  a_op, x_mat, m_inv_mat, y_mat, b_t_mat, c_mat, b_mat, w_mat,
                            f_func, g_func, arg.discreteVelocity(), dest.discreteVelocity(), rhs1 );
        Sk_Solver sk_solver( sk_op, 0.001, 0.01, 2000, 1 );
        sk_solver( f_func, dest.discreteVelocity() );

//
//        logInfo << "- build global matrices - " << std::endl;
//        typedef SparseRowMatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteVelocityFunctionSpaceType >
//            AmatrixType;
//        AmatrixType Amatrix( velocitySpace_, velocitySpace_ );
//        Amatrix.reserve();
//
//        XmatrixType neg_X_Minv_mat( velocitySpace_, sigmaSpace_ );
//        neg_X_Minv_mat.reserve();
//        Xmatrix.matrix().multiply( Mmatrix.matrix(), neg_X_Minv_mat.matrix() );
//        logInfo << "-    1te feritg- " << std::endl;
//        neg_X_Minv_mat.matrix().scale( -1 );
//
//        neg_X_Minv_mat.matrix().multiply( Wmatrix.matrix(), Amatrix.matrix() );
//        Ymatrix.matrix().add( Amatrix.matrix() );
////            Amatrix = Ymatrix;
//
//
//        RHSType Fmat ( velocitySpace_.size(), 1, 1 );
//        neg_X_Minv_mat.matrix().scale ( mu );
//        neg_X_Minv_mat.matrix().multiply( H1rhs, Fmat );
////            H2rhs.add( Fmat );
//        //Fmat = H2rhs;
//
//        H3rhs.scale( -1 );
//        logInfo << "- build global matrices - done" << std::endl;





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
    // reference to operator which should be inverted
//    const StokesPassType & pass_;

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
    const StokesPassType& pass_;

  };

} // end namespace Dune

#endif

