#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */

// OEMBICGSQOp will NOT compile
#ifndef INNER_CG_SOLVERTYPE
    #define INNER_CG_SOLVERTYPE OEMCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
    #define OUTER_CG_SOLVERTYPE OEMCGOp
#endif


#ifdef USE_BFG_CG_SCHEME
    #include <utility>
    //< iteration no , < absLimit, residuum > >
    typedef std::pair<int,std::pair<double,double> >
        IterationInfo;

#endif

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stokes/cghelper.hh>

#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>

#include <cmath>



namespace Dune {

struct SaddlepointInverseOperatorInfo {
    int iterations_inner_avg;
    int iterations_inner_min;
    int iterations_inner_max;
    int iterations_outer_total;
    double max_inner_accuracy;
};

  //!CG Verfahren fuer Sattelpunkt Problem
  //
  /**   \brief
        \tparam StokesPassImp discrete function types etc. get extracted from this
   */

template < class StokesPassImp >
class NestedCgSaddlepointInverseOperator
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
    NestedCgSaddlepointInverseOperator()
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
    SaddlepointInverseOperatorInfo solve( const DomainType& arg,
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

        if ( Parameters().getParam( "disableSolver", false ) ) {
            logInfo.Resume();
            logInfo << "solving disabled via parameter file" << std::endl;
            return SaddlepointInverseOperatorInfo();
        }

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


        typedef A_SolverCaller< WmatrixType,
                                MmatrixType,
                                XmatrixType,
                                YmatrixType,
                                DiscreteSigmaFunctionType,
                                DiscreteVelocityFunctionType >
            A_Solver;
        typedef typename A_Solver::ReturnValueType
                A_SolverReturnType;
        A_Solver a_solver( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space(), relLimit, absLimit, solverVerbosity );

        DiscreteVelocityFunctionType tmp_f ( "tmp_f", f_func.space() );
		DiscretePressureFunctionType new_f ( "new_f", g_func.space() );
        tmp_f.clear();
		new_f.clear();

		// new_f := ( B * A^-1 * f_func ) + g_func
		A_SolverReturnType a_ret;
		a_solver.apply( f_func, tmp_f, a_ret );
#ifdef ADAPTIVE_SOLVER
		if ( isnan(a_ret.second) ) {
		    logInfo << "\n\t\t NaNs detected, lowering error tolerance" << std::endl;
            int max_adaptions = Parameters().getParam( "max_adaptions", 8 );
            int adapt_step = 1;
            double a_relLimit = relLimit;
            double a_absLimit = absLimit;
            while ( true ) {
                if ( adapt_step > max_adaptions ) {
                    logInfo << "\n\t\t max adaption depth reached, aborting" << std::endl;
                    break;
                }
                a_relLimit /= 10.0;
                a_absLimit /= 10.0;
                logInfo << "\n\t\t\t trying with relLimit " << a_relLimit
                        << " and absLimit " << a_absLimit << std::endl;

                A_Solver a_solver_adapt( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space(), relLimit, absLimit, solverVerbosity );
                a_solver_adapt.apply( f_func, tmp_f, a_ret );

                if ( !isnan(a_ret.second) ) {
                    logInfo << "\n\t\t adaption produced NaN-free solution" << std::endl;
                    break;
                }
                adapt_step++;
            }
		}
#endif
		b_t_mat.apply( tmp_f, new_f );
        new_f -= g_func;

        logInfo << " \n\tend calc new_f,f_func " << std::endl;

        typedef SchurkomplementOperator<    A_Solver,
                                            B_t_matrixType,
                                            CmatrixType,
                                            BmatrixType,
                                            MmatrixType,
                                            DiscreteVelocityFunctionType,
                                            DiscretePressureFunctionType >
                Sk_Operator;

        typedef OUTER_CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
                Sk_Solver;

        typedef typename Sk_Solver::ReturnValueType
                SolverReturnType;

        logInfo << " \n\tbegin S*p=new_f " << std::endl;
        Sk_Operator sk_op(  a_solver, b_t_mat, c_mat, b_mat, m_inv_mat,
                            velocity.space(), pressure.space() );
        Sk_Solver sk_solver( sk_op, relLimit, absLimit, 2000, solverVerbosity );
        pressure.clear();

		// p = S^-1 * new_f = ( B_t * A^-1 * B + rhs3 )^-1 * new_f
		SolverReturnType ret;
		sk_solver.apply( new_f, pressure, ret );

		long total_inner = sk_op.getTotalInnerIterations();
		logInfo << "\n\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;
#ifdef ADAPTIVE_SOLVER
        if ( isnan(ret.second) ) {
            logInfo << "\n\t\t NaNs detected, lowering error tolerance" << std::endl;
            int max_adaptions = Parameters().getParam( "max_adaptions", 8 );
            int adapt_step = 1;
            double a_relLimit = relLimit;
            double a_absLimit = absLimit;
            while ( true ) {
                if ( adapt_step > max_adaptions ) {
                    logError << "\n\t\t max adaption depth reached, aborting" << std::endl;
                    break;
                }
                a_relLimit /= 10.0;
                a_absLimit /= 10.0;
                logInfo << "\n\t\t\t trying with relLimit " << a_relLimit
                        << " and absLimit " << a_absLimit << std::endl;

                A_Solver a_solver_adapt( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space(), relLimit, absLimit, solverVerbosity );
                Sk_Operator sk_op_adapt(  a_solver_adapt, b_t_mat, c_mat, b_mat, m_inv_mat,
                            velocity.space(), pressure.space() );
                Sk_Solver sk_solver_adapt( sk_op_adapt, a_relLimit, a_absLimit, 2000, solverVerbosity );
                pressure.clear();
                sk_solver_adapt.apply( new_f, pressure, ret );
                long total_inner = sk_op_adapt.getTotalInnerIterations();
                logInfo << "\n\t\t\t #avg inner iter | #outer iter: " << total_inner / (double)ret.first << " | " << ret.first << std::endl;
                if ( !isnan(ret.second) ) {
                    logInfo << "\n\t\t adaption produced NaN-free solution" << std::endl;
                    break;
                }
                adapt_step++;
            }
        }
#endif
		//
		logInfo << "\n\tend  S*p=new_f" << std::endl;

        DiscreteVelocityFunctionType Bp_temp ( "Bp_temp", velocity.space() );
        Bp_temp.clear();
		// velocity = A^-1 * ( ( -1 * ( B * pressure ) ) + f_func )
		b_mat.apply( pressure, Bp_temp );
        Bp_temp *= ( -1 );
        Bp_temp += f_func;
        a_solver.apply ( Bp_temp, velocity );

        logInfo << "\nEnd NestedCgSaddlePointInverseOperator " << std::endl;

        return SaddlepointInverseOperatorInfo();
    }

  };

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
    SaddlepointInverseOperatorInfo solve( const DomainType& arg,
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

        if ( Parameters().getParam( "disableSolver", false ) ) {
            logInfo.Resume();
            logInfo << "solving disabled via parameter file" << std::endl;
            return SaddlepointInverseOperatorInfo();
        }

		// relative min. error at which cg-solvers will abort
        const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
		// aboslute min. error at which cg-solvers will abort
        double outer_absLimit = Parameters().getParam( "absLimit", 1e-8 );
        const double inner_absLimit = Parameters().getParam( "inner_absLimit", 1e-8 );
        const int solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
        const int maxIter = Parameters().getParam( "maxIter", 500 );

#ifdef USE_BFG_CG_SCHEME
        const double tau = Parameters().getParam( "bfg-tau", 0.1 );
        const bool do_bfg = Parameters().getParam( "do-bfg", true );
#endif
        logInfo.Resume();
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

/*** making our matrices kuhnibert compatible ****/
        b_t_mat.scale( -1 ); //since B_t = -E
        w_mat.scale( m_inv_mat(0,0) );
        rhs1 *=  m_inv_mat(0,0);
        m_inv_mat.scale( 1 / m_inv_mat(0,0) );

        //transformation from StokesPass::buildMatrix
        VelocityDiscreteFunctionType v_tmp ( "v_tmp", velocity.space() );
        x_mat.apply( rhs1, v_tmp );
        rhs2 -= v_tmp;
/***********/

        typedef A_SolverCaller< WmatrixType,
                                MmatrixType,
                                XmatrixType,
                                YmatrixType,
                                DiscreteSigmaFunctionType,
                                DiscreteVelocityFunctionType >
            A_Solver;
        typedef typename A_Solver::ReturnValueType
            ReturnValueType;
#ifdef USE_BFG_CG_SCHEME
        //the bfg scheme uses the outer acc. as a base
        double current_inner_accuracy = do_bfg ? tau * outer_absLimit : inner_absLimit;
        double max_inner_accuracy = current_inner_accuracy;
#else
        double current_inner_accuracy = inner_absLimit;
#endif
        A_Solver a_solver( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space(), relLimit, current_inner_accuracy, solverVerbosity > 3 );
        ReturnValueType a_solver_info;
/*****************************************************************************************/

        int iteration = 0;
        int total_inner_iterations = 0;
        int min_inner_iterations = std::numeric_limits<int>::max();
        int max_inner_iterations = 0;
        const int max_adaptions = Parameters().getParam( "max_adaptions", 2 ) ;
        int current_adaption = 0;
        double delta; //norm of residuum
        double gamma=0;
        double rho;

        VelocityDiscreteFunctionType F( "f", velocity.space() );
        F.assign(rhs2);
        VelocityDiscreteFunctionType tmp1( "tmp1", velocity.space() );
        tmp1.clear();
        VelocityDiscreteFunctionType xi( "xi", velocity.space() );
        xi.clear();
        PressureDiscreteFunctionType tmp2( "tmp2", pressure.space() );
        PressureDiscreteFunctionType residuum( "residuum", pressure.space() );
        residuum.assign(rhs3);

        PressureDiscreteFunctionType d( "d", pressure.space() );
        PressureDiscreteFunctionType h( "h", pressure.space() );

        // u^0 = A^{-1} ( F - B * p^0 )
        b_mat.apply( pressure, tmp1 );
        F-=tmp1; // F ^= rhs2 - B * p
        a_solver.apply(F,velocity);

        // r^0 = G - B_t * u^0 + C * p^0
        b_t_mat.apply( velocity, tmp2 );
        residuum -= tmp2;
        tmp2.clear();
        c_mat.apply( pressure, tmp2 );
        residuum += tmp2;

        // d^0 = r^0
        d.assign( residuum );

        delta = residuum.scalarProductDofs( residuum );
        while( (delta > outer_absLimit ) && (iteration++ < maxIter ) ) {
            if ( iteration > 1 ) {
                // gamma_{m+1} = delta_{m+1} / delta_m
                gamma = delta / gamma;
                // d_{m+1} = r_{m+1} + gamma_m * d_m
                d *= gamma;
                d += residuum;
            }
            if ( iteration >=  maxIter && current_adaption < max_adaptions ) {
                current_adaption++;
                iteration = 2;//do not execute first step in next iter again
                outer_absLimit /= 0.01;
#ifdef USE_BFG_CG_SCHEME
                current_inner_accuracy /= 0.01;
                a_solver.setAbsoluteLimit( current_inner_accuracy );
#endif
                logInfo << "\n\t\t Outer CG solver reset, tolerance lowered" << std::endl;

            }

#ifdef USE_BFG_CG_SCHEME
                if ( do_bfg ) {
                    //the form from the precond. paper (does not work properly)
//                    current_inner_accuracy = tau * std::min( 1. , absLimit / std::min ( std::pow( delta, int(iteration) ), 1.0 ) );
                    //my form, works
                    current_inner_accuracy = tau * std::min( 1. , outer_absLimit / std::min ( delta , 1.0 ) );
                    a_solver.setAbsoluteLimit( current_inner_accuracy );
                    max_inner_accuracy = std::max( max_inner_accuracy, current_inner_accuracy );
                    if( solverVerbosity > 1 )
                        logInfo << "\t\t\t set inner limit to: " << current_inner_accuracy << "\n";
                }
#endif
            // xi = A^{-1} ( B * d )
            tmp1.clear();
            b_mat.apply( d, tmp1 );
            a_solver.apply( tmp1, xi, a_solver_info );

#ifdef USE_BFG_CG_SCHEME
            if( solverVerbosity > 1 )
                logInfo << "\t\t inner iterations: " << a_solver_info.first << std::endl;
            total_inner_iterations += a_solver_info.first;
            min_inner_iterations = std::min( min_inner_iterations, a_solver_info.first );
            max_inner_iterations = std::max( max_inner_iterations, a_solver_info.first );
#endif

            // h = B_t * xi  + C * d
            b_t_mat.apply( xi, h );
            tmp2.clear();
            c_mat.apply( d, tmp2 );
            h += tmp2;

            rho = delta / d.scalarProductDofs( h );

            // p_{m+1} = p_m - ( rho_m * d_m )
            pressure.addScaled( d, -rho );

            // u_{m+1} = u_m + ( rho_m * xi_m )
            velocity.addScaled( xi, +rho );

            // r_{m+1} = r_m - rho_m * h_m
            residuum.addScaled( h, -rho );

            //save old delta for new gamma calc in next iter
            gamma = delta;

            // d_{m+1} = < r_{m+1} ,r_{m+1} >
            delta = residuum.scalarProductDofs( residuum );

            if( solverVerbosity > 2 )
                logInfo << "\t" << iteration << " SPcg-Iterationen  " << iteration << " Residuum:" << delta << std::endl;
        }

        logInfo << "End SaddlePointInverseOperator " << std::endl;

        SaddlepointInverseOperatorInfo info; //left blank in case of no bfg
#ifdef USE_BFG_CG_SCHEME
        const double avg_inner_iterations = total_inner_iterations / (double)iteration;
        if( solverVerbosity > 0 )
            logInfo << "\n #avg inner iter | #outer iter: "
                    <<  avg_inner_iterations << " | " << iteration << std::endl;

        info.iterations_inner_avg = avg_inner_iterations;
        info.iterations_inner_min = min_inner_iterations;
        info.iterations_inner_max = max_inner_iterations;
        info.iterations_outer_total = iteration;
        info.max_inner_accuracy = max_inner_accuracy;
#endif
        return info;
    }

  };

} // end namespace Dune

#endif

