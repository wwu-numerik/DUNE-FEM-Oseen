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

        if ( Parameters().getParam( "disableSolver", false ) ) {
            logInfo.Resume();
            logInfo << "solving disabled via parameter file" << std::endl;
            return;
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
		//
//#ifndef NLOG
//        Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
//        Stuff::printDoubleVectorMatlabStyle( f_func.leakPointer(), f_func.size(), "f_func", matlabLogStream );
//#endif

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

//		typedef CG_SOLVERTYPE< DiscreteVelocityFunctionType, A_OperatorType >
//		F_Solver;
//        F_Solver f_solver( a_op, relLimit, absLimit, 2000, solverVerbosity );
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
//#ifndef NLOG
//        Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
//        Stuff::printDoubleVectorMatlabStyle( new_f.leakPointer(), new_f.size(), "new_f", matlabLogStream );
//#endif
        if ( Parameters().getParam( "solution-print", true ) ) {
            Stuff::oneLinePrint( logDebug, new_f );
        }
        logInfo << " \n\tend calc new_f,f_func " << std::endl;

        typedef SchurkomplementOperator<    A_Solver,
                                            B_t_matrixType,
                                            CmatrixType,
                                            BmatrixType,
                                            MmatrixType,
                                            DiscreteVelocityFunctionType,
                                            DiscretePressureFunctionType >
                Sk_Operator;

        #ifdef USE_BFG_CG_SCHEME
            typedef OUTER_CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator, true >
                    Sk_Solver;
        #else
            typedef OUTER_CG_SOLVERTYPE< DiscretePressureFunctionType, Sk_Operator >
                    Sk_Solver;
        #endif

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
//#ifndef NLOG
//        Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
//        Stuff::printDoubleVectorMatlabStyle( pressure.leakPointer(), pressure.size(), "pressure", matlabLogStream );
//#endif

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
                    logInfo << "\n\t\t max adaption depth reached, aborting" << std::endl;
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

        logInfo << "\nEnd SaddlePointInverseOperator " << std::endl;
    }

  };

template < class StokesPassImp >
class MirkoSaddlepointInverseOperator
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
    MirkoSaddlepointInverseOperator()
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
            return;
        }

		// relative min. error at which cg-solvers will abort
        const double relLimit = Parameters().getParam( "relLimit", 1e-4 );
		// aboslute min. error at which cg-solvers will abort
        const double absLimit = Parameters().getParam( "absLimit", 1e-8 );
        const bool solverVerbosity = Parameters().getParam( "solverVerbosity", 0 );
        const unsigned int maxIter = Parameters().getParam( "maxIter", 500 );

        logInfo.Resume();
        logInfo << "Begin MirkoSaddlePointInverseOperator " << std::endl;

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

/*** making our matrices mirko compatible ****/
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

        A_Solver a_solver( w_mat, m_inv_mat, x_mat, y_mat, rhs1.space(), relLimit, absLimit*0.001, solverVerbosity );

/*****************************************************************************************/

        unsigned int count = 0;
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
        while( (delta > absLimit ) && (count++ < maxIter ) ) {
            if ( count > 1 ) {
                // gamma_{m+1} = delta_{m+1} / delta_m
                gamma = delta / gamma;
                // d_{m+1} = r_{m+1} + gamma_m * d_m
                d *= gamma;
                d += residuum;
            }

            // xi = A^{-1} ( B * d )
            tmp1.clear();
            b_mat.apply( d, tmp1 );
            a_solver.apply( tmp1, xi );

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

            if( solverVerbosity > 0 )
                std::cerr << count << " SPcg-Iterationen  " << count << " Residuum:" << delta << "\n";
        }

    }

  };

} // end namespace Dune

#endif

