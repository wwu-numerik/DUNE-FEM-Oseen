#ifndef DUNE_OSEEN_SOLVERS_SADDLE_POINT_HH
#define DUNE_OSEEN_SOLVERS_SADDLE_POINT_HH

#include <dune/fem/oseen/solver/solver_interface.hh>

namespace Dune {

	/**
		\brief Saddlepoint Solver
		The inner CG iteration is implemented by passing a custom Matrix operator to a given\n
		Dune solver. The outer iteration is a implementation of the CG algorithm as described in\n
		Kuhnibert
		Optionally the BFG scheme as described in YADDA is uesd to control the inner solver tolerance.
		/todo get references in doxygen right
	**/
	template < class OseenPassImp >
	class SaddlepointInverseOperator
	{
	  private:

		typedef OseenPassImp
			OseenPassType;

		typedef typename OseenPassType::Traits::DiscreteOseenFunctionWrapperType
			DiscreteOseenFunctionWrapperType;

		typedef typename OseenPassType::DomainType
			DomainType;

		typedef typename OseenPassType::RangeType
			RangeType;

		typedef typename DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType
			PressureDiscreteFunctionType;
		typedef typename DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType
			VelocityDiscreteFunctionType;


	  public:
		SaddlepointInverseOperator()
		{}

		/** takes raw matrices and right hand sides from pass as input, executes nested cg algorithm and outputs solution
		*/
		template <  class X_MatrixType,
					class M_invers_matrixType,
					class Y_MatrixType,
					class O_MatrixType,
					class E_MatrixType,
					class R_MatrixType,
					class Z_MatrixType,
					class W_MatrixType,
					class DiscreteSigmaFunctionType,
					class DiscreteVelocityFunctionType,
					class DiscretePressureFunctionType  >
		SaddlepointInverseOperatorInfo solve( const DomainType& /*arg*/,
					RangeType& dest,
					X_MatrixType& Xmatrix,
					M_invers_matrixType& Mmatrix,
					Y_MatrixType& Ymatrix,
					O_MatrixType& Omatrix,
					E_MatrixType& Ematrix,
					R_MatrixType& Rmatrix,
					Z_MatrixType& Zmatrix,
					W_MatrixType& Wmatrix,
					const DiscreteSigmaFunctionType& rhs1_orig,
					const DiscreteVelocityFunctionType& rhs2_orig,
					const DiscretePressureFunctionType& rhs3 ) const
		{

            auto& logDebug = DSC_LOG_DEBUG;
            auto& logInfo = DSC_LOG_INFO;


			// relative min. error at which cg-solvers will abort
			const double relLimit = DSC_CONFIG_GET( "relLimit", 1e-4 );
			// aboslute min. error at which cg-solvers will abort
			double outer_absLimit = DSC_CONFIG_GET( "absLimit", 1e-8 );
			const double inner_absLimit = DSC_CONFIG_GET( "inner_absLimit", 1e-8 );
			const int solverVerbosity = DSC_CONFIG_GET( "solverVerbosity", 0 );
			const int maxIter = DSC_CONFIG_GET( "maxIter", 500 );
			const bool use_velocity_reconstruct = DSC_CONFIG_GET( "use_velocity_reconstruct", true );

			const double tau = DSC_CONFIG_GET( "bfg-tau", 0.1 );
			const bool do_bfg = DSC_CONFIG_GET( "do-bfg", true );
			logInfo.resume();
			logInfo << "Begin SaddlePointInverseOperator " << std::endl;

			logDebug.resume();
			//get some refs for more readability
            auto& pressure = dest.discretePressure();
            auto& velocity = dest.discreteVelocity();

            const auto& x_mat      = Xmatrix;
            const auto& m_inv_mat  = Mmatrix;
            const auto& y_mat      = Ymatrix;
            const auto& o_mat      = Omatrix;
            const auto& b_t_mat = Ematrix; //! renamed
            const auto& c_mat      = Rmatrix; //! renamed
            const auto& b_mat      = Zmatrix; //! renamed
            const auto& w_mat      = Wmatrix;

	/*** making our matrices kuhnibert compatible ****/
			//rhs1 = M^{-1} * rhs1
			//B_t = -B_t
            const double m_scale = m_inv_mat.matrix()(0,0);
            b_t_mat.matrix().scale( -1 ); //since B_t = -E
			DiscreteSigmaFunctionType rhs1 = rhs1_orig;
			rhs1 *=  m_scale;

			//rhs2 = rhs2 - X * M^{-1} * rhs1
			VelocityDiscreteFunctionType v_tmp ( "v_tmp", velocity.space() );
			x_mat.apply( rhs1, v_tmp );
			DiscreteVelocityFunctionType rhs2 = rhs2_orig;
			rhs2 -= v_tmp;
	/***********/

			typedef InnerCGSolverWrapper< W_MatrixType,
									M_invers_matrixType,
									X_MatrixType,
									Y_MatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType >
				InnerCGSolverWrapperType;

			typedef typename InnerCGSolverWrapperType::ReturnValueType
				ReturnValueType;
			ReturnValueType a_solver_info;

			//the bfg scheme uses the outer acc. as a base
            double current_inner_accuracy = do_bfg ? tau * outer_absLimit : inner_absLimit;
            double max_inner_accuracy = current_inner_accuracy;

			InnerCGSolverWrapperType innerCGSolverWrapper( w_mat, m_inv_mat, x_mat, y_mat,
														   o_mat, rhs1.space(),rhs2.space(), relLimit,
														   current_inner_accuracy, solverVerbosity > 3 );

	/*****************************************************************************************/

			int iteration = 0;
			int total_inner_iterations = 0;
			int min_inner_iterations = std::numeric_limits<int>::max();
			int max_inner_iterations = 0;

			const int max_adaptions = DSC_CONFIG_GET( "max_adaptions", 2 ) ;
			int current_adaption = 0;
			double delta; //norm of residuum
			double gamma=0;
			double rho;

			VelocityDiscreteFunctionType F( "f", velocity.space() );
			F.assign(rhs2);//should be rhs2_orig ?!?
			VelocityDiscreteFunctionType tmp1( "tmp1", velocity.space() );
			tmp1.clear();
			VelocityDiscreteFunctionType xi( "xi", velocity.space() );
			xi.clear();
			PressureDiscreteFunctionType tmp2( "tmp2", pressure.space() );
			PressureDiscreteFunctionType residuum( "residuum", pressure.space() );
			residuum.assign(rhs3);

			PressureDiscreteFunctionType d( "d", pressure.space() );
			PressureDiscreteFunctionType h( "h", pressure.space() );

			// u^0 = A^{-1} ( F - B * p^0 ) (3.95a)
			b_mat.apply( pressure, tmp1 );
			F-=tmp1; // F = rhs2 - X * M^{-1} * rhs1 - B * p
			innerCGSolverWrapper.apply(F,velocity);

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
					current_inner_accuracy /= 0.01;
					innerCGSolverWrapper.setAbsoluteLimit( current_inner_accuracy );
					logInfo << "\n\t\t Outer CG solver reset, tolerance lowered" << std::endl;
				}

                if ( do_bfg ) {
                    //the form from the precond. paper
                    current_inner_accuracy = tau * std::min( 1. , outer_absLimit / std::min ( delta , 1.0 ) );
                    innerCGSolverWrapper.setAbsoluteLimit( current_inner_accuracy );
                    max_inner_accuracy = std::max( max_inner_accuracy, current_inner_accuracy );
                    if( solverVerbosity > 1 )
                        logInfo << "\t\t\t set inner limit to: " << current_inner_accuracy << "\n";
                }
				// xi = A^{-1} ( B * d )
				tmp1.clear();
				b_mat.apply( d, tmp1 );

				innerCGSolverWrapper.apply( tmp1, xi, a_solver_info );

				if( solverVerbosity > 1 )
					logInfo << "\t\t inner iterations: " << a_solver_info.first << std::endl;
				total_inner_iterations += a_solver_info.first;
				min_inner_iterations = std::min( min_inner_iterations, a_solver_info.first );
				max_inner_iterations = std::max( max_inner_iterations, a_solver_info.first );

				// h = B_t * xi  + C * d
				b_t_mat.apply( xi, h );
				tmp2.clear();
				c_mat.apply( d, tmp2 );
				h += tmp2;

				rho = delta / d.scalarProductDofs( h );

				// p_{m+1} = p_m - ( rho_m * d_m )
				pressure.axpy( -rho, d );
				if ( !use_velocity_reconstruct ) {
					// u_{m+1} = u_m + ( rho_m * xi_m )
					velocity.axpy( +rho, xi );
				}
				// r_{m+1} = r_m - rho_m * h_m
				residuum.axpy( -rho, h );

				//save old delta for new gamma calc in next iter
				gamma = delta;

				// d_{m+1} = < r_{m+1} ,r_{m+1} >
				delta = residuum.scalarProductDofs( residuum );

				if( solverVerbosity > 2 )
					logInfo << "\t" << iteration << " SPcg-Iterationen  " << iteration << " Residuum:" << delta << std::endl;
			}

			if ( use_velocity_reconstruct ) {
				// u^0 = A^{-1} ( F - B * p^0 )
				F.assign(rhs2);
				tmp1.clear();
				b_mat.apply( pressure, tmp1 );
				F-=tmp1; // F ^= rhs2 - B * p
				innerCGSolverWrapper.apply(F,velocity);
			}

			logInfo << "End SaddlePointInverseOperator " << std::endl;
			SaddlepointInverseOperatorInfo info; //left blank in case of no bfg
			const double avg_inner_iterations = total_inner_iterations / (double)iteration;
			if( solverVerbosity > 0 )
				logInfo << "\n #avg inner iter | #outer iter: "
						<<  avg_inner_iterations << " | " << iteration << std::endl;
			info.iterations_inner_avg = avg_inner_iterations;
			info.iterations_inner_min = min_inner_iterations;
			info.iterations_inner_max = max_inner_iterations;
			info.iterations_outer_total = iteration;
			info.max_inner_accuracy = max_inner_accuracy;
			return info;
		} //end SaddlepointInverseOperator::solve

	  };//end class SaddlepointInverseOperator


} //end namespace Dune

#endif // DUNE_OSEEN_SOLVERS_SADDLE_POINT_HH

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

