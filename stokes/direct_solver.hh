#ifndef DUNE_STOKES_DIRECT_SOLVER_HH
#define DUNE_STOKES_DIRECT_SOLVER_HH

#include <dune/stokes/solver_interface.hh>

namespace Dune {

	/** \brief Operator wrapping Matrix vector multiplication for
				matrix \f$ S :=  B_t * A^-1 * B + rhs3 \f$
				**/
	template <  class A_OperatorType,
				class B_t_matrixType,
				class CmatrixType,
				class BmatrixType,
				class DiscreteVelocityFunctionType ,
				class DiscretePressureFunctionType>
	class FullSytemOperator //: public OEMSolver::PreconditionInterface
	{
		public:

			typedef FullSytemOperator <   A_OperatorType,
												B_t_matrixType,
												CmatrixType,
												BmatrixType,
												DiscreteVelocityFunctionType,
												DiscretePressureFunctionType>
					ThisType;

			FullSytemOperator( A_OperatorType& a_operator,
									const B_t_matrixType& b_t_mat,
									const CmatrixType& c_mat,
									const BmatrixType& b_mat,
									const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& velocity_space,
									const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space )
				: a_operator_(a_operator),
				b_t_mat_(b_t_mat),
				c_mat_(c_mat),
				b_mat_(b_mat),
				tmp_velocity ( "tmp1", velocity_space ),
				tmp_pressure ( "tmp2", pressure_space ),
				do_bfg( Parameters().getParam( "do-bfg", true ) ),
				total_inner_iterations( 0 ),
				pressure_space_(pressure_space),
				velocity_space_(velocity_space),
				precond_( c_mat_.rows() + b_mat_.rows() )
			{
			}

//			double ddotOEM(const double*v, const double* w) const
//			{
//				DiscretePressureFunctionType V( "ddot V", pressure_space_, v );
//				DiscretePressureFunctionType W( "ddot W", pressure_space_, w );
//				return V.scalarProductDofs( W );
//			}

			template <class VECtype>
			void multOEM(const VECtype *x, VECtype * ret) const
			{
				const size_t numDofs_velocity = velocity_space_.size();
				a_operator_.multOEM( x, ret );
				b_mat_.multOEMAdd( x + numDofs_velocity, ret );
				b_t_mat_.multOEM( x, ret + numDofs_velocity );
				c_mat_.multOEMAdd( x + numDofs_velocity, ret + numDofs_velocity );
			}

		#ifdef USE_BFG_CG_SCHEME
			template <class VECtype>
			void multOEM(const VECtype *x, VECtype * ret, const IterationInfo& info ) const
			{
				multOEM(x,ret);
			}
		#endif

			ThisType& systemMatrix ()
			{
				return *this;
			}

			IdentityMatrix<CmatrixType>& preconditionMatrix()
			{
				return precond_;
			}

			bool hasPreconditionMatrix () const
			{
				return false;
			}

			bool rightPrecondition() const
			{
				return false;
			}

			template <class VecType>
			void precondition( const VecType* tmp, VecType* dest ) const
			{
				assert( false );
				precond_.multOEM( tmp, dest );
			}

			long getTotalInnerIterations()
			{
				return total_inner_iterations;
			}

		private:
			mutable A_OperatorType& a_operator_;
			const B_t_matrixType& b_t_mat_;
			const CmatrixType& c_mat_;
			const BmatrixType& b_mat_;
			mutable DiscreteVelocityFunctionType tmp_velocity;
			mutable DiscretePressureFunctionType tmp_pressure;
			bool do_bfg;
			mutable long total_inner_iterations;
			const typename DiscretePressureFunctionType::DiscreteFunctionSpaceType& pressure_space_;
			const typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType& velocity_space_;
			IdentityMatrix<CmatrixType> precond_;
	};


template < class StokesPassImp >
class DirectKrylovSolver
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
	DirectKrylovSolver()
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
				YmatrixObjectType& Omatrix,
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

		logInfo << "Begin DirectKrylovSolver " << std::endl;

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
		YmatrixType& o_mat      = Omatrix.matrix();
		B_t_matrixType& b_t_mat = Ematrix.matrix(); //! renamed
		CmatrixType& c_mat      = Rmatrix.matrix(); //! renamed
		BmatrixType& b_mat      = Zmatrix.matrix(); //! renamed
		WmatrixType& w_mat      = Wmatrix.matrix();

		c_mat.scale( -1 ); //since B_t = -E

		DiscretePressureFunctionType& g_func = rhs3;
		g_func *= ( -1 ); //since G = -H_3

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

		DomainType rhs_wrapper( arg.space(), f_func, g_func );

		typedef CombinedDiscreteFunction< DomainType >
			CombinedDiscreteFunctionType;
		CombinedDiscreteFunctionType combined_dest( dest );
		CombinedDiscreteFunctionType combined_rhs( rhs_wrapper );

		typedef MatrixA_Operator<   WmatrixType,
									MmatrixType,
									XmatrixType,
									YmatrixType,
									DiscreteSigmaFunctionType,
									DiscreteVelocityFunctionType>
				A_OperatorType;
		A_OperatorType a_operator( w_mat, m_inv_mat, x_mat, y_mat, o_mat, rhs1.space() , velocity.space() );

		typedef FullSytemOperator< A_OperatorType,
									B_t_matrixType,
									CmatrixType,
									BmatrixType,
									DiscreteVelocityFunctionType,
									DiscretePressureFunctionType >
			FullSytemOperatorType;

		FullSytemOperatorType fullsystem_operator( a_operator, b_t_mat, c_mat, b_mat,
												   velocity.space(), pressure.space() );

		typedef SOLVER_NAMESPACE::OUTER_CG_SOLVERTYPE< CombinedDiscreteFunctionType, FullSytemOperatorType >
				KrylovSolverType;
		KrylovSolverType kr( fullsystem_operator , relLimit, absLimit, 2000, solverVerbosity );

		kr.apply( combined_rhs, combined_dest );
		combined_dest.copyBack( dest );

		return SaddlepointInverseOperatorInfo();

	}

  };

} //namespace Dune

#endif // DUNE_STOKES_DIRECT_SOLVER_HH
