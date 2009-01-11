#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_original.hh
 */

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include<dune/fem/operator/matrix/spmatrix.hh>

//using istl results in compile errors atm /rene
//#define USE_ISTL

#if defined(USE_ISTL)
    #include<dune/fem/operator/matrix/istlmatrix.hh>
    #include "matrixoperator.hh"
#endif

#include "../src/logging.hh"

namespace Dune {
  //!CG Verfahren fuer Sattelpunkt Problem
  //!siehe NavSt Skript Siebert,Seite 36, Alg 3.34
  /**   \brief Inversion operator using CG algorithm
        \tparam EllipticInverseOperatorType Operator type used to do the inner A^-1 inversion
   */

  template < class StokesPassImp >
            //class EllipticInverseOperatorType>
    class SaddlepointInverseOperator
//        : public Operator< typename StokesPassImp :: DiscreteStokesFunctionWrapperType::DomainFieldType,
//                           typename  StokesPassImp ::DiscreteStokesFunctionWrapperType::RangeFieldType,
//                            typename StokesPassImp ::DiscreteStokesFunctionWrapperType, typename StokesPassImp :: DiscreteStokesFunctionWrapperType>
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
//    typedef typename StokesPassType::PressureDiscreteFunctionType
//        PressureDiscreteFunctionType;
//    typedef typename StokesPassType::VelocityDiscreteFunctionType
//        VelocityDiscreteFunctionType;


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


    template <  class AMatrixType,
                class B1MatrixType,
                class B2MatrixType,
                class CMatrixType,
                class FMatrixType,
                class GMatrixType >
    void solve( const DomainType& arg,
                RangeType& dest,
                AMatrixType& mat_A,
                B1MatrixType& mat_B1,
                B2MatrixType& mat_B2,
                CMatrixType& mat_C,
                FMatrixType& mat_F,
                GMatrixType& mat_G ) const
    {
        Logging::LogStream& logDebug = Logger().Dbg();
        Logging::LogStream& logError = Logger().Err();
        Logging::LogStream& logInfo = Logger().Info();
        logInfo << "Begin SaddlePointInverseOperator " << std::endl;

        //schurkomplement Operator ist: B2 * A_inv * B1 + C

        //schritt 1
        //löse S * u_0 = F - B1 * p_0

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

