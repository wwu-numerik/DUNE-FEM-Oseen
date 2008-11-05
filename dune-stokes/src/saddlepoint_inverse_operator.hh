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

#include "logging.hh"

namespace Dune {
  //!CG Verfahren fuer Sattelpunkt Problem
  //!siehe NavSt Skript Siebert,Seite 36, Alg 3.34
  /**   \brief Inversion operator using CG algorithm
        \tparam EllipticInverseOperatorType Operator type used to do the inner A^-1 inversion
   */

  template <class VelocityDiscreteFunctionType,
            class PressureDiscreteFunctionType,
            class StokesPassImp,
            class EllipticInverseOperatorType>
    class SaddlepointInverseOperator
        : public Operator< typename PressureDiscreteFunctionType::DomainFieldType,
                            typename PressureDiscreteFunctionType::RangeFieldType,
                            PressureDiscreteFunctionType,PressureDiscreteFunctionType>
{
  private:

    typedef typename PressureDiscreteFunctionType::FunctionSpaceType PressureSpaceType;
    typedef typename VelocityDiscreteFunctionType::FunctionSpaceType VelocitySpaceType;

    typedef StokesPassImp StokesPassType;

    typedef typename StokesPassType::MatrixType MatrixType;

    typedef typename StokesPassType::B_OperatorType B_OperatorType;
    typedef typename StokesPassType::B_Transposed_OperatorType B_Transposed_OperatorType;
    typedef typename StokesPassType::C_OperatorType C_OperatorType;

  public:
    /** \todo Please doc me!
     * \brief Constructor:
    * aufSolver is the InverseOperator for Solving the elliptic Problem A^-1
    * rhs1 is  stored as member,no better idea
	  **/
    SaddlepointInverseOperator( const StokesPassType& op,
			 double redEps,
			 double absLimit,
			 int maxIter,
			 int verbose,
			 const EllipticInverseOperatorType& aufSolver,
			 PressureSpaceType& pressurespc,
			 VelocitySpaceType& spc
			 //const PressureDiscreteFunctionType rhs2
			 )
      : pass_(op),
        error_reduction_per_step_ ( redEps ),
        epsilon_ ( absLimit ),
        max_iterations_ (maxIter ),
        verbosity_ ( verbose ),
        aufSolver_(aufSolver),
        b_op_(pass_.Get_B_Operator()),
        bT_op_(pass_.Get_B_Transposed_Operator()),
        c_op_(pass_.Get_C_Operator()),
        rhs1_(pass_.rhs1()),
        pressure_space_(pressurespc),
        velocity_space_(spc),
        velocity_("velocity",velocity_space_)
    {

    }

    /** \todo Please doc me! */
    virtual void operator()(const PressureDiscreteFunctionType& arg,
                            PressureDiscreteFunctionType& dest ) const
    {
        Logging::LogStream& logDebug = Logger().Dbg();
        Logging::LogStream& logError = Logger().Err();
        Logging::LogStream& logInfo = Logger().Info();
        logInfo << "Begin SaddlePointInverseOperator " << std::endl;

        typedef typename VelocityDiscreteFunctionType::DiscreteFunctionSpaceType
            FunctionSpaceType;
        typedef typename FunctionSpaceType::RangeFieldType
            RangeFieldType;
        typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType
            PressureFunctionSpaceType;

        int iterations_count = 0;
        RangeFieldType spa=0, residuum, q, quad;

        VelocityDiscreteFunctionType f("f",velocity_space_);
        f.assign(rhs1_);
        VelocityDiscreteFunctionType u("u",velocity_space_);
        u.clear();
        VelocityDiscreteFunctionType tmp1("tmp1",velocity_space_);//part one of our pseedup loop ;-)
        tmp1.clear();
        VelocityDiscreteFunctionType xi("xi",velocity_space_);
        xi.clear();

        PressureDiscreteFunctionType tmp2("tmp2",pressure_space_);//part two
        //  tmp2.assign(arg);
        PressureDiscreteFunctionType r("r",pressure_space_);
        r.assign(arg);
        //p<->d
        PressureDiscreteFunctionType p("p",pressure_space_);
        PressureDiscreteFunctionType h("h",pressure_space_);
        PressureDiscreteFunctionType g("g",pressure_space_);

        //   double bla=r.scalarProductDofs(r);
        //  std::cout<<"arg dot arg="<< bla <<"\n";


        //(3.95a)
        //  std::cout<<"SPoperator dest (3.95a)\n";
        //       dest.print(std::cout);

        b_op_.apply(dest,tmp1);
        //bT_op_.apply_t(dest,tmp1);
        //  std::cout<<"SPoperator bop(dest) (3.95a)\n";
        //       tmp1.print(std::cout);

        f-=tmp1; // f = - ( B * dest )
        //  std::cout<<"*********************************************************\n";
        // std::cout<<"SPoperator f (3.95a)\n";
        //  f.print(std::cout);
        //      std::cout<<"*********************************************************\n";

        aufSolver_(f,u);
        //  std::cout<<"*********************************************************\n";
        //    std::cout<<"SPoperator Result auf solver (3.95a)\n";
        //   u.print(std::cout);
        // std::cout<<"*********************************************************\n";
        //(3.95b)

        //  bT_op_(u,tmp2);
        bT_op_.apply(u,tmp2);
        //b_op_.apply_t(u,tmp2);
        r-=tmp2;


        tmp2.clear();

        // c_op_(dest,tmp2);
        c_op_.apply(dest,tmp2);


        r+=tmp2;

        //   std::cout<<"r=\n";

        //(3.95c)
        p.assign(r);


        //(3.95d)
        residuum = r.scalarProductDofs( r );
        // std::cout<<"rho ="<<residuum<<"\n";
        //	       int foo;
        //   std::cin>>foo;
        //   assert(foo==1);
        while((residuum > epsilon_) && (iterations_count++ < max_iterations_))
        {
            // fall ab der zweiten iteration *************

            if(iterations_count > 1)
            {  //(3.95l)
                const RangeFieldType e = residuum / spa;

                //(3.95m)
                p *= e;
                p += r;
            }

            // grund - iterations - schritt **************
            //(3.95e)
            tmp1.clear();
            // b_op_(p,tmp1);
            b_op_.apply(p,tmp1);
            //bT_op_.apply_t(p,tmp1);

            aufSolver_(tmp1,xi);


            //(3.95f)

            //b_op_.apply_t(xi,h);
            bT_op_.apply(xi,h);
            tmp2.clear();


            //   c_op_(p,tmp2);

            c_op_.apply(p,tmp2);

            h+=tmp2;


            //(3.95g)
            quad = p.scalarProductDofs( h );
            //   std::cout<<"quad ="<<quad<<"\n";
            q    = residuum / quad;

            //  std::cout<<"q ="<<q<<"\n";
            //(3.95h)
            dest.addScaled( p, -q );
            //(3.95i)
            u.addScaled(xi,+q);
            //(3.95j)
            r.addScaled( h, -q );
            // std::cout<<"r=\n";
            // r.print(std::cout);
            //\ro_m merken
            spa = residuum;

            // residuum neu berechnen *********************
            //(3.95K)
            //\ro_{m+1} berechnen
            residuum = r.scalarProductDofs( r );


            logInfo << " SPcg-Iteration Nr. " << std::setw(3) << iterations_count
                    << " Residuum:" << residuum << "        \r";
        }


        velocity_.assign(u);
        logInfo << "End SaddlePointInverseOperator " << std::endl;
    }

    VelocityDiscreteFunctionType& velocity() { return velocity_; }

  private:
    // reference to operator which should be inverted
    const StokesPassType & pass_;

    // reduce error each step by
    double error_reduction_per_step_;

    // minial error to reach
    typename VelocityDiscreteFunctionType::RangeFieldType epsilon_;

    // number of maximal iterations
    int max_iterations_;

    // level of output
    int verbosity_ ;

    //the CGSolver for A^-1
    const EllipticInverseOperatorType& aufSolver_;

    B_OperatorType& b_op_;
    B_Transposed_OperatorType& bT_op_;
    C_OperatorType& c_op_;

    mutable VelocityDiscreteFunctionType& rhs1_;
    PressureSpaceType& pressure_space_;
    VelocitySpaceType& velocity_space_;
    mutable VelocityDiscreteFunctionType velocity_;
    // PressureDiscreteFunctionType rhs2_;
  };

} // end namespace Dune

#endif

