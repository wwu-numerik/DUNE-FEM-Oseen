#ifndef DUNE_SPINVERSE_OPERATORS_HH
#define DUNE_SPINVERSE_OPERATORS_HH

/** \file
    \brief SPinverseoperator_symm.hh
 */

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include<dune/fem/operator/matrix/spmatrix.hh>
#include "matrixoperator.hh"

namespace Dune {
  //!CG Verfahren fuer Sattelpunkt Problem 
  //!siehe NavSt Skript Siebert,Seite 36, Alg 3.34
  /** \brief Inversion operator using CG algorithm
   */
 
//   typedef SPDGoperator < ScalDiscreteFunctionType ,VecDiscreteFunctionSpaceType > OperatorType; 
  



  template <class DiscreteFunctionType,class PressureDiscreteFunctionType,class OperatorType,class InverseOperatorType>
  class SPCGInverseOperator : public Operator<
    typename PressureDiscreteFunctionType::DomainFieldType,
    typename PressureDiscreteFunctionType::RangeFieldType,
    PressureDiscreteFunctionType,PressureDiscreteFunctionType> 
  {


    typedef SparseRowMatrix<double> MatrixType;
    //    typedef OperatorType SPOPType;
    
    typedef typename PressureDiscreteFunctionType::FunctionSpaceType PressureSpaceType;	   
    typedef typename DiscreteFunctionType::FunctionSpaceType VeloSpaceType;

    typedef OperatorType MappingType;
    /*
    typedef MatrixOperator<MatrixType,PressureDiscreteFunctionType,DiscreteFunctionType> BOPType;
    typedef MatrixOperator<MatrixType,DiscreteFunctionType ,PressureDiscreteFunctionType> BTOPType;
    typedef MatrixOperator<MatrixType,PressureDiscreteFunctionType,PressureDiscreteFunctionType> COPType;
    */
  
//     typedef typename MappingType::PressureGMType BOPType;
//     typedef typename MappingType::PressureDMType BTOPType;
//     typedef typename MappingType::PressureSMType COPType;

    /**************** for non istl version*/
    typedef MatrixType BOPType;
    //  typedef MatrixType BTOPType;
    typedef MatrixType COPType;
    /*****************************************/

    // typedef typename OperatorType::BOPType BOPType;
//     typedef typename OperatorType::BTOPType BTOPType;
//     typedef typename OperatorType::COPType COPType;
  public:
    /** \todo Please doc me! */
    //!Constructor:
    //!aufSolver is the InverseOperator for Solving the elliptic Problem A^-1
    //!rhs1 is  stored as member,no better idea 
    SPCGInverseOperator( const MappingType& op,
			 double redEps,
			 double absLimit,
			 int maxIter,
			 int verbose,
			 const InverseOperatorType& aufSolver,
			 PressureSpaceType& pressurespc,
			 VeloSpaceType& spc
			 //const PressureDiscreteFunctionType rhs2
			 ) 
      : op_(op), _redEps ( redEps ), epsilon_ ( absLimit ) , 
        maxIter_ (maxIter ) , _verbose ( verbose ),aufSolver_(aufSolver), bop_(op_.getBOP()),
	//	btop_(op_.getBTOP()),
	cop_(op_.getCOP()),
	rhs1_(op_.rhs1()),
	pressurespc_(pressurespc),
	spc_(spc),
	velocity_("velocity",spc_)
    {
     
    }
    /** \todo Please doc me! */
    virtual void operator()(const PressureDiscreteFunctionType& arg, 
                            PressureDiscreteFunctionType& dest ) const 
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType Field;
      typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureFunctionSpaceType;
      int count = 0;
      Field spa=0, spn, q, quad;
     
      DiscreteFunctionType f("f",spc_);
      f.assign(rhs1_);
      DiscreteFunctionType u("u",spc_);
      u.clear();
      DiscreteFunctionType tmp1("tmp1",spc_);
      tmp1.clear();
      DiscreteFunctionType xi("xi",spc_);
      xi.clear();
      PressureDiscreteFunctionType tmp2("tmp2",pressurespc_);
      //  tmp2.assign(arg);
      PressureDiscreteFunctionType r("r",pressurespc_);
      r.assign(arg);

      
      //p<->d
      PressureDiscreteFunctionType p("p",pressurespc_);
      PressureDiscreteFunctionType h("h",pressurespc_);
      PressureDiscreteFunctionType g("g",pressurespc_);  
      
   


      //(3.95a)
     //  std::cout<<"SPoperator dest (3.95a)\n";
      //       dest.print(std::cout);

      //bop_(dest,tmp1);
      bop_.apply(dest,tmp1);
   
      //  std::cout<<"SPoperator bop(dest) (3.95a)\n";
    //       tmp1.print(std::cout);
      
      f-=tmp1;
     //  std::cout<<"*********************************************************\n";
     // std::cout<<"SPoperator f (3.95a)\n";
    //  f.print(std::cout);
//      std::cout<<"*********************************************************\n";
   
      aufSolver_(f,u);
  //     int foo;
//       std::cin>>foo;
//       assert(foo==1);
      //  std::cout<<"*********************************************************\n";
 //    std::cout<<"SPoperator Result auf solver (3.95a)\n";
  //   u.print(std::cout);
  // std::cout<<"*********************************************************\n";
       //(3.95b)
   
      // btop_(u,tmp2);
      //        btop_.apply(u,tmp2); 
      bop_.apply_t(u,tmp2);
      r-=tmp2;
	  
      
      tmp2.clear();

      // cop_(dest,tmp2);
      cop_.apply(dest,tmp2);
	  

      r+=tmp2; 
     
	  //   std::cout<<"r=\n";

     
     
	  //(3.95c) 
	  p.assign(r);
    
      
	  //(3.95d)
	  spn = r.scalarProductDofs( r );
	  // std::cout<<"rho ="<<spn<<"\n"; 
	  //	       int foo;
	  //   std::cin>>foo;
	  //   assert(foo==1);
	  while((spn > epsilon_) && (count++ < maxIter_)  ) 
	    {
	      // fall ab der zweiten iteration *************
      
	      if(count > 1)
		{  //(3.95l)
		  const Field e = spn / spa;

		  //(3.95m)
		  p *= e;
		  p += r;
		}

	      // grund - iterations - schritt **************
	      //(3.95e)
	      tmp1.clear();
	      // bop_(p,tmp1);
	      bop_.apply(p,tmp1);
	      aufSolver_(tmp1,xi);
	  
// 	      int foo;
//       std::cin>>foo;
//       assert(foo==1);
	      //(3.95f)
	      // btop_(xi,h);
	      bop_.apply_t(xi,h);
	      //  btop_.apply(xi,h);
	      tmp2.clear();

	  
	      //   cop_(p,tmp2);

	      cop_.apply(p,tmp2);
	      
	      h+=tmp2;

	  
	      //(3.95g)
	      quad = p.scalarProductDofs( h );
	      //   std::cout<<"quad ="<<quad<<"\n"; 
	      q    = spn / quad; 
	      
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
	      spa = spn;
      
	      // residuum neu berechnen *********************
	      //(3.95K)
	      //\ro_{m+1} berechnen
	      spn = r.scalarProductDofs( r ); 
	    
	 
	      if(_verbose > 0)
		std::cerr << count << " SPcg-Iterationen  " << count << " Residuum:" << spn << "        \r";
	    }
      if(_verbose > 0)
	std::cerr << "\n";
    
           velocity_.assign(u);
  }
    
    DiscreteFunctionType& velocity()     {
      return velocity_;
    }

  private:
    // reference to operator which should be inverted 
    const OperatorType & op_;
  
    // reduce error each step by 
    double _redEps; 

    // minial error to reach 
    typename DiscreteFunctionType::RangeFieldType epsilon_;

    // number of maximal iterations
    int maxIter_;

    // level of output 
    int _verbose ;
    
    //the CGSolver for A^-1
    const InverseOperatorType& aufSolver_;
    
    BOPType& bop_;
    //  BTOPType& btop_;
    COPType& cop_;

    mutable DiscreteFunctionType& rhs1_;
    PressureSpaceType& pressurespc_;
    VeloSpaceType& spc_;
    mutable DiscreteFunctionType velocity_;
    // PressureDiscreteFunctionType rhs2_;
  };

//   /** \todo Please doc me! */
//   template <class DiscreteFunctionType, class OperatorType>
//   class CGInverseOp : public Operator<
//     typename DiscreteFunctionType::DomainFieldType,
//     typename DiscreteFunctionType::RangeFieldType,
//     DiscreteFunctionType,DiscreteFunctionType> 
//   {
//   public:
//     /** \todo Please doc me! */
//     CGInverseOp( OperatorType & op , double  redEps , double absLimit , int maxIter , int verbose ) : 
//       op_(op),
//       _redEps ( redEps ),
//       epsilon_ ( absLimit*absLimit ) , 
//       maxIter_ (maxIter ) ,
//       _verbose ( verbose ) , 
//       r_(0),
//       p_(0), 
//       h_(0) {} 

//     /** \todo Please doc me! */             
//     ~CGInverseOp()
//     {
//       if(p_) delete p_;
//       if(r_) delete r_;
//       if(h_) delete h_;
//     }

//     /** \todo Please doc me! */      
//     virtual void operator() (const DiscreteFunctionType& arg,
//                              DiscreteFunctionType& dest ) const 
//     {
//       typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
//       typedef typename FunctionSpaceType::RangeFieldType Field;

//       int count = 0;
//       Field spa=0, spn, q, quad;

//       if(!p_) p_ = new DiscreteFunctionType ( arg );
//       if(!r_) r_ = new DiscreteFunctionType ( arg );
//       if(!h_) h_ = new DiscreteFunctionType ( arg );
    
//       DiscreteFunctionType & r = *r_;
//       DiscreteFunctionType & p = *p_;
//       DiscreteFunctionType & h = *h_;

//       //op_.prepareGlobal(arg,dest);

//       op_( dest, h );

//       r.assign(h) ;
//       r -= arg;

//       p.assign(arg);
//       p -= h;

//       spn = r.scalarProductDofs( r );
   
//       while((spn > epsilon_) && (count++ < maxIter_)) 
//         {
//           // fall ab der zweiten iteration *************
      
//           if(count > 1)
//             { 
//               const Field e = spn / spa;
//               p *= e;
//               p -= r;
//             }

//           // grund - iterations - schritt **************
      
//           op_( p, h );
      
//           quad = p.scalarProductDofs( h );
//           q    = spn / quad;

//           dest.add( p, q );
//           r.add( h, q );

//           spa = spn;
      
//           // residuum neu berechnen *********************
      
//           spn = r.scalarProductDofs( r ); 
//           if(_verbose > 0)
//             std::cerr << count << " cg-Iterationen  " << count << " Residuum:" << spn << "        \r";
//         }
//       if(_verbose > 0)
//         std::cerr << "\n";
//       op_.finalizeGlobal();
//     }
  
//   private:
//     // no const reference, we make const later 
//     MatrixType &B_;
//     MatrixType &BT_;
//     MatrixType &C_;

//     // reduce error each step by 
//     double _redEps; 

//     // minial error to reach 
//     typename DiscreteFunctionType::RangeFieldType epsilon_;

//     // number of maximal iterations
//     int maxIter_;

//     // level of output 
//     int _verbose ;

//     //the Solver for the elliptic Problem
    
//     ElliptInverseOperatorType aufSolver_; 
    
//     // tmp variables 
//     mutable DiscreteFunctionType* r_;
//     mutable DiscreteFunctionType* p_;
//     mutable DiscreteFunctionType* h_;
//   };


} // end namespace Dune

#endif

