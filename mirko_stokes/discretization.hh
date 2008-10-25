/***********************************************************************************************
 Discretization is the interface class used to define the discretization of the problem
***********************************************************************************************/
#ifndef __STOKES_DISCRETIZATION_HH__
#define __STOKES_DISCRETIZATION_HH__

/** \file
    \brief discretization.hh is the interface class used to define the discretization of the problem
 */

#define GLATT 1			
#define SNG 0
#define STEP 0

/* include definition of the physical problem (coefficient functions, etc.) */
#include "mystokesmodell.hh"

#include "fluxes.hh"

// choose discrete function type
//#include <dune/fem/space/lagrangespace.hh>
//#include <dune/fem/discretefunction/dfadapt.hh>
//#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/gridpart/gridpart.hh>

// include L2 projection into discrete funcxtion spaces
//#include "stuff.cc"
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/misc/l2error.hh>
// ascii parser for parameter files 
#include <dune/fem/io/file/asciiparser.hh>

// include rouines to generate  the grid 
//#include "makegrid.cc"

// grape data io 
#include <dune/fem/io/file/grapedataio.hh>

// if grape was configured then include headers 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>

#endif

#include "discretemodels.hh"

#include <dune/fem/solver/inverseoperators.hh>
// #include <dune/fem/operator/SPinverseoperator_istl.hh>
#include "SPinverseoperator_original.hh"
//#include "SPinverseoperator.hh"
#include <dune/fem/solver/oemsolver/oemsolver.hh>
//#include "spmatrix.hh"
//#include <dune/fem/solver/istlsolver.hh>
//#include<dune/fem/operator/matrix/spmatrix.hh>
//#include<dune/fem/operator/matrix/istlmatrix.hh>

#include "Stokespass.hh"

//#include "Stokespasslocal2.hh"


using namespace Dune;


namespace LDGExample { 

  template <class ModelImpType,class GradientSpaceImp , int polOrd=0 >
struct DiscrParam
{
  typedef ModelImpType  ModelType;

  enum { dimRange = ModelType::dimRange };
 
  typedef typename ModelType::GridType             GridType;

  enum { dim      = GridType::dimension };
  enum { dimworld = GridType::dimensionworld };
  enum { polyOrder = polOrd };

  typedef typename ModelType::FieldType            FieldType;
  typedef typename ModelType::FuncSpaceType        FuncSpaceType;
 
  typedef LeafGridPart < GridType > GridPartType ;
};

  

// The actual operator
// StokesModelType is StokesDiscreteModel siehe discretemodel.hh
template <class StokesModelType> 
class MySpaceOperator :
 public Operator<
    double,
    double,
    typename StokesModelType::Traits::DestinationType,
    typename StokesModelType::Traits::DestinationType>
{
  typedef typename StokesModelType :: Traits Traits;
public:
  typedef MySpaceOperator<StokesModelType> ThisType;
  enum { polOrd = StokesModelType :: polynomialOrder };
  
  typedef typename Traits:: DestinationType DestinationType;
  typedef typename Traits:: DestinationType DiscreteFunctionType;
  typedef typename Traits:: DestinationType SolutionType;
  typedef typename Traits::GridPartType GridPartType ;
  typedef typename GridPartType :: Traits :: GridType GridType;

  typedef DofManager<GridType> DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  //typedef typename Traits::ContainedSpaceType ContainedSpaceType;

  //  typedef TimeDiscrParam TimeDiscParamType;
  enum PassIdType{ u , pass1 };  
  typedef StartPass<DiscreteFunctionType,(int)u> Pass0Type;
  // note, the destination type of the pass 0 is the argument type of pass 1
  typedef StokesFEPass<StokesModelType , Pass0Type,(int) pass1> Pass1Type;
  
  typedef typename Pass1Type::PressureSpaceType PressureSpaceType;
  typedef typename Pass1Type::DiscretePressureSpaceType DiscretePressureSpaceType;
  typedef typename Pass1Type::DiscretePressureType DiscretePressureFunctionType;
 


  typedef OEMCGOp <DestinationType, Pass1Type> InverseOperatorType;
     typedef DiscreteFunctionSpaceType Space1Type;
 
  typedef typename Space1Type:: FunctionSpaceType FuncSpaceType;
  
 
  typedef typename FuncSpaceType::RangeType RangeType;

  typedef SPCGInverseOperator<DestinationType,DiscretePressureFunctionType,Pass1Type,InverseOperatorType> SPInverseOperatorType;
  //! interface methods for the initial data
  
  //  #include "testproblem.cc"w
#if 1
  class f_RHSData : public Function < FuncSpaceType , f_RHSData > {
  public:
    typedef typename FuncSpaceType :: DomainType DomainType; 
    typedef typename FuncSpaceType :: RangeType  RangeType; 
    
    f_RHSData (FuncSpaceType& f)
      : Function < FuncSpaceType , f_RHSData > ( f ) { } ;
#if GLATT    
    template <class FaceDomainType>
    RangeType frhs(const FaceDomainType & x) const
    {
      enum { dim = FaceDomainType::dimension };
      RangeType ret;
      
      double factor=0;
      
      ret[0]=sin(x[1]);
      ret[0]*=2;
      ret[0]*=exp(x[0]);
      ret[0]*=factor;
      ret[1]=cos(x[1]);
      ret[1]*=2;
      ret[1]*=exp(x[0]);
      ret[1]*=factor;
      //ret[1]+=100;
      return ret;
 }
#endif
#if SNG
   template <class FaceDomainType>
    RangeType frhs(const FaceDomainType & x) const
    {
      enum { dim = FaceDomainType::dimension };
      RangeType ret;
      
       ret[0]=sin(x[1]);
       ret[0]*=2;
       ret[0]*=exp(x[0]);
   
   
     
      ret[1]=cos(x[1]);
      ret[1]*=2;
      ret[1]*=exp(x[0]);
   
     
      
      return ret;
 }
#endif
#if STEP
   template <class FaceDomainType>
    RangeType frhs(const FaceDomainType & x) const
    {
      enum { dim = FaceDomainType::dimension };
      RangeType ret;
      
       ret[0]=0;
       ret[1]=-0.02;
       return ret;
 }
#endif





    void evaluate (const DomainType & x , RangeType & ret) const
    {
      enum { dimR = RangeType :: dimension }; 
      ret = 0.0;
      ret = frhs( x );
      return;
      }
  
    void evaluate (const DomainType & x , double time, RangeType & ret) const
    { evaluate (x,ret); return;};
  };


  //! the exact solution to the problem for EOC calculation 
  class ExactSolution : public Function < FuncSpaceType , ExactSolution >
  {
    typedef typename FuncSpaceType::RangeType RangeType;
    typedef typename FuncSpaceType::RangeFieldType RangeFieldType;
    typedef typename FuncSpaceType::DomainType DomainType;
  public:
    ExactSolution (FuncSpaceType &f) : Function < FuncSpaceType , ExactSolution > ( f ) {}

#if GLATT
    //! u1(x,y) = -exp(x)(y cos y+sin y)
    //! u2(x,y) = exp(x)(y sin y)
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      /////////////////////////////////////////////////
      //paper-beispiel      
      ret[0]=cos(x[1]);
      ret[0]*=x[1];
      ret[0]+=sin(x[1]);
      ret[0]*=exp(x[0]);
      ret[0]*=-1.0;


      ret[1]=sin(x[1]);
      ret[1]*=x[1];
      ret[1]*=exp(x[0]);



    }
#endif
#if SNG
  void evaluate (const DomainType & x , RangeType & ret) const
    {
     
     
      DomainType tmp=x;
      tmp[0]-=1.01;
      tmp[1]-=1.01;
      double abs=tmp*tmp  ;
      
      ret[0]=tmp[0];
      ret[0]/=abs;
      ret[1]=tmp[1];
      ret[1]/=abs;


    }
#endif
#if STEP
 void evaluate (const DomainType & x , RangeType & ret) const
    {
      ret[0]=0;
      ret[1]=0;
    }
#endif

    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };
 

 //! the exact solution to the problem for EOC calculation 
  class ExactPressure : public Function < PressureSpaceType , ExactPressure >
  {
    typedef typename PressureSpaceType::RangeType RangeType;
    typedef typename PressureSpaceType::RangeFieldType RangeFieldType;
    typedef typename PressureSpaceType::DomainType DomainType;
  public:
    ExactPressure (PressureSpaceType &f) : Function<PressureSpaceType , ExactPressure > ( f ) {}

   
    //! p(x,y) = 2exp(x)(sin y)
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      RangeType tmp(0.0);
    //   ret=sin(x[1])*cos(x[0]);
      //    ret=0.0;
      
      ret = sin(x[1]);
      ret *=exp(x[0]);
      ret*=2.0;
      
      //for step
      //   ret*=0;
   
    }
    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };
  
#endif



  //! the exact solution to the problem for EOC calculation 
  template <class FuncSPCType>
  class ExactGradient : public Function < FuncSPCType , ExactGradient<
                        FuncSPCType > >
  {
    typedef typename FuncSPCType::RangeType RangeType;
    typedef typename FuncSPCType::RangeFieldType RangeFieldType;
    typedef typename FuncSPCType::DomainType DomainType;
  public:
    ExactGradient (const FuncSPCType &f) 
      : Function < FuncSPCType , ExactGradient< FuncSPCType > > ( f ) {}

    //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      ret = 1.0;

      enum { dimR = RangeType :: dimension }; 
      enum { dimD = DomainType::dimension  };

      for(int j=0; j<dimR; ++j)
      { 
        for(int i=0; i<dimD; ++i)
        {
          if(i == j) ret[j] *= (1.0-2.0*x[i]);
          else       ret[j] *= ( x[i] - SQR(x[i]) );
        }
      }
    }

    void evaluate (const DomainType & x , double time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };


  MySpaceOperator ( GridType & grid , StokesModelType & lpm,double c11,double c12,
		    double d11,double d12) 
    : grid_(grid)
    , dm_(DofManagerFactoryType::getDofManager(grid_))
    , model_ (lpm) 
    , gridPart_(grid_)
    , space1_(gridPart_)
    , pressurespc_(gridPart_)  
    , pass0_()
    , pass1_(lpm , pass0_,space1_/*,c11,c12,d11,d12*/)
    , pressure_("press",pressurespc_)
    , rhs_("rhs",space1_)
    , sol_("sol",space1_)
    , errorout_(0.0)
    , presserrorout_(0.0)
  {
    std::cout << "Created SpaceOperator \n";

    std::cout.flush();
  }

  void operator()(const DestinationType& arg, DestinationType& dest) const 
  { 
    const_cast<ThisType&> (*this).apply(arg,dest);
  }
  
  // apply space discretisation 
  void apply(const DestinationType& arg, DestinationType& dest)
  {
    
    const int steps = 2;
    double error[steps];  
    //RangeType error[steps];

    double press[steps];
    
    //GradRangeType gradError[steps];
    
    for(int i=0; i<steps; ++i)
    {
      if(i > 0)
      {
        // refineGlobal is defined in description.hh
	grid_.globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
        dm_.resize();   
	data[0]=grid_.size(0);
      }

      
      DestinationType rhs(dest);
      //DestinationType tmp(dest);
      DiscretePressureFunctionType pressuretest("press_sol",pressurespc_);
      
      pressure_.clear();
      dest.clear();
      FuncSpaceType sp; 
      //      f_RHSData rhsfunc(sp); 
      //ExactPressure testf(pressurespc_);
      DestinationType exactsol("exact",space1_) ;
      
//        StokesSingExactSolution<FuncSpaceType> exact(sp); 
//           StokesSingExactPressure< DiscretePressureSpaceType> testf(pressurespc_);
//       L2Projection <double,double,f_RHSData,DestinationType> l2pro; 
//        L2Projection <double,double,ExactPressure,DiscretePressureFunctionType> l2pro_pr;
      StokesExactSolution<FuncSpaceType> exact(sp); 
      StokesExactPressure< PressureSpaceType> testf(pressurespc_);
      //      L2Projection <double,double,f_RHSData,DestinationType> l2pro;   
     //  L2Projection <double,double,StokesExactSolution<FuncSpaceType>,DestinationType> l2pro; 
//       L2Projection <double,double,StokesExactPressure<PressureSpaceType>,DiscretePressureFunctionType> l2pro_pr; 
//       //    l2pro(exact,exactsol);
       
      /*
      {
       rhs.clear();
       typedef f_RHSData RHSDataType;
       L2Projection < DestinationType> l2pro;
       l2pro.project( rhsfunc , rhs );
       l2pro.project( exact,dest);

	L2Projection<DiscretePressureFunctionType> pl2pro;
	pl2pro.project(testf,pressure_);
	}*/
      
      ///exactsol.assign(dest);
    //   if(0)
    //    GrapeDataDisplay < GridType > grape(grid_); 
//        grape.dataDisplay(pressure_ , false );
//       grape.addData(dest);
//       grape.addData(pressure_);
//       grape.display( );
      
     
      pass1_.buildMatrix(rhs);
      InverseOperatorType aufSolver(pass1_,1e-10,1e-10,5000,false);
      SPInverseOperatorType invOp(pass1_,1e-8,1e-8,5000,true,aufSolver,pressurespc_,space1_);
      DiscretePressureFunctionType start(pass1_.rhs2());

      DestinationType rhs1(pass1_.rhs1());

      //
       pressure_.clear();
      dest.clear();
      //  pressure_.assign(pressuretest);
      
      invOp(start,pressure_); 
      sol_.assign(invOp.velocity());
     //  pressuretest-=pressure_;     
       
	 
// 	 GrapeDataDisplay < GridType > grape( grid_ ); 
// 	 grape.dataDisplay(pressuretest , false );
	 
       
       // sol_.assign(invOp.velocity());
       pass1_.finalizeGlobal();

       L2Error < DiscretePressureFunctionType > l2perr;
       L2Error < DestinationType > l2err;
      
       // pol ord for calculation the error chould by higher than 
       // pol for evaluation the basefunctions 
       //exact.clear();
       //   testf.clear();
       error[i] = 0.0;/*l2err.norm3(exact , sol_);*/
       press[i] = l2perr.norm(testf,pressure_);
       
       //for(int k=0; k<GradRangeType::dimension; k++)
       //  std::cout << "\nGradError : " << gradError[i][k] << "\n";
       std::cout << "\npressure Error : " << press[i] << "\n";
   //     for(int k=0; k<RangeType::dimension; k++)
       std::cout << "\nError : " << error[i] << "\n";
       if(i==steps-1)
	 {
	   data[1]= error[i];
	   data[3]=press[i];
	 }
       if( i > 0 )
	 { 
	   // for(int k=0; k<RangeType::dimension; k++)
// 	    {
      	      double eoc = log( error[i-1]/error[i]) / M_LN2;
	      std::cout << "EOC["<<"] = " << eoc << " \n";
	      
// 	    }
	   double pressureeoc = log( press[i-1]/press[i]) / M_LN2;
	   std::cout << "EOC pressure = " << pressureeoc << " \n";
	   
	   if(i==steps-1)
	       {
		 data[2]= eoc;
		 data[4]=pressureeoc;
	       }

	   
	 }
    }
  }
  
  template <class TimeProviderType>
  void timeProvider(TimeProviderType & tp )
  {
    pass1_.timeProvider(&tp);
  }
  
  SolutionType& solution () 
  {
    //SolutionType& dest = const_cast<SolutionType&> (pass1_.destination());
    return sol_;
  } 
  
  DiscretePressureFunctionType& pressure(){
    return pressure_;
  }



  GridPartType & gridPart () { return gridPart_; }
  DestinationType * createDestinationFct (std::string name) 
    {
    return new DestinationType (name , space1_ );
  }

  double errorout(){return errorout_;}
  double presserrorout(){return presserrorout_;}


  double* dataptr(){return data;}


private:
  GridType & grid_;
  DofManagerType & dm_;
  const StokesModelType & model_;

  // we use the same index set and grid part for all spaces
  GridPartType gridPart_;
  
  mutable Space1Type space1_;
  mutable DiscretePressureSpaceType pressurespc_;
  
  mutable Pass0Type pass0_;
  mutable Pass1Type pass1_;
  mutable DiscretePressureFunctionType pressure_;
  mutable DestinationType rhs_;
  mutable DestinationType sol_; 
  mutable double errorout_;
  mutable double presserrorout_;
  mutable double data[5];
};



/*****************************************************************************************************/
/*                                                                                                   */
/*                                      SIMULATION ROUTINE                                           */
/*                                                                                                   */
/*****************************************************************************************************/

template <class DiscrType> 
void simul(typename DiscrType::ModelType & model, std::string paramFile) 
{
try{
   typedef typename DiscrType::GridType              GridType;
   typedef typename DiscrType::ModelType             ModelType;
   enum { polOrd = DiscrType::polyOrder };
   
   std::ofstream myfile;
   std::ofstream plotfile;
   //  myfile.open("dataglatt10d11.txt");
   //  plotfile.open("plotglatt10d11.txt");
   //myfile<<"level: "<<level<<"\n";
   //  myfile<<"\\begin{tabular}{l|l|l|l|l}\n";
   //    myfile<<"d11 & $||u-u_{h}||_2$ & EOC & $||p-p_{h}||_2$ & EOC  \\\\ \n";
   //   myfile.open("dataSING.txt");
   //    plotfile.open("plotSING.txt");
  
   // choice of fluxes 
   typedef StokesFlux <ModelType> NumericalFluxType;
   
   //typedef GradientFlux GradientFluxType;
    
   typedef StokesDiscreteModel < ModelType, NumericalFluxType, polOrd > StokesModelType;
   
    
    typedef MySpaceOperator < StokesModelType> 
      SpaceOperatorType; 
    typedef typename SpaceOperatorType :: DestinationType DestinationType;
  




    
    // initialize grid
    char dummyfile [4096];
    const char * paramfile = paramFile.c_str();
    readParameter(paramfile,"Grid",dummyfile);
    std::string macroGridName(dummyfile);

    int level;
    readParameter(paramfile,"StartLevel",level);
    
    readParameter(paramfile,"Output",dummyfile);
    std::string myfileName(dummyfile);
    myfile.open(myfileName.c_str());
    myfile<<"level: "<<level<<"\n";
    myfile<<"\\begin{tabular}{l|l|l|l|l}\n";
    // myfile<<"d11 & $||u-u_{h}||_2$ & EOC & $||p-p_{h}||_2$ & EOC  \\\\ \n";
    readParameter(paramfile,"Plot",dummyfile);
    std::string myplotName(dummyfile);
    plotfile.open(myplotName.c_str());   
  
    myfile<<"level & $||u-u_{h}||_2$ & EOC & $||p-p_{h}||_2$ & EOC  \\\\ \n";
 
    
    
   
    //int bla;
    //fscanf(stdin,"%d",&bla);
  
    

    double c11;
    readParameter(paramfile,"c11",c11);
    double c12; 
    readParameter(paramfile,"c12",c12);
    double d11;
    readParameter(paramfile,"d11",d11);
    double d12; 
    readParameter(paramfile,"d12",d12);
    int runs;
    readParameter(paramfile,"Runs",runs);
      
    double nu;
    readParameter(paramfile,"nu",nu);
    //    nu=1.0;

    int display = 0;
    readParameter(paramfile,"display",display);
        
    
    for(int i=0;i<runs;i++)
      {
	//nu*=10.0;
	GridPtr<GridType> gridptr(macroGridName,MPIHelper::getCommunicator()); 
	GridType & grid = *gridptr;

	grid.globalRefine(DGFGridInfo<GridType>::refineStepsForHalf() * level+i);
	grid.loadBalance();
	std::cout << "Grid size = " << grid.size(0) << "\n";
	std::cout << "using nu = "<< nu << "\n";
	model.setnu(nu);
	
	
	NumericalFluxType numericalFlux(model,c11,c12,d11,d12);
	
	StokesModelType lpm(model, numericalFlux,nu);
	
	SpaceOperatorType spaceOp( grid,lpm,c11,c12,d11,d12 );
    
  
	//! storage for the discrete solution and its update
	DestinationType *solution = spaceOp.createDestinationFct("solution");
	DestinationType *tmpRhs = spaceOp.createDestinationFct("tmpRhs");

	// initial data != 0
	solution->clear();
	spaceOp(*tmpRhs,*solution);
	double entry=spaceOp.errorout();
	double entry2=spaceOp.presserrorout();
	double* data=spaceOp.dataptr();
	//myfile<<"d11="<<d11<<" velo error= " << entry <<" pressure error= "<<entry2<<"\n";
 	myfile<< "$" << data[0]  <<"$ & $ "<<data[1]<<"$ & $"<<data[2]<<"$ & $"<<data[3] <<"$ & $"<<data[4]<<"$  \\\\ \n";
	plotfile << data[0]  <<" "<<data[1]<<" "<<data[2]<<" "<<data[3] <<" "<<data[4]<<" \n"; 

	std::string name;
	name="velocity";
#if HAVE_GRAPE
	if( display )
	  {
	    GrapeDataDisplay < GridType > grape( spaceOp.gridPart() ); 
	    grape.addData(spaceOp.solution());
	    grape.addData(spaceOp.pressure());
	    grape.dataDisplay( spaceOp.solution() );
	    //	vtkout(spaceOp.gridPart,spaceOp.solution,"velocity",1);
	    // 	VTKIO<GridPartType> vtkio( spaceOp.gridPart() );
	    // 	vtkio.addVertexData(spaceOp.solution() ); 
	    // 	vtkio.write(name.c_str(),Dune::VTKOptions::ascii);
	    
	  }
#endif
       
    delete solution; 
    delete tmpRhs;
      }
    myfile<<"\\end{tabular}";
    myfile.close(); 
    std::cout<<myfileName<<"\n";
    plotfile.close();
}
	catch ( Dune::Exception & e )
	{
		std::cout << "DUNE ERROR: " << e.what() << std::endl;
	}
	catch ( std::exception f )
	{
		std::cout << "Unknown ERROR: " << std::endl;
	}     
  }


} // end namespace LDGExample
#endif
