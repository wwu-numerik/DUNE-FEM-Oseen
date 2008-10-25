#ifndef DUNE_STOKESPASS_HH
#define DUNE_STOKESPASS_HH

/** \file
    \brief Stokespasslocal2.hh
 */

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/modelcaller.hh>
#include <dune/fem/operator/matrixoperator.hh>
#include <dune/fem/operator/2order/dgmatrixsetup.hh>

// * needs to move
// #include "../misc/timenew.hh"
#include <dune/fem/misc/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  //
  template <class DiscreteModelImp, class PreviousPassImp>
  class StokesFEPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;
    typedef StokesFEPass<DiscreteModelImp,PreviousPassImp> ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType; //velocity function

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType VelocityDiscreteFunctionSpaceType;
    
    
    // Types extracted from the discrete function space type
    typedef typename VelocityDiscreteFunctionSpaceType::GridType GridType;
    typedef typename VelocityDiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename VelocityDiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename VelocityDiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename VelocityDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
        
    // Types extracted from the underlying grid
    typedef typename GridType::Traits::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalVeloType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;
   
    typedef SparseRowMatrix<double> MatrixType;
    
    




    // Range of the destination
    enum { dimDomain = VelocityDiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = VelocityDiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //meine Typedefs

    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd =VelocityDiscreteFunctionSpaceType::polOrd};

    typedef typename VelocityDiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename VelocityDiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    
    typedef MatrixFunctionSpace<DomainFieldType,RangeFieldType,dimDomain,dimRange,dimRange> GradientSpaceType;
    //typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,dimGradRange> GradientSpaceType;
    typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd> DiscreteGradientSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd> DiscreteGradientSpaceType;
    typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,1> PressureSpaceType;
    typedef DiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd> DiscretePressureSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd> DiscretePressureSpaceType;
    
    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradJacobianRangeType;
    typedef AdaptiveDiscreteFunction<DiscreteGradientSpaceType> DiscreteStressType;
    typedef typename DiscreteStressType::LocalFunctionType LocalStressType;
    
    typedef typename DiscretePressureSpaceType::RangeType PressureRangeType;
    typedef typename DiscretePressureSpaceType::JacobianRangeType PressureJacobianRangeType;
    typedef AdaptiveDiscreteFunction<DiscretePressureSpaceType> DiscretePressureType;
    typedef typename DiscretePressureType::LocalFunctionType LocalPressureType;
    
    
    typedef MatrixOperator<MatrixType,PressureRangeType,RangeType> BOPType;
    typedef MatrixOperator<MatrixType,RangeType,PressureRangeType> BTOPType;
    typedef MatrixOperator<MatrixType,PressureRangeType,PressureRangeType> COPType;
    
    typedef SparseRowMatrixObject<DiscreteGradientSpaceType,VelocityDiscreteFunctionSpaceType> GradMatType;
    typedef typename GradMatType::LocalMatrixType LocalGradMatType;
    typedef typename GradMatType::MatrixType  GMType;
   
   
    typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> StabMatType;
    typedef typename StabMatType::LocalMatrixType LocalStabMatType;
    typedef typename StabMatType::MatrixType SMType;

   
    typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,DiscretePressureSpaceType> PressureGradMatType;
    typedef typename PressureGradMatType::LocalMatrixType LocalPressureGradMatType;
    typedef typename PressureGradMatType::MatrixType  PGMType;
   
   
    typedef SparseRowMatrixObject<DiscretePressureSpaceType,DiscretePressureSpaceType> PressureStabMatType;
    typedef typename PressureStabMatType::LocalMatrixType LocalPressureStabMatType;
    typedef typename PressureStabMatType::MatrixType PSMType;
    

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param quadOrd0 defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param quadOrd1 defines the order of the face quadrature which is by default 2* space polynomial order 
    StokesFEPass(DiscreteModelType& problem, 
		 PreviousPassType& pass, 
		 VelocityDiscreteFunctionSpaceType& spc,
		 int volumeQuadOrd =-1,int faceQuadOrd=-1) :
      BaseType(pass, spc),
      problem_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gradientspc_(spc_.gridPart()),
      pressurespc_(spc_.gridPart()),
      tmp_("FEPass::tmp",gradientspc_),
      tmp2_("FEPass::tmp2",gradientspc_),
      tmp1_("FEPass::tmp1",spc_),
      dtMin_(std::numeric_limits<double>::max()),
      VeloRhs_("VeloRhs",spc_),
      StressRhs_("StressRhs",gradientspc_),
      PressureRhs_("PressureRhs",pressurespc_),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0), 
      phi_(0.0),
      phiNeigh_(0.0),
      psi_(0.0),
      tau_(0.0),
      tauNeigh_(0.0),
      grads_(0.0),
      time_(0),
      diffVar_(),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (5*spc_.order()+1) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (5*spc_.order()+1) : faceQuadOrd ),
      maxNumberUnknowns_(5* (gradientspc_.baseFunctionSet(*(spc_.begin())).numBaseFunctions())),
      maxNum_(5*(spc_.baseFunctionSet(*(spc_.begin())).numBaseFunctions())),
      gradMatrix_(gradientspc_,spc_,""),//VDM
      stabMatrix_(spc_,spc_,""),//VSM
      pressureGradMatrix_(spc_,pressurespc_,""),//PGM
      pressureStabMatrix_(pressurespc_,pressurespc_,""),//PSM
      matrixAssembled_(false),
      direction_(0.0)
    {
      
      std::cout<<"unknowms:"<<maxNumberUnknowns_<<"\n";
      direction_[0]=1;
      direction_[1]=(1.0+sqrt(5.0))/2.0;
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );
    }
    
    //! Destructor
    virtual ~StokesFEPass() {
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void processTimeProvider(TimeProvider* time) {
      time_ = time;
    }

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return dtMin_;
    }


    template <class DFType> 
    void buildMatrix(DFType & rhs )
    {
      int size = spc_.size(); 
      int singleSize = spc_.size();
      int gradSize = gradientspc_.size();
      int pressureSize = pressurespc_.size();
      std::cout << "Resize Matrix with " << size << "\n";
      
      
      
      tmp1_.clear();
      VeloRhs_.clear();
      VeloRhs_.assign(rhs);
      rhs.clear();
      dest_=&rhs;
      StressRhs_.clear();
      PressureRhs_.clear();
      
      ElementAndNeighbors stencil;

      gradMatrix_.reserve(stencil);
      gradMatrix_.clear();
      
      stabMatrix_.reserve(stencil);
      stabMatrix_.clear();
      
      
      pressureGradMatrix_.reserve(stencil);
      pressureGradMatrix_.clear();
      
      pressureStabMatrix_.reserve(stencil);
      pressureStabMatrix_.clear();


     

      this->operator()(rhs,rhs);
      matrixAssembled_ = true;
    
  
      double* stressptr=StressRhs_.leakPointer();
      double* veloptr=VeloRhs_.leakPointer();
      double* tmp1ptr=tmp1_.leakPointer();
 //     std::cout<<"************************************************\n";
//              std::cout<<"Stokespass:pressuregradmatrix=\n"; 
//              pressureGradMatrix_.matrix().print(std::cout);
//              std::cout<<"************************************************\n";
//   int foo;
//        std::cin>> foo;
//       assert(foo==1);
   //    gradMatrix_.matrix().multOEM_t(stressptr,tmp1ptr);
//       for(register int i=0;i<spc_.size();++i)
// 	{
// 	  veloptr[i]-=tmp1ptr[i];
// 	}
      gradMatrix_.matrix().multOEM_t(stressptr,tmp1ptr);
      (*dest_)-=tmp1_;
    
    }
    
    
    MatrixType& getBOP()const {return pressureGradMatrix_.matrix();}
    MatrixType& getCOP() const {return pressureStabMatrix_.matrix();}

    DiscretePressureSpaceType& getScalSpc() const {return pressurespc_;}
    VelocityDiscreteFunctionSpaceType getVecSpc() const {return spc_;}
    


    void multOEM(const double * arg, double * dest) const
    {
      //  if(!matrixAssembled_) 
      {

	const int size = spc_.size();

	tmp_.clear();
	tmp1_.clear();
	tmp2_.clear();
      
	double * tmpP  = tmp_.leakPointer();
	double * tmpP1 = tmp1_.leakPointer();
      
	gradMatrix_.matrix().multOEM(arg,tmpP);
	gradMatrix_.matrix().multOEM_t(tmpP,tmpP1);
	
	stabMatrix_.matrix().multOEM(arg,dest);

	for(register int i=0; i<size; ++i) dest[i] += tmpP1[i];
	
	return ;
      }
    }
    const ThisType & systemMatrix () const { return *this; }
    
  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      //dest_->clear();
 
      //  caller_.setArgument(*arg_);

      // time initialisation
      dtMin_ = std::numeric_limits<double>::max();
      if (time_) {
	// caller_.setTime(time_->time());
      }
      else {
	//    caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      if (time_) {
        time_->provideTimeStepEstimate(dtMin_);
      }

      // caller_.finalize();
    }

  public:
    const VelocityDiscreteFunctionSpaceType & getSpace() const {
      return spc_;
    }

  public:
    void prepareGlobal(const ArgumentType& arg, DestinationType& dest) const{
      prepare(arg,dest);
    }

    void finalizeGlobal() {
      //finalize(*arg_,*dest_); 
      matrixAssembled_ =false;
    }

    void applyLocal(EntityType& en) const
    {
      //- typedefs
     
      typedef typename VelocityDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      
      
      typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
      
     
      typedef typename DiscretePressureSpaceType::BaseFunctionSetType PressureBaseFunctionSetType;
     
      
      typedef typename DestinationType::LocalFunctionType LocalFuncType;
     
      LocalFuncType localrhs=dest_->localFunction(en);
  
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      double massVolElInv;
      double vol = volumeElement(en, volQuad,massVolElInv);
      
      const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      
      const GradientBaseFunctionSetType& gradbsetEn = gradientspc_.baseFunctionSet(en);
      const int gradientNumDofs = gradbsetEn.numBaseFunctions();
      
      const PressureBaseFunctionSetType& pressurebsetEn = pressurespc_.baseFunctionSet(en);
      const int pressureNumDofs = pressurebsetEn.numBaseFunctions();
     
      LocalGradMatType en_grad = gradMatrix_.localMatrix(en,en);
      LocalStabMatType en_stab=stabMatrix_.localMatrix(en,en);
      LocalPressureGradMatType en_pressGrad=pressureGradMatrix_.localMatrix(en,en);
      LocalPressureStabMatType en_pressStab=pressureStabMatrix_.localMatrix(en,en);
      
      double nu=problem_.get_nu();
      
      double nu_rt=sqrt(nu);
      
      // Volumetric integral part
      const int quadNop = volQuad.nop();
      
      GradJacobianRangeType helpmatr(0.0);
      JacobianRangeType helpmatr2(0.0); 
    
      //if(0)     
      for (int l = 0; l < quadNop ; ++l) 
	{
	  GradJacobianRangeType helpmatr(0.0);
	  JacobianRangeType helpmatr2(0.0); 
	  
	  double intel = volQuad.weight(l);    
	   
	  for(int k = 0;k < gradientNumDofs ;++k)
	    {	      
	      // eval tau_k
	      gradbsetEn.evaluate(k, volQuad, l,tau_);
	      
	      for (int j = 0; j < numDofs; ++j) 
		{  
		  bsetEn.evaluate(j, volQuad, l, psi_[0]);
		  if(k==0)
		    {
		      //Rhs Assembling
		      RangeType fval(0.0);
		      problem_.rightHandSide(en,volQuad,l,fval); 
		      localrhs[j]+=bsetEn.evaluateSingle(j,volQuad,l,fval)*intel;
		    }
	
		  
		  // needed for calculating phi_j*div tau_k with the evaluateGradientSingle method
		  helpmatr*=0.0;
		  
		  helpmatr[0][0]=psi_[0][0];
		  helpmatr[1][1]=psi_[0][0];
		  helpmatr[2][0]=psi_[0][1];
		  helpmatr[3][1]=psi_[0][1];
		  			  
		  
		  // eval -phi_j*div tau_k
		  // -VGM
		  double VGM =-1.0* gradbsetEn.evaluateGradientSingle(k,en,volQuad,l,helpmatr)*intel*nu_rt;
		  en_grad.add(k,j,VGM);
		  
		  ////////////////////////////////////////////////
		  //pressureGradientMatrix , pressureDivMatrix  //
		  ////////////////////////////////////////////////
		  
		  
		  if(k==0)
		    {
		     

		      for(int n = 0; n < pressureNumDofs ; ++n)
			{
			  //eval q_n
			  pressurebsetEn.evaluate(n,volQuad,l,qu_);
			  
			  helpmatr2=0.0;
			  for(int m = 0; m < dimDomain ; ++m)
			    {
			      helpmatr2[m][m]=qu_[0];
			    }
			  int PGM_cols  =  pressurespc_.mapToGlobal(en,n);  
			  int PDM_rows  =  PGM_cols;
			
			  //eval (-q_n*div phi_j)
			  double PGM =-1.0*bsetEn.evaluateGradientSingle(j,en, volQuad,l,helpmatr2)*intel; 
			 
			  en_pressGrad.add(j,n,PGM);

			}
		    }    
		}
	      
	    }
	}
	
      // Surface integral part
      // We use the unitOuternormal for consistency with the flux definitions
      // The integrationelement is calculated seperately

      IntersectionIterator endnit = en.ileafend();

      double dtLocal = 0.0;
      double minvol = vol; 

      // if(0)
      for (IntersectionIterator nit = en.ileafbegin(); nit != endnit; ++nit) 
	{ 

	  FaceQuadratureType faceQuadInner(spc_.gridPart(),nit, faceQuadOrd_,  
					   FaceQuadratureType::INSIDE);
	  
	  if (nit.neighbor()) 
	    { 
	      //Access to the neighbor Element
	      EntityPointerType neighEp=nit.outside();
	      EntityType& nb=  *neighEp;
				
	      { 
	
		FaceQuadratureType faceQuadOuter(spc_.gridPart(),nit, faceQuadOrd_, 
						 FaceQuadratureType::OUTSIDE);
	

		//the Basefunctiosets on the Neighbor element
		const BaseFunctionSetType& bsetNeigh = spc_.baseFunctionSet(nb);
		const GradientBaseFunctionSetType& gradbsetNeigh = gradientspc_.baseFunctionSet(nb);
		const PressureBaseFunctionSetType& pressurebsetNeigh = pressurespc_.baseFunctionSet(nb);

		//local Matrices on the Neighbor element
 		LocalGradMatType nb_grad = gradMatrix_.localMatrix(en,nb);
		LocalStabMatType nb_stab=stabMatrix_.localMatrix(en,nb);
		LocalPressureGradMatType nb_pressGrad=pressureGradMatrix_.localMatrix(en,nb);
		LocalPressureStabMatType nb_pressStab=pressureStabMatrix_.localMatrix(en,nb);


		double massVolNbInv;
		double nbvol = volumeElement(nb, volQuad, massVolNbInv);
		if (nbvol<minvol) minvol=nbvol;

		const int quadNop = faceQuadInner.nop();
		for (int l = 0; l < quadNop ; ++l) 
		  {
		   
		    DomainType normal=nit.unitOuterNormal(faceQuadInner.localPoint(l));
		    		      
		    double int_element = nit.intersectionGlobal().integrationElement(faceQuadInner.localPoint(l));
	
	      	
		    for(int j=0; j< numDofs; ++j)
		      { 
		
			  
			//calculatin phi_j on the element- and neighbor side of the current intersection
			bsetEn.evaluate(j, faceQuadInner, l, phi_); 
		
			bsetNeigh.evaluate(j, faceQuadOuter, l, phiNeigh_ );
		


			//this stuff should be made member of the class
			RangeType zerophi(0.0); 
			GradientRangeType zerotau(0.0);
			GradientRangeType enflux(0.0);
			GradientRangeType neighflux(0.0);
			RangeType stab_en(0.0);
			RangeType stab_nb(0.0);
		
			  
			PressureRangeType zeroqu(0.0);
			PressureRangeType pressure_en(0.0);
			PressureRangeType pressure_neigh(0.0);
		
			

			//calculating the fluxes
			//\hat{u}_{\sigma}
			problem_.numericalFlux_u_sigma(nit,faceQuadInner,l,phi_,phiNeigh_,enflux,neighflux);
			
			//\check{\sigma}
			problem_.numericalFlux_u_stab(nit,faceQuadInner,l,phi_,phiNeigh_,stab_en,stab_nb);
			
			problem_.numericalFlux_u_p(nit,faceQuadInner,l,phi_,phiNeigh_,pressure_en,pressure_neigh);
			
			//Loop over Grad Dofs calculation  VGM and VDM 
			for(int k=0;k < gradientNumDofs;++k)
			  {
			  
			      
			    gradbsetEn.evaluate(k,faceQuadInner,l,tau_);
			    gradbsetNeigh.evaluate(k,faceQuadOuter,l,tauNeigh_);
			      
			    RangeType sflux_en;
			    RangeType sflux_nb;

			    
			   problem_.numericalFlux_sigma(nit,faceQuadInner,l,tau_,tauNeigh_,sflux_en,sflux_nb);
			    
			    double VGM_en =  gradbsetEn.evaluateSingle(k,faceQuadInner,l,enflux)
			      *faceQuadInner.weight(l)* massVolElInv*int_element*nu_rt;
			    double VGM_nb = gradbsetEn.evaluateSingle(k,faceQuadInner,l,neighflux)
			      *faceQuadInner.weight(l)* massVolElInv*int_element*nu_rt;
			     
		
			    en_grad.add(k,j,VGM_en);
			    nb_grad.add(k,j,VGM_nb);
			    
		
			      
			  }
			 
			for(int n=0;n<pressureNumDofs;++n)
			  {
		
			      
			     
			     
			      
			    pressurebsetEn.evaluate(n,faceQuadInner,l,qu_);
			    pressurebsetNeigh.evaluate(n,faceQuadOuter,l,quneigh_);
			      
			    RangeType prflux_en(0.0);
			    RangeType prflux_nb(0.0);
			    PressureRangeType zeroq(0.0);
			 
			  
			    problem_.numericalFlux_p(nit,faceQuadInner,l,qu_,quneigh_,prflux_en,prflux_nb);
			    /////////////////////////////////////////////////////////////
			    //PGM
			    double PGM_en=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_en)
			      *faceQuadInner.weight(l)* massVolElInv*int_element;
			    double PGM_nb=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_nb)
			      *faceQuadInner.weight(l)* massVolElInv*int_element;
			    //////////////////////////////////////////////////////////////////


			    
			    en_pressGrad.add(j,n,PGM_en);
			    nb_pressGrad.add(j,n,PGM_nb);
			    
			    
			    
			    
			      
			    if(j==0)
			      {
			
				PressureRangeType pressValen(0.0);
				PressureRangeType pressValNeigh(0.0);
				  
				//\check{u}_p
			
				problem_.numericalFlux_p_stab(nit,faceQuadInner,l,qu_,quneigh_,pressValen,pressValNeigh);
				
				for(int m=0;m<pressureNumDofs;++m)
				  { 
		
				    ////////////////////////////////////////////////////////////////////////////
				    // -PSM
				    double PSM_en =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner,l,pressValen)
				      *faceQuadInner.weight(l)* massVolElInv*int_element;
				    double PSM_nb =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner,l,pressValNeigh)
				      *faceQuadOuter.weight(l)* massVolElInv*int_element;
				      	 
			
				    en_pressStab.add(m,n,PSM_en);
				    nb_pressStab.add(m,n,PSM_nb);
				  }
			      }
			  }
			 
			for(int i=0;i < numDofs;++i)
			  {
			  
			    RangeType tau(0.0);
			    bsetEn.evaluate(i,faceQuadInner,l,tau);
			      
			    //Stabiliztion
			    double VSM_en =stab_en*tau;
			    VSM_en*=-1.0;
			    VSM_en*=faceQuadInner.weight(l)* massVolElInv*int_element;
			
			    double VSM_nb =stab_nb*tau;
			    VSM_nb*=-1.0;
			    VSM_nb*=faceQuadInner.weight(l)* massVolElInv*int_element;
			  
			      
			    en_stab.add(i,j,VSM_en);
			    nb_stab.add(i,j,VSM_nb);
			  }		
		      }
		  }
	      }
	    }

	  //if(0)	 //fuer stabmatrix sigmaflux(0,0,phi,0) auswerten
	  if (nit.boundary()) 
	    {

	      RangeType zerophi(0.0); 
	      LocalVeloType localVelo=VeloRhs_.localFunction(en);
	      LocalStressType localStress =StressRhs_.localFunction(en);
	      LocalPressureType localPressure=PressureRhs_.localFunction(en);
	      
	      const int quadNop = faceQuadInner.nop();
	      for (int l = 0; l < quadNop ; ++l) 
		{
		  DomainType normal=nit.unitOuterNormal(faceQuadInner.localPoint(l));
		  double int_element=nit.intersectionGlobal().integrationElement(faceQuadInner.localPoint(l));
	
		  RangeType DirichletValue;
		  RangeType zero(0.0),dummy(0.0),stabflux(0.0);
		  
		  problem_.dirichletValues(nit,faceQuadInner,l,DirichletValue);
		  
		  //calculate C11*g_D
		  problem_.numericalFlux_u_stab(nit,faceQuadInner,l,DirichletValue,zero,stabflux,dummy);
		  // is to be added on the right hand side
		  stabflux*=-1.0;
		  
		  
		  
		  for(int j=0; j<numDofs; ++j)
		    { 
		      
		      bsetEn.evaluate(j,faceQuadInner,l,phi_);

		      //Boundary contribution for velo-RHS
		      localrhs[j]+=bsetEn.evaluateSingle(j,faceQuadInner,l,stabflux)
			*faceQuadInner.weight(l)*massVolElInv*int_element;
		      
		      for(int k=0;k<gradientNumDofs;++k)
			{
			  
			  gradbsetEn.evaluate(k,faceQuadInner,l,tau_);

			  RangeType sflux_bnd(0.0);
			  RangeType dummy(0.0);
			  GradientRangeType zero(0.0);

			  FMatrixHelp::multAssign<double,2,2>(tau_,normal,sflux_bnd);
			 
			  //VDM
			  double VDM_bnd=-1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,sflux_bnd)
			    *faceQuadInner.weight(l)*massVolElInv*int_element*nu_rt;
 			 
			 
			  GradientRangeType bnd_stress(0.0);
			  tensorProduct2d(DirichletValue,normal,bnd_stress);
			  
 
			  //Boundary contribution for stress-RHS
			  if(j==0)
			    {		
			      localStress[k]+=gradbsetEn.evaluateSingle(k,faceQuadInner,l,bnd_stress)
				  *faceQuadInner.weight(l)*massVolElInv*int_element*nu_rt;
			    }
			
			   
			}

		      for(int i=0;i<numDofs;++i)
			{
			 
			  bsetEn.evaluate(j,faceQuadInner,l,phi_); 	  
			  
			  RangeType stab_bnd(0.0);
			  problem_.numericalFlux_u_stab(nit,faceQuadInner,l,phi_,zerophi,stab_bnd,dummy);
			  

			 
			  double VSM_bnd=-1.0*bsetEn.evaluateSingle(i,faceQuadInner,l,stab_bnd)
			    *faceQuadInner.weight(l)*massVolElInv*int_element;
			  
			  en_stab.add(i,j,VSM_bnd);
			
		
			}
		      for(int n=0;n<pressureNumDofs;++n)
			{
			
			  pressurebsetEn.evaluate(n,faceQuadInner,l,qu_);
			  
			  RangeType prflux_bnd(0.0);
			
			  prflux_bnd=normal;
			  prflux_bnd*=qu_;
			  
			  
			  double PGM_bnd=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_bnd)
			    *faceQuadInner.weight(l)*massVolElInv*int_element;
		
			      en_pressGrad.add(j,n,PGM_bnd);
		
			  PressureRangeType bnd_pressure=DirichletValue*normal;
			
			  //Boundary contribution for pressure-RHS
			  //-H implementieren siehe Uzawa-Algorithmus 
			  if(j==0)
			    { 
			      if(nit.boundaryId()==1 || nit.boundaryId()==2)
				{	 
				  localPressure[n]+=1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner,l,bnd_pressure)
				    *faceQuadInner.weight(l)*massVolElInv*int_element;
				}
			    }
			}
		    }
		}	
	    } // end if boundary
      	}      
      
    
    }

    
    DestinationType& rhs1() const
    {
      return *dest_;
    }
    
    DiscretePressureType& rhs2()
    {
      return PressureRhs_;
    }
  
    DiscreteStressType& stress() const
    {
      return StressRhs_;
    }
    


  private:
    StokesFEPass();
    StokesFEPass(const StokesFEPass&);

  private:
    double volumeElement(const EntityType& en,
                         const VolumeQuadratureType& quad, 
                         double & massVolInv ) const {
      double result = 0.0;
      massVolInv = 0.0;
      const int quadNop = quad.nop();
      for (int qp = 0; qp < quadNop; ++qp) {
        massVolInv += quad.weight(qp);//volumen referenzelement
        result += 
          quad.weight(qp) * en.geometry().integrationElement(quad.point(qp));
      }
      massVolInv /= result;
      //   massVolInv = 1.0;
      return result;    }





    inline  void tensorProduct2d(const RangeType& v,const DomainType& n, GradientRangeType& ret) const
    {
      ret(0,0)= v[0]*n[0];
      ret(0,1)= v[0]*n[1];
      ret(1,0)= v[1]*n[0];
      ret(1,1)= v[1]*n[1];
    }

    void tensorProduct3d(const RangeType& v,const DomainType& n, GradientRangeType& ret)
    {
      ret[0]= v[0]*n[0];
      ret[1]= v[0]*n[1];
      ret[2]= v[0]*n[2];
      ret[3]= v[1]*n[0];
      ret[4]= v[1]*n[1];
      ret[5]= v[1]*n[2];
      ret[6]= v[2]*n[0];
      ret[7]= v[2]*n[1];
      ret[8]= v[2]*n[2];
    }
   

    
    
    
  private:
    // mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    VelocityDiscreteFunctionSpaceType& spc_;
    mutable DiscreteGradientSpaceType gradientspc_;
    mutable DiscretePressureSpaceType pressurespc_;
    mutable DiscreteStressType tmp_; 
    mutable DiscreteStressType tmp2_; 
    mutable DestinationType tmp1_;
    mutable double dtMin_;

    //Rhs functions for Dirichlet Valus

    mutable DestinationType VeloRhs_;
    mutable DiscreteStressType StressRhs_;
    mutable DiscretePressureType PressureRhs_;

    //! Some helper variables
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType phi_;
    mutable RangeType phiNeigh_;
    mutable PressureJacobianRangeType psi_;
    mutable GradientRangeType tau_;
    mutable GradientRangeType tauNeigh_;
    mutable RangeType gradEval_;
    mutable PressureRangeType qu_;
    mutable PressureRangeType quneigh_;
    mutable DomainType grads_;
    TimeProvider* time_;
    FieldVector<int, 0> diffVar_;
  

    //neue Variablen
    // mutable GridPartType& gridPrt_; 
  

   
    int volumeQuadOrd_,faceQuadOrd_;

    int maxNumberUnknowns_;
    int maxNum_;
    
  
    mutable GradMatType gradMatrix_;
    mutable StabMatType stabMatrix_;
    mutable PressureGradMatType pressureGradMatrix_;
    mutable PressureStabMatType pressureStabMatrix_;

    mutable bool matrixAssembled_;

    mutable DomainType direction_;
    
   
  };
  
} // end namespace Dune

#endif
