#ifndef DUNE_STOKESPASS_HH
#define DUNE_STOKESPASS_HH

/** \file
    \brief Stokespass.hh
 */

#include <dune/fem/pass/pass.hh>
// #include <dune/fem/pass/selection.hh>
// #include <dune/fem/pass/discretemodel.hh>
// #include <dune/fem/pass/modelcaller.hh>

#include "matrixoperator.hh"
#include <dune/fem/operator/2order/dgmatrixsetup.hh>

// * needs to move
// #include "../misc/timenew.hh"
#include <dune/fem/solver/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  //
  template <class DiscreteModelImp, class PreviousPassImp,int passIdImp>
  class StokesFEPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp,passIdImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp,passIdImp> BaseType;
    typedef StokesFEPass<DiscreteModelImp,PreviousPassImp,passIdImp> ThisType;

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
   //  typedef typename DiscreteModelType::SelectorType SelectorType;
    // typedef DiscreteModelCaller<
//       DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;
   
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
    // typedef FieldMatrix<RangeFieldType,dimDomain,dimDomain> HelperMatrixType;
    
    typedef MatrixOperator<MatrixType,PressureRangeType,RangeType> BOPType;
    typedef MatrixOperator<MatrixType,RangeType,PressureRangeType> BTOPType;
    typedef MatrixOperator<MatrixType,PressureRangeType,PressureRangeType> COPType;
   
    typedef SparseRowMatrixObject<DiscreteGradientSpaceType,VelocityDiscreteFunctionSpaceType> GradMatType;
    typedef typename GradMatType::LocalMatrixType LocalGradMatType;
    typedef typename GradMatType::MatrixType  GMType;
   
    typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,DiscreteGradientSpaceType> DivMatType;
    typedef typename DivMatType::LocalMatrixType LocalDivMatType;
    typedef typename DivMatType::MatrixType  DMType;

    typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> StabMatType;
    typedef typename StabMatType::LocalMatrixType LocalStabMatType;
    typedef typename StabMatType::MatrixType SMType;

   
    typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,DiscretePressureSpaceType> PressureGradMatType;
    typedef typename PressureGradMatType::LocalMatrixType LocalPressureGradMatType;
    typedef typename PressureGradMatType::MatrixType  PGMType;
   
    typedef SparseRowMatrixObject<DiscretePressureSpaceType,VelocityDiscreteFunctionSpaceType> PressureDivMatType;
    typedef typename PressureDivMatType::LocalMatrixType LocalPressureDivMatType;
    typedef typename PressureDivMatType::MatrixType  PDMType;
   
   
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
	// 	 double c11,double c12,double d11,double d12,
		 int volumeQuadOrd =-1,int faceQuadOrd=-1) :
      BaseType(pass, spc),
      //      caller_(problem),
      problem_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gradientspc_(spc_.gridPart()),
      pressurespc_(spc_.gridPart()),
      tmp_("FEPass::tmp",gradientspc_),
      tmp2_("FEPass::tmp2",spc_),
      tmp1_("FEPass::tmp1",spc_),
      dtMin_(std::numeric_limits<double>::max()),
      VeloRhs_("VeloRhs",spc_),
      StressRhs_("StressRhs",gradientspc_),
      PressureRhs_("PressureRhs",pressurespc_),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0), 
      phi_(0.0),
      phiNeigh_(0.0),
      psi_(0.0),
      taujac_(0.0),
      tau_(0.0),
      tauNeigh_(0.0),
      grads_(0.0),
      time_(0),
      diffVar_(),
      twistUtil_(spc.grid()),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (5*spc_.order()+1) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (5*spc_.order()+1) : faceQuadOrd ),
      maxNumberUnknowns_(5* (gradientspc_.baseFunctionSet(*(spc_.begin())).numBaseFunctions())),
      maxNum_(5*(spc_.baseFunctionSet(*(spc_.begin())).numBaseFunctions())),
      //  massMatrix_(gradientspc_.size(),gradientspc_.size(),maxNumberUnknowns_,0.0),
      //       massMatrixInv_(gradientspc_.size(),gradientspc_.size(),maxNumberUnknowns_,0.0),
      gradMatrix_(gradientspc_,spc_,""),//VDM
      divMatrix_(spc_,gradientspc_,""),//VGM
      stabMatrix_(spc_,spc_,""),//VSM
      pressureGradMatrix_(spc_,pressurespc_,""),//PGM
      pressureDivMatrix_(pressurespc_,spc_,""),//PDM
      pressureStabMatrix_(pressurespc_,pressurespc_,""),//PSM
    //   c11_(c11),
//       c12_(c12),
//       d11_(d11),
//       d12_(d12),
      matrixAssembled_(false)
    //   direction_(0.0)
    {
      
      std::cout<<"unknowms:"<<maxNumberUnknowns_<<"\n";
    //   direction_[0]=1;
//       direction_[1]=(1.0+sqrt(5.0))/2.0;
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );
    }
    
    //! Destructor
    virtual ~StokesFEPass() {
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void setTime(const double time) {
      time_ = time;
    }

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return 0.0;
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
      StressRhs_.clear();
      PressureRhs_.clear();
     
      ElementAndNeighbors stencil;

      gradMatrix_.reserve();
      gradMatrix_.clear();
      
      divMatrix_.reserve();
      divMatrix_.clear();
      
      stabMatrix_.reserve();
      stabMatrix_.clear();
      
      
      pressureGradMatrix_.reserve();
      pressureGradMatrix_.clear();
      pressureDivMatrix_.reserve();
      pressureDivMatrix_.clear();
      
      pressureStabMatrix_.reserve();
      pressureStabMatrix_.clear();
     
      this->operator()(rhs,rhs);
      matrixAssembled_ = true;
     
      double* stressptr=StressRhs_.leakPointer();
      double* veloptr=VeloRhs_.leakPointer();
      double* tmp1ptr=tmp1_.leakPointer();
      double* tmp2ptr=tmp2_.leakPointer();

      divMatrix_.multOEM(stressptr,tmp1ptr);
      //gradMatrix_.multOEM_t(stressptr,tmp1ptr);
      
      for(register int i=0;i<spc_.size();++i)
	{
	  veloptr[i]-=tmp1ptr[i];
	}
      //  gradMatrix_.multOEM_t(stressptr,tmp1ptr);
  
    }

    
    MatrixType& getBOP()const {return pressureGradMatrix_.matrix();}
    MatrixType& getBTOP() const {return pressureDivMatrix_.matrix();}
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
	
 	gradMatrix_.multOEM(arg,tmpP);
	divMatrix_.multOEM(tmpP,tmpP1);
	
       
	stabMatrix_.multOEM(arg,dest);
	

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
   //    if (time_) {
//         time_->provideTimeStepEstimate(dtMin_);
      //   }

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
      typedef typename VelocityDiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename VelocityDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      
      typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
      typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
      
      typedef typename DiscretePressureSpaceType::IndexSetType PressureIndexSetType;
      typedef typename DiscretePressureSpaceType::BaseFunctionSetType PressureBaseFunctionSetType;
      //- statements

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
      LocalDivMatType en_div= divMatrix_.localMatrix(en,en);
      LocalStabMatType en_stab=stabMatrix_.localMatrix(en,en);
      LocalPressureGradMatType en_pressGrad=pressureGradMatrix_.localMatrix(en,en);
      LocalPressureDivMatType en_pressDiv=pressureDivMatrix_.localMatrix(en,en);
      LocalPressureStabMatType en_pressStab=pressureStabMatrix_.localMatrix(en,en);
       double nu=problem_.get_nu();
      
      nu=1;      // Volumetric integral part
      const int quadNop = volQuad.nop();
      GradJacobianRangeType helpmatr(0.0);
  
      JacobianRangeType helpmatr2(0.0); 
      
      //if(0)     
      for (int l = 0; l < quadNop ; ++l) 
	{
	  GradJacobianRangeType helpmatr(0.0);
	  JacobianRangeType helpmatr2(0.0); 
	  
	  double intel = volQuad.weight(l);    
	  double intelement =en.geometry().integrationElement(volQuad.point(l));  

	  for(int k = 0;k < gradientNumDofs ;++k)
	    {
	      int VGM_rows  = gradientspc_.mapToGlobal(en,k);
	      int VDM_cols = VGM_rows;
	      
	      // eval tau_k

	      gradbsetEn.evaluate(k, volQuad[l],tau_);
	      
	      for (int j = 0; j < numDofs; ++j) 
		{  
		  bsetEn.evaluate(j, volQuad[l], psi_[0]);
		   

		  if(k==0)
		    {
		      //Rhs Assembling
		      RangeType fval(0.0);
		      problem_.rightHandSide(en,volQuad,l,fval); 
		      localrhs[j]+=bsetEn.evaluateSingle(j,volQuad[l],fval)*intel*intelement;
		    }
		  
		  int VGM_cols = spc_.mapToGlobal(en,j);
		  int VDM_rows = VGM_cols;
		  
		  int PGM_rows = spc_.mapToGlobal(en,j);
		  int PDM_cols = PGM_rows;
		  
		  // needed for calculating phi_j*div tau_k with the evaluateGradientSingle method
		  helpmatr*=0.0;
		  
		  helpmatr[0][0]=psi_[0][0];
		  helpmatr[1][1]=psi_[0][0];
		  helpmatr[2][0]=psi_[0][1];
		  helpmatr[3][1]=psi_[0][1];
		  
		  
		  // eval -phi_j*div tau_k
		  // -VGM
		  double VGM =-1.0* gradbsetEn.evaluateGradientSingle(k,en,volQuad[l],helpmatr)*intel*nu;
		  
		  
		  // eval tau_k:grad phi_j
		  // VDM
		  //double VDM = bsetEn.evaluateGradientSingle(j,en,volQuad,l,tau_mat)*intel;
		  double VDM = bsetEn.evaluateGradientSingle(j,en,volQuad[l],tau_)*intel*intelement;
		  

		  en_div.add(j,k,VDM);
		  en_grad.add(k,j,VGM);

		  ////////////////////////////////////////////////
		  //pressureGradientMatrix , pressureDivMatrix und Convectionsanteil der Stabmatrix  //
		  ////////////////////////////////////////////////
		  
		  
		  if(k==0)
		    {
		   

		      for(int n = 0; n < pressureNumDofs ; ++n)
			{ //eval q_n
			  pressurebsetEn.evaluate(n,volQuad[l],qu_);
			
			  helpmatr2*=0.0;
			  for(int m = 0; m < dimDomain ; ++m)
			    {
			      helpmatr2[m][m]=qu_[0];
			    }
			  int PGM_cols  =  pressurespc_.mapToGlobal(en,n);  
			  int PDM_rows  =  PGM_cols;
			
			  //eval (-q_n*div phi_j)
			  double PGM =-1.0*bsetEn.evaluateGradientSingle(j,en, volQuad[l],helpmatr2)*intel*intelement; 
			  // std::cout<<"PGM="<< PGM<<"\n";
			  //eval -(-phi_j*grad q_n)
			  double PDM =1.0*pressurebsetEn.evaluateGradientSingle(n,en,volQuad[l],psi_)*intel*intelement;
			  //	    std::cout<<"PDM="<< PDM<<"\n"; 
			  en_pressGrad.add(j,n,PGM);
			  en_pressDiv.add(n,j,PDM);
		
			}
		    }    
		}
	      
	    }
	}
	
      // Surface integral part
     
      IntersectionIterator endnit = en.ileafend();

      double dtLocal = 0.0;
      double minvol = vol; 
      //   int counter=0;
      //if(0)
      for (IntersectionIterator nit = en.ileafbegin(); nit != endnit; ++nit) 
	{ 
	  //	  int twistSelf = twistUtil_.twistInSelf(nit); 
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
		
		 LocalGradMatType nb_grad = gradMatrix_.localMatrix(en,nb);
		 LocalDivMatType nb_div= divMatrix_.localMatrix(en,nb);
		 LocalStabMatType nb_stab=stabMatrix_.localMatrix(en,nb);
		 LocalPressureGradMatType nb_pressGrad=pressureGradMatrix_.localMatrix(en,nb);
		 LocalPressureDivMatType nb_pressDiv=pressureDivMatrix_.localMatrix(en,nb);
		 LocalPressureStabMatType nb_pressStab=pressureStabMatrix_.localMatrix(en,nb);





  
		double massVolNbInv;
		double nbvol = volumeElement(nb, volQuad, massVolNbInv);
		if (nbvol<minvol) minvol=nbvol;

		const int quadNop = faceQuadInner.nop();
		for (int l = 0; l < quadNop ; ++l) 
		  {
		    //We use int_elemt seperately,which makes the flux-implemantion more readable and flexible
		    //This should be changed in serious applications for efficiency
		    DomainType normal=nit.unitOuterNormal(faceQuadInner.localPoint(l));
		    double int_element = nit.intersectionGlobal().integrationElement(faceQuadInner.localPoint(l));
		    
		    double hinv=1.0;
		    hinv/=int_element;
	      	
		    for(int j=0; j< numDofs; ++j)
		      { 
			//Note:The cols always depend on neighbour elements
			int VGM_cols_en=spc_.mapToGlobal(en,j);
			int VGM_cols_nb=spc_.mapToGlobal(nb,j);
			  
			int VSM_cols_en=VGM_cols_en;
			int VSM_cols_nb=VGM_cols_nb;
			  
			int PDM_cols_en=VGM_cols_en;
			int PDM_cols_nb=VGM_cols_nb;
			  
			int VDM_rows=VGM_cols_en;
			// int VSM_rows=VDM_rows;
			//			  int PGM_rows=VDM_rows;
			  
			//calculatin phi_j on the element- and neighbor side of the current intersection
			bsetEn.evaluate(j, faceQuadInner[l], phi_); 
		// 	DomainType x=nit.intersectionGlobal().global(faceQuadInner.localPoint(l));
			bsetNeigh.evaluate(j, faceQuadOuter[l], phiNeigh_ );
	


			//this stuff should be made member of the class
			RangeType zerophi(0.0); 
			GradientRangeType zerotau(0.0);
			GradientRangeType enflux(0.0);
			GradientRangeType neighflux(0.0);
			RangeType stab_en(0.0);
			RangeType stab_nb(0.0);
			RangeType conv_en(0.0);
			RangeType conv_nb(0.0);
			  
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
			    int VGM_rows = gradientspc_.mapToGlobal(en,k);

			    int VDM_cols_en =VGM_rows;
			      
			    int VDM_cols_nb =gradientspc_.mapToGlobal(nb,k);
			      
			    gradbsetEn.evaluate(k,faceQuadInner[l],tau_);
			    gradbsetNeigh.evaluate(k,faceQuadOuter[l],tauNeigh_);
			      
			    RangeType sflux_en(0.0);
			    RangeType sflux_nb(0.0);

			    
				
			    problem_.numericalFlux_sigma(nit,faceQuadInner,l,tau_,tauNeigh_,sflux_en,sflux_nb);
			   
			    double VGM_en =  gradbsetEn.evaluateSingle(k,faceQuadInner[l],enflux)
			      *faceQuadInner.weight(l)* massVolElInv*int_element*nu;
			    double VGM_nb = gradbsetEn.evaluateSingle(k,faceQuadInner[l],neighflux)
			      *faceQuadInner.weight(l)* massVolElInv*int_element*nu;
			    
;
			    en_grad.add(k,j,VGM_en);
			    nb_grad.add(k,j,VGM_nb);
			    
			    double VDM_en=-1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],sflux_en)
			      *faceQuadInner.weight(l)*/* massVolElInv*/int_element;
			    double VDM_nb=-1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],sflux_nb)
			      *faceQuadInner.weight(l)*/* massVolElInv*/int_element;
			    en_div.add(j,k,VDM_en);
			    nb_div.add(j,k,VDM_nb);
			    
			 
			  }
			 
			for(int n=0;n<pressureNumDofs;++n)
			  {
			    int PGM_cols_en = pressurespc_.mapToGlobal(en,n);
			    int PGM_cols_nb = pressurespc_.mapToGlobal(nb,n);
			    int PDM_rows=PGM_cols_en;
			    int PGM_rows = VDM_rows;
			      
			     
			     
			      
			    pressurebsetEn.evaluate(n,faceQuadInner[l],qu_);
			    pressurebsetNeigh.evaluate(n,faceQuadOuter[l],quneigh_);
			      
			    RangeType prflux_en(0.0);
			    RangeType prflux_nb(0.0);
			    PressureRangeType zeroq(0.0);
			 
			    // \hat{p}
			 
			    problem_.numericalFlux_p(nit,faceQuadInner,l,qu_,quneigh_,prflux_en,prflux_nb);
			    /////////////////////////////////////////////////////////////
			    //PGM
			    double PGM_en=1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],prflux_en)
			      *faceQuadInner.weight(l)/* massVolElInv*/*int_element;
			    double PGM_nb=1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],prflux_nb)
			      *faceQuadInner.weight(l)/* massVolElInv*/*int_element;
			    //////////////////////////////////////////////////////////////////
 
			    en_pressGrad.add(j,n,PGM_en);
			    nb_pressGrad.add(j,n,PGM_nb);
			    
		
			      
			    ///////////////////////////////////////////////////////////////////
			    //-PDM 
			    double PDM_en =-1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner[l],pressure_en)
			      *faceQuadInner.weight(l)/* massVolElInv*/*int_element;
			    double PDM_nb =-1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner[l],pressure_neigh)
			      *faceQuadInner.weight(l)/* massVolElInv*/*int_element;

			    ///////////////////////////////////////////////////////////////////////////
			  //   pressureDivMatrix_.add(PDM_rows,PDM_cols_en,PDM_en);
// 			    pressureDivMatrix_.add(PDM_rows,PDM_cols_nb,PDM_nb);
			    
			    en_pressDiv.add(n,j,PDM_en);
			    nb_pressDiv.add(n,j,PDM_nb);
			      
			    if(j==0)
			      {
				// int PGM_en = pressurespc_.mapToGlobal(en,j);
				PressureRangeType pressValen(0.0);
				PressureRangeType pressValNeigh(0.0);
				  
				//\check{u}_p
				problem_.numericalFlux_p_stab(nit,faceQuadInner,l,qu_,quneigh_,pressValen,pressValNeigh);
				
				for(int m=0;m<pressureNumDofs;++m)
				  { 
				    int PSM_rows = pressurespc_.mapToGlobal(en,m);
				    int PSM_cols_en= PGM_cols_en;
				    int PSM_cols_nb= PGM_cols_nb;
				    ////////////////////////////////////////////////////////////////////////////
				    // -PSM
				    double PSM_en =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner[l],pressValen)
				      *faceQuadInner.weight(l)/* massVolElInv*/*int_element;
				    double PSM_nb =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner[l],pressValNeigh)
				      *faceQuadOuter.weight(l)/* massVolElInv*/*int_element;
				    en_pressStab.add(m,n,PSM_en);
				    nb_pressStab.add(m,n,PSM_nb);
				   
				  }
			      }
			  }
			 
			for(int i=0;i < numDofs;++i)
			  {
			    int VSM_rows=spc_.mapToGlobal(en,i);
			      
			    RangeType tau(0.0);
			    bsetEn.evaluate(i,faceQuadInner[l],tau);
			      
			    //Stabiliztion
			    double VSM_en =stab_en*tau;
			    VSM_en*=-1.0;
			    VSM_en*=faceQuadInner.weight(l)/* massVolElInv*/*int_element;
			
			    double VSM_nb =stab_nb*tau;
			    VSM_nb*=-1.0;
			    VSM_nb*=faceQuadInner.weight(l)/* massVolElInv*/*int_element;
			  
			       
			    en_stab.add(i,j,VSM_en);
			    nb_stab.add(i,j,VSM_nb);
			      
// 			   
  
			  }		
		  

		      }




		  }
		
	      }
	    
	    }
	  
	  if (nit.boundary()) 
	    {

	      RangeType zerophi(0.0); 
	      LocalVeloType localVelo=VeloRhs_.localFunction(en);
	      LocalStressType localStress =StressRhs_.localFunction(en);
	      LocalPressureType localPressure=PressureRhs_.localFunction(en);
	      double* Veloptr =VeloRhs_.leakPointer();
	      double* Stressptr =StressRhs_.leakPointer();
	      double* Pressptr =PressureRhs_.leakPointer();
	      const int quadNop = faceQuadInner.nop();
	      for (int l = 0; l < quadNop ; ++l) 
		{
		  DomainType normal=nit.unitOuterNormal(faceQuadInner.localPoint(l));
		  double int_element=nit.intersectionGlobal().integrationElement(faceQuadInner.localPoint(l));
		 
		  double hinv=1;
		  hinv/=int_element;
		  RangeType DirichletValue;
		  RangeType zero(0.0),dummy(0.0),stabflux(0.0);
	
		  problem_.dirichletValues(nit,faceQuadInner,l,DirichletValue);
		  problem_.numericalFlux_u_stab(nit,faceQuadInner,l,DirichletValue,zero,stabflux,dummy);
// 		  // is to be added on the right hand side
 		  stabflux*=-1.0;
		  
		  
		  for(int j=0; j<numDofs; ++j)
		    { 
		      int VDM_rows=spc_.mapToGlobal(en,j);
		      int VSM_rows=VDM_rows;
		      int VGM_cols_bnd=VDM_rows;
		      //int VDM_cols_bnd = StressRhs_rows;

		      bsetEn.evaluate(j,faceQuadInner[l],phi_);

		     
			  Veloptr[VDM_rows]+=bsetEn.evaluateSingle(j,faceQuadInner[l],stabflux)
			    *faceQuadInner.weight(l)/*massVolElInv*/*int_element;
			 
		      for(int k=0;k<gradientNumDofs;++k)
			{
			  int VDM_cols_bnd=gradientspc_.mapToGlobal(en,k);
			  int VGM_rows_bnd=VDM_cols_bnd;
			  int Rhsrow=VDM_cols_bnd;
			  gradbsetEn.evaluate(k,faceQuadInner[l],tau_);

			  RangeType sflux_bnd(0.0);RangeType dummy(0.0);
			  GradientRangeType zero(0.0);

			     
			  multi2d(tau_,normal,sflux_bnd);
			 
			     double VDM_bnd=-1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],sflux_bnd)
			       *faceQuadInner.weight(l)/*massVolElInv*/*int_element;
			     en_div.add(j,k,VDM_bnd);
			     
			  GradientRangeType bnd_stress(0.0);
			  tensorProduct2d(DirichletValue,normal,bnd_stress);
			  GradientRangeType natval(0.0);
			  tensorProduct2d(phi_,normal,natval);

			  //Boundary contribution for stress-RHS
			  if(j==0)
			  
				Stressptr[Rhsrow]+=gradbsetEn.evaluateSingle(k,faceQuadInner[l],bnd_stress)
				  *faceQuadInner.weight(l)*massVolElInv*int_element*nu;
			  
			}

		      for(int i=0;i<numDofs;++i)
			{
			  int VSM_cols_bnd=spc_.mapToGlobal(en,i);

			  bsetEn.evaluate(j,faceQuadInner[l],phi_); 	  
			  
			  RangeType stab_bnd(0.0);
			  RangeType conv_bnd(0.0);
		
			  problem_.numericalFlux_u_stab(nit,faceQuadInner,l,phi_,zerophi,stab_bnd,dummy);

			  double VSM_bnd=-1.0*bsetEn.evaluateSingle(i,faceQuadInner[l],stab_bnd)
			    *faceQuadInner.weight(l)/*massVolElInv*/*int_element;
				  
			      
		
			      
			      en_stab.add(i,j,VSM_bnd);
			      
		
			  
			}
		      for(int n=0;n<pressureNumDofs;++n)
			{
			  int PGM_rows=spc_.mapToGlobal(en,j);
			  int PGM_cols_bnd=pressurespc_.mapToGlobal(en,n);
			  int PSM_row= PGM_cols_bnd;
			  pressurebsetEn.evaluate(n,faceQuadInner[l],qu_);
			  
			  RangeType prflux_bnd(0.0);
			
			  prflux_bnd=normal;
			  prflux_bnd*=qu_;
			  
			 
			  
			  double PGM_bnd=1.0*bsetEn.evaluateSingle(j,faceQuadInner[l],prflux_bnd)
			    *faceQuadInner.weight(l)/*massVolElInv*/*int_element;
			  
			  //PGM_bnd*=-1.0;
			  
			  
			  en_pressGrad.add(j,n,PGM_bnd);
			     

			  PressureRangeType bnd_pressure=DirichletValue*normal;
			
			  //Boundary contribution for pressure-RHS
			  //-H implementieren siehe Uzawa-Algorithmus 
			  if(j==0)
			    { 
			      if(nit.boundaryId()==1 || nit.boundaryId()==2)
				{
				  Pressptr[PSM_row]+=1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner[l],bnd_pressure)
				    *faceQuadInner.weight(l)/*massVolElInv*/*int_element;
				}
			    }
			}
		    }
		}	
	    } // end if boundary
      	}      

      if (dtLocal>2.*std::numeric_limits<double>::min()) 
	{
	  dtLocal = minvol/dtLocal;
	  if (dtLocal < dtMin_) dtMin_ = dtLocal;
	}
    }

    
    DestinationType& rhs1() const
    {
      return VeloRhs_;
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
      //    massVolInv = 1.0;
      return result;    }
  





    void multi2d(const GradientRangeType& tau,const DomainType& n,RangeType& ret) const
    { 
      ret[0]=tau[0]*n[0]+tau[1]*n[1];
      ret[1]=tau[2]*n[0]+tau[3]*n[1];
    }

    //converts a GradientRangeType to a JaobianRangeType:FielVector<n*n>->FieldMatrix<n,n>
    void convertutil2d(const GradientRangeType& tau,JacobianRangeType& mat)const
    { assert(false);
      mat[0][0]=tau[0];
      mat[0][1]=tau[1]; 
      mat[1][0]=tau[2];
      mat[1][1]=tau[3];
    }
    
   
    void tensorProduct2d(const RangeType& v,const DomainType& n, GradientRangeType& ret) const
    { 
      ret(0,0)= v[0]*n[0];
      ret(0,1)= v[0]*n[1];
      ret(1,0)= v[1]*n[0];
      ret(1,1)= v[1]*n[1];
     
  
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
    mutable DestinationType tmp2_;
    mutable DestinationType tmp1_;
    mutable double dtMin_;

    //Rhs functions from Dirichlet Valus

    mutable DestinationType VeloRhs_;
    mutable DiscreteStressType StressRhs_;
  mutable DiscretePressureType PressureRhs_;

    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType phi_;
    mutable RangeType phiNeigh_;
    mutable PressureJacobianRangeType psi_;
    mutable JacobianRangeType taujac_;
    mutable GradientRangeType tau_;
    mutable GradientRangeType tauNeigh_;
    mutable RangeType gradEval_;
    mutable PressureRangeType qu_;
    mutable PressureRangeType quneigh_;
    mutable DomainType grads_;
    double time_;
    FieldVector<int, 0> diffVar_;
  

    //neue Variablen
    // mutable GridPartType& gridPrt_; 
  

    TwistUtility<GridType> twistUtil_;

    int volumeQuadOrd_,faceQuadOrd_;

    int maxNumberUnknowns_;
    int maxNum_;
    
    //    mutable MatrixType massMatrix_ ,massMatrixInv_;

    mutable GradMatType gradMatrix_;
    mutable DivMatType divMatrix_;
    mutable StabMatType stabMatrix_;
    mutable PressureGradMatType pressureGradMatrix_;
    mutable PressureDivMatType pressureDivMatrix_;
    mutable PressureStabMatType pressureStabMatrix_;

//     mutable MatrixType gradMatrix_,divMatrix_,stabMatrix_;
//     mutable MatrixType pressureGradMatrix_,pressureDivMatrix_,pressureStabMatrix_;
//     mutable double c11_;
//     mutable double c12_;
//     mutable double d11_;
//     mutable double d12_;
    mutable bool matrixAssembled_;

    mutable DomainType direction_;
    
    //     mutable BOPType BOP_;
    //     mutable BTOPType BTOP_;
    //     mutable COPType COP_;

};
  
} // end namespace Dune

#endif
