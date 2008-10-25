#ifndef DUNE_STOKESPASS_HH
#define DUNE_STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/modelcaller.hh>
#include <dune/fem/operator/matrixoperator.hh>








// * needs to move
// #include "../misc/timenew.hh"
#include <dune/fem/misc/timeutility.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  //
  template <class DiscreteModelImp, class PreviousPassImp,class MatroxObjectImp>
  class StokesLDGoperator :
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
    
    
    typedef MatrixObjectImp MatrixHandelType;

    typedef typename MatrixHandlerType::MatrixAddHandleType MatrixAddHandleType;
    typedef typename MatrixHandlerType::MatrixType MatrixType;
    typedef typename MatrixHandlerType::PreconditionMatrixType PreconditionMatrixType;


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
    
    typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,dimGradRange> GradientSpaceType;
    typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd> DiscreteGradientSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd> DiscreteGradientSpaceType;
    typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,1> PressureSpaceType;
    typedef DiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd> DiscretePressureSpaceType;
    // typedef LegendreDiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd> DiscretePressureSpaceType;
    
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
		 double c11,double c12,double d11,double d12,
		int volumeQuadOrd =-1,int faceQuadOrd=-1) :
      BaseType(pass, spc),
      caller_(problem),
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
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (2*spc_.order()+1) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (2*spc_.order()+1) : faceQuadOrd ),
      maxNumberUnknowns_(20* (gradientspc_.baseFunctionSet(*(spc_.begin())).numBaseFunctions())),
      massMatrix_(gradientspc_.size(),gradientspc_.size(),maxNumberUnknowns_,0.0),
      massMatrixInv_(gradientspc_.size(),gradientspc_.size(),maxNumberUnknowns_,0.0),
      gradMatrix_(gradientspc_.size(),spc_.size(),maxNumberUnknowns_,0.0),//VDM
      divMatrix_(spc_.size(),gradientspc_.size(),maxNumberUnknowns_,0.0),//VGM
      stabMatrix_(spc_.size(),spc_.size(),maxNumberUnknowns_,0.0),//VSM
      pressureGradMatrix_(spc_.size(),pressurespc_.size(),maxNumberUnknowns_,0.0),//PGM
      pressureDivMatrix_(pressurespc_.size(),spc_.size(),maxNumberUnknowns_,0.0),//PDM
      pressureStabMatrix_(pressurespc_.size(),pressurespc_.size(),maxNumberUnknowns_,0.0),//PSM
      c11_(c11),
      c12_(c12),
      d11_(d11),
      d12_(d12),
      matrixAssembled_(false)
    {
      

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
      StressRhs_.clear();
      PressureRhs_.clear();


      gradMatrix_.resize(gradSize,singleSize);
      gradMatrix_.clear();
        
      massMatrix_.resize(gradSize,gradSize);
      massMatrix_.clear();
	
      divMatrix_.resize(singleSize,gradSize);
      divMatrix_.clear();
	
      stabMatrix_.resize(singleSize,singleSize);
      stabMatrix_.clear();

   //    massMatrixInv_.resize(gradientspc_.size());
//       massMatrixInv_.clear();

      

      pressureGradMatrix_.resize(singleSize,pressureSize);
      pressureGradMatrix_.clear();


      pressureDivMatrix_.resize(pressureSize,singleSize);
      pressureDivMatrix_.clear();



    pressureStabMatrix_.resize(pressureSize,pressureSize);
    pressureStabMatrix_.clear();



  this->operator()(rhs,rhs);
      matrixAssembled_ = true;
     
      double* stressptr=StressRhs_.leakPointer();
      double* veloptr=VeloRhs_.leakPointer();
      double* tmp1ptr=tmp1_.leakPointer();
     
      divMatrix_.multOEM(stressptr,tmp1ptr);
      for(register int i=0;i<spc_.size();++i)
	{
	  veloptr[i]-=tmp1ptr[i];
	}
	//   VeloRhs_-=tmp1_;
       
    
      
//       std::cout<<"Stress\n";
//       StressRhs_.print(std::cout);
//       std::cout<<"Pressure\n";
//       PressureRhs_.print(std::cout);
   //    BOP_(pressureGradMatrix_);
//       BTOP_(pressureDivMatrix_);
//       COP_(pressureStabMatrix_);
  //     std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:pressuregradmatrix=\n"; 
//       pressureGradMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";


//       std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:pressuredivmatrix=\n"; 
//       pressureDivMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";


  //     std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:pressurestabmatrix=\n"; 
//       pressureStabMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";
//       std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:stabmatrix=\n"; 
//       stabMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";
    //   std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:gradmatrix=\n"; 
//       gradMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";

//       std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:divmatrix=\n"; 
//       divMatrix_.print(std::cout);
//       std::cout<<"************************************************\n";
      
//       int foo;
//       std::cin>> foo;
//       assert(foo==1);




      // for(int i=0; i<gradSize; ++i) 
//       {
//         double val = massMatrix_(i,i);
//         if( std::abs(val) > 0.0 ) massMatrixInv_.add(i,i,1.0/val);
//       }

      

      
      //std::cout<<"massmatrix\n";
      //massMatrixInv_.print(std::cout);
      //assert(false);
    }

    
    MatrixType& getBOP()const {return pressureGradMatrix_;}
    MatrixType& getBTOP() const {return pressureDivMatrix_;}
    MatrixType& getCOP() const {return pressureStabMatrix_;}

    DiscretePressureSpaceType& getScalSpc() const {return pressurespc_;}
    VelocityDiscreteFunctionSpaceType getVecSpc() const {return spc_;}
    


    void multOEM(const double * arg, double * dest) const
    {
      //  if(!matrixAssembled_) 
	{

	  const int size = spc_.size();
// 	  std::cout << "Start multOEM \n"; 
 	//   std::cout<<"************************************************\n";
// 	  std::cout<<"Stokespass:ARGUMENT of multoem\n";
// 	  for(register int i=0; i<size; ++i)  
// 	    {
// 	      std::cout<<arg[i]<<"\n";}
// 	  std::cout<<"************************************************\n";
	  tmp_.clear();
      tmp1_.clear();
      tmp2_.clear();
      
      double * tmpP  = tmp_.leakPointer();
      //      double * tmp2P = tmp2_.leakPointer();
      double * tmpP1 = tmp1_.leakPointer();
      
      gradMatrix_.multOEM(arg,tmpP);

      divMatrix_.multOEM(tmpP,tmpP1);
      // std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:graddiv of multoem\n";
//       for(register int i=0; i<size; ++i)  
// 	{
// 	  std::cout<<tmpP1[i]<<"\n";}
//       std::cout<<"************************************************\n";
     
      stabMatrix_.multOEM(arg,dest);
     //  std::cout<<"************************************************\n";
//       std::cout<<"Stokespass:Result of multoem\n";
//       for(register int i=0; i<size; ++i)  
// 	{
// 	  std::cout<<dest[i]<<"\n";}
//       std::cout<<"************************************************\n";
      for(register int i=0; i<size; ++i) dest[i] += tmpP1[i];
      
      //   std::cout << "End multOEM \n";  
     
   //    int foo;
//       std::cin>> foo;
//       assert(foo==1);
      
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
 
      caller_.setArgument(*arg_);

      // time initialisation
      dtMin_ = std::numeric_limits<double>::max();
      if (time_) {
        caller_.setTime(time_->time());
      }
      else {
        caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      if (time_) {
        time_->provideTimeStepEstimate(dtMin_);
      }

      caller_.finalize();
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
      caller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      double massVolElInv;
      double vol = volumeElement(en, volQuad,massVolElInv);
      
      const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      
      const GradientBaseFunctionSetType& gradbsetEn = gradientspc_.baseFunctionSet(en);
      const int gradientNumDofs = gradbsetEn.numBaseFunctions();
      
      const PressureBaseFunctionSetType& pressurebsetEn = pressurespc_.baseFunctionSet(en);
      const int pressureNumDofs = pressurebsetEn.numBaseFunctions();


      

      // Volumetric integral part
      const int quadNop = volQuad.nop();
      GradJacobianRangeType helpmatr(0.0);
      JacobianRangeType helpmatr2(0.0); 
      //    if(0)     
      for (int l = 0; l < quadNop ; ++l) 
	{
	  GradJacobianRangeType helpmatr(0.0);
	  JacobianRangeType helpmatr2(0.0); 
	  
	  
	  double intel = volQuad.weight(l);    
	   
	  for(int k = 0;k < gradientNumDofs ;++k)
	    {
	      int VGM_rows  = gradientspc_.mapToGlobal(en,k);
	      int VDM_cols = VGM_rows;
	      
	      // eval tau_k
	      gradbsetEn.eval(k, volQuad, l,tau_);
	      JacobianRangeType tau_mat(0.0);
	      // 	     //convert to Jacobianrangetype
 	      convertutil2d(tau_,tau_mat);

	      for (int j = 0; j < numDofs; ++j) 
		{  
		  bsetEn.eval(j, volQuad, l, psi_[0]);
		 
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
		  
		
		//   std::cout<<"helpermat=\n"<< helpmatr<<"\n";
		  


		  // eval -phi_j*div tau_k
		  // -VGM
		  double VGM =-1.0* gradbsetEn.evaluateGradientSingle(k,en,volQuad,l,helpmatr)*intel;

		  // eval tau_k:grad phi_j
		  // VDM
		  double VDM = bsetEn.evaluateGradientSingle(j,en,volQuad,l,tau_mat)*intel;
		    
		  divMatrix_.add(VDM_rows,VDM_cols,VDM);
		  
		  gradMatrix_.add(VGM_rows,VGM_cols,VGM);

		    ////////////////////////////////////////////////
		    //pressureGradientMatrix und pressureDivMatrix//
		  ////////////////////////////////////////////////
		  
	
		  if(k==0)
		  for(int n = 0; n < pressureNumDofs ; ++n)
		    { //eval q_n
			    pressurebsetEn.eval(n,volQuad,l,qu_);
			    
			    helpmatr2*=0.0;
			    for(int m = 0; m < dimDomain ; ++m)
			      {
				helpmatr2[m][m]=qu_[0];
			      }
			    int PGM_cols  =  pressurespc_.mapToGlobal(en,n);  
			    int PDM_rows  =  PGM_cols;
			    
			    
			   
			   
			
			    
			    //eval (-q_n*div phi_j)
			    double PGM =-1.0*bsetEn.evaluateGradientSingle(j,en, volQuad,l,helpmatr2)*intel; 
			  //   std::cout<<"PGM="<< PGM<<"\n";
			    //eval -(-phi_j*grad q_n)
			    double PDM =1.0*pressurebsetEn.evaluateGradientSingle(n,en,volQuad,l,psi_)*intel;
		// 	    std::cout<<"PDM="<< PDM<<"\n";
			    pressureDivMatrix_.add(PDM_rows,PDM_cols,PDM);
			    pressureGradMatrix_.add(PGM_rows,PGM_cols,PGM);
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
      //   int counter=0;
      // if(0)
      for (IntersectionIterator nit = en.ileafbegin(); nit != endnit; ++nit) 
	  { 
	    int twistSelf = twistUtil_.twistInSelf(nit); 
	    FaceQuadratureType faceQuadInner(nit, faceQuadOrd_, twistSelf, 
					     FaceQuadratureType::INSIDE);
	    if (nit.neighbor()) 
	      { 
		//Access to the neighbor Element
		EntityPointerType neighEp=nit.outside();
		EntityType& nb=  *neighEp;
				
		{ 
		   int twistNeighbor = twistUtil_.twistInNeighbor(nit);
		   FaceQuadratureType faceQuadOuter(nit, faceQuadOrd_, twistNeighbor,
						   FaceQuadratureType::OUTSIDE);
	// 	   FaceQuadratureType faceQuadOuter(nit, faceQuadOrd_, twistSelf, 
// 						   FaceQuadratureType::INSIDE);
		  caller_.setNeighbor(nb);

		  //the Basefunctiosets on the Neighbor element
		  const BaseFunctionSetType& bsetNeigh = spc_.baseFunctionSet(nb);
		  const GradientBaseFunctionSetType& gradbsetNeigh = gradientspc_.baseFunctionSet(nb);
		  const PressureBaseFunctionSetType& pressurebsetNeigh = pressurespc_.baseFunctionSet(nb);
		  
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
			  //	  int VSM_rows=VDM_rows;
			  //			  int PGM_rows=VDM_rows;
			  
			  //calculatin phi_j on the element- and neighbor side of the current intersection
			  bsetEn.eval(j, faceQuadInner, l, phi_); 
			  DomainType x=nit.intersectionGlobal().global(faceQuadInner.localPoint(l));
			  bsetNeigh.eval(j, faceQuadOuter, l, phiNeigh_ );
			  DomainType y=nit.intersectionGlobal().global(faceQuadOuter.localPoint(l));
		
			 //  std::cout<<j<<"\n"; 
// 			  std::cout<<" basefunc evaluation\n";
// 			  //std::cout<<x<<"\n";
// 			  std::cout<<phi_<<"\n";
// 			  //std::cout<<y<<"\n";
// 			  std::cout<<phiNeigh_<<"\n";
// 			 std::cout<<"*************************\n";


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
		
			  u_sigma_flux(normal,int_element,phi_,zerophi,enflux);
			  u_sigma_flux(normal,int_element,zerophi,phiNeigh_,neighflux);
		
			  //\check{\sigma}
			  stab_sigmaflux(normal,hinv,phi_,zerophi,stab_en);
			  stab_sigmaflux(normal,hinv,zerophi,phiNeigh_,stab_nb);
			  //if(j%2==0)
			  // std::cout<<"fluxes\n";
// 			  std::cout<<stab_en<<"\n";
// 			  std::cout<<stab_nb<<"\n";
			   
// 			  int foo;std::cout<<"*********\n";
// 			  std::cin>> foo;
// 			  std::cout<<"*********\n";
// 			  assert(foo==1);
			  /////HACK!!!!!!!!!
			//   if(j>1)
// 			    {//stab_en*=-1.0;
//  			      stab_nb*=-1.0;}
// 			  //\tilde{u}_p

			  uflux(normal,int_element,phi_,zerophi,pressure_en);
			  uflux(normal,int_element,zerophi,phiNeigh_,pressure_neigh);
			  //Loop over Grad Dofs calculation  VGM and VDM 
			  for(int k=0;k < gradientNumDofs;++k)
			    {
			      int VGM_rows = gradientspc_.mapToGlobal(en,k);

			      int VDM_cols_en =VGM_rows;
			      
			      int VDM_cols_nb =gradientspc_.mapToGlobal(nb,k);
			      
			      gradbsetEn.eval(k,faceQuadInner,l,tau_);
			      gradbsetNeigh.eval(k,faceQuadOuter,l,tauNeigh_);
			      
			      RangeType sflux_en;
			      RangeType sflux_nb;

			      sigmaflux(normal,int_element,tau_,zerotau,sflux_en);
			      sigmaflux(normal,int_element,zerotau,tauNeigh_,sflux_nb);
				
			      
			      double VGM_en =  gradbsetEn.evaluateSingle(k,faceQuadInner,l,enflux)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      double VGM_nb = gradbsetEn.evaluateSingle(k,faceQuadInner,l,neighflux)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			     
			      gradMatrix_.add(VGM_rows,VGM_cols_en,VGM_en);
			      gradMatrix_.add(VGM_rows,VGM_cols_nb,VGM_nb);

			      double VDM_en=-1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,sflux_en)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      double VDM_nb=-1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,sflux_nb)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			     				   
			      divMatrix_.add(VDM_rows,VDM_cols_en,VDM_en);
			      divMatrix_.add(VDM_rows,VDM_cols_nb,VDM_nb);
			      
			    }
			 
	 	  for(int n=0;n<pressureNumDofs;++n)
			    {
			      int PGM_cols_en = pressurespc_.mapToGlobal(en,n);
			      int PGM_cols_nb = pressurespc_.mapToGlobal(nb,n);
			      int PDM_rows=PGM_cols_en;
			      int PGM_rows = VDM_rows;
			      
			     
			     
			      
			      pressurebsetEn.eval(n,faceQuadInner,l,qu_);
			      pressurebsetNeigh.eval(n,faceQuadOuter,l,quneigh_);
			      
			      RangeType prflux_en(0.0);
			      RangeType prflux_nb(0.0);
			      PressureRangeType zeroq(0.0);
			 
			      // \hat{p}
			      pressureflux(normal,int_element,qu_,zeroq,prflux_en);
			      pressureflux(normal,int_element,zeroq,quneigh_,prflux_nb);
			      /////////////////////////////////////////////////////////////
			      //PGM
			      double PGM_en=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_en)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      double PGM_nb=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_nb)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      //////////////////////////////////////////////////////////////////

			      pressureGradMatrix_.add(PGM_rows,PGM_cols_en,PGM_en);
			      pressureGradMatrix_.add(PGM_rows,PGM_cols_nb,PGM_nb);
			      
			      ///////////////////////////////////////////////////////////////////
			      //-PDM 
			      double PDM_en =-1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner,l,pressure_en)
				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      double PDM_nb =-1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner,l,pressure_neigh)
				*faceQuadInner.weight(l)* massVolElInv*int_element;

			      ///////////////////////////////////////////////////////////////////////////
			      pressureDivMatrix_.add(PDM_rows,PDM_cols_en,PDM_en);
			      pressureDivMatrix_.add(PDM_rows,PDM_cols_nb,PDM_nb);
			    
			    
			      
			      if(j==0)
				{
				  // int PGM_en = pressurespc_.mapToGlobal(en,j);
				  PressureRangeType pressValen(0.0);
				  PressureRangeType pressValNeigh(0.0);
				  
				  //\check{u}_p
				  stab_uflux(normal,int_element,qu_,zeroqu,pressValen);
				  stab_uflux(normal,int_element,zeroqu,quneigh_,pressValNeigh);
				  
				  for(int m=0;m<pressureNumDofs;++m)
				    { 
				      int PSM_rows = pressurespc_.mapToGlobal(en,m);
				      int PSM_cols_en= PGM_cols_en;
				      int PSM_cols_nb= PGM_cols_nb;
				      ////////////////////////////////////////////////////////////////////////////
				      // -PSM
				      double PSM_en =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner,l,pressValen)
					*faceQuadInner.weight(l)* massVolElInv*int_element;
				      double PSM_nb =1.0*pressurebsetEn.evaluateSingle(m,faceQuadInner,l,pressValNeigh)
					*faceQuadOuter.weight(l)* massVolElInv*int_element;
				      	 
				      pressureStabMatrix_.add(PSM_rows,PSM_cols_en,PSM_en);
				      pressureStabMatrix_.add(PSM_rows,PSM_cols_nb,PSM_nb);
				    }
				}
			    }
			 
			  for(int i=0;i < numDofs;++i)
			    {
			      int VSM_rows=spc_.mapToGlobal(en,i);
			      // int VSM_cols_en=spc_.mapToGlobal(en,i);
			      // int VSM_cols_nb=spc_.mapToGlobal(nb,i);
			      RangeType tau(0.0);
			      bsetEn.eval(i,faceQuadInner,l,tau);
			    //   if(i>1)
// 				tau*=-1.0;
			      
			      double VSM_en =stab_en*tau;
			      VSM_en*=-1.0;
			      VSM_en*=faceQuadInner.weight(l)* massVolElInv*int_element;
			
			      double VSM_nb =stab_nb*tau;
			      VSM_nb*=-1.0;
			      VSM_nb*=faceQuadInner.weight(l)* massVolElInv*int_element;
			  
			      
			     
			// 	{std::cout<<"**Matrixvalues\n";
// 				  std::cout<<VSM_en<<"\n";
// 				  std::cout<<VSM_nb<<"\n";
// 				  int foo;
// 				  std::cout<<"**************\n";
// 				  std::cin>> foo;
// 				  assert(foo==1);}
			      // 	  if(j>1)
			     //  if(i>1)
// 				VSM_nb*=-1.0;
			     //  double VSM_en =-1.0*bsetEn.evaluateSingle(i,faceQuadInner,l,stab_en)
// 				*faceQuadInner.weight(l)* massVolElInv*int_element;
// 			      double VSM_nb =-1.0*bsetEn.evaluateSingle(i,faceQuadInner,l,stab_nb)
// 				*faceQuadInner.weight(l)* massVolElInv*int_element;
			      stabMatrix_.add(VSM_rows,VSM_cols_en,VSM_en);
// 			      std::cout<<"entity row="<<VSM_rows<<" cols="<<VSM_cols_en<<" value="<<VSM_en<<"\n";
			      stabMatrix_.add(VSM_rows,VSM_cols_nb,VSM_nb);
			     //  if(i%2==0 && j%2==0)
// 				{ std::cout<<"entity row="<<VSM_rows<<" cols="<<VSM_cols_en<<" value="<<VSM_en<<"\n";
// 				  std::cout<<"neigh row="<<VSM_rows<<" cols= "<<VSM_cols_nb<<" value="<<stab_nb<<"\n";}
			    }		


		    }




		    }
		
	      }
  }
	    // if(0)	 //fuer stabmatrix sigmaflux(0,0,phi,0) auswerten
	      if (nit.boundary()) 
		{ RangeType zerophi(0.0); 
	// 	  LocalVeloType localVelo=VeloRhs_.localFunction(en);
// 		  LocalStressType localStress =StressRhs_.localFunction(en);
// 		  LocalPressureType localPressure=PressureRhs_.localFunction(en);
		  double* Veloptr =VeloRhs_.leakPointer();
		  double* Stressptr =StressRhs_.leakPointer();
		  double* Pressptr =PressureRhs_.leakPointer();
		  const int quadNop = faceQuadInner.nop();
		  for (int l = 0; l < quadNop ; ++l) 
		  {
		    DomainType normal=nit.unitOuterNormal(faceQuadInner.localPoint(l));
		    // double  int_element=nit.integrationOuterNormal(faceQuadInner.localPoint(l)).two_norm();
		    double int_element=nit.intersectionGlobal().integrationElement(faceQuadInner.localPoint(l));
		    double hinv=1;
		    hinv/=int_element;
		    RangeType DirichletValue;


		    
		    boundaryvalue(nit,faceQuadInner,l,DirichletValue);
		   //  std::cout<<"dirichletval \n";
// 		    std::cout<<DirichletValue[0]<<"\n";
// 		    std::cout<<DirichletValue[1]<<"\n";
		  
		    
		    for(int j=0; j<numDofs; ++j)
		      { 
			int VDM_rows=spc_.mapToGlobal(en,j);
			int VSM_rows=VDM_rows;
		
			//int VDM_cols_bnd = StressRhs_rows;
			    
			bsetEn.eval(j,faceQuadInner,l,phi_);
		
			//Boundary contribution for velo-RHS
		// 	localVelo[j]+=bsetEn.evaluateSingle(j,faceQuadInner,l,DirichletValue)
// 			  *faceQuadInner.weight(l)*massVolElInv*int_element*int_element;

			Veloptr[VDM_rows]-=c11_*bsetEn.evaluateSingle(j,faceQuadInner,l,DirichletValue)
   			  *faceQuadInner.weight(l)*massVolElInv*int_element*hinv;
			
			for(int k=0;k<gradientNumDofs;++k)
			  {
			    int VDM_cols_bnd=gradientspc_.mapToGlobal(en,k);
			    int Rhsrow=VDM_cols_bnd;
			    gradbsetEn.eval(k,faceQuadInner,l,tau_);
			    
			    RangeType sflux_bnd(0.0);
			    
			    multi2d(tau_,normal,sflux_bnd);
			    sigmaflux(normal,int_element,tau_,tau_,sflux_bnd);
			    //VDM
			    double VDM_bnd=-1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,sflux_bnd)
			      *faceQuadInner.weight(l)*massVolElInv*int_element;
			    
			    divMatrix_.add(VDM_rows,VDM_cols_bnd,VDM_bnd);
			    GradientRangeType bnd_stress(0.0);
			    tensorProduct2d(DirichletValue,normal,bnd_stress);
			   
			    //Boundary contribution for stress-RHS
			    if(j==0)
			      // localStress[k]+=gradbsetEn.evaluateSingle(k,faceQuadInner,l,bnd_stress)
// 				*faceQuadInner.weight(l)*massVolElInv*int_element;
			      Stressptr[Rhsrow]+=gradbsetEn.evaluateSingle(k,faceQuadInner,l,bnd_stress)
				*faceQuadInner.weight(l)*massVolElInv*int_element;
			  }


    
			for(int i=0;i<numDofs;++i)
			    {
			      int VSM_cols_bnd=spc_.mapToGlobal(en,i);
			      
			      bsetEn.eval(j,faceQuadInner,l,phi_); 	  
			      
			      RangeType stab_bnd(0.0);
			      stab_sigmaflux(normal,hinv,phi_,zerophi,stab_bnd);
			        
			      double VSM_bnd=c11_*bsetEn.evaluateSingle(i,faceQuadInner,l,stab_bnd)
				*faceQuadInner.weight(l)*massVolElInv*int_element;
				  
			      
			      stabMatrix_.add(VSM_rows,VSM_cols_bnd,VSM_bnd);
			    }
			for(int n=0;n<pressureNumDofs;++n)
			  
			  { int PGM_rows=spc_.mapToGlobal(en,j);
			    int PGM_cols_bnd=pressurespc_.mapToGlobal(en,n);
			    int PSM_row= PGM_cols_bnd;
			    pressurebsetEn.eval(n,faceQuadInner,l,qu_);
				  
			    RangeType prflux_bnd(0.0);
			    
			   //  prflux_bnd=normal;
// 			    prflux_bnd*=qu_;
			    
			   pressureflux(normal,int_element,qu_,qu_,prflux_bnd);

			   double PGM_bnd=1.0*bsetEn.evaluateSingle(j,faceQuadInner,l,prflux_bnd)
			     *faceQuadInner.weight(l)*massVolElInv*int_element;
			    
			   pressureGradMatrix_.add(PGM_rows,PGM_cols_bnd,PGM_bnd);
			    
		 	    PressureRangeType bnd_pressure=DirichletValue*normal;
			    //Boundary contribution for pressure-RHS
			    //-H implementieren siehe Uzawa-Algorithmus 
			   //  if(j==0)
		  	    //   localPressure[n]+=1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner,l,bnd_pressure)
			      // 				*faceQuadInner.weight(l)*massVolElInv*int_element;	
			    if(j==0)
			     { Pressptr[PSM_row]+=1.0*pressurebsetEn.evaluateSingle(n,faceQuadInner,l,bnd_pressure)
				 *faceQuadInner.weight(l)*massVolElInv*int_element;}
			    
				  
			  }
		      
		  }
			
		
		    
		      
		  
		  
	      
	      
		}	
		  //std::cout<<"divmatrix=\n";
		  //divMatrix_.print(std::cout);
	  } // end if boundary
	      
    /*std::cout<<"gradmatrix=\n";
      gradMatrix_.print(std::cout);*/
// 	std::cout<<"helpermatrix=\n";
// 	std::cout<<helpmatr<<"\n";
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
      //   massVolInv = 1.0;
      return result;    }
    
    
#if GLATT
    void boundaryvalue(const IntersectionIterator& nit, 
		       const FaceQuadratureType& FaceQuad,
		       const int l,
		       RangeType& ret) const
    {
      DomainType x=nit.intersectionGlobal().global(FaceQuad.localPoint(l));
      //dirichletValue=problem_.dirichlet_velo(point_l);
    //   std::cout<<"x1: "<<x[0]<<"\n";
//       std::cout<<"x2: "<<x[1]<<"\n";
//paper 
    //u1
      ret[0]=cos(x[1]);
      ret[0]*=x[1];
      ret[0]+=sin(x[1]);
      ret[0]*=exp(x[0]);
      ret[0]*=-1.0;

      //u2
      ret[1]=sin(x[1]);
      ret[1]*=x[1];
      ret[1]*=exp(x[0]);
    }
    
#endif

#if SNG
    void boundaryvalue(const IntersectionIterator& nit, 
		    const FaceQuadratureType& FaceQuad,
		       const int l,
		       RangeType& ret) const
    {
      DomainType x=nit.intersectionGlobal().global(FaceQuad.localPoint(l));
      //dirichletValue=problem_.dirichlet_velo(point_l);
    //   std::cout<<"x1: "<<x[0]<<"\n";
//       std::cout<<"x2: "<<x[1]<<"\n";


     
      

      DomainType tmp=x;
      tmp[0]-=1.1;
      tmp[1]-=1.1;
      double abs=tmp*tmp  ;
      
      ret[0]=tmp[0];
      ret[0]/=abs;
      ret[1]=tmp[1];
      ret[1]/=abs;

    }

#endif


    // \hat{u}_{\sigma}
    void u_sigma_flux(const DomainType& normal,
		      const double & scaling,
		      const RangeType& phileft,
		      const RangeType& phiright, 
		      GradientRangeType & ret)const
    {//DomainType C12(0.0);
       DomainType C12=normal;
       //double factor=normal*C12;
        C12*=scaling;
       C12*=c12_;
      
      RangeType tmp,tmp2;
      GradientRangeType tmp3;
      tmp2=phileft;
      tmp2-=phiright;
      // tmp2*=factor;
      tensorProduct2d(tmp2,normal,tmp3);
      multi2d(tmp3,C12,tmp2);                                                                                                                                                                                                                                                                                                                                                                                                                

      tmp=phileft+phiright;
      tmp *= 0.5;
      tmp+=tmp2;
      
      
      
      tensorProduct2d(tmp,normal,ret);
 
   }

    
    //\tilde{sigma}
    void sigmaflux(const DomainType& normal,
		   const double& scaling,
		   const GradientRangeType& tauleft,
		   const GradientRangeType& tauright,
		   RangeType& ret)const
    {
      //  DomainType C12(0.0); 
      DomainType C12=normal;
            C12*=scaling; 
      C12*=c12_; 
      GradientRangeType tmp1,tmp2;
      RangeType tmp3;
      tmp1=tauleft;
      tmp1+=tauright;
      tmp1*=0.5;

      tmp2=tauleft;
      tmp2-=tauright;
      
      multi2d(tmp2,normal,tmp3);

      tensorProduct2d(tmp3,C12,tmp2);
      //   tmp2*=scaling;
      tmp1-=tmp2;

      multi2d(tmp1,normal,ret);
      
     }
    
    //\check{\sigma}
    void stab_sigmaflux(const DomainType& normal,
			const double& scaling,
			const RangeType& phileft,
			const RangeType& phiright,
			RangeType & ret)const
      {
	double C11=c11_;
	C11*=scaling;
	// std::cout<<"scaling"<<c11 <<"\n";
	RangeType tmp(0.0);
	GradientRangeType tmp2;
	
	tmp=phileft;
	tmp-=phiright;
	tmp*=C11;
	
// 	std::cout<< tmp<<"\n";

	ret=tmp;
// 	tensorProduct2d(tmp,normal,tmp2);

//  	multi2d(tmp2,normal,ret);
// 	std::cout<< ret<<"\n";

// 	int foo;
// 	  std::cin>>foo;
// 	assert(foo==1);
      }
    
    //\tilde{u}_p
    void uflux(const DomainType& normal,
	       const double& scaling,
	       const RangeType& phileft,
	       const RangeType& phiright,
	       PressureRangeType & ret) const
    { //RangeType D12(0.0);
      RangeType D12=normal;
       D12*=scaling;
       D12*=d12_;
      RangeType tmp1,tmp2;
      
      tmp1=phileft;
      tmp1+=phiright;
      tmp1*=0.5;

      tmp2=phileft;
      tmp2-=phiright;

      double factor=tmp2*normal;
      tmp2=D12;
      tmp2*=factor;
      tmp2+=tmp1;
      ret=tmp2*normal;
    }
    
    
    
    
    
    //\hat{p}
    void pressureflux(const DomainType& normal,	
		      const double& scaling,
		      const PressureRangeType& phileft,
		      const  PressureRangeType& phiright, 
		      RangeType & ret)const
    {
      PressureRangeType tmp;
      PressureRangeType tmp2;
      RangeType tmp3;
      //  RangeType D12(0.0);
       RangeType D12=normal;
        D12*=scaling;
      //
        D12*=d12_;;
      
      
      tmp=phileft;
      tmp+=phiright;
      tmp *= 0.5;
      
      tmp2=phileft;
      tmp2-=phiright;
      
      tmp3=normal;
      
      tmp3*=tmp2;
      tmp2=D12*tmp3;
     
      
      tmp-=tmp2;
      
      ret=normal;
      ret*=tmp;
     
     
    }


  //\check{u}_p
    // multiplication with normal twice yields 1!!!
    void stab_uflux(const DomainType& normal,	
		    const double& scaling,
		    const PressureRangeType& qleft,
		    const  PressureRangeType& qright, 
		    PressureRangeType & ret)const
    {    
      double D11=d11_;
      D11*=scaling;
      ret=qleft-qright;

      ret*=D11;
    }

    void multi2d(const GradientRangeType& tau,const DomainType& n,RangeType& ret) const
    {
      ret[0]=tau[0]*n[0]+tau[1]*n[1];
      ret[1]=tau[2]*n[0]+tau[3]*n[1];
    }
	
    void multi3d(const GradientRangeType& tau,const DomainType& n,RangeType& ret)
    {
      ret[0]=tau[0]*n[0]+tau[1]*n[1]+tau[2]*n[2];
      ret[1]=tau[3]*n[0]+tau[4]*n[1]+tau[5]*n[2];
      ret[3]=tau[6]*n[0]+tau[7]*n[1]+tau[8]*n[2];
    }
    //converts a GradientRangeType to a JaobianRangeType:FielVector<n*n>->FieldMatrix<n,n>
    void convertutil2d(const GradientRangeType& tau,JacobianRangeType& mat)const
    {
      mat[0][0]=tau[0];
      mat[0][1]=tau[1]; 
      mat[1][0]=tau[2];
      mat[1][1]=tau[3];
    }
    
    //converts a GradientRangeType to a JaobianRangeType:FielVector<n*n>->FieldMatrix<n,n>
    void convertutil3d(const GradientRangeType& tau,JacobianRangeType& mat)const
    {
      mat[0][0]=tau[0];
      mat[0][1]=tau[1]; 
      mat[0][2]=tau[2];
      mat[1][0]=tau[3];
      mat[1][1]=tau[4];
      mat[1][2]=tau[5]; 
      mat[2][0]=tau[6];
      mat[2][1]=tau[7];
      mat[2][2]=tau[8];
     }
    void tensorProduct2d(const RangeType& v,const DomainType& n, GradientRangeType& ret) const
    {
      ret[0]= v[0]*n[0];
      ret[1]= v[0]*n[1];
      ret[2]= v[1]*n[0];
      ret[3]= v[1]*n[1];
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
    mutable DiscreteModelCallerType caller_;
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
    TimeProvider* time_;
    FieldVector<int, 0> diffVar_;
  

    //neue Variablen
    // mutable GridPartType& gridPrt_; 
  

    TwistUtility<GridType> twistUtil_;

    int volumeQuadOrd_,faceQuadOrd_;

    int maxNumberUnknowns_;
 
    
    mutable MatrixType massMatrix_ ,massMatrixInv_;
    mutable MatrixType gradMatrix_,divMatrix_,stabMatrix_;
    mutable MatrixType pressureGradMatrix_,pressureDivMatrix_,pressureStabMatrix_;
    mutable double c11_;
    mutable double c12_;
    mutable double d11_;
    mutable double d12_;
    mutable bool matrixAssembled_;

//     mutable BOPType BOP_;
//     mutable BTOPType BTOP_;
//     mutable COPType COP_;

  };
  
} // end namespace Dune

#endif
