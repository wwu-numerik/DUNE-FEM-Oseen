#ifndef __LDGEX_FLUX_CC__
#define __LDGEX_FLUX_CC__



template<class ModelImp>
class StokesFlux
{
public:
  typedef ModelImp ModelType;
  typedef typename ModelType::Traits::DomainType DomainType;
  typedef typename ModelType::Traits ModelTraits;
  enum{dimRange =  ModelTraits::dimRange};
    
  typedef typename ModelType::RangeType RangeType;
  typedef typename ModelType::GradRangeType GradientRangeType;
  typedef typename ModelType::PressureRangeType PressureRangeType;
  




public:
  StokesFlux(const ModelType& mod,
	     const double c11,
	     const double c12,
	     const double d11, 
	     const double d12):
    mod_(mod),
    c11_(c11),
    c12_(c12),
    d11_(d11),
    d12_(d12)
  {
    
     upw_[0]=M_PI;
     upw_[1]=exp(1.0);
  }

// \hat{u}_{\sigma}
  void u_sigma_flux(const DomainType& normal,
		    const double & scaling,
		    const RangeType& phileft,
		    const RangeType& phiright, 
		    GradientRangeType & left,
		    GradientRangeType & right)const
{
  double swtch=0.5;
  
  swtch=c12_;
  
  swtch*=scaling;
  if(normal*upw_<0)
    swtch*=-1.0;
  
  RangeType tmp,tmp2(0.0);
  //GradientRangeType tmp3;
      
  ///////left value ///////
    
    tmp=phileft;
    tmp *= 0.5;
    
    tmp2=phileft;
    tmp2*=swtch;
    // tmp+=tmp2;
    
    tensorProduct2d(tmp,normal,left); 
   
    ///////right value //////
      tmp=phiright;
      tmp *= 0.5;
      tmp2=phiright;
      tmp2*=-1.0;
      tmp2*=swtch;
    //   tmp+=tmp2;
      
      tensorProduct2d(tmp,normal,right);


     
 //  DomainType C12=normal;
//   if((normal[0]+normal[1])<0)
//     C12*=-1.0;
      
      
//   C12*=scaling;
//   C12*=c12_;
     
    
//   RangeType tmp,tmp2(0.0);
//   GradientRangeType tmp3;
//   /////////////////leftvalue/////////////////    
//     tmp=phileft;
//     tmp *= 0.5;
//     tmp2=phileft;
//     tensorProduct2d(tmp2,normal,tmp3);
//     tmp2=0.0;
//     tmp3.umv(C12,tmp2);         
//     tmp+=tmp2;
      
//     tensorProduct2d(tmp,normal,left); 
//     ////////////////rightvalue//////////////////  
//     tmp=phiright;
//     tmp *= 0.5;
//     tmp2=phiright;
//     tmp2*=-1.0;
//     tensorProduct2d(tmp2,normal,tmp3);
//     tmp2=0.0;
//     tmp3.umv(C12,tmp2);         
//     tmp+=tmp2;
        
//     tensorProduct2d(tmp,normal,right);
 
} 







//\tilde{sigma}
void sigmaflux(const DomainType& normal,
	       const double& scaling,
	       const GradientRangeType& tauleft,
	       const GradientRangeType& tauright,
	       RangeType& left,
	       RangeType& right)const
{
      double swtch=0.5;
     
      swtch=c12_;
      swtch*=scaling;
      if(normal*upw_<0)
	swtch*=-1.0;
      GradientRangeType tmp1(0.0),tmp2; 
      
      /////////////////leftvalue/////////////////

	tmp1=tauleft;
	tmp2=tmp1;
	tmp1*=0.5;
	tmp2*=swtch;
	tmp1-=tmp2;
      
      
	tmp1.umv(normal,left);
      
	//////////////rightvalue///////////////////////////  
	tmp1=0.0;

	tmp1=tauright;
	tmp2=tmp1;	
	tmp1*=0.5; 
	tmp2*=-1.0;
	tmp2*=swtch;
// 	tmp1-=tmp2;
	tmp1.umv(normal,right);

}



//\tilde{u}_p
void uflux(const DomainType& normal,
	   const double& scaling,
	   const RangeType& phileft,
	   const RangeType& phiright,
	   PressureRangeType & left,
	   PressureRangeType & right ) const
{ //RangeType D12(0.0);
  RangeType D12=normal;
     
  D12[0]=1.0;
  D12[1]=1.0;
  //   if((normal[0]+normal[1])<0)
  // 	d12_*=-1.0;
      
  D12*=scaling;
  D12*=d12_;
  RangeType tmp1,tmp2;
  ////////leftvalue///////////////
  tmp1=phileft;
  //  tmp1+=phiright;
  tmp1*=0.5;

  tmp2=phileft;
  // tmp2-=phiright;

  double factor=tmp2*normal;
      
  tmp2=D12;
  tmp2*=factor;
  tmp2+=tmp1;
  left=tmp2*normal;
    
  //////////////rightvalue//////////////////
     
  tmp1=phiright;
      
  tmp1*=0.5;

  tmp2=phiright;
  tmp2*=-1.0;
      
  factor=tmp2*normal;
      
  tmp2=D12;
      
  tmp2*=factor;
  tmp2+=tmp1;
  right=tmp2*normal;
} 


//\hat{p}
void pressureflux(const DomainType& normal,	
		  const double& scaling,
		  const PressureRangeType& phileft,
		  const  PressureRangeType& phiright, 
		  RangeType & left,
		  RangeType & right )const
{
  PressureRangeType tmp;
  PressureRangeType tmp2;
  RangeType tmp3;
  //  RangeType D12(0.0);
  RangeType D12=normal;
  D12[0]=1.0;
  D12[1]=1.0;
  //if((normal[0]+normal[1])<0)
  //	d12_*=-1.0;
  D12*=scaling;
      
  D12*=d12_;
      
  //////////////////left value///////////////
           
  tmp=phileft;
  tmp *= 0.5;
      
  tmp2=phileft;
            
  tmp3=normal;
      
  tmp3*=tmp2;
   
  tmp2=D12*tmp3;
     
      
  // tmp-=tmp2;
      
  left=normal;
  left*=tmp;
  ///////////////////////right value/////////////////////
	
	
    tmp=phiright;
    tmp *= 0.5;
      
    tmp2=phiright;
    tmp2*=-1.0;
      
    tmp3=normal;
      
    tmp3*=tmp2;
      
    tmp2=D12*tmp3;
     
      
//     tmp-=tmp2;
      
    right=normal;
    right*=tmp; 
}

 //\check{\sigma}
    void stab_sigmaflux(const DomainType& normal,
			const double& scaling,
			const RangeType& phileft,
			const RangeType& phiright,
			RangeType & gLeft,
			RangeType & gRight)const
      {
	double C11=c11_;
	C11*=scaling;
	C11*=-1.0;
	//	RangeType tmp(0.0);
	//GradientRangeType tmp2;
	
	gLeft =phileft;
	gLeft*=C11;

	gRight=phiright;
	gRight*=-1.0;
	gRight*=C11;

       }



//\check{u}_p
// multiplication with normal twice yields 1!!!
void stab_uflux(const DomainType& normal,	
		const double& scaling,
		const PressureRangeType& qleft,
		const  PressureRangeType& qright, 
		PressureRangeType & gLeft,
		PressureRangeType & gRight)const
    {    
      double D11=d11_;
      D11*=scaling;
      gLeft=qleft;
      gLeft*=D11;

      gRight=qright;
      gRight*=-1.0;
      gRight*=D11;
    }

inline  void tensorProduct2d(const RangeType& v,const DomainType& n, GradientRangeType& ret) const
    {
      ret(0,0)= v[0]*n[0];
      ret(0,1)= v[0]*n[1];
      ret(1,0)= v[1]*n[0];
      ret(1,1)= v[1]*n[1];


//  ret(0,0)= v[0]*n[0];
//  ret(0,1)= v[1]*n[0];
//  ret(1,0)= v[0]*n[1];
//  ret(1,1)= v[1]*n[1];


    }

private:
const ModelType& mod_;
DomainType upw_;
const double c11_;
const double c12_;
const double d11_;
const double d12_;
};

#endif
