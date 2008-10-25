#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

// FR passes 
// #include <dune/fem/pass/dgpass.hh>
// #include <dune/fem/pass/discretemodel.hh>
// #include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/dgspace.hh>

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/istl/bvector.hh>
//#include <dune/fem/function/blockvectorfunction.hh>

#include <dune/fem/quadrature/cachequad.hh>
// #include <dune/grid/common/referenceelements.cc>
#include <dune/fem/quadrature/elementquadrature.hh>

#include <dune/istl/bvector.hh>
#include <dune/fem/function/blockvectorfunction.hh>

//*************************************************************
namespace LDGExample {  

  using namespace Dune;

  // MethodOrderTraits
  template <class Model,int polOrd,int dimRange>
  class PassTraits {
  public:
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };
    
    typedef LeafGridPart<GridType> GridPartType;
    //typedef HierarchicGridPart<GridType> GridPartType;
    

    
    /* Velocity Space */
    typedef FunctionSpace<double, double, dimDomain, dimRange> FunctionSpaceType; //Velocity Space
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    /* Gradient Space */
    typedef MatrixFunctionSpace<DomainFieldType,RangeFieldType,dimDomain,dimRange,dimRange> GradientSpaceType;
    /*Pressure Space */
    typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,1> PressureSpaceType;
    
    
        
    
    // typedef LegendreDiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    /* The Discrete Function Space for the Velocity*/
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<DiscreteFunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    /* The Discrete Function Space for the Gradient*/
    typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd-1,CachingStorage> DiscreteGradientSpaceType;   
  
    //typedef LegendreDiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, polOrd> DiscreteGradientSpaceType;
    /* The Discrete Function Space for the prerssure*/
    typedef DiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd-1,CachingStorage> DiscretePressureSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<PressureSpaceType, GridPartType, polOrd> DiscretePressureSpaceType;


    //  typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
    //typedef typename DiscreteFunctionType::LeakPointerType VeloLeakPointerType;
    typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
   
  //typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
    typedef double* VeloLeakPointerType; 
  };

  //**********************************************************
  //
  //  StokesTraits 
  //
  //**********************************************************

  template <class Model,class NumFlux,int polOrd>
  class StokesDiscreteModel;

  template <class Model,class NumFlux,int polOrd>
  struct StokesTraits
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimDomain = ModelTraits::dimDomain };
    enum { dimRange = ModelTraits::dimRange };

    typedef PassTraits<Model,polOrd,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename Traits::GradientSpaceType GradientSpaceType;
    typedef typename GradientSpaceType::RangeType GradRangeType;

    typedef typename Traits::PressureSpaceType PressureSpaceType;
    typedef typename PressureSpaceType::RangeType PressureRangeType;
    
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef DiscreteFunctionType DestinationType;
    
    typedef StokesDiscreteModel<Model,NumFlux,polOrd> DiscreteModelType;
  };

  template <class Model,class NumFlux,int polOrd>
  class StokesDiscreteModel
  { 
  public:
    enum { polynomialOrder = polOrd };

    typedef StokesTraits<Model,NumFlux,polOrd> Traits;
    //typedef Dune::Selector<0,1> SelectorType;
 //    typedef Dune::Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GradRangeType GradRangeType;
  
    typedef typename Traits::PressureRangeType PressureRangeType;
    
    typedef typename GridType::Traits::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename EntityType::EntityPointer EntityPointerType;
 
  public:
    StokesDiscreteModel(const Model& mod,const NumFlux& numf,const double nu) :
      model_(mod),
      numflux_(numf),
      nu_(nu)
    {}

    // const Model & data () const { return model_; }

    bool hasSource() const { return true; }
    bool hasFlux() const { return true; }


    template<class QuadratureType>
    inline void dirichletValues(const IntersectionIterator & it, 
			   const QuadratureType& quad,
			   int p,
			   RangeType &ret)
      {
      	model_.dirichletValues(it,quad,p,ret);
      }

    



    /* The numerical Flux for the variable u for the first equation */
    template <class QuadratureType>
    void numericalFlux_u_sigma(IntersectionIterator& it,
			       const QuadratureType quad,
			       int p,
			       const RangeType& uLeft,
			       const RangeType& uRight,
			       GradRangeType& gLeft,
			       GradRangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.u_sigma_flux(normal,
			    h,
			    uLeft,
			    uRight,
			    gLeft,
			    gRight);
      
    }
  


    /* The numerical Flux for the variable sigma */
    template <class QuadratureType>
    void numericalFlux_sigma (IntersectionIterator& it,
			      const QuadratureType quad,
			      int p,
			      //const ArgumentTuple& uLeft,
			      //const ArgumentTuple& uRight,
			      const GradRangeType& uLeft,
			      const GradRangeType& uRight,
			      RangeType& gLeft,
			      RangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.sigmaflux(normal,
			 h,
			 uLeft,
			 uRight,
			 gLeft,
			 gRight);

      
    }

    /* The numerical Flux for the variable u for the divergence constraint */
    template <class QuadratureType>
    void numericalFlux_u_p(IntersectionIterator& it,
			   const QuadratureType quad,
			   int p,
			   const RangeType& uLeft,
			   const RangeType& uRight,
			   PressureRangeType& gLeft,
			   PressureRangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.uflux(normal,
			    h,
			    uLeft,
			    uRight,
			    gLeft,
			    gRight);
      
    }
  
    /* The numerical Flux for the variable p */
    template <class QuadratureType>
    void numericalFlux_p (IntersectionIterator& it,
			      const QuadratureType quad,
			      int p,
			      const PressureRangeType& uLeft,
			      const PressureRangeType& uRight,
			      RangeType& gLeft,
			      RangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.pressureflux(normal,
			    h,
			    uLeft,
			    uRight,
			    gLeft,
			    gRight);

      
    }


    /* Stablization of u */
    template <class QuadratureType>
    void numericalFlux_u_stab (IntersectionIterator& it,
			     const QuadratureType quad,
			     int p,
			     const RangeType& uLeft,
			     const RangeType& uRight,
			     RangeType& gLeft,
			     RangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      //double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      double hinv=1.0;				
      hinv/=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.stab_sigmaflux(normal,
			hinv,
			uLeft,
			uRight,
 			gLeft,
			gRight);
      
    }
    


    /* Stablization of p */
    template <class QuadratureType>
    void numericalFlux_p_stab (IntersectionIterator& it,
			       const QuadratureType quad,
			       int p,
			       const PressureRangeType& uLeft,
			       const PressureRangeType& uRight,
			       PressureRangeType& gLeft,
			       PressureRangeType& gRight) const 
    {
      const DomainType normal = it.unitOuterNormal(quad.localPoint(p));
      //double h=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      double hinv=1.0;				
      hinv/=it.intersectionGlobal().integrationElement(quad.localPoint(p));
      numflux_.stab_uflux(normal,
			hinv,
			uLeft,
			uRight,
 			gLeft,
			gRight);
      
    }


    template<class QuadratureType>
    void mass(const EntityType & en, 
	      const double time ,
	      const QuadratureType& quad,
	      int p,
	      const RangeType& phi,
	      RangeType &ret ) const
    {

      
      model_.mass(en,time,quad,p,phi,ret);
      return;
    }


    template<class QuadratureType>
    void convectiveFlux(const EntityType & en, 
		       const double time ,
		       const QuadratureType& quad,
		       int p,
		       const RangeType& phi,
		       JacobianRangeType &ret)
    {
      assert(false);
      model_.convectiveFlux(en,time,quad,p,phi,ret);
    }
    
    double get_nu(){return nu_;}
    
    template<class QuadratureType>
      inline void rightHandSide(const EntityType& en,
			 const QuadratureType quad,
			 int p,
			 RangeType& ret)const
      {
	//  ret=0.0;
	//  DomainType x;
	//       x=en.geometry().global(quad.localPoint(p));
      
	model_.source(en,quad,p,ret);
	//     std::cout<<ret<<std::endl;
      }



   
  private:
    const Model& model_;
    const NumFlux& numflux_;
    const double nu_;
  };

}
#endif
