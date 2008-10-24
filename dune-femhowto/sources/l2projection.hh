#include <dune/fem/operator/projection/l2projection.hh>
//***********************************************************************
class ProjectionExample {
 //! define the function space our unkown belongs to 
  typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
 public:
  //! define the type of discrete function we are using 
  typedef DFAdapt < DiscreteFunctionSpaceType > DiscreteFunctionType;
    //! perform l2-projection
  template <class Function>
  void solve (const Function f, DiscreteFunctionType& solution ) {
    L2Projection<double,double,Function,DiscreteFunctionType> l2proj;
    l2proj(f, solution);
  }
};
typedef ProjectionExample ExampleType;

