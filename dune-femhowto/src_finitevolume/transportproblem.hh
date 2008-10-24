#ifndef TRANSPORTPROBLEM_HH 
#define TRANSPORTPROBLEM_HH 

//===============================================================
// We want to solve the PDE
//
//   dc/dt + div( uc ) = 0   in Omega
//                   c = b   on Gamma_in
//                   c = c0  for t<0
//
// We define functions for the parameters c0, u and b.
//===============================================================
template <class GridType>
struct ProblemImp 
{

  // **** type of analytical function space, here f: R^dim --> R
  typedef Dune::FunctionSpace<double, double, GridType::dimension, 1>  FunctionSpaceType;
  
  // the exact solution 
  template<int dimworld, class ct>
  static inline double exact(const Dune::FieldVector<ct,dimworld>& x, double t)
  {
   
    if(x[0]<=0.25+t)
      return 2.0;
    else
      return 1.0;
   

    /******************************
    if (x.infinity_norm()<t+0.25)
      return 1.0;
    else
      return 0.0;
    ***************************/
    

    /****************************************************
    Dune::FieldVector<ct,dimworld> velocity;
    velocity(t,x,velocity);
    velocity*=t;
    double twoNorm = (x-velocity).two_norm();
    bool outside = true;
    for (unsigned int i = 0;i<x.dim();++i)
      if (x[i]<velocity[i]*t) outside=false;
    

    if ( (x.infinity_norm()<t+0.25) && ( (!outside   ) ||  ((outside) && (twoNorm<t+0.25)    ) )  )
      return 1;
    else
      return 0;
    ************************************************/

  }
  
  // the initial condition c0
  template<int dimworld, class ct>
  static inline double initial(const Dune::FieldVector<ct,dimworld>& x)
  {
    return exact(x,0.0);
  }
  
  // the boundary condition b on inflow boundary
  template<int dimworld, class ct>
  static inline double boundary (const Dune::FieldVector<ct,dimworld>& x, double t)
  {
    return exact(x,t);
  }
  
  // the vector field u set to velocity 
  template<int dimworld, class ct>
  static inline void velocity(const double time,
			      const Dune::FieldVector<ct,dimworld>& x, 
			      Dune::FieldVector<ct,dimworld>& velocity)
  {
    velocity = 0.0;
    velocity[0]=1.0;
  }
  
}; // end class Problem 



//===============================================================
//
// The numerical Flux
//
//===============================================================
template <class DomainType, class ctype> 
struct LinearUpwindFlux
{
  //problem type
  typedef ProblemImp<GridType> Problem;
  
  static inline double numericalFlux(const DomainType& normal,
                                     const DomainType& velocity,
                                     const ctype& uLeft,
                                     const ctype& uRight,
                                     ctype& gLeft,
                                     ctype& gRight)
  {
    const ctype upwind = normal * velocity;
    if(upwind > 0)
    {
      gLeft = uLeft * upwind;
    }
    else 
    {
      gLeft = uRight * upwind;
    }
    gRight = gLeft;

    return std::abs(upwind);
  }


  static inline double boundaryFlux(const double time,
                                    const DomainType& xGlobal,
                                    const DomainType& normal,
                                    const DomainType& velocity,
                                    const ctype& uLeft,
                                    ctype& gLeft)
  {
    const ctype upwind = normal * velocity;
    if(upwind > 0)
    {
      gLeft = uLeft * upwind;
    }
    else 
    {
      gLeft = Problem::boundary(xGlobal,time) * upwind;
    }
    return std::abs(upwind);
  }

}; // end class LinearUpwindFlux 


#endif
