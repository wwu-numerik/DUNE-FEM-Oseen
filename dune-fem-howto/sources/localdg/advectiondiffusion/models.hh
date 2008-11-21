#ifndef ADVDIFF_MODELS_HH
#define ADVDIFF_MODELS_HH

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/solver/odesolver.hh>

#include "advectdiff.hh"
// determine polynomial order 
enum { 
  order=POLORDER, // order of function spaces 
  rksteps=POLORDER+1 // order of ode solver 
}; 

using namespace Dune;
using namespace std;

// Select grid part 
//typedef LeafGridPart<GridType> GridPartType;
//typedef HierarchicGridPart<GridType> GridPartType;
typedef DGAdaptiveLeafGridPart<GridType> GridPartType;

//////////////////////////////////
// scalar model and flux choice 
////////////////////////////////// 
#if PROBLEM == 0 // Rotating Pulse 
   #include "scalarmodels.hh"
   #include "initpulse.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   //typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   //typedef DGLimitedAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   //typedef DuneODE::ExplRungeKutta<DgType> ODEType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
   //typedef DuneODE::ImplTimeStepper<DgType> ODEType;
   
#elif PROBLEM == 1 // Advection Diffusion 
   #include "scalarmodels.hh"
   #include "initadvectdiff.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   //typedef UpwindFlux<ModelType> FluxType;

   typedef LLFFlux<ModelType> FluxType;
   //typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgType;

   // choose ode solver 
   typedef DuneODE::ExplRungeKutta<DgType> ODEType;
   //typedef DuneODE::ExplTimeStepper<DgType> ODEType;
   //typedef DuneODE::ImplTimeStepper<DgType> ODEType;

#elif PROBLEM == 2 // Burgers equation 
   #include "scalarmodels.hh"
   #include "initburgers.cc"
   typedef U0<GridType> InitialDataType;
   typedef BurgersModel<GridPartType,InitialDataType > ModelType;
   typedef LLFFlux<ModelType> FluxType;
   typedef DGLimitedAdvectionOperator<ModelType,LLFFlux,order> DgType;
   //typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType> ODEType;

#elif PROBLEM == 3 // discontinuous Advection Diffusion
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   //typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   //typedef DuneODE::ExplTimeStepper<DgType> ODEType;
   typedef DuneODE::ImplTimeStepper<DgType> ODEType;

#elif PROBLEM == 4 // discontinuous Advection  
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
#endif

#endif // endif ADVDIFF_MODELS_HH
