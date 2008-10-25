/** \file
    \brief stokes.cc
 */

/***********************************************************************************************
 
   Sourcefile:  bucklev.cc

   Titel:       main for Buckley-Leverett problem

   Decription:  the parameters for the simulation of a convection diffusion problem are 
                read from a parameter file and the calculation of the numerical solution 
                is started.

                - The problem description is done in the files convdiff.hh and convdiff.cc
                - The time-explicit solver is defined in ../../src/solver/thetascheme.cc
                - The physical parameters for the convection diffusion problem are to be 
                  find in ./ConvDiffData


***********************************************************************************************/

#include <iostream>
#include <fstream>

#include <config.h>
//#include "vtkout.hh"
#include <dune/common/stdstreams.cc>
#include <dune/common/misc.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

//#include <dune/grid/utility/gridtype.hh>

using namespace Dune;
const int dim = dimworld;
const int ncomp = dimworld;

#define SMOOTH 1
#define SNG 0


#include "discretization.hh"

using namespace LDGExample;

//namespace LDGExample {

template <typename Field, class GridImp, int dimR, int polOrd=0 >
struct DescriptionTraits
{
  enum { dim      = GridImp::dimension };
  enum { dimworld = GridImp::dimensionworld };
  
  enum { dimRange = dimR }; 

  typedef Field   FieldType;
  typedef GridImp GridType;

  // the physical datas 
  //  typedef Data DataType;

  // the model 
  typedef ModelParam  <FieldType,GridType,dimRange>  ModelParamType;

  //  typedef Model       <ModelParamType>  ModelType;
  typedef FunctionSpace<FieldType,FieldType,dimworld,dimR> FunctionSpaceType;
  typedef FunctionSpace<FieldType,FieldType,dimworld,1> PressureSpaceType;
  
  typedef MatrixFunctionSpace<double,double,dimworld,dimR,dimworld> GradientSpaceType;

#if SNG  
  typedef StokesSingModel<ModelParamType,FunctionSpaceType,PressureSpaceType> ModelType;  
#endif
#if SMOOTH  
  typedef StokesModel<ModelParamType,FunctionSpaceType,PressureSpaceType> ModelType; 
#endif  
// the discretisation 
  typedef DiscrParam  <ModelType,GradientSpaceType, polOrd>  DiscrParamType;
;

};




// main routine for simulation
  
int main (int argc, char **argv)
{ 
  // error message if called without parameter file
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <parameter file>\n",argv[0]);
  }
  
  // read data from the parameter file
  const char * paramname = "parameter";
  if(argc == 2) paramname = argv[1];

  std::string paramfile ( paramname );
  





  //initialize the problem
 
  typedef DescriptionTraits <double,GridType,ncomp,1> DescrType;
  typedef DescrType::ModelType ModelType;

  typedef DescrType::DiscrParamType DiscrParamType;
  //typedef Discretization <DiscrParamType>  DiscrType;
  
  ModelType model;
  //Timer timer;
  simul<DiscrParamType>(model,paramfile);
  //std::cout << "CPUtime: " << timer.elapsed() << " sec!\n";
  
  //DiscrType discr(model, paramfile);
  
  
  return 0;
}

