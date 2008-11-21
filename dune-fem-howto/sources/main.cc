#include <iostream>
#include <config.h>

// polynom approximation order of approximation (fix in Makefile.am)
const int polOrd = POLORDER;

// grid and grid part
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/fem/space/common/adaptiveleafgridpart.hh> 
//! the grid part to use (the grid type is fixed in gridtype.hh 
//! and influenced by the setting in the Makefile.am
typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;

// disctrete function spaces
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
//! define the function space
typedef Dune::FunctionSpace < double , double, dimworld , 2 > FuncSpace;

// discrete function 
#include <dune/fem/discretefunction/dfadapt.hh>

// we need the l2-error for the EOC
#include <dune/fem/misc/l2error.hh>

using namespace Dune;

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

// the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > {
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
   //! f(x1,x2,x3,..,xn) = sin(2pi x1)sinh(x2-0.5)x3...xn
  void evaluate (const DomainType & x , RangeType & ret) const {
    ret = sin(2.*M_PI*x[0])*sinh(2.*M_PI*(x[1]-0.5)); 
    for(int i=2; i<DomainType::dimension; i++)
      ret *= x[i];
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const {
    evaluate ( x , ret );
  }
};
// *************************************************
#include "l2projection.hh"
//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  ExampleType example;

  if(argc != 3)
  {
    fprintf(stderr,"usage: %s grid-file <maxlevel> \n",argv[0]);
    exit(1);
  }
  GridPtr<GridType> gridptr(argv[1]);
  GridType& grid=*gridptr;
  int ml = atoi( argv[2] );
  double* error = new double[ml];
  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();

  typedef ExampleType::DiscreteFunctionType DiscreteFunctionType;
  typedef DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType::GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType uh ( "sol", linFuncSpace );
  ExactSolution u(linFuncSpace);
  L2Error < DiscreteFunctionType > l2err;

  for(int i=0; i<ml; i+=step)
  {
    grid.globalRefine(step);
    DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );
    dm.resize();
    example.solve(u,uh);
    DiscreteFunctionType :: RangeType err = l2err.norm(u , uh);
    error[i] = err.two_norm();
    std::cout << "\n L2 Error: " << error[i] << " : (" 
	      << err << ")\n";
    if (i>0) {
      double eoc = log( error[i-step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  delete [] error;
  return 0;
}

