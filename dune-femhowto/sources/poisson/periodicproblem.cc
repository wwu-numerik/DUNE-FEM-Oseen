/**************************************************************************
**       Title: periodicproblem.cc
**    $RCSfile$
**   $Revision: 1649 $$Name$
**       $Date: 2007-05-31 09:58:52 +0200 (Thu, 31 May 2007) $
**   Copyright: GPL Author: robertk/patrickh
** Description: File demonstrating a simple numerics problem on arbitrary
**              grids: poisson-problem with known solution is given
**              and compared with numerical solution (EOC)
**              Dune grid parser is used.
**
**              For changing grid-types, compile with
**
**              make clean
**
**              and one of the following
**           
**              make GRIDTYPE=ALBERTAGRID       (default)
**                    -> compiles and works correctly
**              make
**                    -> compiles and works correctly
**              make GRIDTYPE=YASPGRID
**                    -> does not compile
**              make GRIDTYPE=SGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=ALUGRID_SIMPLEX
**                    -> compiles and works correctly
**
**************************************************************************/

// uncomment the following line to use grape
#define USE_GRAPE HAVE_GRAPE

#define VERBOSE false


//- system includes
#include <iostream>
#include <config.h>

//- Dune includes 
#include <dune/common/stdstreams.cc>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/space/common/periodicgridpart.hh>
#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/misc/l2error.hh>

//- local includes 
#include "operators/ellipticfeop.hh"

#include "problems/problem_2.hh"

#ifndef POLORDER
  #define POLORDER 1
#endif


using namespace Dune;

//***********************************************************************
/*! Essentially we are interested in solving the following problem:
      
            \Omega = (0,1)^dimworld

     -\triangle u  =  f    in   \Omega
               
 	  u is a periodic function
	       
         f(x) = sin(2 Pi x_1) * sin (2 Pi x_2) * ... * sin (2 Pi x_N)
                  

 Unfortunately since we have a periodic boundary condition the solution is no more unique. Generally we have to fix one solution by demanding that the solution has mean value zero. But this seems to be quite awkward for implementation. (In our case) one possiblity for avoiding such difficulties is to solve the following problem instead:

  \epsilon u -\triangle u  =  f    in   \Omega
             
 	     u is a periodic function
  
      f(x) = (1 + \epsilon (1 / ( 4 Pi^2 N ) ) ) ( sin(2 Pi x_1) * sin (2 Pi x_2) * ... * sin (2 Pi x_N) )
	   
 where \epsilon is a very small parameter. Now this problem has a unique solution! 
	     
 If \epsilon is sufficiently small, the influence of '\epsilon u' tends to zero and we have got one possible solution of our orginal problem. 
	       
 
 An exact solution of both problems is given by 
         
	 u(x) = 1 / ( 4 Pi^2 N ) * sin(2 Pi x_1) * sin (2 Pi x_2) * ... * sin (2 Pi x_N)

NOTE: This construction only works because we already know the exact solution. In general both solutions are (faintly) different since we do not know how to adapt f.
In such cases we need to solve problem 2 with the same right hand side (f) like problem 1. 
	 
*/

//***********************************************************************

// For more detailed descriptions see: poisson.cc. The following program was created on the basis of the code of poisson.cc. 


//! GridPartType will be the index set we are using. This index set describes the total number of simplices after a certain number of refinement-steps.
//! We can say: LeafGridPart<GridType> selects the leaf level of the grid. Now we have some periodic grid
typedef PeriodicLeafGridPart<GridType> GridPartType;


typedef FunctionSpace < double , double , dimworld , 1 > FuncSpace; 

typedef Tensor< FuncSpace > TensorType;
typedef MassTerm< FuncSpace > MassTermType;
typedef RHSFunctionPeriodic< FuncSpace > RHSFunctionType;
typedef ExactSolutionPeriodic< FuncSpace > ExactSolutionType;

typedef FuncSpace::DomainType DomainType; 


//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FuncSpace::RangeType RangeType;


//! defines the function space to which the numerical solution belongs to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace
        < FuncSpace, GridPartType, POLORDER, CachingStorage >
  DiscreteFunctionSpaceType;

  
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;


//! defines the discrete laplace operator, see "dune-femhowto/sources/poisson/operators"
typedef EllipticFEOp< DiscreteFunctionType, TensorType, MassTermType >
  EllipticOperatorType;


//! finding the inverse operator that we are using to solve the system 
typedef OEMCGOp<DiscreteFunctionType,EllipticOperatorType> InverseOperatorType; 


typedef DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

double algorithm ( std :: string &filename, int maxlevel, int turn )
{
  GridPtr< GridType > gridptr( filename ); 
  
  gridptr->globalRefine( maxlevel );

  GridPartType gridPart( *gridptr );

  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  std::cout << std :: endl << "Solving for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << DiscreteFunctionSpaceType :: polynomialOrder << "." 
            << std :: endl << std :: endl;
  
  //! define the discrete function 'solution' (of the type 'DiscreteFunctionType') we are using, see	    
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();
  
  DiscreteFunctionType discrete_f( "discrete_f", discreteFunctionSpace );
  discrete_f.clear();
  
  const RangeFieldType epsilon = 0.01;
  //we solve  \epsilon u - \triangle u = f
      
  TensorType stiff( discreteFunctionSpace );
  MassTermType mass( discreteFunctionSpace, epsilon );
  
  RHSFunctionType f( discreteFunctionSpace, epsilon );
    
  EllipticOperatorType elliptic(stiff, mass, discreteFunctionSpace,
                               EllipticOperatorType :: ASSEMBLED, epsilon );
  //In comparison to poisson.cc we have got a stiff tensor and a mass term!
			       
			       
   //! build right hand side, does not allocate b!
  RightHandSideAssembler< DiscreteFunctionType >
    :: assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder >( f , discrete_f );
   
  // solve the linear system (with CG)
  double dummy = 12345.67890;
  InverseOperatorType cg( elliptic, dummy, 1e-8, 20000, VERBOSE); //accurarcy of cg? = 1 e-8 //opMode (operator Mode) = Verbose
  cg( discrete_f, solution );

  // calculation of L2 error
  // polynomial order for this calculation should be higher than the polynomial
  // order of the base functions
  ExactSolutionType u( discreteFunctionSpace ); 
  L2Error< DiscreteFunctionType > l2error;
  DiscreteFunctionSpaceType :: RangeType error = l2error.norm( u, solution );
  std :: cout << "L2 Error: " << error << std :: endl << std :: endl;
  
  DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
    u_discrete( "exact solution", u, discreteFunctionSpace.gridPart() );
   
  #if USE_GRAPE
  // if grape was found then display solution 
  if( turn > 0 ) {
    GrapeDataDisplay < GridType > grape( *gridptr );
    grape.addData( solution );
    grape.addData( u_discrete );
    grape.display();
  }
  #endif

  return error[ 0 ];
}



// main programm, run algorithm twice to calc EOC 
int main( int argc, char **argv )
{
  if( argc != 2 ) {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    return 1;
  }
  
  int level = atoi( argv[ 1 ] );
  double error[ 2 ];

  std :: string macroGridName( "grids/square.dgf" );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;
  
  const int step = DGFGridInfo< GridType > :: refineStepsForHalf();
  level = (level > step ? level - step : 0);
  
  for( int i = 0; i < 2; ++i )
    error[ i ] = algorithm( macroGridName, level + i*step, i );

  const double eoc = log( error[ 0 ] / error[ 1 ] ) / M_LN2;
  std :: cout << "EOC = " << eoc << std :: endl;
  
  return 0;
}
