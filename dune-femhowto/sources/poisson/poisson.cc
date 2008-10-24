/**************************************************************************
**       Title: poisson.cc
**    $RCSfile$
**   $Revision: 1649 $$Name$
**       $Date: 2007-05-31 09:58:52 +0200 (Thu, 31 May 2007) $
**   Copyright: DUNE community
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
**                    -> compiles and works correctly
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

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/misc/l2error.hh>

//- local includes 

#include "operators/ellipticfeop.hh" 

#ifndef POLORDER
  #define POLORDER 1
#endif


using namespace Dune;

//***********************************************************************
/*! Poisson problem: 

  This is an example how to solve the equation on 
      
           \Omega = (0,1)^dimworld

    -\triangle u  =  f    in   \Omega
               u  =  0    on   \partial\Omega
	       
         f(x,y,z) = 2 (x-x^2) (y-y^2) +
                    2 (z-z^2) (y-y^2) +    
                    2 (x-x^2) (z-z^2)
  

  An exact solution to this problem is given by 
         
	 u(x,y,z) = ( x - x^2 ) ( y - y^2 ) ( z - z^2 )

  with the finite element method using lagrangian elements of polynom order 1.
*/

//***********************************************************************

// First of all we need some 'forward declaration' of the class 'Tensor'.
// That means we tell the program that there WILL BE a class 'Tensor' (declared in line 137).
class Tensor; // (1)
// It is important that 'Tensor' is a known class to define for instance:
//  'typedef LaplaceFEOp< DiscreteFunctionType, Tensor > LaplaceOperatorType;',
// but it is not yet important HOW 'Tensor' is defined. 
// Only when we want to use it explictly it is getting important and there should be a real declaration!

class MassTerm; // see 'class Tensor'

//! GridPartType will be the index set we are using. This index set describes the total number of simplices after a certain number of refinement-steps.
//! We can say: LeafGridPart<GridType> selects the leaf level of the grid.
typedef LeafGridPart<GridType> GridPartType;


//! define the function space of elements ' v :  \Omega  \rightarrow  \R '
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double , dimworld , 1 > FuncSpace; // (2)
//Let 'v : domain --> codomain' be an element of FuncSpace.
//'FunctionSpace <double,double,dimworld,1>' declares what kind of Function Space (FuncSpace) we are using.
//domain and codomain must be sets of 'vectors' (Fieldvectors!). (Typically: domain contains 'vectors' with two or three entries, codomain 'vectors' with one entry)
//Now for each 'vector' (or each set of 'vectors'), the type of a vector entry must be declared.
//Therfore FunctionSpace recieves the arguments: < double , double , dimworld , 1 >, which stand for:
//< double (= type of vector entries of a domain-vector), double (= type of vector entries of a codomain-vector), dimworld (= dimension of domain = number of entries of a domain-vector), 1 (= dimension of codomain = number of entries of a codomain-vector)>


//! define the type of elements of the domain \Omega
typedef FuncSpace::DomainType DomainType; // (3)
//'FuncSpace::DomainType' gets the type of elements of the domain, i.e. all possible fieldvectors (of dimension 1, or 2, or 3, ...) and the type of the vector entries (double, int, ...).


//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FuncSpace::RangeType RangeType; // (4)
//'FuncSpace::RangeType' gets the type of elements of the codomain, i.e. all possible fieldvectors (of dimension 1, or 2, or...) and the type of the vector entries (double, int, ...).


//! defines the function space to which the numerical solution belongs to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace
        < FuncSpace, GridPartType, POLORDER, CachingStorage >
  DiscreteFunctionSpaceType;
//< FuncSpace (= the Function Space that is discretisized), GridPartType, POLORDER (= polynomial order of the discretization), CachingStorage >


//! defines the type of discrete function we are using, see
//! dune/fem/function/adaptivefunction/adaptivefunction.hh
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//Discrete Function is made dependendent on the refinments

//! defines the discrete laplace operator, see "dune-femhowto/sources/poisson/operators"
typedef EllipticFEOp< DiscreteFunctionType, Tensor, MassTerm> LaplaceOperatorType;
// Tensor is just a dummy! To implement the Laplace-operator you do not need no tensor/matrix at all. Not even the unity-matrix is used, since the Laplacian can be instantiated directly.

//! finding the inverse operator that we are using to solve the system 
//the following alternativs are possible to define/find the inverseoperator:

//see dune/fem/operator/inverseoperators.hh for the following alternative: 
//typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType; //the operator is found by using a CG method
/****************************************/
// or see dune/fem/solver/oemsolver/oemsolvers.hh for the following alternatives:
typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType; //the operator is also found by using a CG method. it is equal to the first alternative
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType; //The operator is found by using a Bi-CG method. It can be used if the matrix is not symmetric. See oemsolvers.hh for more details.
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;

#include "problems/problem_1.hh" //includes the definition of the 'example-problem' that we want to solve
//We need to include problem_1.hh at this position because it uses some of the "typedef's" that we just declared. If we included problem_1.hh earlier, the required "typedef's" would have been unknowns and therefor would lead to errors.

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
  
  // by 'discrete_f' we mean: the right hand side of the linear system of equations that we have to solve in order to compute 'u_h' (the discrete solution of the poisson problem)
  DiscreteFunctionType discrete_f( "discrete_f", discreteFunctionSpace );
  discrete_f.clear();
      
  RHSFunction f( discreteFunctionSpace );
   
  LaplaceOperatorType laplace( discreteFunctionSpace, 
                               LaplaceOperatorType :: ASSEMBLED );
  //we initialize the LaplaceOperatorType element (laplace) without tensors! To initialize it with tensors there need to be a stiff tensor and a mass term instantiated separately. For instance: LaplaceOperatorType laplace( stifftensor, massterm, discreteFunctionSpace, LaplaceOperatorType :: ASSEMBLED);
   
   //! build right hand side, does not allocate b!
  RightHandSideAssembler< DiscreteFunctionType >
    :: assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder >( f , discrete_f );
    
  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType; 
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    boundaryTreatment( *it , discrete_f );
    
  //laplace.print();
  //discrete_f.print( std::cout );
   
  // solve the linear system (with CG)
  double dummy = 12345.67890;
  InverseOperatorType cg( laplace, dummy, 1e-8, 20000, VERBOSE );
  cg( discrete_f, solution );

  // calculation of L2 error
  // polynomial order for this calculation should be higher than the polynomial
  // order of the base functions
  ExactSolution u( discreteFunctionSpace ); 
  L2Error< DiscreteFunctionType > l2error;
  DiscreteFunctionSpaceType :: RangeType error = l2error.norm( u, solution );
  std :: cout << "L2 Error: " << error << std :: endl << std :: endl;
   
  #if USE_GRAPE
  // if grape was found then display solution 
  if( turn > 0 ) {
    GrapeDataDisplay < GridType > grape( *gridptr );
    grape.dataDisplay( solution );
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
