/**************************************************************************
**       Title: poisson.cc
**    $RCSfile$
**   $Revision: 1988 $$Name$
**       $Date: 2007-08-01 22:53:29 +0200 (Wed, 01 Aug 2007) $
**   Copyright: GPL Author: robertk 
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
**              make GRIDTYPE=YASPGRID       (default)
**                    -> compiles and works correctly
**              make
**                    -> compiles and works correctly
**              make GRIDTYPE=ALBERTAGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=SGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=ALUGRID_SIMPLEX
**                    -> compiles and works correctly
**
**
**              Similarly, the polynomial order can be specified by
**
**              make POLORDER=2
**
**************************************************************************/

#include <config.h>

//#define USE_DUNE_ISTL HAVE_DUNE_ISTL
#define USE_DUNE_ISTL 0

//- system includes
#include <iostream>
#include <sstream>

//- Dune includes 
// grid parts 
//#include <dune/grid/common/gridpart.hh>
#include <dune/fem/gridpart/gridpart.hh>
// mpi helper 
#include <dune/common/mpihelper.hh>

// dgf parser 
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

// graphical output
// uncomment the following line to use grape
//#define USE_GRAPE HAVE_GRAPE
#if HAVE_GRAPE
// Grape
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
// VTK 
#include <dune/fem/io/file/vtkio.hh>

// discrete function space 
#include <dune/fem/space/lagrangespace.hh>
// discrete function 
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>

// solvers 
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>

#if USE_DUNE_ISTL
#warning Using DUNE-Istl 
#include <dune/fem/solver/istlsolver.hh>
#endif

// L2 Error 
#include <dune/fem/misc/l2error.hh>

//- local includes
// laplace operator 
#include "laplace.hh"
// mass matrix 
#include "massmatrix.hh"
// problem description 
#include "problem.hh"
// boundary treatment 
#include "boundary.hh"


#ifndef POLORDER
  #define POLORDER 1
#endif

//***********************************************************************
/*! Poisson problem: 

  This is an example solving the poisson problem
  \f{displaymath}
  \begin{array}{rcll}
  -\triangle u &=& f & \quad \mbox{in}\ \Omega\\
             u &=& 0 & \quad \mbox{on}\ \partial\Omega
  \end{array}
  \f}
  with the finite element method using Lagrangian elements. The polynomial
  order is given by POLORDER.
  
  In this example, $\Omega = ]0,1[^{dimworld}$ and
  \f[
  f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
  \f]

  The exact solution to the poisson problem is then given by
  \f[
  u( x, y, z ) = \prod_{i=1}^{dimworld} (x_i - x_i^2).
  \f]
*/
//***********************************************************************

using namespace Dune;


template < class GridImp,
           template <class> class ProblemImp,
           template <class> class OperatorImp >
struct Algorithm
{
	
  //******** typedefs **********

  // type of grid 
  typedef GridImp GridType;

  // GridPart is build up by all leaf elements
  typedef LeafGridPart< GridType > GridPartType;

  // type of analytical function space, here f: R^dim --> R
  typedef FunctionSpace<double, double, GridType::dimension, 1>  FunctionSpaceType;

  // define the discrete function space our unkown belongs to, here we use a lagrange space
  typedef LagrangeDiscreteFunctionSpace
     < FunctionSpaceType, GridPartType, POLORDER >  DiscreteFunctionSpaceType;

  // define the type of discrete function we are using
  #if USE_DUNE_ISTL
    typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  #else
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  #endif

  // problem type
  typedef ProblemImp<FunctionSpaceType> ProblemType;

  // The data functions (as defined in problem.hh)
  typedef typename ProblemType :: RHSFunctionType RHSFunctionType;
  typedef typename ProblemType :: ExactSolutionType ExactSolutionType;

  // define the discrete operator, see laplace.hh and massmatrix.hh
  typedef OperatorImp< DiscreteFunctionType > OperatorType;

  // define the inverse operator we are using to solve the system
  // see also dune/fem/solver/oemsolver/oemsolver.hh or dune/fem/solver/inverseoperators.hh
  #if USE_DUNE_ISTL
    //typedef ISTLCGOp <DiscreteFunctionType, OperatorType > InverseOperatorType;
    typedef ISTLBICGSTABOp <DiscreteFunctionType, OperatorType > InverseOperatorType;
  #else
    //typedef CGInverseOperator<DiscreteFunctionType> InverseOperatorType;
    //typedef OEMCGOp<DiscreteFunctionType,OperatorType> InverseOperatorType;
  typedef OEMBICGSTABOp<DiscreteFunctionType,OperatorType> InverseOperatorType;
    //typedef OEMBICGSQOp<DiscreteFunctionType,OperatorType> InverseOperatorType;
    //typedef OEMGMRESOp<DiscreteFunctionType,OperatorType> InverseOperatorType;
  #endif


  //******** algorithm **********
  static double calc (GridType& grid, 
                      const std::string name,
                      bool verbose = false )
  {
    // **** Initialisation stuff
    // create grid part
    GridPartType gridPart( grid );
    // create discrete function space
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    // print info
    std::cout << std :: endl << "Solving for " << discreteFunctionSpace.size()
              << " unkowns and polynomial order "
              << DiscreteFunctionSpaceType :: polynomialOrder << "." 
              << std :: endl << std :: endl;


    // create discrete function for solution
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();

    // create discrete function for right hand side
    DiscreteFunctionType rhs( "rhs", discreteFunctionSpace );
    rhs.clear();

    // create analytical right hand side
    RHSFunctionType f( discreteFunctionSpace ); 

    // **** build right hand side (see boundary.hh)
    assembleRHS<2*DiscreteFunctionSpaceType::polynomialOrder> ( f , rhs );

    // **** set the dirichlet boundary points to the corresponding values of the
    //      exact solution u (see boundary.hh, problem.hh)
    ExactSolutionType u( discreteFunctionSpace ); 
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType; 
    const IteratorType endit = discreteFunctionSpace.end();
    for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    {
      boundaryTreatment( *it , u , rhs );
    }


    // create problem operator (Laplace or MassMatrix Operator)
    OperatorType problemOperator( discreteFunctionSpace );

    // **** solve the linear system (with CG)
    InverseOperatorType cg( problemOperator, 1e-8, 1e-8, 20000, verbose );
    cg( rhs, solution );


    // **** calculation of L2 error (see dune/fem/misc/l2error.hh)
    // polynomial order for this calculation should be higher than the polynomial
    // order of the base functions
    L2Error< DiscreteFunctionType > l2error;
    typename DiscreteFunctionSpaceType :: RangeType error = l2error.norm( u, solution );
    std :: cout << "L2 Error: " << error << std :: endl << std :: endl;


    // **** graphical output
#if USE_GRAPE
    // if grape was found then display solution 
    if( verbose ) 
    {
      GrapeDataDisplay < GridType > grape( grid );
      grape.dataDisplay( solution );
    }
#endif
    VTKIO<GridPartType> vtkio( gridPart );
    vtkio.addVertexData( solution ); 
    vtkio.write(name.c_str(),Dune::VTKOptions::ascii);


    // **** return calculated error
    return error[ 0 ];
  }
};



// main programm, run algorithm twice to calc EOC 
int main( int argc, char **argv )
{
  // initialize MPI 
  MPIHelper::instance(argc,argv);

  // try-catch block to catch exceptions 
  try 
  {
    if( argc != 2 ) 
    {
      std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
      return 1;
    }

    int level = atoi( argv[ 1 ] );
    double poissonError[ 2 ];
    double projectionError[ 2 ];


    // create macro grid name 
    std::stringstream dgfFileName;
    dgfFileName << "../macrogrids/unitcube" << GridType :: dimension << ".dgf";
    std :: cout << "loading dgf: " << dgfFileName << std :: endl;

    // create grid pointer 
    GridPtr<GridType> gridPtr( dgfFileName.str() );

    // get grid reference 
    GridType& grid = *gridPtr;

    // make start refinement 
    grid.globalRefine( level * DGFGridInfo<GridType>::refineStepsForHalf());

    for( int i = 0; i < 2; ++i )
    {
      // calculate L2-projection 
      projectionError[ i ] = Algorithm<GridType,L2ProjectionProblem, MassOperator>::calc(grid,"l2-projection", false);

      // calculate poisson problem 
      poissonError[ i ] = Algorithm<GridType,PoissonProblem, LaplaceOperator>::calc(grid,"poisson",true);

      // make refinement for next step 
      grid.globalRefine( DGFGridInfo<GridType>::refineStepsForHalf() ); 
    }

    // calculate EOC
    const double poissonEOC    = log( poissonError[ 0 ] / poissonError[ 1 ] ) / M_LN2;
    const double projectionEOC = log( projectionError[ 0 ] / projectionError[ 1 ] ) / M_LN2;

    std :: cout << "Poisson EOC = " << poissonEOC << std :: endl;
    std :: cout << "L2-Projection EOC = " << projectionEOC << std :: endl;

    return 0;
  } 
  catch( Exception& exception ) 
  {
    std :: cerr << exception << std :: endl;
  }
}
