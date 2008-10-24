/**************************************************************************
**  gaussintegration.cc
**  demo program for use of element quadratures und intersection iterators
**  program computes an integral 
**  int_Omega div*sigma  by simple element decomposition and gauss-formula
**   
**  Bernard Haasdonk 19.6.2006
**************************************************************************/


#include <iostream>
#include <config.h>
#include <dune/common/stdstreams.cc>
#include <dune/common/exceptions.hh>

// include discontinuous lagrange functions
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/referenceelements.hh>

// #if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
// #endif

using namespace Dune;

// polynom approximation order of quadratures, 
// at least polynom order of basis functions 
const int polOrd = POLORDER;

//! The GridPart provides the part of the grid we are working on
//! by choosing the right index set 
typedef LeafGridPart < GridType > GridPartType; 

//! define the abstract function space, \f[ \R^2 \rightarrow \R^2 \f]
//! the function spacces defines the RangeType,and the Types of the Jacobian and Hessian
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, GridType::dimensionworld, GridType::dimensionworld> FunctionSpaceType;

//! define the function space our unkown belongs to 
//! Template parameters are: 
//! the abstrcat functionspace
//! the GridPart
//! the Polynomialorder of the DiscreteSpace
//! the storage method
//

typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType,	polOrd,CachingStorage> DiscreteFunctionSpaceType;
//! define the type of discrete function belonging to DiscreteFunctionSpaceType
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//! the exact function to the problem for projection 
//! Second template parameter is for Barton-Neckman 
class VectorFunction : public Function < FunctionSpaceType , VectorFunction > 
{
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;
public:
  VectorFunction (FunctionSpaceType &f) : 
    Function < FunctionSpaceType , VectorFunction > ( f ) {}
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret[0] = 1.0; // normalization factor
    for(int i=0; i<DomainType::dimension; i++)
      ret[0] *= x[i];
    ret[1] = 1.0; // normalization factor
    for(int i=0; i<DomainType::dimension; i++)
	      ret[1] *= (0.5-x[i]);
  }
  void evaluate (const DomainType & x , RangeFieldType time , 
		 RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};

/**************************************************************************/
// class L2Projection
//
// the class performing an L2-projection of the Vectorfunction to a 
// vector-valued Discretefunction (minor modification of preceding
// exercise)
/**************************************************************************/
template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

 public:
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
      FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    //! The IteratorType of the DiscreteFunctionSpace is a leafIterator,
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;

    const FunctionSpaceType& space =  discFunc.space();

    //! set discfunc to zero
    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename FunctionSpaceType::RangeType fval (0.0);
    double phi_times_f (0.0);
    
    
    Iterator it = space.begin();
    Iterator endit = space.end();

    // Get quadrature rule
    const int codim = 0;
    CachingQuadrature<GridPartType, codim> quad(*it, 2*polOrd);
    
    for( ; it != endit ; ++it) {
      LocalFuncType lf = discFunc.localFunction(*it);
      const typename FunctionSpaceType::BaseFunctionSetType & baseset =
        lf.baseFunctionSet();
      const typename GridType::template Codim<0>::Entity::Geometry& 
	itGeom = (*it).geometry();
      for(int qP = 0; qP < quad.nop(); qP++) {
	f.evaluate(itGeom.global(quad.point(qP)), fval);
	for(int i=0; i<lf.numDofs(); i++) {
          phi_times_f = baseset.evaluateSingle(i,quad,qP,fval);
          lf[i] += quad.weight(qP) * phi_times_f ;
        }
      }
    }
  }
};

/**************************************************************************/
// class StandardDivergenceIntegrator 
//
// computing the integral of the divergence by
// sum_e  int_e  div f
/**************************************************************************/
template <class DiscreteFunctionType>
class StandardDivergenceIntegrator
{
public:
  template <int polOrd>
  double integrate (DiscreteFunctionType &discFunc, double time)
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::RangeType RangeType;

    const typename DiscreteFunctionType::FunctionSpaceType
      & space = discFunc.space();
    
    RangeType ret (0.0);
    JacobianRangeType jacobian (0.0);
        
    double sum = 0.0;
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();
    
    for(; it != endit ; ++it)
      {
	const int codim = 0;
	CachingQuadrature<GridPartType,0> quad(*it, polOrd);
	LocalFuncType lf = discFunc.localFunction(*it);
	GeometryType geo = (*it).geometry().type();
	for(int qP = 0; qP < quad.nop(); qP++)
	  {
	    double det = (*it).geometry().integrationElement(quad.point(qP));
	    lf.jacobian(quad.point(qP),jacobian);
	    //lf.jacobianLocal((*it),quad.point(qP),jacobian);
	    double trace = 0.0;
	    for (int i = 0;i < RangeType::dimension; ++i)
	      trace += jacobian[i][i];
	    sum += det * quad.weight(qP) * trace;
	  }
      }
    return sum;
  }
};


/**************************************************************************/
// class GaussDivergenceIntegrator
//
// class computing the integral of the divergence by
// sum_e  int_boundary_e  f * n
/**************************************************************************/
template <class DiscreteFunctionType>
class GaussDivergenceIntegrator
{
public:
  template <int polOrd>
  double integrate (DiscreteFunctionType &discFunc, double time)
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename GridType::template Codim<1>::LocalGeometry LocalGeometryType;    
    typedef typename EntityType::LeafIntersectionIterator IntersectionIteratorType;
    typedef typename GridType::ctype CoordType;

    // Elementquadrature with codim 1
    typedef CachingQuadrature<GridPartType, 1> QuadratureType;
        const typename DiscreteFunctionType::FunctionSpaceType
      & space = discFunc.space();
    
    RangeType fval (0.0);

    double sum = 0.0;
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();
    
    for(; it != endit ; ++it) {
	IntersectionIteratorType endnit=it->ileafend();
	for (IntersectionIteratorType nit = it->ileafbegin(); nit!=endnit; ++nit) {
	    QuadratureType quad(space.gridPart(),nit, polOrd, QuadratureType::INSIDE);
	    LocalFuncType lf = discFunc.localFunction(*it);
	    FieldVector<CoordType, DomainType::dimension-1> local = 0.0;
	    DomainType integrationNormal = nit.integrationOuterNormal(local);
	    for(int qP = 0; qP < quad.nop(); qP++) {
                DomainType global = it->geometry().global(quad.point(qP));                    
	        lf.evaluate(quad,qP,fval) ;		
		sum += fval * integrationNormal * quad.weight(qP);
	    }
	}
    }
    return sum;
  }
};

//**************************************************
//  main programm
//**************************************************
int main (int argc, char **argv)
{
  // initialize MPI if available 
  MPIHelper::instance(argc,argv);

// try and catch exceptions 
try {
    if(argc != 3)
    {
      fprintf(stderr,"usage: %s <gridfilename> <maxlevel> \n",argv[0]);
      exit(1);
    }
    int ml = atoi( argv[1] );
    GridPtr<GridType> grid(argv[1]);  
    
    double error[2];
    const int dim = GridType::dimension; 
    int n[dim];
    double h[dim];
    for(int i=0; i<dim; i++)  { n[i] = 2; h[i] = 1.0; }
    grid->globalRefine(ml);
    bool visualize = false;

    GridPartType part ( *grid);
    DiscreteFunctionSpaceType discFuncSpace ( part );
    DiscreteFunctionType discfunc ( "DGfunc", discFuncSpace );
    discfunc.clear();
    
    VectorFunction f ( discFuncSpace ); 
    
    //! perform l2-projection componentwise for getting a DG-vector-function
    L2Projection<DiscreteFunctionType, VectorFunction, polOrd>::
      project(f, discfunc);
    
    // calculation of divergence with standard decomposition
    std::cout << "integrating standard divergence integral : " << std::flush;
    StandardDivergenceIntegrator < DiscreteFunctionType > stdDivIntegrator;
    double stddivint = stdDivIntegrator.integrate<polOrd>(discfunc,0.0);
    std::cout << stddivint << "\n";

    // calculation of divergence with gauss-formula
    std::cout << "integrating gauss divergence integral : " << std::flush;
    GaussDivergenceIntegrator < DiscreteFunctionType > gaussDivIntegrator;
    double gaussdivint = gaussDivIntegrator.integrate<polOrd>(discfunc,0.0);
    std::cout << gaussdivint << "\n";
    
#if HAVE_GRAPE
    // if Grape was found, then display last discretefunction 
    if(visualize)
    {
      GrapeDataDisplay < GridType > grape(*grid); 
      grape.dataDisplay( discfunc );
    }
#endif
  }
  catch (Dune::Exception &e) 
  {
    std::cerr << e << std::endl;
    return 1;
  } 
  catch (...) 
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 2;
  }
  
  return 0;
}

