#ifndef DUNE_STOKES_SOLVER_INTERFACE_HH
#define DUNE_STOKES_SOLVER_INTERFACE_HH

// OEMBICGSQOp will NOT compile
#ifndef INNER_CG_SOLVERTYPE
	#define INNER_CG_SOLVERTYPE OEMCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
	#define OUTER_CG_SOLVERTYPE OEMCGOp
#endif


#ifdef USE_BFG_CG_SCHEME
	#include <utility>
	//< iteration no , < absLimit, residuum > >
	typedef std::pair<int,std::pair<double,double> >
		IterationInfo;
	#include <dune/stokes/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE DuneStokes
#else
	#include <dune/fem/solver/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE Dune
#endif

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/stokes/cghelper.hh>

#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>

#include <cmath>



namespace Dune {
//! utility struct used to expose runtime statistics
struct SaddlepointInverseOperatorInfo {
	double iterations_inner_avg;
	int iterations_inner_min;
	int iterations_inner_max;
	int iterations_outer_total;
	double max_inner_accuracy;
};

} //end namespace Dune

#endif // DUNE_STOKES_SOLVER_INTERFACE_HH
