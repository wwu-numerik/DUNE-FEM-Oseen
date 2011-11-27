#ifndef DUNE_STOKES_SOLVER_DEFINES_HH
#define DUNE_STOKES_SOLVER_DEFINES_HH

#ifdef HAVE_CMAKE_CONFIG
	#include "cmake_config.h"
#endif

#ifndef STOKES_USE_ISTL

// OEMBICGSQOp will NOT compile
#ifndef INNER_CG_SOLVERTYPE
	#define INNER_CG_SOLVERTYPE OEMCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
	#define OUTER_CG_SOLVERTYPE OEMCGOp
#endif

#if defined(USE_BFG_CG_SCHEME) || defined(FORCE_CUSTOM_SOLVER)
	#include <utility>
	//< iteration no , < absLimit, residuum > >
	typedef std::pair<int,std::pair<double,double> >
		IterationInfo;
	#include <dune/stokes/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE DuneStokes
	#define SOLVER_INTERFACE_NAMESPACE StokesOEMSolver
#else
	#include <dune/fem/solver/oemsolver/oemsolver.hh>
	#define SOLVER_NAMESPACE Dune
	#define SOLVER_INTERFACE_NAMESPACE OEMSolver
#endif

#else //ifndef STOKES_USE_ISTL

#include <dune/fem/solver/istlsolver.hh>

#ifdef USE_BFG_CG_SCHEME
#   warning "undefining bfg for istl solver backends"
#   undef USE_BFG_CG_SCHEME
#endif

#ifndef INNER_CG_SOLVERTYPE
#   define INNER_CG_SOLVERTYPE ISTLCGOp
#endif

#ifndef OUTER_CG_SOLVERTYPE
#   define OUTER_CG_SOLVERTYPE ISTLCGOp
#endif

#define SOLVER_NAMESPACE Dune
#define SOLVER_INTERFACE_NAMESPACE OEMSolver

#endif

#endif // DUNE_STOKES_SOLVER_DEFINES_HH
