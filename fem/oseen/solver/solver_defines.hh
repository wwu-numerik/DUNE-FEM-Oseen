#ifndef DUNE_OSEEN_SOLVER_DEFINES_HH
#define DUNE_OSEEN_SOLVER_DEFINES_HH

#ifdef HAVE_CMAKE_CONFIG
	#include "cmake_config.h"
#endif

#if ! STOKES_USE_ISTL

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
	#include <dune/fem/oseen/oemsolver/oemsolver.hh>
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

#endif // DUNE_OSEEN_SOLVER_DEFINES_HH

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

