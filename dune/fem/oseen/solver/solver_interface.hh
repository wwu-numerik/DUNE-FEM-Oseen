#ifndef DUNE_OSEEN_SOLVER_INTERFACE_HH
#define DUNE_OSEEN_SOLVER_INTERFACE_HH

#include <cmake_config.h>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/oseen/assembler/ported_matrixobject.hh>
#include <dune/fem/oseen/solver/cghelper.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/common/logging.hh>

#include <cmath>
#include <boost/utility.hpp>

namespace Dune {
//! utility struct used to expose runtime statistics
struct SaddlepointInverseOperatorInfo {
	double iterations_inner_avg;
	int iterations_inner_min;
	int iterations_inner_max;
	int iterations_outer_total;
	double max_inner_accuracy;

	SaddlepointInverseOperatorInfo()
		:iterations_inner_avg(-1.0f),iterations_inner_min(-1),
		iterations_inner_max(-1),iterations_outer_total(-1),
		max_inner_accuracy(-1.0f)
	{}
};


} //end namespace Dune

#endif // DUNE_OSEEN_SOLVER_INTERFACE_HH

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

