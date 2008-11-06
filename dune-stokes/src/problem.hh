/** \file problem.hh
    \brief
 **/

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include "logging.hh"
#include "velocity.hh"
#include "pressure.hh"
#include "force.hh"
#include "dirichletdata.hh"

template < int grid_dim >
class Problem
{
    public:
        typedef Velocity< grid_dim >
            VelocityType;
        typedef Pressure< grid_dim >
            PressureType;
        typedef Force< grid_dim >
            ForceType;
        typedef DirichletData< grid_dim >
            DirichletDataType;

    private:
};

#endif  // end of problem.hh
