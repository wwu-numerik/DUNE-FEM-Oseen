/** \file problem.hh
    \brief contains a class Problem with traitsclass ProblemTraits
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

/**
 *  \brief containing typedefs needed by Problem
 *  \tparam int gridDim dimension of the grid
 **/
template < int gridDim >
class ProblemTraits
{
    public:
        typedef Velocity< gridDim >
            VelocityType;
        typedef Pressure< gridDim >
            PressureType;
        typedef Force< gridDim >
            ForceType;
        typedef DirichletData< gridDim >
            DirichletDataType;
};


template < int gridDim >
class Problem
{
    public:
        typedef ProblemTraits< gridDim >
            Traits;
        typedef typename Traits::VelocityType
            VelocityType;
        typedef typename Traits::PressureType
            PressureType;
        typedef typename Traits::ForceType
            ForceType;
        typedef typename Traits::DirichletData
            DirichletDataType;

    private:
};

#endif  // end of problem.hh
