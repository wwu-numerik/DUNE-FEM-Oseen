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

/**
 *  \brief  a collection of some analytical functions solving a stokes problem
 *  namely velocity, pressure, force term and dirichlet boundary data
 **/
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
        typedef typename Traits::DirichletDataType
            DirichletDataType;

    /**
     *  \brief  constructor
     *  doing nothing
     **/
    Problem()
    {
    }

    /**
     *  \brief destructor
     *  doing nothing
     **/
    ~Problem()
    {
    }

    /**
     *  \brief to get the velocity
     *  \return VelocityType velocity
     **/
    VelocityType velocity()
    {
        return velocity_;
    }

    /**
     *  \brief to get the pressure
     *  \return PressureType pressure
     **/
    PressureType pressure()
    {
        return pressure_;
    }
    /**
     *  \brief to get the force term
     *  \return ForceType force
     **/
    ForceType force()
    {
        return force_;
    }
    /**
     *  \brief to get the dirichlet boundary data
     *  \return DirichletDataType dirichlet boundary data
     **/
    DirichletDataType dirichletData()
    {
        return dirichletData_;
    }

    private:
        VelocityType velocity_;
        PressureType pressure_;
        ForceType force_;
        DirichletDataType dirichletData_;
};

#endif  // end of problem.hh
