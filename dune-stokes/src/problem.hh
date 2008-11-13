/**
 *  \file   problem.hh
 *
 *  \brief  contains a class Problem with traitsclass ProblemTraits
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
 *  \brief  containing typedefs needed by Problem
 *
 *  \tparam gridDim
 *          dimension of the grid
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
 *
 *  namely velocity, pressure, force term and dirichlet boundary data
 *
 *  \tparam gridDim
 *          dimension of the grid
 *
 *  \todo   extensive docu with latex
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
     *
     *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
     **/
    Problem( const double viscosity )
        :force_( viscosity )
    {
    }

    /**
     *  \brief  destructor
     *
     *  doing nothing
     **/
    ~Problem()
    {
    }

    /**
     *  \brief  to get the velocity
     *
     *  \return velocity
     **/
    VelocityType& velocity()
    {
        return velocity_;
    }

    /**
     *  \brief  to get the pressure
     *
     *  \return pressure
     **/
    PressureType& pressure()
    {
        return pressure_;
    }
    /**
     *  \brief  to get the force term
     *
     *  \return force
     **/
    ForceType& force()
    {
        return force_;
    }
    /**
     *  \brief  to get the dirichlet boundary data
     *
     *  \return dirichlet boundary data
     **/
    DirichletDataType& dirichletData()
    {
        return dirichletData_;
    }

    /**
     *  \brief  a simple test of all class' functionalities
     **/
    void testMe()
    {
        // some logstreams
        Logging::LogStream& infoStream = Logger().Info();
        infoStream << "testing class Problem..." << std::endl;
        //tests
        velocity_.testMe();
        pressure_.testMe();
        force_.testMe();
        dirichletData_.testMe();
        // happy
        infoStream << "...test passed!" << std::endl;
    }

    private:
        VelocityType velocity_;
        PressureType pressure_;
        ForceType force_;
        VelocityType velocity_;
        PressureType pressure_;
        ForceType force_;
        DirichletDataType dirichletData_;

};

#endif // end of problem.hh
