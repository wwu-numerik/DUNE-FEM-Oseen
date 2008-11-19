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
template < int griddim, class VelocityFunctionSpaceImp, class PressureFunctionSpaceImp >
class ProblemTraits
{
    public:
        static const unsigned int gridDim = griddim;
        typedef VelocityFunctionSpaceImp
            VelocityFunctionSpaceType;
        typedef Velocity< VelocityTraits < gridDim, VelocityFunctionSpaceType > >
            VelocityType;
        typedef PressureFunctionSpaceImp
            PressureFunctionSpaceType;
        typedef Pressure< PressureTraits< gridDim, PressureFunctionSpaceType > >
            PressureType;
        typedef Force< ForceTraits< gridDim, VelocityFunctionSpaceType > >
            ForceType;
        typedef DirichletData< DirichletDataTraits< gridDim, VelocityFunctionSpaceType > >
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
template < class TraitsImp >
class Problem
{
    public:
        typedef TraitsImp
            Traits;
        typedef typename Traits::VelocityType
            VelocityType;
        typedef typename Traits::PressureType
            PressureType;
        typedef typename Traits::ForceType
            ForceType;
        typedef typename Traits::DirichletDataType
            DirichletDataType;
        typedef typename VelocityType::FunctionSpaceType
            VelocityFunctionSpaceType;
        typedef typename PressureType::FunctionSpaceType
            PressureFunctionSpaceType;

        static const unsigned int gridDim = Traits::gridDim;
    /**
     *  \brief  constructor
     *
     *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
     **/
    Problem( const double viscosity, const VelocityFunctionSpaceType& velocity_space, const PressureFunctionSpaceType& press_space )
        : velocity_( velocity_space ),
          pressure_ ( press_space ),
          force_( viscosity, velocity_space ),
          dirichletData_( velocity_space )

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
    const VelocityType& velocity() const
    {
        return velocity_;
    }

    /**
     *  \brief  to get the pressure
     *
     *  \return pressure
     **/
    PressureType& pressure() const
    {
        return pressure_;
    }

    /**
     *  \brief  to get the force term
     *
     *  \return force
     **/
    ForceType& force() const
    {
        return force_;
    }

    /**
     *  \brief  to get the dirichlet boundary data
     *
     *  \return dirichlet boundary data
     **/
    DirichletDataType& dirichletData() const
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
        DirichletDataType dirichletData_;
};

#endif // end of problem.hh
