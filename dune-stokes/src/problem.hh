/**
 *  \file   problem.hh
 *
 *  \brief  contains a class Problem with traitsclass ProblemTraits
 **/

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dune/stuff/logging.hh>
#include "velocity.hh"
#include "pressure.hh"
#include "analyticaldata.hh"

namespace StokesProblem {
	#ifdef SIMPLE_PROBLEM
		static const std::string ProblemIdentifier = "SIMPLE-PROBLEM";
	#elif defined(CONSTANT_PROBLEM)
		static const std::string ProblemIdentifier = "CONSTANT-PROBLEM";
	#elif defined(ROTATE_PROBLEM)
		static const std::string ProblemIdentifier = "ROTATE-PROBLEM";
	#elif defined(MICRO_PROBLEM)
		static const std::string ProblemIdentifier = "MICRO-PROBLEM";
	#elif defined(MICRO_PROBLEM_WOIDS)
		static const std::string ProblemIdentifier = "MICRO-PROBLEM-WOIDS";
	#elif defined(GENRALIZED_STOKES_PROBLEM)
		static const std::string ProblemIdentifier = "GENRALIZED-STOKES-PROBLEM";
	#elif defined(COCKBURN_PROBLEM)
		static const std::string ProblemIdentifier = "COCKBURN-PROBLEM";
	#else
		static const std::string ProblemIdentifier = "UNKNOWN-PROBLEM";
	#endif
}

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
template < int griddim, class DiscreteFunctionWrapperImp, bool hasAnalyticalSolution = false >
class Problem
{
    public:
        static const unsigned int gridDim = griddim;
        static const bool hasMeaningfulAnalyticalSolution = hasAnalyticalSolution;

        typedef DiscreteFunctionWrapperImp
            DiscreteFunctionWrapperType;

        typedef typename DiscreteFunctionWrapperType::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceWrapperType;

        typedef typename DiscreteFunctionSpaceWrapperType
                ::DiscreteVelocityFunctionSpaceType
                ::FunctionSpaceType
		  VelocityFunctionSpaceType;

        typedef Velocity< VelocityTraits < gridDim, VelocityFunctionSpaceType > >
            VelocityType;

        typedef typename DiscreteFunctionSpaceWrapperType
                ::DiscretePressureFunctionSpaceType
                ::FunctionSpaceType
            PressureFunctionSpaceType;

        typedef Pressure< PressureTraits< gridDim, PressureFunctionSpaceType > >
            PressureType;
        typedef Force< VelocityFunctionSpaceType >
            ForceType;
        typedef DirichletData< VelocityFunctionSpaceType >
            DirichletDataType;

    /**
     *  \brief  constructor
     *
     *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
     **/
    Problem( const double viscosity, const DiscreteFunctionWrapperType& funcWrapper )
        : velocity_( funcWrapper.discreteVelocity().space() ),
          pressure_ ( funcWrapper.discretePressure().space() ),
          force_( viscosity, funcWrapper.discreteVelocity().space() ),
          dirichletData_( funcWrapper.discreteVelocity().space() )

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
    const PressureType& pressure() const
    {
        return pressure_;
    }

    /**
     *  \brief  to get the force term
     *
     *  \return force
     **/
    const ForceType& force() const
    {
        return force_;
    }

    /**
     *  \brief  to get the dirichlet boundary data
     *
     *  \return dirichlet boundary data
     **/
    const DirichletDataType& dirichletData() const
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
//        force_.testMe();
//        dirichletData_.testMe();
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
