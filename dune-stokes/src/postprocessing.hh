/**

 *  \file   postprocessing.hh
 *  \brief  postprocessing.hh
 **/

#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/io/file/vtkio.hh>

#include "logging.hh"
#include "problem.hh"
#include "parametercontainer.hh"

#include <dune/fem/misc/l2norm.hh>


#define VTK_WRITE(z)    vtkWriter_.addVertexData(z); \
                        vtkWriter_.write(( "data/z" ) ); \
                        vtkWriter_.clear();


template <  class ProblemImp,
            class GridPartImp,
            class DiscreteVelocityFunctionImp,
            class DiscretePressureFunctionImp >
class PostProcessor
{
    public:
        typedef ProblemImp
            ProblemType;
        typedef typename ProblemType::VelocityType
            ContinuousVelocityType;
        typedef typename ProblemType::PressureType
            ContinuousPressureType;
        typedef typename ProblemType::ForceType
            ForceType;
        typedef typename ProblemType::DirichletDataType
            DirichletDataType;

        typedef GridPartImp
            GridPartType;
        typedef typename GridPartType::GridType
            GridType;

        typedef Dune::VTKIO<GridPartType>
            VTKWriterType;

        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;



        PostProcessor( const ProblemType& problem, const GridPartType& gridPart,
                        const DiscreteVelocityFunctionSpaceType& velocity_space,
                        const DiscretePressureFunctionSpaceType& press_space)
            : problem_( problem ),
            gridPart_( gridPart ),
            velocitySpace_ ( velocity_space ),
            discreteExactVelocity_( "u_exact", velocity_space ),
            discreteExactForce_( "f_exact", velocity_space ),
            discreteExactDirichlet_( "gd_exact", velocity_space ),
            discreteExactPressure_( "p_exact", press_space ),
            errorFunc_velocity_( "err_velocity", velocity_space ),
            errorFunc_pressure_( "err_pressure", press_space ),
            solutionAssembled_(false),
            l2_error_pressure_( - std::numeric_limits<double>::max() ),
            l2_error_velocity_( - std::numeric_limits<double>::max() ),
            vtkWriter_( gridPart )
        {

        }

        ~PostProcessor()
        {
        }

        void assembleExactSolution()
        {
            typedef Dune::L2Projection< double, double, ContinuousVelocityType, DiscreteVelocityFunctionType > ProjectionV;
                ProjectionV projectionV;
            projectionV( problem_.velocity(), discreteExactVelocity_ );

            typedef Dune::L2Projection< double, double, DirichletDataType, DiscreteVelocityFunctionType > ProjectionD;
                ProjectionD projectionD;
            projectionD( problem_.dirichletData(), discreteExactDirichlet_ );

            typedef Dune::L2Projection< double, double, ForceType, DiscreteVelocityFunctionType > ProjectionF;
                ProjectionF projectionF;
            projectionF( problem_.force(), discreteExactForce_ );

            typedef Dune::L2Projection< double, double, ContinuousPressureType, DiscretePressureFunctionType > ProjectionP;
                ProjectionP projectionP;
            projectionP( problem_.pressure(), discreteExactPressure_ );
        }

        void save( const GridType& grid, const DiscretePressureFunctionType& pressure, const DiscreteVelocityFunctionType& velocity )
        {
            if ( !solutionAssembled_ )
                assembleExactSolution();

            calcError( pressure, velocity );

            VTK_WRITE( discreteExactVelocity_ );
			VTK_WRITE( discreteExactPressure_ );
            VTK_WRITE( discreteExactForce_ );
			VTK_WRITE( discreteExactDirichlet_ );
			VTK_WRITE( errorFunc_pressure_ );
			VTK_WRITE( errorFunc_velocity_ );

        }

        void calcError( const DiscretePressureFunctionType& pressure, const DiscreteVelocityFunctionType& velocity )
        {
            if ( !solutionAssembled_ )
                assembleExactSolution();

            errorFunc_pressure_.assign( discreteExactPressure_ );
            errorFunc_pressure_ -= pressure;
            errorFunc_velocity_.assign( discreteExactVelocity_ );
            errorFunc_velocity_ -= velocity;

            Dune::L2Norm< GridPartType > l2_Error( gridPart_ );
            l2_error_pressure_ =
                l2_Error.norm( errorFunc_pressure_ );
            l2_error_velocity_ =
                l2_Error.norm( errorFunc_velocity_ );

            Logger().Info()  << "L2-Error Pressure: " << std::setw(8) << l2_error_pressure_ << "\n"
                                << "L2-Error Velocity: " << std::setw(8) << l2_error_velocity_ << std::endl;
        }

    private:
        const ProblemType& problem_;
        const GridPartType& gridPart_;
        const DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscreteVelocityFunctionType discreteExactVelocity_;
        DiscreteVelocityFunctionType discreteExactForce_;
        DiscreteVelocityFunctionType discreteExactDirichlet_;
        DiscretePressureFunctionType discreteExactPressure_;
        DiscreteVelocityFunctionType errorFunc_velocity_;
        DiscretePressureFunctionType errorFunc_pressure_;
        bool solutionAssembled_;
        double l2_error_pressure_;
        double l2_error_velocity_;
        VTKWriterType vtkWriter_;

};

#undef VTK_WRITE

#endif // end of postprocessing.hh
