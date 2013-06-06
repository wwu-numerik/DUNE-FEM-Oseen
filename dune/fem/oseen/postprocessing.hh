/**

 *  \file   postprocessing.hh
 *  \brief  postprocessing.hh
 **/

#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <cmake_config.h>

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/customprojection.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/fem/functions/transform.hh>
#include <dune/stuff/fem/functions/integrals.hh>

#include <boost/format.hpp>
#include <cmath>
#include <sstream>
#include <tuple>

//! Error and vtk output wrapper class for Stokes problem/pass
template <  class OseenPassImp, class ProblemImp >
class PostProcessor
{
    public:
        typedef ProblemImp
            ProblemType;

        typedef OseenPassImp
            OseenPassType;
		typedef typename OseenPassType::Traits::DiscreteOseenFunctionSpaceWrapperType
            DiscreteOseenFunctionSpaceWrapperType;
		typedef typename OseenPassType::Traits::DiscreteOseenFunctionWrapperType
            DiscreteOseenFunctionWrapperType;

        typedef typename ProblemType::VelocityType
            ContinuousVelocityType;
        typedef typename ProblemType::PressureType
            ContinuousPressureType;
        typedef typename ProblemType::ForceType
            ForceType;
        typedef typename ProblemType::DirichletDataType
            DirichletDataType;

		typedef typename OseenPassType::Traits::GridPartType
            GridPartType;
        typedef typename GridPartType::GridType
            GridType;

        typedef Dune::SubsamplingVTKIO<GridPartType>
            VTKWriterType;

		typedef typename OseenPassType::Traits::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;
		typedef typename OseenPassType::Traits::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

		typedef typename OseenPassType::Traits::DiscretePressureFunctionType
            DiscretePressureFunctionType;
		typedef typename OseenPassType::Traits::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;


        PostProcessor( const DiscreteOseenFunctionSpaceWrapperType& wrapper, const ProblemType& prob )
            :
            problem_( prob ),
            spaceWrapper_( wrapper ),
            gridPart_( wrapper.gridPart() ),
            velocitySpace_ ( wrapper.discreteVelocitySpace() ),
            discreteExactVelocity_( "u_exact", wrapper.discreteVelocitySpace() ),
            discreteExactForce_( "f_exact", wrapper.discreteVelocitySpace() ),
            discreteExactDirichlet_( "gd_exact", wrapper.discreteVelocitySpace() ),
            discreteExactPressure_( "p_exact", wrapper.discretePressureSpace() ),
            errorFunc_velocity_( "err_velocity", wrapper.discreteVelocitySpace() ),
            errorFunc_pressure_( "err_pressure", wrapper.discretePressureSpace() ),
            solutionAssembled_( false ),
            current_refine_level_( std::numeric_limits<int>::min() ),
            l2_error_pressure_( - std::numeric_limits<double>::max() ),
            l2_error_velocity_( - std::numeric_limits<double>::max() ),
            vtkWriter_( wrapper.gridPart() ),
            datadir_( DSC_CONFIG_GET( "fem.io.datadir", std::string("data") ) + "/" )
        {
            DSC::testCreateDirectory( datadir_ );
        }

		/** \brief analytical data is L2 projected
			\todo only use DSC::CustomProjection when really necessary
			**/
        void assembleExactSolution()
        {
            DSFe::CustomProjection::project( problem_.dirichletData(), discreteExactDirichlet_ );

            typedef Dune::L2Projection< double, double, ContinuousVelocityType, DiscreteVelocityFunctionType > ProjectionV;
                ProjectionV projectionV;
            projectionV( problem_.velocity(), discreteExactVelocity_ );

            typedef Dune::L2Projection< double, double, ForceType, DiscreteVelocityFunctionType > ProjectionF;
                ProjectionF projectionF;
            projectionF( problem_.force(), discreteExactForce_ );

            typedef Dune::L2Projection< double, double, ContinuousPressureType, DiscretePressureFunctionType > ProjectionP;
                ProjectionP projectionP;
            projectionP( problem_.pressure(), discreteExactPressure_ );
			if ( DSC_CONFIG_GET( "save_matrices", false ) ) {
                auto& matlabLogStream = DSC_LOG_ERROR;
				DSC::printDiscreteFunctionMatlabStyle( discreteExactVelocity_, "u_exakt", matlabLogStream );
				DSC::printDiscreteFunctionMatlabStyle( discreteExactPressure_, "p_exakt", matlabLogStream );
			}
        }

		//! output function that 'knows' function output mode; assembles filename
        template <class Function>
        void vtk_write( const Function& f ) {
            if ( Function::FunctionSpaceType::DimRange > 1 ) {
                vtkWriter_.addVectorVertexData( f );
                vtkWriter_.addVectorCellData( f );
            }
            else {
                vtkWriter_.addVertexData( f );
                vtkWriter_.addCellData( f );
            }

            std::stringstream path;
            if ( DSC_CONFIG_GET( "per-run-output", false ) )
                path    << datadir_ << "/ref"
                        <<  current_refine_level_ << "_" << f.name();
            else
                path << datadir_ << "/" << f.name();

            vtkWriter_.write( path.str().c_str() );
            vtkWriter_.clear();
        }

		//! use this function if no reference (ie. coarser/finer) solution is available, or an analytical one is
        void save( const GridType& grid, const DiscreteOseenFunctionWrapperType& wrapper, int refine_level )
        {
            if ( ProblemType:: hasMeaningfulAnalyticalSolution ) {
                if ( !solutionAssembled_ || current_refine_level_ != refine_level ) //re-assemble solution if refine level has changed
                    assembleExactSolution();
                current_refine_level_ = refine_level;

				calcError( wrapper );
                vtk_write( discreteExactVelocity_ );
                vtk_write( discreteExactPressure_ );
                vtk_write( discreteExactForce_ );
                vtk_write( discreteExactDirichlet_ );
                vtk_write( errorFunc_pressure_ );
                vtk_write( errorFunc_velocity_ );
            }

            save_common( grid, wrapper, refine_level );
        }

		//! use this save in eoc runs with no analytical solution available
        void save( const GridType& grid, const DiscreteOseenFunctionWrapperType& wrapper, const DiscreteOseenFunctionWrapperType& reference, int refine_level )
        {
            current_refine_level_ = refine_level;
            calcError( wrapper, reference );
            vtk_write( discreteExactVelocity_ );
            vtk_write( discreteExactPressure_ );
            vtk_write( errorFunc_pressure_ );
            vtk_write( errorFunc_velocity_ );

            save_common( grid, wrapper, refine_level );
        }

		//! used by both PostProcessor::save modes, outputs solutions (in grape/vtk form), but no errors or analytical functions
        void save_common( const GridType& grid, const DiscreteOseenFunctionWrapperType& wrapper, int refine_level )
        {
            current_refine_level_ = refine_level;

            vtk_write( wrapper.discretePressure() );
            vtk_write( wrapper.discreteVelocity() );
#ifndef NLOG
			entityColoration();
#endif
        }

		void calcError( const DiscreteOseenFunctionWrapperType& wrapper )
		{
			calcError( wrapper.discretePressure() , wrapper.discreteVelocity() );
		}

		//! proxy function that is to be used if no analytical solutions are availble to calculate errors against
        void calcError( const DiscreteOseenFunctionWrapperType& computed, const DiscreteOseenFunctionWrapperType& reference )
        {
            discreteExactPressure_.assign( reference.discretePressure() );
            discreteExactVelocity_.assign( reference.discreteVelocity() );
			//set to to true so calcError call does not try to assemble exact solutions again
            solutionAssembled_ = true;
            calcError( computed.discretePressure(), computed.discreteVelocity() );
        }

		//! print and save L2 error(functions)
        void calcError( const DiscretePressureFunctionType& pressure, const DiscreteVelocityFunctionType& velocity )
        {
            if ( !solutionAssembled_ )
                assembleExactSolution();

            errorFunc_pressure_.assign( discreteExactPressure_ );
            errorFunc_pressure_ -= pressure;
            errorFunc_velocity_.assign( discreteExactVelocity_ );
            errorFunc_velocity_ -= velocity;

            Dune::L2Norm< GridPartType > l2_Error( gridPart_ );

            l2_error_pressure_ = l2_Error.norm( errorFunc_pressure_ );
            l2_error_velocity_ = l2_Error.norm( errorFunc_velocity_ );

            const double boundaryInt = DSFe::boundaryIntegral( problem_.dirichletData(), discreteExactVelocity_.space() );
            const double pressureMean = DSFe::integralAndVolume( pressure, pressure.space() ).first;
            const double exactPressureMean = DSFe::integralAndVolume( problem_.pressure(), discreteExactPressure_.space() ).first;

            DSC_LOG_INFO.resume();
            DSC_LOG_INFO << "L2-Error Pressure: " << std::setw(8) << l2_error_pressure_ << "\n"
                            << "L2-Error Velocity: " << std::setw(8) << l2_error_velocity_ << "\n"
							<< boost::format( "Pressure volume integral: %f (discrete), %f (exact)\n") % pressureMean % exactPressureMean
							<< boost::format( "g_D boundary integral: %f\n") % boundaryInt;
        }

		//! used to sore errors in runinfo structure (for eoc latex output)
        std::vector<double> getError()
        {
            std::vector<double> ret;
            ret.push_back( l2_error_velocity_ );
            ret.push_back( l2_error_pressure_ );
            return ret;
        }

		//! assign each entity it's 'id' int and save/(vtk)output it in a discrete function
        void entityColoration()
        {
            DiscretePressureFunctionType cl ( "entitiy-num", spaceWrapper_.discretePressureSpace() );
            unsigned long numberOfEntities = 0;

            typedef typename GridPartType::GridType::template Codim< 0 >::Entity
                EntityType;
            typedef typename GridPartType::template Codim< 0 >::IteratorType
                EntityIteratorType;
            typedef typename GridPartType::IntersectionIteratorType
                IntersectionIteratorType;

            EntityIteratorType entityItEndLog = velocitySpace_.end();
            for (   EntityIteratorType entityItLog = velocitySpace_.begin();
                    entityItLog != entityItEndLog;
                    ++entityItLog, ++numberOfEntities ) {
                const EntityType& entity = *entityItLog;
                typename DiscretePressureFunctionType::LocalFunctionType
                    lf = cl.localFunction( entity );

                for ( int i = 0; i < lf.numDofs(); ++i ){
                    lf[i] = numberOfEntities;
                }
            }
            vtk_write( cl );
        }

    private:

        const ProblemType& problem_;
        const DiscreteOseenFunctionSpaceWrapperType& spaceWrapper_;
        const GridPartType& gridPart_;
        const DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscreteVelocityFunctionType discreteExactVelocity_;
        DiscreteVelocityFunctionType discreteExactForce_;
        DiscreteVelocityFunctionType discreteExactDirichlet_;
        DiscretePressureFunctionType discreteExactPressure_;
        DiscreteVelocityFunctionType errorFunc_velocity_;
        DiscretePressureFunctionType errorFunc_pressure_;
        bool solutionAssembled_;
        int current_refine_level_;
        double l2_error_pressure_;
        double l2_error_velocity_;
        VTKWriterType vtkWriter_;
        std::string datadir_;
};

#undef vtk_write

#endif // end of postprocessing.hh

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


