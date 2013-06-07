/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <cmake_config.h>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/oseen/assembler/ported_matrixobject.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/fem/oseen/defaulttraits.hh>
#include <dune/fem/oseen/datacontainer.hh>
#include <dune/fem/oseen/solver/solvercaller.hh>
#include <dune/fem/oseen/assembler/all.hh>
#include <dune/fem/oseen/runinfo.hh>

#include <dune/stuff/fem/customprojection.hh>
#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/common/tuple.hh>
#ifndef NLOG
#   include <dune/stuff/common/print.hh>
#   include <dune/stuff/common/misc.hh>
#   include <dune/stuff/common/logging.hh>
#endif
#include <dune/stuff/grid/entity.hh>
#include <dune/stuff/common/profiler.hh>

namespace Dune {

/**
 *  \brief  OseenPass
 *  \todo   doc
 **/
template <  class DiscreteModelImp,
		  template <class> class TraitsImp = StokesTraits >
class OseenPass
{
    public:
        typedef OseenPass< DiscreteModelImp, TraitsImp >
            ThisType;
        typedef DiscreteModelImp
            DiscreteModelType;
        typedef TraitsImp< DiscreteModelType >
            Traits;
        typedef typename Traits::DiscreteOseenFunctionWrapperType
            DomainType;
        typedef typename Traits::DiscreteOseenFunctionWrapperType
            RangeType;

        //!
        OseenPass( DiscreteModelType discreteModel,
					typename Traits::GridPartType& gridPart,
                    const typename Traits::DiscreteOseenFunctionSpaceWrapperType& spaceWrapper,
                    const typename Traits::DiscreteVelocityFunctionType beta,
					const bool do_oseen_discretization )//! \todo move to model
            : discreteModel_( discreteModel ),
            gridPart_( gridPart ),
            spaceWrapper_( spaceWrapper ),
            velocitySpace_( spaceWrapper.discreteVelocitySpace() ),
            pressureSpace_( spaceWrapper.discretePressureSpace() ),
			sigmaSpace_( gridPart ),
			beta_( beta ),
			do_oseen_discretization_( do_oseen_discretization ),
	      info_( SaddlepointInverseOperatorInfo() )
        {}

        //! used in Postprocessing to get refs to gridparts, spaces
        const typename Traits::DiscreteOseenFunctionSpaceWrapperType& GetFunctionSpaceWrapper()
        {
            return spaceWrapper_;
        }

        /**
         *  \todo doc
         **/
        template < class OtherTraitsImp = Traits  >
        void apply( const DomainType &arg, RangeType &dest, Dune::Oseen::RhsDatacontainer<OtherTraitsImp>* rhs_datacontainer = nullptr )
        {
            // profiler information
            DSC_PROFILER.startTiming("Pass_init");
            typedef Oseen::Assembler::Factory< Traits >
                Factory;
            // M\in R^{M\times M}
            auto MInversMatrix = Factory::matrix( sigmaSpace_, sigmaSpace_ );
            ASSERT_EQ( MInversMatrix->matrix().rows(), MInversMatrix->matrix().cols() );
            // W\in R^{M\times L}
            auto Wmatrix = Factory::matrix( sigmaSpace_, velocitySpace_ );
            // X\in R^{L\times M}
            auto Xmatrix = Factory::matrix( velocitySpace_, sigmaSpace_ );
            // O,Y\in R^{L\times L}
            auto Ymatrix = Factory::matrix( velocitySpace_, velocitySpace_ );
            auto Omatrix = Factory::matrix( velocitySpace_, velocitySpace_ );
            // Z\in R^{L\times K}
            auto Zmatrix = Factory::matrix( velocitySpace_, pressureSpace_ );
            // E\in R^{K\times L}
            auto Ematrix = Factory::matrix( pressureSpace_, velocitySpace_ );
            // R\in R^{K\times K}
            auto Rmatrix = Factory::matrix( pressureSpace_, pressureSpace_ );
            // H_{1}\in R^{M}
            auto H1rhs = Factory::rhs( "H1", sigmaSpace_ );
            // H_{2}\in R^{L}
            auto H2rhs = Factory::rhs( "H2", velocitySpace_ );
            auto H2_O_rhs = Factory::rhs( "H2_O", velocitySpace_ );
            // H_{3}\in R^{K}
            auto H3rhs = Factory::rhs( "H3", pressureSpace_ );
            auto m_integrator = Factory::integrator( MInversMatrix );
            auto w_integrator = Factory::integrator( Wmatrix );
            auto x_integrator = Factory::integrator( Xmatrix );
            auto y_integrator = Factory::integrator( Ymatrix );
            auto o_integrator = Factory::integratorO( Omatrix, beta_ );
            auto z_integrator = Factory::integrator( Zmatrix );
            auto e_integrator = Factory::integrator( Ematrix );
            auto r_integrator = Factory::integrator( Rmatrix );
            auto h1_integrator = Factory::integrator( H1rhs );
            auto h2_integrator = Factory::integrator( H2rhs );
            auto h2_o_integrator = Factory::integratorO( H2_O_rhs, beta_ );
            H2rhs += H2_O_rhs;
            auto h3_integrator = Factory::integrator( H3rhs );
            DSC_PROFILER.stopTiming("Pass_init");

#ifndef STOKES_CONV_ONLY
            if ( do_oseen_discretization_ )
            {
                Oseen::Assembler::Coordinator< Traits, typename Factory::OseenIntegratorTuple >
                        coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

                typename Factory::OseenIntegratorTuple tuple(	m_integrator, w_integrator, x_integrator, y_integrator,
                                        o_integrator, z_integrator, e_integrator, r_integrator,
                                        h1_integrator, h2_integrator,h2_o_integrator, h3_integrator );
                coordinator.apply( tuple );
            }
            else
            {
                Oseen::Assembler::Coordinator< Traits, typename Factory::StokesIntegratorTuple >
                        coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

                typename Factory::StokesIntegratorTuple tuple(	m_integrator, w_integrator, x_integrator, y_integrator,
                                        z_integrator, e_integrator, r_integrator,
                                        h1_integrator, h2_integrator,h3_integrator );
                coordinator.apply( tuple );
            }
#else
            Oseen::Assembler::Coordinator< Traits, typename Factory::ConvIntegratorTuple >
                    coordinator ( discreteModel_, gridPart_, velocitySpace_, pressureSpace_, sigmaSpace_  );

            typename Factory::ConvIntegratorTuple tuple(	o_integrator, h2_integrator, h2_o_integrator );
            coordinator.apply( tuple );
#endif
            // do the actual lgs solving
            DSC_LOG_INFO << "Solving system with " << dest.discreteVelocity().size() << " + " << dest.discretePressure().size() << " unknowns" << std::endl;
            info_ = Oseen::SolverCallerProxy< ThisType >::call( do_oseen_discretization_, rhs_datacontainer, dest,
                                            arg, Xmatrix, MInversMatrix, Ymatrix, Omatrix, Ematrix,
                                            Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs, beta_ );
        } // end of apply

		void getRuninfo( DSC::RunInfo& info )
        {
			info.iterations_inner_avg = int( info_.iterations_inner_avg );
            info.iterations_inner_min = info_.iterations_inner_min;
            info.iterations_inner_max = info_.iterations_inner_max;
            info.iterations_outer_total = info_.iterations_outer_total;
            info.max_inner_accuracy = info_.max_inner_accuracy;
        }

    private:
        DiscreteModelType discreteModel_;
		const typename Traits::GridPartType& gridPart_;
        const typename Traits::DiscreteOseenFunctionSpaceWrapperType& spaceWrapper_;
		const typename Traits::DiscreteVelocityFunctionSpaceType& velocitySpace_;
		const typename Traits::DiscretePressureFunctionSpaceType& pressureSpace_;
		typename Traits::DiscreteSigmaFunctionSpaceType sigmaSpace_;
        const typename Traits::DiscreteVelocityFunctionType beta_;
		const bool do_oseen_discretization_;
        SaddlepointInverseOperatorInfo info_;

	public:
		void printInfo() const
		{
#ifndef NLOG
            auto& infoStream = DSC_LOG_INFO;
			infoStream << boost::format( "pressure_gradient/convection scaling: %e | %e\npass viscosity: %e\n")
								% discreteModel_.pressure_gradient_scaling()
								% discreteModel_.convection_scaling()
								% discreteModel_.viscosity();
			int numberOfEntities = 0;
			int numberOfIntersections = 0;
			int numberOfBoundaryIntersections = 0;
			int numberOfInnerIntersections = 0;
            infoStream << "this is OseenPass::apply()" << std::endl;

			// do an empty grid walk to get informations
			double maxGridWidth( 0.0 );
			typename Traits::EntityIteratorType entityItEndLog = velocitySpace_.end();
			for (   typename Traits::EntityIteratorType entityItLog = velocitySpace_.begin();
					entityItLog != entityItEndLog;
					++entityItLog ) {
				const typename Traits::EntityType& entity = *entityItLog;
				// count entities
				++numberOfEntities;
				// walk the intersections
				typename Traits::IntersectionIteratorType intItEnd = gridPart_.iend( entity );
				for (   typename Traits::IntersectionIteratorType intIt = gridPart_.ibegin( entity );
						intIt != intItEnd;
						++intIt ) {
					// count intersections
					++numberOfIntersections;
                    maxGridWidth = std::max(intIt->geometry().volume(), maxGridWidth );
					// if we are inside the grid
					if ( intIt->neighbor() && !intIt->boundary() ) {
						// count inner intersections
						++numberOfInnerIntersections;
					}
					// if we are on the boundary of the grid
					if ( !intIt->neighbor() && intIt->boundary() ) {
						// count boundary intersections
						++numberOfBoundaryIntersections;
					}
				}
			}
			if ( numberOfEntities > 19 ) {
				infoStream << "found " << numberOfEntities << " entities," << std::endl;
				infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
				infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
				infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
				infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
				infoStream << "- starting gridwalk" << std::endl;
			} else {
				infoStream << "found " << numberOfEntities << " entities," << std::endl;
				infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
				infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
				infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
				infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
				infoStream << "- starting gridwalk" << std::endl;
			}
			infoStream.suspend();
#endif
		}
};

} // namespace Dune
#endif  // end of stokespass.hh

/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
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

