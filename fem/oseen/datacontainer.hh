#ifndef DUNE_OSEEN_DATACONTAINER_HH
#define DUNE_OSEEN_DATACONTAINER_HH

namespace Dune { namespace Oseen {

//! when requested we store \f$ \varDelta u, \nabla p (u \cdot \nabla ) u\f$ in this struct after the solver
template < class Traits >
struct RhsDatacontainer {
    typename Traits::DiscreteVelocityFunctionType velocity_laplace;
    typename Traits::DiscreteVelocityFunctionType pressure_gradient;
    typename Traits::DiscreteSigmaFunctionType velocity_gradient;
    typename Traits::DiscreteVelocityFunctionType convection;

    RhsDatacontainer( const typename Traits::DiscreteVelocityFunctionSpaceType& space,
                      const typename Traits::DiscreteSigmaFunctionSpaceType& sigma_space)
        : velocity_laplace( "velocity_laplace", space ),
        pressure_gradient( "pressure_gradient", space ),
        velocity_gradient( "velocity_gradient", sigma_space ),
        convection( "convection", space )
    {}
    void scale( double factor ) {
        velocity_laplace	*= factor;
        pressure_gradient	*= factor;
        velocity_gradient	*= factor;
        convection			*= factor;
    }
    void clear() {
        velocity_laplace.clear();
        pressure_gradient.clear();
        velocity_gradient.clear();
        convection.clear();
    }
};

} }  // namespace Dune { namespace Oseen {

#endif // DUNE_OSEEN_DATACONTAINER_HH

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

