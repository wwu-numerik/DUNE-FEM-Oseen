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
