#ifndef DUNE_STOKES_INTEGRATORS_FACTORY_HH
#define DUNE_STOKES_INTEGRATORS_FACTORY_HH

#include <dune/common/shared_ptr.hh>

#ifdef STOKES_USE_ISTL
#   define STOKES_MATRIX_OBJECT ISTLMatrixObject
#else
#   define STOKES_MATRIX_OBJECT SparseRowMatrixObject
#endif

#ifdef STOKES_USE_ISTL
#   include <dune/fem/operator/2order/dgmatrixtraits.hh>
#   define STOKES_MATRIX_OBJECT_TRAITS Dune::Fem::DGMatrixTraits
#else
    template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
    struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
        struct StencilType {
            template < typename T >
            static int nonZerosEstimate( const T& rangeSpace ) {
                return rangeSpace.maxNumLocalDofs() * 1.5f;
            }
        };
    };
#define STOKES_MATRIX_OBJECT_TRAITS MatrixTraits
#endif

#if 0
template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct DiagonalMatrixTraits :
#ifdef STOKES_USE_ISTL
        public Dune::ISTLMatrixTraits<RowSpaceImp,ColSpaceImp> {
#else
        public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
#endif
    struct StencilType {
        static int nonZerosEstimate( const ColSpaceImp& ) {
            return 1;
        }
    };
};
#endif

namespace Dune { namespace Stokes { namespace Integrators {

template < class FactoryType, class MatrixType >
struct IntegratorSelector {};
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::MInversMatrixIntegratorType> { typedef typename FactoryType::MInversMatrixIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::WmatrixType> { typedef typename FactoryType::WmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::XmatrixType> { typedef typename FactoryType::XmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::YmatrixType> { typedef typename FactoryType::YmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::OmatrixType> { typedef typename FactoryType::OmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::ZmatrixType> { typedef typename FactoryType::ZmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::EmatrixType> { typedef typename FactoryType::EmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::RmatrixType> { typedef typename FactoryType::RmatrixTypeIntegratorType Type; };


template < class StokesTraitsType >
class Factory {
    typedef Factory< StokesTraitsType >
        ThisType;
    //necessary so it can extract the MatrixTypes
    template < class G, class T > friend class IntegratorSelector;
public:
    typedef typename StokesTraitsType::DiscreteSigmaFunctionSpaceType
        DiscreteSigmaFunctionSpaceType;
    typedef typename StokesTraitsType::DiscreteVelocityFunctionSpaceType
        DiscreteVelocityFunctionSpaceType;
    typedef typename StokesTraitsType::DiscreteVelocityFunctionType
        DiscreteVelocityFunctionType;
    typedef typename StokesTraitsType::DiscretePressureFunctionSpaceType
        DiscretePressureFunctionSpaceType;

    typedef STOKES_MATRIX_OBJECT<  DiscreteSigmaFunctionSpaceType,
                                    DiscreteSigmaFunctionSpaceType,
                                    STOKES_MATRIX_OBJECT_TRAITS<DiscreteSigmaFunctionSpaceType,DiscreteSigmaFunctionSpaceType> >
        MInversMatrixInternalType;
    typedef Dune::shared_ptr< MInversMatrixInternalType >
        MInversMatrixType;

    typedef STOKES_MATRIX_OBJECT<  DiscreteVelocityFunctionSpaceType,
                                        DiscreteSigmaFunctionSpaceType,
                                        STOKES_MATRIX_OBJECT_TRAITS<DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType> >
            WmatrixInternalType;
    typedef Dune::shared_ptr< WmatrixInternalType >
        WmatrixType;
    typedef STOKES_MATRIX_OBJECT<  DiscreteSigmaFunctionSpaceType,
                                        DiscreteVelocityFunctionSpaceType,
                                        STOKES_MATRIX_OBJECT_TRAITS<DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType> >
            XmatrixInternalType;
    typedef Dune::shared_ptr< XmatrixInternalType >
        XmatrixType;
    typedef STOKES_MATRIX_OBJECT<  DiscreteVelocityFunctionSpaceType,
                                    DiscreteVelocityFunctionSpaceType,
                                    STOKES_MATRIX_OBJECT_TRAITS<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
        YmatrixInternalType;
    typedef Dune::shared_ptr< YmatrixInternalType >
        YmatrixType;
    typedef STOKES_MATRIX_OBJECT<  DiscretePressureFunctionSpaceType,
                                        DiscreteVelocityFunctionSpaceType,
                                        STOKES_MATRIX_OBJECT_TRAITS<DiscretePressureFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
            ZmatrixInternalType;
    typedef Dune::shared_ptr< ZmatrixInternalType >
        ZmatrixType;
    typedef STOKES_MATRIX_OBJECT<  DiscreteVelocityFunctionSpaceType,
                                        DiscretePressureFunctionSpaceType,
                                        STOKES_MATRIX_OBJECT_TRAITS<DiscreteVelocityFunctionSpaceType,DiscretePressureFunctionSpaceType> >
            EmatrixInternalType;
    typedef Dune::shared_ptr< EmatrixInternalType >
        EmatrixType;
    typedef STOKES_MATRIX_OBJECT<  DiscretePressureFunctionSpaceType,
                                    DiscretePressureFunctionSpaceType,
                                    STOKES_MATRIX_OBJECT_TRAITS<DiscretePressureFunctionSpaceType,DiscretePressureFunctionSpaceType> >
        RmatrixInternalType;
    typedef Dune::shared_ptr< RmatrixInternalType >
        RmatrixType;

    static const bool verbose_ = true;

public:
    typedef Stokes::Integrators::M< MInversMatrixType, StokesTraitsType >
        MInversMatrixIntegratorType;
    typedef Stokes::Integrators::W< WmatrixType, StokesTraitsType >
        WmatrixTypeIntegratorType;
    typedef Stokes::Integrators::X< XmatrixType, StokesTraitsType >
        XmatrixTypeIntegratorType;
    typedef Stokes::Integrators::Y< YmatrixType, StokesTraitsType >
        YmatrixTypeIntegratorType;
    typedef Stokes::Integrators::O< YmatrixType, StokesTraitsType, DiscreteVelocityFunctionType >
        OmatrixTypeIntegratorType;
    typedef Stokes::Integrators::Z< ZmatrixType, StokesTraitsType >
        ZmatrixTypeIntegratorType;
    typedef Stokes::Integrators::E< EmatrixType, StokesTraitsType >
        EmatrixTypeIntegratorType;
    typedef Stokes::Integrators::R< RmatrixType, StokesTraitsType >
        RmatrixTypeIntegratorType;

    static MInversMatrixType get( const DiscreteSigmaFunctionSpaceType& a, const DiscreteSigmaFunctionSpaceType& b )
    { return make( a, b); }

    template < class F, class G >
    static Dune::shared_ptr< STOKES_MATRIX_OBJECT<  F, G,
                                    STOKES_MATRIX_OBJECT_TRAITS<F,G> >
                           > matrix( const F& f, const G& g )
    {
        typedef STOKES_MATRIX_OBJECT<  F, G,
                    STOKES_MATRIX_OBJECT_TRAITS<F,G> >
            InternalMatrixType;
        Dune::shared_ptr< InternalMatrixType > m( new InternalMatrixType(f,g) );
        m->reserve( verbose_ );
        return m;
    }

    template < class F >
    static typename IntegratorSelector< ThisType, F >::Type integrator( F& f )
    {
        return typename IntegratorSelector< ThisType, F >::Type( f );
    }
    //more magic so I don't need two func names please
    static OmatrixTypeIntegratorType integratorO( YmatrixType& g,
                                                          const DiscreteVelocityFunctionType& f )
    {
        return OmatrixTypeIntegratorType( g, f );
    }
};

#undef STOKES_MATRIX_OBJECT
#undef STOKES_MATRIX_OBJECT_TRAITS
/*namespace Dune*/ } /*namespace Stokes*/ } /*namespace Integrators*/ }
#endif // DUNE_STOKES_INTEGRATORS_FACTORY_HH
