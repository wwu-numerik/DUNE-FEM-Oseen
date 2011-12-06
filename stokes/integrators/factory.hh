#ifndef DUNE_STOKES_INTEGRATORS_FACTORY_HH
#define DUNE_STOKES_INTEGRATORS_FACTORY_HH

#include <dune/common/shared_ptr.hh>

#ifdef STOKES_USE_ISTL
#   include <dune/stokes/integrators/mod_istlmatrix.hh>
#   define STOKES_MATRIX_OBJECT ModifiedISTLMatrixObject
#   include <dune/fem/operator/2order/dgmatrixtraits.hh>
#   include <dune/stokes/integrators/bcrstraits.hh>
#   define STOKES_MATRIX_OBJECT_TRAITS Dune::Stokes::Integrators::ModifiedDGMatrixTraits
#else
#   define STOKES_MATRIX_OBJECT SparseRowMatrixObject
    template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
    struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
        struct StencilType {
            template < typename T >
            static int nonZerosEstimate( const T& rangeSpace ) {
                return rangeSpace.maxNumLocalDofs() * 1.5f;
            }
        };
    };
#   define STOKES_MATRIX_OBJECT_TRAITS MatrixTraits
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

//! A Static map of Matrix-/DiscreteFunctionType onto IntegratorType
template < class FactoryType, class MatrixType >
struct IntegratorSelector {};
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::MInversMatrixType>
{ typedef typename FactoryType::MInversMatrixIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::WmatrixType>
{ typedef typename FactoryType::WmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::XmatrixType>
{ typedef typename FactoryType::XmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::YmatrixType>
{ typedef typename FactoryType::YmatrixTypeIntegratorType Type; };
//template < class FactoryType > //O mapping doesn't work because of the extra  arg
//struct IntegratorSelector< FactoryType, typename FactoryType::OmatrixType>
//{ typedef typename FactoryType::OmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::ZmatrixType>
{ typedef typename FactoryType::ZmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::EmatrixType>
{ typedef typename FactoryType::EmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::RmatrixType>
{ typedef typename FactoryType::RmatrixTypeIntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::DiscreteSigmaFunctionType>
{ typedef typename FactoryType::H1_IntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::DiscreteVelocityFunctionType>
{ typedef typename FactoryType::H2_IntegratorType Type; };
//template < class FactoryType >//H2_O does not work because of extra arg
//struct IntegratorSelector< FactoryType, typename FactoryType::DiscreteVelocityFunctionType>
//{ typedef typename FactoryType::H2_IntegratorType Type; };
template < class FactoryType >
struct IntegratorSelector< FactoryType, typename FactoryType::DiscretePressureFunctionType>
{ typedef typename FactoryType::H3_IntegratorType Type; };

//! A static map of DiscreteFunctionSpace onto DiscreteFunction
template < class FactoryType, class DiscreteFunctionSpaceType >
struct DiscreteFunctionSelector {};
template < class FactoryType >
struct DiscreteFunctionSelector< FactoryType, typename FactoryType::DiscreteSigmaFunctionSpaceType >
{ typedef typename FactoryType::DiscreteSigmaFunctionType Type; };
template < class FactoryType >
struct DiscreteFunctionSelector< FactoryType, typename FactoryType::DiscreteVelocityFunctionSpaceType >
{ typedef typename FactoryType::DiscreteVelocityFunctionType Type; };
template < class FactoryType >
struct DiscreteFunctionSelector< FactoryType, typename FactoryType::DiscretePressureFunctionSpaceType >
{ typedef typename FactoryType::DiscretePressureFunctionType Type; };

//!
template < class StokesTraitsType >
class Factory {
    typedef Factory< StokesTraitsType >
        ThisType;
public:
    typedef typename StokesTraitsType::DiscreteSigmaFunctionSpaceType
        DiscreteSigmaFunctionSpaceType;
    typedef typename StokesTraitsType::DiscreteSigmaFunctionType
        DiscreteSigmaFunctionType;
    typedef typename StokesTraitsType::DiscreteVelocityFunctionSpaceType
        DiscreteVelocityFunctionSpaceType;
    typedef typename StokesTraitsType::DiscreteVelocityFunctionType
        DiscreteVelocityFunctionType;
    typedef typename StokesTraitsType::DiscretePressureFunctionSpaceType
        DiscretePressureFunctionSpaceType;
    typedef typename StokesTraitsType::DiscretePressureFunctionType
        DiscretePressureFunctionType;

    template < class T, class R >
    struct MatrixObject {
        typedef STOKES_MATRIX_OBJECT_TRAITS<T,R> Traits;
        typedef STOKES_MATRIX_OBJECT<  T, R, Traits >Type;
    };

    typedef typename MatrixObject< DiscreteSigmaFunctionSpaceType, DiscreteSigmaFunctionSpaceType >::Type
        MInversMatrixInternalType;
    typedef Dune::shared_ptr< MInversMatrixInternalType >
        MInversMatrixType;
    typedef typename MatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType >::Type
        WmatrixInternalType;
    typedef Dune::shared_ptr< WmatrixInternalType >
        WmatrixType;
    typedef typename MatrixObject< DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType >::Type
        XmatrixInternalType;
    typedef Dune::shared_ptr< XmatrixInternalType >
        XmatrixType;
    typedef typename MatrixObject< DiscreteVelocityFunctionSpaceType, DiscreteVelocityFunctionSpaceType >::Type
        YmatrixInternalType;
    typedef Dune::shared_ptr< YmatrixInternalType >
        YmatrixType;
    typedef typename MatrixObject< DiscretePressureFunctionSpaceType, DiscreteVelocityFunctionSpaceType >::Type
        ZmatrixInternalType;
    typedef Dune::shared_ptr< ZmatrixInternalType >
        ZmatrixType;
    typedef typename MatrixObject< DiscreteVelocityFunctionSpaceType, DiscretePressureFunctionSpaceType >::Type
        EmatrixInternalType;
    typedef Dune::shared_ptr< EmatrixInternalType >
        EmatrixType;
    typedef typename MatrixObject< DiscretePressureFunctionSpaceType, DiscretePressureFunctionSpaceType >::Type
        RmatrixInternalType;
    typedef Dune::shared_ptr< RmatrixInternalType >
        RmatrixType;

    static const bool verbose_ = true;

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
    typedef Stokes::Integrators::H1< DiscreteSigmaFunctionType, StokesTraitsType >
        H1_IntegratorType;
    typedef Stokes::Integrators::H2< DiscreteVelocityFunctionType , StokesTraitsType >
        H2_IntegratorType;
    typedef Stokes::Integrators::H2_O< DiscreteVelocityFunctionType , StokesTraitsType, DiscreteVelocityFunctionType >
        H2_O_IntegratorType;
    typedef Stokes::Integrators::H3< DiscretePressureFunctionType, StokesTraitsType >
        H3_IntegratorType;
    typedef tuple<	MInversMatrixIntegratorType,
                    WmatrixTypeIntegratorType,
                    XmatrixTypeIntegratorType,
                    YmatrixTypeIntegratorType,
                    OmatrixTypeIntegratorType,
                    ZmatrixTypeIntegratorType,
                    EmatrixTypeIntegratorType,
                    RmatrixTypeIntegratorType,
                    H1_IntegratorType,
                    H2_IntegratorType,
                    H2_O_IntegratorType,
                    H3_IntegratorType >
        OseenIntegratorTuple;
    typedef tuple<	MInversMatrixIntegratorType,
                    WmatrixTypeIntegratorType,
                    XmatrixTypeIntegratorType,
                    YmatrixTypeIntegratorType,
                    ZmatrixTypeIntegratorType,
                    EmatrixTypeIntegratorType,
                    RmatrixTypeIntegratorType,
                    H1_IntegratorType,
                    H2_IntegratorType,
                    H3_IntegratorType >
        StokesIntegratorTuple;

    template < class F, class G >
    static Dune::shared_ptr< typename MatrixObject< F, G >::Type > matrix( const F& f, const G& g )
    {
        typedef typename MatrixObject< F, G >::Type
            InternalMatrixType;
        Dune::shared_ptr< InternalMatrixType > m( new InternalMatrixType(f,g) );
        m->reserve( verbose_ );
        return m;
    }
    template < class DiscreteFunctionSpaceType >
    static typename DiscreteFunctionSelector< ThisType, DiscreteFunctionSpaceType >::Type
        rhs( const std::string name, const DiscreteFunctionSpaceType& space )
    {
        typename DiscreteFunctionSelector< ThisType, DiscreteFunctionSpaceType >::Type f( name, space );
        f.clear();
        return f;
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
    static H2_O_IntegratorType integratorO( DiscreteVelocityFunctionType& g,
                                                          const DiscreteVelocityFunctionType& f )
    {
        return H2_O_IntegratorType( g, f );
    }
};

#undef STOKES_MATRIX_OBJECT
#undef STOKES_MATRIX_OBJECT_TRAITS
/*namespace Dune*/ } /*namespace Stokes*/ } /*namespace Integrators*/ }
#endif // DUNE_STOKES_INTEGRATORS_FACTORY_HH
