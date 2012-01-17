#ifndef DUNE_OSEEN_INTEGRATORS_FACTORY_HH
#define DUNE_OSEEN_INTEGRATORS_FACTORY_HH

#include <dune/common/shared_ptr.hh>
#include <dune/common/static_assert.hh>

#if STOKES_USE_ISTL
#   include <dune/oseen/assembler/mod_istlmatrix.hh>
#   define STOKES_MATRIX_OBJECT ModifiedISTLMatrixObject
#   include <dune/fem/operator/2order/dgmatrixtraits.hh>
#   include <dune/oseen/assembler/bcrstraits.hh>
#   define STOKES_MATRIX_OBJECT_TRAITS Dune::Stokes::Integrators::ModifiedDGMatrixTraits
#else
#   define STOKES_MATRIX_OBJECT SparseRowMatrixObject
    template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
    struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
        struct StencilType {
            template < typename T >
            static int nonZerosEstimate( const T& rangeSpace ) {
                return rangeSpace.mapper().maxNumDofs() * 1.5f;
            }
        };
    };
#   define STOKES_MATRIX_OBJECT_TRAITS MatrixTraits
#endif

#if 0
template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct DiagonalMatrixTraits :
#if STOKES_USE_ISTL
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

#define MK_FUNC_NAME(name) Discrete ## name ## FunctionType
#define TYPEDEF_MATRIX_AND_INTEGRATOR( Name, Row, Col ) \
    typedef typename MatrixObject< MK_FUNC_NAME(Row), MK_FUNC_NAME(Col) >::Type \
        Name ## matrixInternalType; \
    typedef Dune::shared_ptr< Name ## matrixInternalType > \
        Name ## matrixType ; \
    typedef Stokes::Integrators:: Name < Name ## matrixType, StokesTraitsType > \
        Name ## matrixIntegratorType

#define SPECIALIZE_IntegratorSelector(Name) \
    template < class FactoryType > \
    struct IntegratorSelector< FactoryType, typename FactoryType:: Name ## matrixType> \
    { typedef typename FactoryType:: Name ## matrixIntegratorType Type; }




namespace Dune { namespace Stokes { namespace Integrators {

//! A Static map of Matrix-/DiscreteFunctionType onto IntegratorType
template < class FactoryType, class MatrixType >
struct IntegratorSelector {};
SPECIALIZE_IntegratorSelector(M);
SPECIALIZE_IntegratorSelector(W);
SPECIALIZE_IntegratorSelector(X);
SPECIALIZE_IntegratorSelector(Y);
SPECIALIZE_IntegratorSelector(Z);
SPECIALIZE_IntegratorSelector(E);
SPECIALIZE_IntegratorSelector(R);

//template < class FactoryType > //O mapping doesn't work because of the extra  arg
//struct IntegratorSelector< FactoryType, typename FactoryType::OmatrixType>
//{ typedef typename FactoryType::OmatrixTypeIntegratorType Type; };

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
    TYPEDEF_MATRIX_AND_INTEGRATOR( M, Sigma, Sigma );
    TYPEDEF_MATRIX_AND_INTEGRATOR( W, Sigma, Velocity );
    TYPEDEF_MATRIX_AND_INTEGRATOR( X, Velocity, Sigma );
    TYPEDEF_MATRIX_AND_INTEGRATOR( Y, Velocity, Velocity );
    TYPEDEF_MATRIX_AND_INTEGRATOR( Z, Velocity, Pressure );
    TYPEDEF_MATRIX_AND_INTEGRATOR( E, Pressure, Velocity );
    TYPEDEF_MATRIX_AND_INTEGRATOR( R, Pressure, Pressure );
    static const bool verbose_ = true;

    typedef Stokes::Integrators::O< YmatrixType, StokesTraitsType, DiscreteVelocityFunctionType >
        OmatrixIntegratorType;
    typedef Stokes::Integrators::H1< DiscreteSigmaFunctionType, StokesTraitsType >
        H1_IntegratorType;
    typedef Stokes::Integrators::H2< DiscreteVelocityFunctionType , StokesTraitsType >
        H2_IntegratorType;
    typedef Stokes::Integrators::H2_O< DiscreteVelocityFunctionType , StokesTraitsType, DiscreteVelocityFunctionType >
        H2_O_IntegratorType;
    typedef Stokes::Integrators::H3< DiscretePressureFunctionType, StokesTraitsType >
        H3_IntegratorType;
    typedef tuple<	MmatrixIntegratorType,
                    WmatrixIntegratorType,
                    XmatrixIntegratorType,
                    YmatrixIntegratorType,
                    OmatrixIntegratorType,
                    ZmatrixIntegratorType,
                    EmatrixIntegratorType,
                    RmatrixIntegratorType,
                    H1_IntegratorType,
                    H2_IntegratorType,
                    H2_O_IntegratorType,
                    H3_IntegratorType >
        OseenIntegratorTuple;
    typedef tuple<	OmatrixIntegratorType,
                    H2_IntegratorType,
                    H2_O_IntegratorType>
        ConvIntegratorTuple;
    typedef tuple<	MmatrixIntegratorType,
                    WmatrixIntegratorType,
                    XmatrixIntegratorType,
                    YmatrixIntegratorType,
                    ZmatrixIntegratorType,
                    EmatrixIntegratorType,
                    RmatrixIntegratorType,
                    H1_IntegratorType,
                    H2_IntegratorType,
                    H3_IntegratorType >
        StokesIntegratorTuple;

    template < class RowSpace, class ColSpace >
    struct magic {
            typedef typename DiscreteFunctionSelector< ThisType, RowSpace >::Type
                RowType;
            typedef typename DiscreteFunctionSelector< ThisType, ColSpace >::Type
                ColType;
        typedef typename MatrixObject< RowType, ColType >::Type
            InternalMatrixType;
        typedef Dune::shared_ptr<InternalMatrixType>
            PointerType;
    };
    template < class F, class G >
    static auto matrix( const F& f, const G& g ) -> typename magic<F,G>::PointerType
    {
        typedef typename magic<F,G>::InternalMatrixType
            InternalMatrixType;
        typename magic<F,G>::PointerType m( new InternalMatrixType(f,g) );
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
    static OmatrixIntegratorType integratorO( YmatrixType& g,
                                                          const DiscreteVelocityFunctionType& f )
    {
        return OmatrixIntegratorType( g, f );
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
#endif // DUNE_OSEEN_INTEGRATORS_FACTORY_HH
