/**
 *  \file   discretestokesfunctionspacewrapper.hh
 *  \todo   doc
 **/

#ifndef DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED
#define DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED

#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune
{

// forward
template < class DiscreteStokesFunctionSpaceWrapperTraitsImp >
class DiscreteStokesFunctionSpaceWrapper;

/**
 *  \todo   doc
 **/
template < class DiscreteVelocitySpaceImp, class DiscretePressureSpaceImp >
class DiscreteStokesFunctionSpaceWrapperTraits
{
    public:

        //! type of discrete velocity function space
        typedef DiscreteVelocitySpaceImp
            DiscreteVelocityFunctionSpaceType;

        //! type of discrete pressure function space
        typedef DiscretePressureSpaceImp
            DiscretePressureFunctionSpaceType;

        /**
         *  \name for interface compliance
         *  \{
         **/
        //! own type (CRTP) (for interface compliance)
        typedef DiscreteStokesFunctionSpaceWrapper<
                    DiscreteStokesFunctionSpaceWrapperTraits<   DiscreteVelocityFunctionSpaceType,
                                                                DiscretePressureFunctionSpaceType > >
            DiscreteFunctionSpaceType;

        //! type of function space
        typedef typename DiscreteVelocityFunctionSpaceType::FunctionSpaceType
            FunctionSpaceType;

        //! type of base function set
        typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper
        typedef typename DiscreteVelocityFunctionSpaceType::MapperType
            MapperType;

        //! type of block mapper
        typedef typename DiscreteVelocityFunctionSpaceType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;
        /**
         *  \}
         **/
}; // end of DiscreteStokesFunctionSpaceWrapperTraits

/**
 *  \todo   doc
 **/
template < class DiscreteStokesFunctionSpaceWrapperTraitsImp >
class DiscreteStokesFunctionSpaceWrapper
    : public DiscreteFunctionSpaceDefault< DiscreteStokesFunctionSpaceWrapperTraitsImp >
{
    public:

        typedef DiscreteStokesFunctionSpaceWrapperTraitsImp
            Traits;

        //! base type
        typedef DiscreteFunctionSpaceDefault< DiscreteStokesFunctionSpaceWrapperTraitsImp >
            BaseType;

        //! type of discrete velocity function space
        typedef typename Traits::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! type of discrete pressure function space
        typedef typename Traits::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! own type (CRTP)
        typedef typename BaseType::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        //! type of function space (of the discrete velocity function space)
        typedef typename BaseType::FunctionSpaceType
            FunctionSpaceType;

        //! type of base function set (of the discrete velocity function space)
        typedef typename BaseType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper (of the discrete velocity function space)
        typedef typename BaseType::MapperType
            MapperType;

        //! type of block mapper (of the discrete velocity function space)
        typedef typename BaseType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part (of the discrete velocity function space)
        typedef typename BaseType::GridPartType
            GridPartType;

        //! type of underlying grid (of the discrete velocity function space)
        typedef typename GridPartType::GridType
            GridType;

        //! type of used index set (of the discrete velocity function space)
        typedef typename GridPartType::IndexSetType
            IndexSetType;

        //! type of iterator of codim 0 entities for grid traversal (of the discrete velocity function space)
        typedef typename GridPartType::template Codim< 0 >::IteratorType
            IteratorType;

        //! type of codim 0 entity (of the discrete velocity function space)
        typedef typename IteratorType::Entity
            EntityType;

        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        DiscreteStokesFunctionSpaceWrapper( GridPartType& gridPart )
            : BaseType( gridPart ),
            velocitySpace_( gridPart ),
            pressureSpace_( gridPart )
        {}

        /**
         *  \brief  destructor
         **/
        ~DiscreteStokesFunctionSpaceWrapper()
        {}

        /**
         *  \brief  returns the discrete velocity space
         *  \todo   doc
         **/
        DiscreteVelocityFunctionSpaceType& discreteVelocitySpace()
        {
            return velocitySpace_;
        }
        const DiscreteVelocityFunctionSpaceType& discreteVelocitySpace() const
        {
            return velocitySpace_;
        }

        /**
         *  \brief  return the discrete pressure space
         *  \todo   doc
         **/
        DiscretePressureFunctionSpaceType& discretePressureSpace()
        {
            return pressureSpace_;
        }
        const DiscretePressureFunctionSpaceType& discretePressureSpace() const
        {
            return pressureSpace_;
        }

        /**
         *  \brief  return type identifier of the discrete velocity function space
         **/
        DFSpaceIdentifier type() const
        {
            return velocitySpace_.type();
        }

        /**
         *  \brief  get base function set of the discrete velocity function space for given entity
         *  \todo   doc
         **/
        inline const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
        {
            return velocitySpace_.baseFunctionSet( entity );
        }

        /**
         *  \brief  returns true if the discrete velocity function space contains DoFs of given codimension
         *  \todo   doc
         **/
        inline bool contains( const int codim ) const
        {
            return velocitySpace_.contains( codim );
        }

        /**
         *  \brief  returns true if the discrete velocity function space contains only globally continuous functions
         *  \todo   doc
         **/
        inline bool continuous() const
        {
            return velocitySpace_.continuous();
        }

        /**
         *  \brief  get index of the sequence in the discrete velocity function spaces grid sequences
         *  \todo   doc
         **/
        inline int sequence() const
        {
            return velocitySpace_.sequence();
        }

        /**
         *  \brief  get global order of the discrete velocity function space
         *  \todo doc
         **/
        inline int order() const
        {
            return velocitySpace_.order();
        }

        /**
         *  \brief  get a reference to the discrete velocity function spacees DoF mapper
         *  \todo   doc
         **/
        inline MapperType& mapper() const
        {
            return velocitySpace_.mapper();
        }

        /**
         *  \brief  get a reference to the discrete velocity function spaces block mapper
         *  \todo   doc
         **/
        inline BlockMapperType& blockMapper() const
        {
            return velocitySpace_.blockMapper();
        }

        /**
         *  \brief  get reference to grid the discrete velocity function space belongs to
         *  \todo   doc
         **/
        inline const GridType& grid() const
        {
            return velocitySpace_.grid();
        }

        /**
         *  \brief  get reference to grid the discrete velocity function space belongs to
         *  \todo   doc
         **/
        inline GridType& grid()
        {
            return velocitySpace_.grid();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space grid partition
         *  \todo   doc
         **/
        inline const GridPartType& gridPart() const
        {
            return velocitySpace_.gridPart();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space grid partition
         *  \todo   doc
         **/
        inline GridPartType& gridPart()
        {
            return velocitySpace_.gridPart();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space index set
         *  \todo   doc
         **/
        inline const IndexSetType& indexSet() const
        {
            return velocitySpace_.indexSet();
        }

        /**
         *  \brief  get number of DoFs for the discrete velocity function space
         *  \todo   doc
         **/
        inline int size() const
        {
            return velocitySpace_.size();
        }

        /**
         *  \brief  get discrete velocity function space iterator pointing to the first entity of the associated grid partition
         *  \todo   doc
         **/
        inline IteratorType begin() const
        {
            return velocitySpace_.begin();
        }

        /**
         *  \brief  get discrete velocity function space iterator pointing behind the last entity of the associated grid partition
         *  \todo   doc
         **/
        inline IteratorType end() const
        {
            return velocitySpace_.end();
        }

        /**
         *  \brief  apply discrete velocity function space functor to each entity in the associated grid partition
         *  \todo   doc
         **/
        template< class FunctorType >
        inline void forEach( FunctorType& f ) const
        {
            return velocitySpace_.forEach( f );
        }

        /**
         *  \brief  returns true if the discrete velocity function space grid has more than one geometry type
         *  \todo   doc
         **/
        inline bool multipleGeometryTypes() const
        {
            return velocitySpace_.multipleGeometryTypes();
        }

        /**
         *  \brief  returns true if discrete velocity fucntion space base function sets depend on the entity
         *  \todo   doc
         **/
        inline bool multipleBaseFunctionSets() const
        {
            return velocitySpace_.multipleBaseFunctionSets();
        }

        /**
         *  \brief  map local DoF number to global DoF number (discrete velocity function space)
         *  \todo   doc
         **/
        inline int mapToGlobal( const EntityType& entity,
                                const int localDof ) const
        {
            return velocitySpace_.mapToGlobal( entity, localDof );
        }

        /**
         *  \brief  return maximal number of local DoFs in discrete velocity funtion space
         *  \todo   doc
         **/
        inline int maxNumLocalDofs() const
        {
            return velocitySpace_.maxNumLocalDofs();
        }

        /**
         *  \brief  creates DataHandle for given discrete function (from dicrete velocity function space)
         *  \todo   doc
         **/
        template< class DiscreteFunction, class Operation >
        inline typename BaseType::template CommDataHandle< DiscreteFunction, Operation >::Type createDataHandle(    DiscreteFunction& discreteFunction,
                                                                                                                    const Operation* operation ) const
        {
            return velocitySpace_.createDataHandle( discreteFunction,
                                                    operation );
        }

    private:

        DiscreteVelocityFunctionSpaceType velocitySpace_;
        DiscretePressureFunctionSpaceType pressureSpace_;


}; // end of DiscreteStokesFunctionSpaceWrapper

//! forward
template < class DiscreteStokesFunctionWrapperTraitsImp >
class DiscreteStokesFunctionWrapper;

/**
 *  \todo   doc
 **/
template <  class DiscreteStokesFunctionSpaceWrapperImp,
            class DiscreteVelocityFunctionImp,
            class DiscretePressureFunctionImp >
class DiscreteStokesFunctionWrapperTraits
{
    public:

        //! own type (CRTP)
        typedef DiscreteStokesFunctionWrapper<
                    DiscreteStokesFunctionWrapperTraits<
                        DiscreteStokesFunctionSpaceWrapperImp,
                        DiscreteVelocityFunctionImp,
                        DiscretePressureFunctionImp > >
            DiscreteFunctionType;

        //! type of associated discrete function space
        typedef DiscreteStokesFunctionSpaceWrapperImp
            DiscreteFunctionSpaceType;

        //! type of discrete velocity function
        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;

        //! type of discrete pressure function
        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;

}; // end of DiscreteStokesFunctionWrapperTraits

/**
 *  \todo   doc,
 *          should comply with the DiscreteFunctionInterface some time
 **/
template < class DiscreteStokesFunctionWrapperTraitsImp >
class DiscreteStokesFunctionWrapper
//    : public DiscreteFunctionInterface< DiscreteStokesFunctionWrapperTraitsImp >
{
    public:

        typedef DiscreteStokesFunctionWrapperTraitsImp
            Traits;

        typedef typename Traits::DiscreteFunctionType
            DiscreteFunctionType;

        //! DiscreteStokesFunctionSpaceWrapper
        typedef typename Traits::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        //! type of discrete velocity function
        typedef typename Traits::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! type of discrete pressure function
        typedef typename Traits::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! type of range field
        typedef typename DiscreteVelocityFunctionType::RangeFieldType
            RangeFieldType;

        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        DiscreteStokesFunctionWrapper(  const std::string name,
                                        DiscreteFunctionSpaceType& space )
            : space_( space ),
            velocity_( name + std::string("_velocity"), space.discreteVelocitySpace() ),
            pressure_( name + std::string("_pressure"), space.discretePressureSpace() )
        {}

        /**
         *  \brief      constructor
         *
         *              wrappes existing velocity and pressure
         *  \attention  by copying
         **/
        DiscreteStokesFunctionWrapper(  DiscreteFunctionSpaceType& space,
                                        DiscreteVelocityFunctionType& velocity,
                                        DiscretePressureFunctionType& pressure )
            : space_( space ),
            velocity_( velocity ),
            pressure_( pressure )
        {}

        /**
         *  \brief  destructor
         *  \todo   doc
         **/
        ~DiscreteStokesFunctionWrapper()
        {}

        /**
         *  \todo   doc
         **/
        const DiscreteVelocityFunctionType& discreteVelocity() const
        {
            return velocity_;
        }

        /**
         *  \todo   doc
         **/
        DiscreteVelocityFunctionType& discreteVelocity()
        {
            return velocity_;
        }

        /**
         *  \todo   doc
         **/
        const DiscretePressureFunctionType& discretePressure() const
        {
            return pressure_;
        }

        /**
         *  \todo   doc
         **/
        DiscretePressureFunctionType& discretePressure()
        {
            return pressure_;
        }

        /**
         *  \brief  obtain the name of the discrete function
         *  \todo   doc
         **/
        inline const std::string &name() const
        {
            return velocity_.name();
        }

//        //! obtain a local function for an entity (read-only)
//        template< class EntityType >
//        inline const LocalFunctionType localFunction( const EntityType& entity ) const
//        {
//            return velocity_.localFunction( entity );
//        }
//
//        //! obtain a local function for an entity
//        template< class EntityType >
//        inline LocalFunctionType localFunction( const EntityType& entity )
//        {
//            return velocity_.localFunction( entity );
//        }

        /**
         *  \brief  set all degrees of freedom to zero
         *  \todo   doc
         **/
        inline void clear()
        {
            velocity_.clear();
        }

        /**
         *  \brief  obtain total number of DoFs
         *  \todo   doc
         **/
        inline int size() const
        {
            return velocity_.size();
        }

        /**
         *  \brief  add another discrete function to this one
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator+= ( const DiscreteFunctionType& arg )
        {
            velocity_ += arg.discreteVelocity();
            return *this;
        }

        /**
         *  \brief  substract all degrees of freedom from given discrete function using the dof iterators
         *  \todo   doc
         **/
        template < class DFType >
        DiscreteFunctionType& operator-=( const DFType& /*arg*/ )
        {
            assert( false );
            return *this;
        }

        /**
         *  \brief  multiply all DoFs by a scalar factor
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator*=( const RangeFieldType& scalar )
        {
            velocity_ *= scalar;
            return *this;
        }

        /**
         *  \brief  devide all DoFs by a scalar factor
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator/=( const RangeFieldType& scalar )
        {
            velocity_ /= scalar;
            return *this;
        }


    private:

        const DiscreteFunctionSpaceType& space_;
        DiscreteVelocityFunctionType velocity_;
        DiscretePressureFunctionType pressure_;

}; // end of DiscreteStokesFunctionWrapper

}; // end of namespace Dune

#endif // end of discretestokesfunctionspacewrapper.hh
