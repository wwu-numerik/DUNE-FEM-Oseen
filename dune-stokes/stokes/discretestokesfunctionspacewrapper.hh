/**
 *  \file   discretestokesfunctionspacewrapper.hh
 *
 *  \brief  brief
 **/

#ifndef DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED
#define DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED

#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune
{

// forward
template < class DiscreteStokesFunctionSpaceWrapperTraitsImp >
class DiscreteStokesFunctionSpaceWrapper;

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
                    DiscreteStokesFunctionSpaceWrapperTraits<   VelocitySpace,
                                                                PressureSpace > >
            DiscreteFunctionSpaceType;

        //! type of function space
        typedef typename DiscreteVelocityFunctionSpace
            FunctionSpaceType;

        //! type of base function set
        typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper
        typedef typename DiscreteFunctionSpaceType::MapperType
            MapperType;

        //! type of block mapper
        typedef typename DiscreteFunctionSpaceType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part
        typedef typename DiscreteFunctionSpaceType::GridPartType
            GridPartType;
        /**
         *  \}
         **/

}; // end of DiscreteStokesFunctionSpaceWrapperTraits

template < class DiscreteStokesFunctionSpaceWrapperTraitsImp >
class DiscreteStokesFunctionSpaceWrapper
    : public DiscreteFunctionSpaceInterface< DiscreteStokesFunctionSpaceWrapperTraitsImp >
{
    public:

        typedef DiscreteStokesFunctionSpaceWrapperTraitsImp
            Traits;

        //! base type
        typedef DiscreteFunctionSpaceInterface< DiscreteStokesFunctionSpaceWrapperTraitsImp >
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
            Entity;

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
        DiscreteVelocityFunctionSpaceType discreteVelocitySpace() const
        {
            return velocitySpace_;
        }

        /**
         *  \brief  return the discrete pressure space
         *  \todo   doc
         **/
        DiscretePressureFunctionSpaceType discretePressureSpace() const
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
            return velocitySpace_.sequence()
        };

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








    private:

        DiscreteVelocityFunctionSpaceType velocitySpace_;
        DiscretePressureFunctionSpaceType pressureSpace_;


}; // end of DiscreteStokesFunctionSpaceWrapper

}; // end of namespace Dune

#ifndef // end of discretestokesfunctionspacewrapper.hh
