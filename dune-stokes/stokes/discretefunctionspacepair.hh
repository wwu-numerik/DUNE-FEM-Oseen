/**
 *  \file   discretefunctionspacepair.cc
 *
 *  \brief  brief
 **/

#ifndef DISCRETEFUNCTIONSPACEPAIR_HH_INCLUDED
#define DISCRETEFUNCTIONSPACEPAIR_HH_INCLUDED

#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune
{

// forward
template < class DiscreteFunctionPairTraitsImp >
class DiscreteFunctionPair;

template < class T1, class T2, class DiscreteFunctionSpacePairImp >
class DiscreteFunctionPairTraits
{
    public:
        typedef DiscreteFunctionSpacePairImp
            DiscreteFunctionSpacePairType;

        typedef T1
            FirstDiscreteFunctionType;

        typedef T2
            SecondDiscreteFunctionType;
};

template < class DiscreteFunctionPairTraitsImp >
class DiscreteFunctionPair
{
    public:
        typedef DiscreteFunctionPairTraitsImp
            Traits;

        typedef typename Traits::DiscreteFunctionSpacePairType
            DiscreteFunctionSpacePairType;

        typedef typename Traits::FirstDiscreteFunctionType
            FirstDiscreteFunctionType;

        typedef typename Traits::SecondDiscreteFunctionType
            SecondDiscreteFunctionType;

        // names schould be a pair
        DiscreteFunctionPair( const std::string name1, const std::string name2, DiscreteFunctionSpacePairType& spacePair )
            : firstFunction_( name1, spacePair.first() ),
            secondFunction_( name2, spacePair.second() )
        {}

        ~DiscreteFunctionPair()
        {}

    private:
        FirstDiscreteFunctionType firstFunction_;
        SecondDiscreteFunctionType secondFunction_;
};

// forward
template < class DiscreteFunctionSpacePairTraitsImp >
class DiscreteFunctionSpacePair;

template < class T1, class T2 >
class DiscreteFunctionSpacePairTraits
{
    public:
        //! CRTP
        typedef DiscreteFunctionSpacePair< DiscreteFunctionSpacePairTraits < T1, T2 > >
            DiscreteFunctionSpaceType;

        typedef Pair< typename T1::FunctionSpaceType, typename T2::FunctionSpaceType >
            FunctionSpaceType;

        typedef Pair< T1, T2 >
            DiscreteFunctionSpaces;

        typedef Pair < typename T1::BaseFunctionSetType, typename T2::BaseFunctionSetType >
            BaseFunctionSetType;

        typedef Pair < typename T1::BlockMapperType, typename T2::BlockMapperType >
            BlockMapperType;

        typedef Pair < typename T1::GridPartType, typename T2::GridPartType >
            GridPartType;

      typedef T1
            FirstDiscreteFunctionSpaceType;

        typedef T2
            SecondDiscreteFunctionSpaceType;


        typedef Pair <  typename FirstDiscreteFunctionSpaceType::MapperType,
                        typename SecondDiscreteFunctionSpaceType::MapperType >
            MapperType;


        typedef Pair <  typename FirstDiscreteFunctionSpaceType::EntityType,
                        typename SecondDiscreteFunctionSpaceType::EntityType >
            EntityType;

        typedef Pair <  typename FirstDiscreteFunctionSpaceType::GridType,
                        typename SecondDiscreteFunctionSpaceType::GridType >
            GridType;

        typedef Pair <  typename FirstDiscreteFunctionSpaceType::IndexSetType,
                        typename SecondDiscreteFunctionSpaceType::IndexSetType >
            IndexSetType;

        typedef Pair <  typename FirstDiscreteFunctionSpaceType::IteratorType,
                        typename SecondDiscreteFunctionSpaceType::IteratorType >
            IteratorType;
};

template < class DiscreteFunctionSpacePairTraitsImp >
class DiscreteFunctionSpacePair
//    : public DiscreteFunctionSpaceInterface< DiscreteFunctionSpacePairTraitsImp >
{
    public:
        typedef DiscreteFunctionSpaceInterface< DiscreteFunctionSpacePairTraitsImp >
            BaseType;

        typedef DiscreteFunctionSpacePairTraitsImp
            Traits;

        typedef typename Traits::DiscreteFunctionSpaces
            DiscreteFunctionSpaces;

        typedef typename DiscreteFunctionSpaces::Type1
            FirstDiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaces::Type2
            SecondDiscreteFunctionSpaceType;

        typedef typename Traits::BaseFunctionSetType
            BaseFunctionSetType;

        typedef typename Traits::MapperType
            MapperType;

        typedef typename Traits::BlockMapperType
            BlockMapperType;

        typedef typename Traits::EntityType
            EntityType;

        typedef typename Traits::GridType
            GridType;

        typedef typename Traits::GridPartType
            GridPartType;

        typedef typename Traits::IndexSetType
            IndexSetType;

        typedef typename Traits::IteratorType
            IteratorType;


        DiscreteFunctionSpacePair( GridPartType& gridPart )
         //   : BaseType(),
            :firstSpace_( gridPart.first() ),
            secondSpace_( gridPart.second() )
        {}

        ~DiscreteFunctionSpacePair()
        {}

        FirstDiscreteFunctionSpaceType& first()
        {
            return firstSpace_;
        }

        SecondDiscreteFunctionSpaceType& second()
        {
            return secondSpace_;
        }


        /**
         *  \brief return type identifier of discrete function space
         *  \return return type identifier of discrete function space
         **/
        DFSpaceIdentifier type() const
        {}

        /** \brief get base function set for given entity
         *
         *  \param[in]  entity  entity (of codim 0) for which base function is
         *                      requested
         *
         *  \returns BaseFunctionSet for the entity
         */
        inline const BaseFunctionSetType baseFunctionSet( const EntityType &entity ) const
        {}

        /** \brief returns true if the space contains DoFs of given codimension
         *
         *  \param[in]  codim  codimension to check for DoFs
         *
         *  \returns \b true if codimension contains DoFs,
         *           \b false otherwise
         */
        inline bool contains( const int codim ) const
        {
            return false;
        }

        /** \brief returns true if the space contains only globally continuous
         *         functions
         *
         *  For example, a \ref Dune::LagrangeDiscreteFunctionSpace
         *  "Lagrange space" returns \b true while a \ref
         *  Dune::DiscontinuousGalerkinSpace "discontiuous Galerkin space" returns
         *  \b false.
         *
         *  \returns \b true  if the space contians only globally continous
         *                    functions,
         *           \b false otherwise
         */
        inline bool continuous() const
        {
            return false;
        }

        /** \brief get index of the sequence in grid sequences
         *
         *  \return number of current sequence
         **/
        inline int sequence() const
        {
            return 0;
        }

        /** \brief get global order of space
         *
         *  \return order of space, i.e., the maximal polynomial order of base
         *          functions
         **/
        inline int order() const
        {
            return 0;
        }

        /** \brief get a reference to the DoF mapper
         *
         *  \returns reference to mapper
         */
        inline MapperType& mapper() const
        {}

        /** \brief get a reference to the block mapper
         *
         *  \returns refernce to the block mapper
         */
        inline BlockMapperType &blockMapper () const
        {}

        /** \brief get reference to grid this discrete function space belongs to
         *
         *  \returns constant reference to grid
         */
        inline const GridType &grid () const
        {}

        /** \brief get reference to grid this discrete function space belongs to
         *
         *  \returns reference to grid
         */
        inline GridType &grid ()
        {}

        /** \brief get a reference to the associated grid partition
         *
         *  \returns constant reference to the grid partition
         */
        inline const GridPartType &gridPart () const
        {}

        /** \brief get a reference to the associated grid partition
         *
         *  \returns reference to the grid partition
         */
        inline GridPartType &gridPart ()
        {}

        /** \brief Get a reference to the associated index set
         *
         *  \returns const reference to index set
         */
        inline const IndexSetType &indexSet () const
        {}

        /** \brief get number of DoFs for this space
         *
         *  \returns number of DoFs (degrees of freedom)
         */
        inline int size () const
        {}

        /** \brief get iterator pointing to the first entity of the associated grid
         *         partition
         *
         *  \returns iterator pointing to first entity
         */
        inline IteratorType begin () const
        {}

        /** \brief get iterator pointing behind the last entity of the associated
         *         grid partition
         *
         *  \returns iterator pointing behind last entity
         */
        inline IteratorType end () const
        {}

        /** \brief apply a functor to each entity in the associated grid partition
         *
         *  The functor must provide an the following operator
         *  \code
         *  template< class EntityType >
         *  void operator() ( const EntityType & );
         *  \endcode
         *
         *  \param[in]  f  functor to apply
         */
        template< class FunctorType >
        inline void forEach ( FunctorType& f ) const
        {}

        /** \brief returns true if the grid has more than one geometry type
         *
         *  \return \b true if the underlying grid has more than one geometry type
         *          (hybrid grid), \b false otherwise
         */
        inline bool multipleGeometryTypes () const
        {}

        /** \brief returns true if base function sets depend on the entity
         *
         *  \returns \b true if base function set depend on entities, \b false
         *           otherwise
         */
        inline bool multipleBaseFunctionSets () const
        {}

        /** \brief Map local DoF number to global DoF number
         *
         *  Maps an entity and a local DoF number to a global DoF number, i.e.,
         *  the index of the DoF within the DoF vector.
         *
         *  \param[in]  entity    entity (of codim 0) for which the mapping is done
         *  \param[in]  localDof  local dof number
         *
         *  \returns global DoF number
         */
        inline int mapToGlobal ( const EntityType& entity,
                                 const int localDof ) const
        {}

        /** \brief return maximal number of local DoFs
         *
         *  \returns an upper bound for the number of local DoFs
         */
        inline int maxNumLocalDofs () const
        {}


        /** \brief defines type of data handle for communication
         *  \param  DiscreteFunction  type of \ref Dune::DiscreteFunctionInterface
         *                            "discrete function" to communicate
         *  \param  Operation         type of operation to perform on communication
         *                            (defaults to copy)
         */
        template< class DiscreteFunction,
                  class Operation = DFCommunicationOperation :: Copy >
        struct CommDataHandle
        {
          //! type of communication data handle
          typedef typename Traits
            :: template CommDataHandle< DiscreteFunction, Operation > :: Type
            Type;

          //! type of operation to perform on scatter
          typedef typename Traits
            :: template CommDataHandle< DiscreteFunction, Operation > :: OperationType
            OperationType;
        };

        /** \brief Creates DataHandle for given discrete function
         *
         *  \param[in]  discreteFunction  \ref DiscreteFunctionInterface
         *                                "discrete function" to create the data
         *                                handle for
         *  \param[in]  operation         operation to perform on scatter
         */
        template< class DiscreteFunction, class Operation >
        inline typename CommDataHandle< DiscreteFunction, Operation > :: Type
        createDataHandle ( DiscreteFunction& discreteFunction,
                          const Operation* operation ) const
        {}

    private:
        FirstDiscreteFunctionSpaceType firstSpace_;
        SecondDiscreteFunctionSpaceType secondSpace_;

};
};

#endif // end of discretefunctionspacepair.cc
