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
// should be derived from interface some day
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

        typedef Pair<   typename FirstDiscreteFunctionType::RangeFieldType,
                        typename SecondDiscreteFunctionType::RangeFieldType >
            RangeFieldType;

        DiscreteFunctionPair(   const Dune::Pair< const std::string, const std::string > namePair,
                                DiscreteFunctionSpacePairType& spacePair )
            : firstFunction_( namePair.second(), spacePair.first() ),
            secondFunction_( namePair.second(), spacePair.second() )
        {}

        ~DiscreteFunctionPair()
        {}

        FirstDiscreteFunctionType& first()
        {
            return firstFunction_;
        }

        const FirstDiscreteFunctionType& first() const
        {
            return firstFunction_;
        }

        SecondDiscreteFunctionType& second()
        {
            return secondFunction_;
        }

        const SecondDiscreteFunctionType& second() const
        {
            return secondFunction_;
        }

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

        typedef typename Traits::GridPartType
            GridPartType;

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

        const FirstDiscreteFunctionSpaceType& first() const
        {
            return firstSpace_;
        }

        SecondDiscreteFunctionSpaceType& second()
        {
            return secondSpace_;
        }

        const SecondDiscreteFunctionSpaceType& second() const
        {
            return secondSpace_;
        }

    private:
        FirstDiscreteFunctionSpaceType firstSpace_;
        SecondDiscreteFunctionSpaceType secondSpace_;

};
};

#endif // end of discretefunctionspacepair.cc
