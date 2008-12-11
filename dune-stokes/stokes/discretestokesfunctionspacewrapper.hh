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

template < class VelocitySpace, class PressureSpace >
class DiscreteStokesFunctionSpaceWrapperTraits
{
    public:

        //! type of discrete velocity function space
        typedef VelocitySpace
            DiscreteVelocityFunctionSpace;

        //! type of discrete pressure function space
        typedef PressureSpace
            DiscretePressureFunctionSpace;

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
        typedef typename FunctionSpaceType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper
        typedef typename FunctionSpaceType::MapperType
            MapperType;

        //! type of block mapper
        typedef typename FunctionSpaceType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part
        typedef typename FunctionSpaceType::GridPartType
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

        //! own type (CRTP) (for interface compliance)
        typedef typename Traits::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        //! type of discrete velocity function space
        typedef typename Traits::DiscreteVelocityFunctionSpace
            DiscreteVelocityFunctionSpace;

        //! type of discrete pressure function space
        typedef typename Traits::DiscretePressureFunctionSpace
            DiscretePressureFunctionSpace;



}; // end of DiscreteStokesFunctionSpaceWrapper

}; // end of namespace Dune

#ifndef // end of discretestokesfunctionspacewrapper.hh
