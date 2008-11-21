/**
 *  \file   model.hh
 *
 *  \brief  contains a class StokesModel with traits class StokesModelTraits
 **/

#ifndef MODEL_HH
#define MODEL_HH

#include "logging.hh"
#include "problem.hh"

/**
 *  \brief  containing typedefs needed by StokesModel
 *
 *  \tparam gridDim
 *          dimension of the grid
 **/
template < int gridDim >
class StokesModelTraits
{
    public:
        typedef typename ProblemTraits< gridDim >::ForceType
            RightHandSideType;
        typedef typename ProblemTraits< gridDim >::DirichletDataType
            DirichletDataType;
};


/**
 *  \brief
 *
 *  \tparam gridDim
 *          dimension of the grid
 *
 *  \todo   extensive docu with latex
 **/
template < int gridDim >
class StokesModel
{
    public:
        typedef StokesModelTraits< gridDim >
            Traits;
        typedef typename Traits::RightHandSideType
            RightHandSideType,


        /**
         *  \brief constructor
         *
         *  doing nothing
         **/
        StokesModel()
        {
        }

        /**
         *  \brief destructor
         *
         *  doing nothing
         **/
        ~StokesModel()
        {
        }

    private:

};

#endif  // end of model.hh
