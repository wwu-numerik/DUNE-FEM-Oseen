/**
 *  \file   discretegradientpass.hh
 *  \brief  contains a class DiscreteGradientPass
 *
 *          see R. Eymard, T. Gallouet, R. Herbin, 2006, "A cell-centered
 *          finite-volume approximation for anisotropic diffusion operators on
 *          unstructured meshes in any space dimension" for details
 **/

#ifndef DISCRETEGRADIENTPASS_HH
#define DISCRETEGRADIENTPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/space/dgspace.hh>

#ifndef NLOG
    #include "../src/stuff.hh"
    #include "../src/logging.hh"
#endif

namespace Dune{

// forward decraration
template < class DiscreteGradientModelTraitsImp >
class DiscreteGradientModel;

/**
 *  \brief DiscreteGradientModelTraits
 *  \todo   doc
 **/
template < class GridPartImp, class DiscreteFunctionImp >
class DiscreteGradientModelTraits
{
    public:

        //! CRTP trick
        typedef DiscreteGradientModel< DiscreteGradientModelTraits >
            DiscreteModelType;

        //! grid part type
        typedef GridPartImp
            GridPartType;

        /**
         *  \name   types needed for the pass
         *  \{
         **/
        //! return type of the pass
        typedef DiscreteFunctionImp
            DestinationType;
        /**
         *  \}
         **/

};

/**
 *  \brief  DiscreteGradientModel
 *  \todo   doc
 **/
template < class DiscreteGradientModelTraitsImp >
class DiscreteGradientModel
{
    public:

        //! traits
        typedef DiscreteGradientModelTraitsImp
            Traits;

        //! grid part type
        typedef typename Traits::GridPartType
            GridPartType;

};

/**
 *  \brief  DiscreteGradientPass
 *  \todo   doc
 **/
template <  class DiscreteModelImp,
            class PreviousPassImp,
            int PassID = 0 >
class DiscreteGradientPass
    : public Pass < DiscreteModelImp, PreviousPassImp, PassID >
{
    public:

        //! base type
        typedef Pass < DiscreteModelImp, PreviousPassImp, PassID >
            BaseType;

        //! previous pass type
        typedef PreviousPassImp
            PreviousPassType;

        //! discrete model type
        typedef DiscreteModelImp
            DiscreteModelType;

        //! grd part type
        typedef typename DiscreteModelType::GridPartType
            GridPartType;

        /**
         *  \name typedefs needed for interface compliance
         *  \{
         **/
        typedef typename BaseType::DestinationType
            DestinationType;

        typedef typename BaseType::DomainType
            DomainType;

        typedef typename BaseType::RangeType
            RangeType;

        typedef typename BaseType::TotalArgumentType
            TotalArgumentType;
        /**
         *  \}
         **/

        typedef DomainType
            DiscreteFunctionType;

        typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        DiscreteGradientPass(   PreviousPassType& prevPass,
                                DiscreteModelType& discreteModel,
                                GridPartType& gridPart,
                                DiscreteFunctionSpaceType& space )
            : BaseType( prevPass ),
            discreteModel_( discreteModel ),
            gridPart_( gridPart ),
            space_( space ){}

        /**
         *  \brief destructor
         *  \todo   doc
         **/
        ~DiscreteGradientPass(){}

        /**
         *  \brief main method
         *  \todo   doc
         **/
        virtual void apply( const DomainType &arg, RangeType &dest) const
        {

#ifndef NLOG
            // logging stuff
            Logging::LogStream& infoStream = Logger().Info();
            Logging::LogStream& debugStream = Logger().Dbg();
#endif

            //! type of the grid
            typedef typename GridPartType::GridType
                GridType;

            //! entity iterator of the gridpart
            typedef typename GridPartType::template Codim< 0 >::IteratorType
                EntityIteratorType;

            //! type of codim 0 entity
            typedef typename GridType::template Codim< 0 >::Entity
                EntityType;

            //! Intersection iterator of the gridpart
            typedef typename GridPartType::IntersectionIteratorType
                IntersectionIteratorType;

#ifndef NLOG
            infoStream << "  \nthis is DiscreteGradientPass::apply()" << std::endl;

            // walk the grid
            EntityIteratorType entityItEnd = space_.end();
            for (   EntityIteratorType entityIt = space_.begin();
                    entityIt != entityItEnd;
                    ++entityIt ) {
                const EntityType& entity = *entityIt;
                // walk the intersections
                IntersectionIteratorType intItEnd = gridPart_.iend( entity );
                for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
                        intIt != intItEnd;
                        ++intIt ) {
                    // if we are inside the grid
                    if ( intIt.neighbor() && !intIt.boundary() ) {
                    }
                    // if we are on the boundary of the grid
                    if ( !intIt.neighbor() && intIt.boundary() ) {
                    }
                }
            }

        } // end of apply()

        /**
         *  \name methods needed for interface compliance
         *  \{
         **/
        virtual void compute( const TotalArgumentType &arg, DestinationType &dest) const
        {}

        virtual void allocateLocalMemory()
        {}
        /**
         *  \}
         **/


    private:

        DiscreteModelType& discreteModel_;
        GridPartType& gridPart_;
        DiscreteFunctionSpaceType& space_;
}; // end of DiscreteGradientPass

}

#endif // end of dicretegradientpass.hh
