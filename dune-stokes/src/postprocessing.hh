/**
 *  \file   postprocessing.hh
 *  \brief  postprocessing.hh
 **/


template <  class GridPartImp,
            class DiscreteVelocityFunctionImp,
            class DiscretePressureFunctionImp >
class PostProcessor
{
    public:
        typedef GridPartImp
            GridPartType;
        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;
        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;
        PostProcessor( const GridPartType& gridPart )
            : gridPart_( gridPart ),
            exactVelocity_( "u_exact", gridPart ),
            exactPressure_( "p_exact", gridPart )
        {
        }

        ~PostProcessor()
        {
        }

        void assembleExactSolution()
        {
            typedef typename GridPartType::Traits::template Codim< 0 >::IteratorType
                EntityIteratorType;
            EntityIteratorType itEnd = gridPart_.end();
            for ( EntityIteratorType it = gridPart_.begin(); it != itEnd; ++it )
            {
                typedef typename GridPartType::EntityCodim0Type
                    EntityType;
                EntityType& entity = *it;
                typedef typename GridPartType::template Codim< 0 >::Geometry
                    GeometryType;
                const GeometryType& geometry = entity->geometry();
                typedef typename DiscreteVelocityFunctionType::LocalFunctionType
                    LocalVelocityFunctionType;
                LocalVelocityFunctionType localVelocity( entity );
                const int numDofs = localVelocity.numDofs();
                typedef typename LocalVelocityFunctionType::RangeFieldType
                    RangeFieldType;

                for ( int dof = 0; dof < numDofs; ++dof )
                {
                    RangeFieldType localCoord = localVelocity[dof];
                }

                typedef typename DiscretePressureFunctionType::LocalFunctionType
                    LocalPressureFunctionType;

            }
        }

    private:
        GridPartType& gridPart_;
        DiscreteVelocityFunctionType exactVelocity_;
        DiscretePressureFunctionType exactPressure_;
};
