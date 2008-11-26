#ifndef DUNE_PROBLEM_HH
#define DUNE_PROBLEM_HH

namespace Dune
{

// **** Problem data: right hand side f of governing problem
template <class FunctionSpaceType>
class myRHSFunction
    : public Function< FunctionSpaceType, myRHSFunction<FunctionSpaceType> >
{
private:
    typedef myRHSFunction<FunctionSpaceType> ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

public:
    inline myRHSFunction ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate( const DomainType &x , RangeType &phi ) const
    {
    enum { dimension = DomainType :: dimension };
    #if 0
    // f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
    phi = 0;
    for( int i = 0; i < dimension; ++i ) {
        RangeType tmp = 2;
        for( int j = 0; j < dimension; ++j ) {
            if( i == j )
                continue;
            const DomainFieldType &x_j = x[ j ];
            tmp *= x_j - SQR( x_j );
        }
        phi += tmp;
    }
    #else
    phi = 8.0 * M_PI * M_PI * cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);
    #endif
    }
}; // end class myRHSFunction



// ****** Problem data: the exact solution u to the problem (for EOC calculation and boundary values)
template <class FunctionSpaceType>
class myExactSolution
    : public Function < FunctionSpaceType, myExactSolution<FunctionSpaceType> >
{
private:
    typedef myExactSolution<FunctionSpaceType> ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

public:
    inline myExactSolution ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate ( const DomainType &x , RangeType &phi ) const
    {
    enum { dimension = DomainType :: dimension };

    #if 0
    // u( x, y, z ) = \prod_{i=1}^{dimworld} (x_i - x_i^2).
    phi = 1;
    for( int i = 0; i < dimension; ++i )
    {
        const DomainFieldType &x_i = x[ i ];
        phi *= x_i - SQR( x_i );
    }
    #else
    phi = cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);
    #endif
    }

    inline void evaluate ( const DomainType &x , RangeFieldType time, RangeType &phi ) const
    {
        evaluate( x, phi );
    }
}; // end class myExactSolution




//! right hand side and exact solution for poisson problem
template <class FunctionSpaceType>
struct PoissonProblem
{
    // type of right hand side
    typedef myRHSFunction<FunctionSpaceType> RHSFunctionType;

    // type of exact solution
    typedef myExactSolution<FunctionSpaceType> ExactSolutionType;
};



//! for the L2-Projection use the exact solution also for the right hand side
template <class FunctionSpaceType>
struct L2ProjectionProblem
{
    // type of right hand side, use exact solution(!)
    typedef myExactSolution<FunctionSpaceType> RHSFunctionType;

    // type of exact solution 
    typedef myExactSolution<FunctionSpaceType> ExactSolutionType;
};

}
#endif
