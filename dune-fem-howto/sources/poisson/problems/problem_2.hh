#ifndef DUNE_PROBLEM_2_HH
#define DUNE_PROBLEM_2_HH

#include <dune/fem/function/common/function.hh>

// For more detailed descriptions see: poisson_1.hh. The following program was created similarly. Some new classes were added such as Tensor and Mass-Term treatment.

namespace Dune
{

  template< class FunctionSpaceImp >
  class RHSFunctionPeriodic
  : public Function< FunctionSpaceImp, RHSFunctionPeriodic< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    
  private:
    typedef RHSFunctionPeriodic< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  protected:
    const RangeFieldType epsilon_;
   
  public:
    inline explicit RHSFunctionPeriodic ( FunctionSpaceType &fSpace,
                                          const RangeFieldType epsilon = 0.001 )
    : BaseType( fSpace ),
      epsilon_( epsilon )
    {
    }
    
    inline void evaluate( const DomainType &x,
                          RangeType &y ) const
    {
      enum { dim = DomainType :: dimension };
      
      // By starting with 'y = 1' we would solve problem 2 with the right hand
      // side of problem 1. The solutions are almost equal if \epsilon is small
      // enough.
      y = epsilon_ * (1.0 / (4.0 * M_PI * M_PI * dim)) + 1.0;
      for( unsigned int i = 0; i < dim; ++i )
        y *= sin( 2.0 * M_PI * x[ i ] );
    }
  };



  //! The exact solution of the problem is defined. This enables us to compare (later on) the numerical solution with exact solution in order to calculate the EOC (Estimated Order of Convergence).
  template< class FunctionSpaceImp >
  class ExactSolutionPeriodic
  : public Function< FunctionSpaceImp, ExactSolutionPeriodic< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef ExactSolutionPeriodic< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef DomainFieldType TimeType;
    //essentially: 'DomainFieldType' is the type of an entry of a domain-element.
    //But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
    //This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an entry of a domain-element.
    
  public:
    inline explicit ExactSolutionPeriodic ( FunctionSpaceType &fSpace )
    : BaseType ( fSpace )
    {
    }
    
    // in case 'u' has NO time-dependency use the following method: 
    inline void evaluate ( const DomainType &x,
                           RangeType &y ) const
    {
      //  x = (x1,x2,x3), y = u(x1,x2,x3)
      enum { dim = DomainType :: dimension };
      
      y = 1.0 / (4.0 * M_PI * M_PI * dim);
      for( unsigned int i = 0; i < dim; ++i ) 
        y *= sin( 2.0 * M_PI * x[ i ] );
    }
    
    // in case 'u' HAS a time-dependency use the following method: 
    // unfortunately GRAPE requires both cases of the method 'evaluate' to be
    // instantiated
    inline void evaluate ( const DomainType &x,
                           const TimeType &time,
                           RangeType &y ) const
    {
      evaluate( x, y );
    }
  };



  template< class FunctionSpaceImp >
  class Tensor
  : public Function< FunctionSpaceImp, Tensor< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef Tensor< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef DomainFieldType TimeType;

  public:
    // Constructor for Tensor
    inline explicit Tensor ( FunctionSpaceType &fSpace )
    : BaseType( fSpace )  
    {
    }
    
    // instantiate all possible cases of the evaluate-method: 
    inline void evaluate ( const int i,
                           const int j,
                           const DomainType &x,
                           RangeType &y ) const
    {
      return evaluate( x, y );
    }
    
    inline void evaluate ( const DomainType &x,
                           RangeType &y ) const
    {
      y[ 0 ] = 1;
    }
    
    inline void evaluate ( const DomainType &x,
                           const TimeType time,
                           RangeType &y ) const
    {
      return evaluate( x, y );
    }
  };

  

  template< class FunctionSpaceImp >
  class MassTerm
  : public Function< FunctionSpaceImp, MassTerm< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef MassTerm< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  protected:
    const RangeFieldType epsilon_;

  public:
    // Constructor
    inline explicit MassTerm ( FunctionSpaceType &fSpace,
                               const RangeFieldType epsilon = 0.001 )
    : BaseType( fSpace ),
      epsilon_( epsilon )
    //epsilon itself was only valid within RHSFunctionPeriodic, therefore we needed to save its value in epsilon_()
    {
    };
    
    inline void evaluate ( const DomainType &x,
                           RangeType &y ) const
    {
      y[ 0 ] = epsilon_; 
    }
  };

}

#endif
