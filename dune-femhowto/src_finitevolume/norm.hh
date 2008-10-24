#ifndef Norm_HH 
#define Norm_HH





struct Norm
{
  template <class DiscreteFunctionImp> 
  static double L1Norm(DiscreteFunctionImp& errFunc)
  {
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    typedef typename FunctionSpaceType::RangeType RangeType;
    
    double sum = 0.0;

    const typename DiscreteFunctionType::FunctionSpaceType
      & space = errFunc.space();

    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty
    assert( it != endit );
    for(; it != endit ; ++it)
    {
      //CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = errFunc.localFunction(*it);
      sum += fabs(lf[0])* (*it).geometry().volume();

      /*******************************************************************
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate(quad[qP],phi);
        sum += det * quad.weight(qP) * SQR(ret[0] - phi[0]);
      }
      ********************************************************************/
    
      

    }
    //return sqrt(sum);    
    return sum;
  }


};  //end of struct


#endif
