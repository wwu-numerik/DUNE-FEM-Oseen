#ifndef DUNE_BINARYFUNCTIONS_HH
#define DUNE_BINARYFUNCTIONS_HH

#include<functional>
#include<algorithm>

namespace Dune
{
  template<typename Type>
  struct Min
    : std::binary_function<Type,Type,Type>
  {
    Type operator()(const Type& t1, const Type& t2) const
    {
      return std::min(t1,t2);
    }
  };
  
  template<typename Type>
  struct Max
    : std::binary_function<Type,Type,Type>
  {
    Type operator()(const Type& t1, const Type& t2) const
    {
      return std::max(t1,t2);
    }
  };
}

#endif
