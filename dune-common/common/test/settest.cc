#include "config.h"

#include<dune/common/enumset.hh>
#include<iostream>
int main()
{
  using namespace Dune;
  std::cout<<Combine<EnumItem<int,1>,EnumItem<int,2>,int>::contains(1)<<
    " "<<Combine<EnumItem<int,1>,EnumItem<int,2>,int>::contains(2)<<
    " "<<Combine<Combine<EnumItem<int,1>,EnumItem<int,2>,int>,EnumItem<int,0>,int>::contains(3)<<
    " "<<EnumRange<int,1,3>::contains(3)<<std::endl;
}
