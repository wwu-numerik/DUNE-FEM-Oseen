/** \file traits.hh
    \brief  brief
 **/

#ifndef TRAITS_HH
#define TRAITS_HH

#include "dune/common/fvector.hh"

struct Traits
{
    typedef Dune::FieldVector< double, 2 > TwoDVector;
    typedef Dune::FieldVector< TwoDVector, 2 > TwoTimesTwoGradient;

};

#endif  // end of TRAITS_HH
