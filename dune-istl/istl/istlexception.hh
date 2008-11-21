#ifndef DUNE_ISTLEXC_HH
#define DUNE_ISTLEXC_HH

#include <stdlib.h>

#include <dune/common/exceptions.hh>

namespace Dune {
   
    /** 
		@addtogroup ISTL
		@{
     */

  //! derive error class from the base class in common
  class ISTLError : public Exception {};

 
  /** @} end documentation */

} // end namespace

#endif
