// $Id: solvercategory.hh 751 2007-04-20 13:59:25Z mblatt $
#ifndef DUNE_SOLVERCATEGORY_HH
#define DUNE_SOLVERCATEGORY_HH


namespace Dune {
   
  /**
     @addtogroup ISTL_Solvers
     @{
  */
  
  /**
   * @brief Categories for the solvers.
   */
  struct SolverCategory
  { 
    enum { 
      //! \brief Category for sequential solvers.
      sequential,
      //! \brief Category for on overlapping solvers.
      nonoverlapping,
      //! \brief Category for ovelapping solvers.
      overlapping
    };
  };
 
  /** @} end documentation */

} // end namespace

#endif
