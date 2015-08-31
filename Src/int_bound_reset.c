/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Reset to zero the right hand side in the internal boundary.

  Set to zero the right hand side of the conservative equations
  in those cells that have been flagged with ::FLAG_INTERNAL_BOUNDARY.

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 22, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if INTERNAL_BOUNDARY == YES
/* *********************************************************************** */
void InternalBoundaryReset (const State_1D *state, Time_Step *Dts, 
                            int beg, int end, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      grid  pointer to Grid structure
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int i,nv;
  
  for (i = beg; i <= end; i++){
    if (state->flag[i] & FLAG_INTERNAL_BOUNDARY){
      NVAR_LOOP(nv) state->rhs[i][nv] = 0.0;
/*      Dts->cmax[i] = 0.0;   */
    }
  }
}
#endif
