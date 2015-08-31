#include "pluto.h"

/* ************************************************************* */
void States (const State_1D *state, int beg, int end, Grid *grid)
/* 
 *  PURPOSE
 *    
 *    provide 1st order flat reconstruction inside each 
 *    cell. 
 *    Here vL and vR are left and right states with respect 
 *    to the cell interface, while vm and vp refer to the cell 
 *    center, that is:
 *
 *                    VL-> <-VR
 *      |--------*--------|--------*--------|
 *       <-am   (i)   ap->       (i+1)
 *    
 *
 **************************************************************** */
{
  int nv, i;

  #if TIME_STEPPING != EULER
   #error FLAT Reconstruction must be used with EULER integration only
  #endif
  
  for (i = beg; i <= end; i++) {
  for (nv = 0; nv < NVAR; nv++) {
    state->vm[i][nv] = state->vp[i][nv] = state->v[i][nv];
  }}
  PrimToCons(state->vm, state->um, beg, end);
  PrimToCons(state->vp, state->up, beg, end);
  
/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg; i <= end-1; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
