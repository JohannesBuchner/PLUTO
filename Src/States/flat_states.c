#include "pluto.h"

/* ************************************************************* */
void States (const Sweep *sweep, int beg, int end, Grid *grid)
/* 
 *  PURPOSE
 *    
 *    provide 1st order flat reconstruction inside each 
 *    cell. 
 *    Here vL and vR are left and right sweeps with respect 
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
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double **v  = stateC->v;
  double **vp = stateL->v;
  double **vm = stateR->v-1;
  double **up = stateL->u;
  double **um = stateR->u-1;

#if TIME_STEPPING != EULER
  #error FLAT Reconstruction must be used with EULER integration only
#endif
  
  for (i = beg; i <= end; i++) {
    NVAR_LOOP(nv) vm[i][nv] = vp[i][nv] = v[i][nv];
  }
 
/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg; i <= end-1; i++) {
     stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
   }
  #endif

/* -------------------------------------------
    compute sweeps in conservative variables
   ------------------------------------------- */

  PrimToCons (vp, up, beg, end);
  PrimToCons (vm, um, beg, end);
}
