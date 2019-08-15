#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const Sweep *sweep, int beg, int end, 
                 real *cmax, Grid *grid)
/*
 *
 *
 **************************************************************************** */
{
  int    nv, i;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double  scrh;
  double *uL, *uR, *SR, *SL;
  double **fL = stateL->flux, **fR = stateR->flux;
  double  *pL = stateL->prs,   *pR = stateR->prs;

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

/* --------------------------------------------------
        compute HLL  flux
   -------------------------------------------------- */

    if (SL[i] >= 0.0){
    
      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fL[i][nv];
      sweep->press[i] = pL[i];

    }else if (SR[i] <= 0.0){

      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fR[i][nv];
      sweep->press[i] = pR[i];

    }else{

      uL = stateL->u[i];
      uR = stateR->u[i];

      scrh = 1.0/(SR[i] - SL[i]);
      for (nv = NFLX; nv--; ){  
        sweep->flux[i][nv]  =   SL[i]*SR[i]*(uR[nv] - uL[nv])
                              + SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        sweep->flux[i][nv] *= scrh;
      }
      sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
    }
  }
}
