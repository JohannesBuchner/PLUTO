#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const State_1D *state, int beg, int end, 
                 real *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   HLL_SOLVER
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the Euler equations 
 *   using th HLLE solver;
 * 
 *   Reference:    --
 *   
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 **************************************************************************** */
{
  int    nv, i;
  real   scrh;
  double *uL, *uR;
  static real **fL, **fR;
  static real *SL, *SR, *pL, *pR;
  static double *a2L, *a2R, *hL, *hR;

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);
    SR = ARRAY_1D(NMAX_POINT, double);
    SL = ARRAY_1D(NMAX_POINT, double);

    hR = ARRAY_1D(NMAX_POINT, double);
    hL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
  }
/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

  HLL_Speed (state->vL, state->vR, a2L, a2R, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

/* --------------------------------------------------
        compute HLL  flux
   -------------------------------------------------- */

    if (SL[i] >= 0.0){
    
      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = pL[i];

    }else if (SR[i] <= 0.0){

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = pR[i];

    }else{

      uL = state->uL[i];
      uR = state->uR[i];

      scrh = 1.0/(SR[i] - SL[i]);
      for (nv = NFLX; nv--; ){  
        state->flux[i][nv]  =   SL[i]*SR[i]*(uR[nv] - uL[nv])
                              + SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
    }
  }
}
