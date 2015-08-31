#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const State_1D *state, int beg, int end, 
                double *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   LF_SOLVER
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the relativistic MHD equations 
 *   using th Lax-Friedrichs Rusanov solver. 
 * 
 * LAST_MODIFIED
 *
 *   Feb 15/2010 by Andrea Mignone  (mignone@ph.unito.it)
 *              
 *
 **************************************************************************** */
{
  int    nv, i;
  double *uL, *uR, *flux, cL, cR;
  static double **fL, *pL, *a2L, *hL, *cmin_L, *cmax_L;
  static double **fR, *pR, *a2R, *hR, *cmin_R, *cmax_R;
  static double **VL, **VR, **UL, **UR;

  if (fR == NULL){
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    pR      = ARRAY_1D(NMAX_POINT, double);
    pL      = ARRAY_1D(NMAX_POINT, double);
    cmin_R  = ARRAY_1D(NMAX_POINT, double);
    cmin_L  = ARRAY_1D(NMAX_POINT, double);
    cmax_R  = ARRAY_1D(NMAX_POINT, double);
    cmax_L  = ARRAY_1D(NMAX_POINT, double);
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif

    a2R  = ARRAY_1D(NMAX_POINT, double);
    a2L  = ARRAY_1D(NMAX_POINT, double);
    hR   = ARRAY_1D(NMAX_POINT, double);
    hL   = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
        redefine states if GLM_MHD is used
   ---------------------------------------------------- */

  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, hL, fL, pL, beg, end);
  Flux (UR, VR, hR, fR, pR, beg, end);

/* -------------------------------------------------------------------
       compute Max and min eigenvalues for the left and right states
   ------------------------------------------------------------------- */

  MaxSignalSpeed (VL, a2L, hL, cmin_L, cmax_L, beg, end);
  MaxSignalSpeed (VR, a2R, hR, cmin_R, cmax_R, beg, end);

  for (i = beg; i <= end; i++) {

  /* -- compute local max eigenvalue -- */

    cL = MAX(fabs(cmax_L[i]), fabs(cmin_L[i]));
    cR = MAX(fabs(cmax_R[i]), fabs(cmin_R[i]));
    cmax[i] = MAX(cL, cR);
    state->SL[i] = -cmax[i];
    state->SR[i] =  cmax[i];

  /* -- compute Rusanov flux -- */

    uL   = UL[i];
    uR   = UR[i];
    flux = state->flux[i];
    for (nv = NFLX; nv--; ) {
      flux[nv] = 0.5*(fL[i][nv] + fR[i][nv] - cmax[i]*(uR[nv] - uL[nv]));
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
 
  #if DIVB_CONTROL == EIGHT_WAVES
   POWELL_DIVB_SOURCE (state, beg, end, grid);
  #endif

}
