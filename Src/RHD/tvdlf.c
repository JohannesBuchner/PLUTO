#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const State_1D *state, int beg, int end, 
                real *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   LF_SOLVER
 *
 *
 * PURPOSE
 *
 *   - Solve Riemann problem using the Lax-Friedrichs 
 *     Rusanov solver:
 *
 *     F(i+1/2) = [FL + FR - c*(UR - UL)]*0.5
 *
 *     where c = max_speed[(UR + UL)*0.5]
 *
 *   - On input, it takes left and right primitive state
 *     vectors state->vL and state->vR at zone edge i+1/2;
 *     On output, return flux and pressure vectors at the
 *     same interface.
 *
 *   - Also, compute maximum wave propagation speed (cmax) 
 *     for  explicit time step computation
 *  
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
  double *uL, *uR, cL, cR;
  static double **fL, **fR;
  static double *pR, *pL, *a2L, *a2R, *hL, *hR;
  static double *cminL, *cmaxL;
  static double *cminR, *cmaxR;

  if (fR == NULL){
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    pL  = ARRAY_1D(NMAX_POINT, double);
    cminL = ARRAY_1D(NMAX_POINT, double);
    cmaxL = ARRAY_1D(NMAX_POINT, double);
    cminR = ARRAY_1D(NMAX_POINT, double);
    cmaxR = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    hR  = ARRAY_1D(NMAX_POINT, double);
    hL  = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

  MaxSignalSpeed (state->vL, a2L, cminL, cmaxL, beg, end);
  MaxSignalSpeed (state->vR, a2R, cminR, cmaxR, beg, end);

  for (i = beg; i <= end; i++) {
    cL = MAX(fabs(cminL[i]), fabs(cmaxL[i]));
    cR = MAX(fabs(cminR[i]), fabs(cmaxR[i]));
    cmax[i] = MAX(cL, cR);

    uL = state->uL[i];
    uR = state->uR[i];
    for (nv = NFLX; nv--;   ) {
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - cmax[i]*(uR[nv] - uL[nv]));
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
