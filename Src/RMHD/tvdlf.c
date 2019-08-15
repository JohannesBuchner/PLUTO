#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const Sweep *sweep, int beg, int end, 
                double *cmax, Grid *grid)
/*!
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
  static double *cmin_L, *cmax_L;
  static double *cmin_R, *cmax_R;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  if (cmin_L == NULL){
    cmin_R  = ARRAY_1D(NMAX_POINT, double);
    cmin_L  = ARRAY_1D(NMAX_POINT, double);
    cmax_R  = ARRAY_1D(NMAX_POINT, double);
    cmax_L  = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
        redefine sweeps if GLM_MHD is used
   ---------------------------------------------------- */

#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* -------------------------------------------------------------------
       compute Max and min eigenvalues for the left and right sweeps
   ------------------------------------------------------------------- */

  MaxSignalSpeed (stateL, cmin_L, cmax_L, beg, end);
  MaxSignalSpeed (stateR, cmin_R, cmax_R, beg, end);

  for (i = beg; i <= end; i++) {

  /* -- compute local max eigenvalue -- */

    cL = MAX(fabs(cmax_L[i]), fabs(cmin_L[i]));
    cR = MAX(fabs(cmax_R[i]), fabs(cmin_R[i]));
    cmax[i] = MAX(cL, cR);
    sweep->SL[i] = -cmax[i];
    sweep->SR[i] =  cmax[i];

  /* -- compute Rusanov flux -- */

    uL   = stateL->u[i];
    uR   = stateR->u[i];
    flux = sweep->flux[i];
    for (nv = NFLX; nv--; ) {
      flux[nv] = 0.5*(fL[i][nv] + fR[i][nv] - cmax[i]*(uR[nv] - uL[nv]));
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
 
  #if DIVB_CONTROL == EIGHT_WAVES
   POWELL_DIVB_SOURCE (sweep, beg, end, grid);
  #endif

}
