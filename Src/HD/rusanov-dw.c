#include "pluto.h"

/* ********************************************************************* */
void RusanovDW_Solver (const Sweep *sweep, int beg, int end, 
                double *cmax, Grid *grid)
/*!
 * Solve Riemann problem using the Lax-Friedrichs Rusanov solver:
 *
 *   F(i+1/2) = [FL + FR - c*(UR - UL)]*0.5
 *
 * where c = max_speed[(UR + UL)*0.5]
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *********************************************************************** */
{
  int    nv, i;
  static double **fL, **fR, **vRL;
  static double *cRL_min, *cRL_max, *pL, *pR, *a2L, *a2R;
  double *uR, *uL, *flux;

double num, den, lambda, dU, dF, vn;
  
  if (fR == NULL){
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    vRL = ARRAY_2D(NMAX_POINT, NFLX, double);
    cRL_min = ARRAY_1D(NMAX_POINT, double);
    cRL_max = ARRAY_1D(NMAX_POINT, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    pL  = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (sweep->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (sweep->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (sweep->uL, sweep->vL, a2L, fL, pL, beg, end);
  Flux (sweep->uR, sweep->vR, a2R, fR, pR, beg, end);

/* ---------------------------------------------------------------------
    Compute average sweep in order to get the local max characteristic
    velocity.
    In order to avoid underestimating this speed at strong reflecting
    boundaries (or when vxL = -vxR in general) , we average the
    absolute value of the normal velocity.
   ------------------------------------------------------------------- */

  for (i = beg; i <= end; i++)           {
    for (nv = NFLX; nv--;  ) {
      vRL[i][nv] = 0.5*(sweep->vL[i][nv] + sweep->vR[i][nv]);
    }
    vRL[i][VXn] = 0.5*(fabs(sweep->vL[i][VXn]) + fabs(sweep->vR[i][VXn]));
  }

/* ---------------------------------------------------------------------
         Compute sound speed, max and min signal speeds
   --------------------------------------------------------------------- */

  SoundSpeed2    (vRL, a2R, NULL, beg, end, FACE_CENTER, grid);
  MaxSignalSpeed (vRL, a2R, cRL_min, cRL_max, beg, end);

  for (i = beg; i <= end; i++) {
    uR   = sweep->uR[i];
    uL   = sweep->uL[i];
    flux = sweep->flux[i];
  
  /* -- compute max eigenvalue -- */

    cmax[i] = MAX(fabs(cRL_max[i]), fabs(cRL_min[i]));
    g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/sqrt(a2R[i]));

  /* -- compute fluxes -- */
num = den = 0.0;
for (nv = 0; nv < NFLX; nv++) {
  dU = uR[nv] - uL[nv];
  dF = fR[i][nv] - fL[i][nv] + (nv == MXn ? (pR[i] - pL[i]):0.0);

  num += dU*dF;
  den += dU*dU;
}
lambda = fabs(num)/(fabs(den) + 1.e-12);

if (lambda > cmax[i])  lambda = cmax[i];

vn = 0.5*( fabs(sweep->vR[i][VXn]) + fabs(sweep->vL[i][VXn]) );
vn = MAX( fabs(sweep->vR[i][VXn]), fabs(sweep->vL[i][VXn]) );
if (lambda < vn) lambda = vn;

/*
cL = MIN(fabs(cmax_L[i]), fabs(cmin_L[i]));
cR = MIN(fabs(cmax_R[i]), fabs(cmin_R[i]));
*/

    for (nv = NFLX; nv--;  ) {
      flux[nv] = 0.5*(fL[i][nv] + fR[i][nv] - lambda*(uR[nv] - uL[nv]));
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
