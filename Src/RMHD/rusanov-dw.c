#include "pluto.h"

/* ********************************************************************* */
void RusanovDW_Solver (const Sweep *sweep, int beg, int end, 
                       double *cmax, Grid *grid)
/*!
 * Solve the Riemann problem for the relativistic MHD equations 
 * using the dominant-wave Riemann solver. 
 * 
 * LAST_MODIFIED
 *
 *   Feb 04/2013 by Andrea Mignone  (mignone@ph.unito.it)
 *              
 **************************************************************************** */
{
  int    nv, i;
  double *uL, *uR, *flux, cL, cR, lambda;
  static double **fL, *pL, *a2L, *hL, *cmin_L, *cmax_L;
  static double **fR, *pR, *a2R, *hR, *cmin_R, *cmax_R;
  static double **VL, **VR, **UL, **UR;
double num, den, dU, dF, vn;

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
        redefine sweeps if GLM_MHD is used
   ---------------------------------------------------- */

  #ifdef GLM_MHD
   GLM_Solve (sweep, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = sweep->vL; UL = sweep->uL;
   VR = sweep->vR; UR = sweep->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, hL, fL, pL, beg, end);
  Flux (UR, VR, hR, fR, pR, beg, end);

/* -------------------------------------------------------------------
       compute Max and min eigenvalues for the left and right sweeps
   ------------------------------------------------------------------- */

  MaxSignalSpeed (VL, a2L, hL, cmin_L, cmax_L, beg, end);
  MaxSignalSpeed (VR, a2R, hR, cmin_R, cmax_R, beg, end);

  for (i = beg; i <= end; i++) {

    uL   = UL[i];
    uR   = UR[i];
    flux = sweep->flux[i];

  /* -- compute local max eigenvalue -- */

    cL = MAX(fabs(cmax_L[i]), fabs(cmin_L[i]));
    cR = MAX(fabs(cmax_R[i]), fabs(cmin_R[i]));
    cmax[i] = MAX(cL, cR);
    sweep->SL[i] = -cmax[i];
    sweep->SR[i] =  cmax[i];

num = den = 0.0;
lambda = 1.0;
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

  /* -- compute Rusanov flux -- */

    for (nv = NFLX; nv--; ) {
      flux[nv] = 0.5*(fL[i][nv] + fR[i][nv] - lambda*(uR[nv] - uL[nv]));
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
