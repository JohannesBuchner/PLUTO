#include "pluto.h"

/* ********************************************************************* */
void GMUSTA1_Solver (const Sweep *sweep, int beg, int end, 
                     double *cmax, Grid *grid)
/*!
 *  Implementation of the GMUSTA1 solver as in the Appendix of
 *
 *  "MUSTA Fluxes for systems of conservation laws"
 *   Toro & Titarev, JCP (2006) 216, 403
 *
 *********************************************************************** */
{
  int      nv, i, k;
  double dtdx, dxdt, om, fG, dtG, dx, cfl;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double *a2L = stateL->a2, *a2R = stateR->a2;
  double **uL = stateL->u, **uR  = stateR->u;
  double **vL = stateL->v, **vR  = stateR->v;
  double **fL = stateL->flux, **fR  = stateR->flux;

  static State stateLW;
  static unsigned char *flag;
  static double **f_LF;
  static double *cmin_L, *cmax_L;
  static double *cmin_R, *cmax_R;
/*
  #if GEOMETRY != CARTESIAN
   print1 ("! GMUSTA1 only works in Cartesian coordinates\n");
   QUIT_PLUTO(1);
  #endif
*/
  if (f_LF == NULL){
    print (" >> GMUSTA1 Solver\n");

    f_LF = ARRAY_2D(NMAX_POINT, NVAR, double);

    cmin_L = ARRAY_1D(NMAX_POINT, double);
    cmax_L = ARRAY_1D(NMAX_POINT, double);
    cmin_R = ARRAY_1D(NMAX_POINT, double);
    cmax_R = ARRAY_1D(NMAX_POINT, double);

    flag = ARRAY_1D(NMAX_POINT, unsigned char);

    StateStructAllocate(&stateLW);
  }

/* ------------------------------------------------------
            Copy main arrays
   ------------------------------------------------------ */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

//  HLL_Speed (stateL, stateR, sweep->SL, sweep->SR, beg, end);
  for (i = beg; i <= end; i++) {
//    cmax[i] = MAX(fabs(sweep->SL[i]), fabs(sweep->SR[i]));
    cmax[i] = 1.0;
  }
      
  dx  = 1.0;
  cfl = 0.9;
  om  = 1.0/(1.0 + cfl);

/* -------------------------------------------------------
    1. Compute fluxes on the original data
   ------------------------------------------------------- */

  Flux (stateL, beg, end);  /* !! For RRMHD, pressure is incorporated */
  Flux (stateR, beg, end);  /* !! into momentum flux */
  #if RESISTIVITY == NO
  for (i = beg; i <= end; i++) {
    stateL->flux[i][MXn] += stateL->prs[i];
    stateL->prs[i] = 0.0;
    stateR->flux[i][MXn] += stateR->prs[i];
    stateR->prs[i] = 0.0;
  }
  #endif

/* ------------------------------------------------------- 
    2. Compute predictor flux using GFORCE
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    dtG  = cfl*dx/cmax[i];
    dtdx = dtG/dx;
    dxdt = 1.0/dtdx;
// dxdt = cmax[i]; dtdx = 1.0/dxdt;

    NVAR_LOOP(nv){
      f_LF[i][nv] =   0.5*(fL[i][nv] + fR[i][nv]) 
                    - 0.5*dxdt*(uR[i][nv] - uL[i][nv]);
      stateLW.u[i][nv] =   0.5*(uL[i][nv] + uR[i][nv]) 
                         - 0.5*dtdx*(fR[i][nv] - fL[i][nv]);
    }
  }

  ConsToPrim (stateLW.u, stateLW.v, beg, end, flag);
  SoundSpeed2 (&stateLW, beg, end, FACE_CENTER, NULL);
  Flux (&stateLW, beg, end);
  #if RESISTIVITY == NO
  for (i = beg; i <= end; i++) {
    stateLW.flux[i][MXn] += stateLW.prs[i];
    stateLW.prs[i] = 0.0;
  }
  #endif

/* -------------------------------------------------------
    3. Perform one-stage evolution of data on MUSTA mesh
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    dtG  = cfl*dx/cmax[i];
    dtdx = dtG/dx;
    dxdt = 1.0/dtdx;
    VAR_LOOP(nv){
      fG = om*stateLW.flux[i][nv] + (1.0 - om)*f_LF[i][nv];
      uL[i][nv] -= dtdx*(fG        - fL[i][nv]);
      uR[i][nv] -= dtdx*(fR[i][nv] - fG);
    }
  }
  ConsToPrim (uL, vL, beg, end, flag);
  ConsToPrim (uR, vR, beg, end, flag);

/* -------------------------------------------------------
    4. Re-compute fluxes fL, fR on the evolved data
   ------------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);
  #if RESISTIVITY == NO
  for (i = beg; i <= end; i++) {
    stateL->flux[i][MXn] += stateL->prs[i];
    stateL->prs[i] = 0.0;
    stateR->flux[i][MXn] += stateR->prs[i];
    stateR->prs[i] = 0.0;
  }
  #endif

/* -------------------------------------------------------
    5. Compute corrector flux using GFORCE
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    dtG  = cfl*dx/cmax[i];
    dtdx = dtG/dx;
    dxdt = 1.0/dtdx;
// dxdt = cmax[i]; dtdx = 1.0/dxdt;

    VAR_LOOP(nv){
      f_LF[i][nv] =   0.5*(fL[i][nv] + fR[i][nv]) 
                    - 0.5*dxdt*(uR[i][nv] - uL[i][nv]);
      stateLW.u[i][nv] =   0.5*(uL[i][nv] + uR[i][nv]) 
                         - 0.5*dtdx*(fR[i][nv] - fL[i][nv]);
    }
  }
  ConsToPrim (stateLW.u, stateLW.v, beg, end, flag);
  SoundSpeed2 (&stateLW, beg, end, FACE_CENTER, NULL);
  Flux (&stateLW, beg, end);
  #if RESISTIVITY == NO
  for (i = beg; i <= end; i++) {
    stateLW.flux[i][MXn] += stateLW.prs[i];
    stateLW.prs[i] = 0.0;
  }
  #endif

  for (i = beg; i <= end; i++) {
    VAR_LOOP(nv) {
      sweep->flux[i][nv] = om*stateLW.flux[i][nv] + (1.0 - om)*f_LF[i][nv];
    }
    sweep->press[i] = 0.0;
  }
}

