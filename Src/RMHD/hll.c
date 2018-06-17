/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLL Riemann solver for RMHD.

  Solve the Riemann problem for the RMHD equations using the
  single-sweep HLL solver by Toro.

  On input, this function takes left and right primitive sweep vectors
  \c sweep->vL and \c sweep->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface
  \c i+1/2 (note that the \c i refers to \c i+1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
   - "Riemann Solver and Numerical Methods for Fluid Dynamics"
      by E.F. Toro (Chapter 10)
   - "On Godunov-Type Method near Low Densities"
     by B. Einfeldt, C.D. Munz, P.L. Roe, JCP 92, 273-295 (1991)

  \authors A. Mignone (mignone@ph.unito.it)
  \date    April 7, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const Sweep *sweep, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic/isothermal MHD equations 
 * using the HLL Riemann solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;
  real   scrh, bmin, bmax;
  double *uL, *uR, *SL, *SR;
  static double **Uhll;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *hL = stateL->h,     *hR = stateR->h;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  if (Uhll == NULL){
    Uhll = ARRAY_2D(NMAX_POINT, NVAR, double); /* NVAR and not NFLX since it needs  */
                                               /* to be converted in HLL_DIVB_SOUCE */
  }

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

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */

  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    uL = stateL->u[i];
    uR = stateR->u[i];

/* ---- Compute fluxes ---- */
    
    bmin = MIN(0.0, SL[i]);
    bmax = MAX(0.0, SR[i]);
    scrh = 1.0/(bmax - bmin);
    for (nv = NFLX; nv--; ){
      sweep->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                         +  bmax*fL[i][nv] - bmin*fR[i][nv];
      sweep->flux[i][nv] *= scrh;
    }
    sweep->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh;
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
 
  #if DIVB_CONTROL == EIGHT_WAVES
/*
   POWELL_DIVB_SOURCE (sweep, beg, end, grid);
*/

  /* ----------------------------------------------------
       to avoid conversion problems in HLL_DIVB_SOURCE, 
       we use the HLL average provided by SR = -SL = 1 
     ---------------------------------------------------- */
/*
   for (i = beg; i <= end; i++) {
     uL = sweep->uL[i]; uR = sweep->uR[i];
     for (nv = 0; nv < NFLX; nv++) {
       Uhll[i][nv] = 0.5*(uR[nv] + uL[nv] + fL[i][nv] - fR[i][nv]);
     }
     Uhll[i][MXn] += (pL[i] - pR[i])*0.5;
     for (nv = NFLX; nv < NVAR; nv++) Uhll[i][nv] = 0.0;
   }
*/
   HLL_DIVB_SOURCE (sweep, Uhll, beg + 1, end, grid);
  #endif

}
