/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLL Riemann solver for HD.

  Solve the Riemann problem for the HD equations using the
  single-state HLL solver by Toro.

  On input, this function takes left and right primitive state vectors
  \c stateL->v and \c stateR->v at zone edge \c i+1/2;
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
  \date    Oct 12, 2016
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

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double scrh;
  double **fL = stateL->flux, **fR = stateR->flux;
  double  *pL = stateL->prs,   *pR = stateR->prs;
  double *uR, *uL, *SL, *SR;

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

    scrh = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i]  = scrh;

    if (SL[i] > 0.0){
    
      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fL[i][nv];
      sweep->press[i] = pL[i];

    }else if (SR[i] < 0.0){

      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fR[i][nv];
      sweep->press[i] = pR[i];

    }else{

      uR = stateR->u[i];
      uL = stateL->u[i];

      scrh = 1.0 / (SR[i] - SL[i]);
      for (nv = NFLX; nv--;  ) {
        sweep->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        sweep->flux[i][nv] *= scrh;
      }
      sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;

    }
  } /* end loops on points */
}
