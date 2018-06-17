/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLL Riemann solver for MHD.

  Solve the Riemann problem for the MHD equations using the
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
  \date   April 16, 2017
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
  int    nv, i, xdface;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double scrh;
  double *uL, *uR, *SR, *SL;
  double **fL = stateL->flux, **fR = stateR->flux;
  double  *pL = stateL->prs,   *pR = stateR->prs;
  static double **Uhll;

/* --------------------------------------------------------
   0. Allocate memory / initialize arrays
   -------------------------------------------------------- */

  if (Uhll == NULL){
    Uhll = ARRAY_2D(NMAX_POINT, NFLX, double);
  }

#if BACKGROUND_FIELD == YES
  GetBackgroundField (stateL, beg, end, FACE_CENTER, grid);
#endif

/* --------------------------------------------------------
   1. Solve 2x2 Riemann problem with GLM cleaning
   -------------------------------------------------------- */

#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   2. Compute sound speed & fluxes at zone interfaces
   -------------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

#if HALL_MHD == EXPLICIT
  HallMHD_WhistlerSpeed (stateL, beg, end, grid);
  HallMHD_WhistlerSpeed (stateR, beg, end, grid);
#endif

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* --------------------------------------------------------
   3. Get max and min signal velocities
   -------------------------------------------------------- */
             
  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

/* --------------------------------------------------------
   4. Compute HLL flux
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    if (SL[i] > 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        sweep->flux[i][nv] = fL[i][nv];
      }
      sweep->press[i] = pL[i];
      
    }else if (SR[i] < 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        sweep->flux[i][nv] = fR[i][nv];
      }
      sweep->press[i] = pR[i];
      
    }else{
    
      uL   = stateL->u[i];
      uR   = stateR->u[i];
      scrh = 1.0/(SR[i] - SL[i]);
    
      for (nv = 0; nv < NFLX; nv++) {
        sweep->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        sweep->flux[i][nv] *= scrh;
      }
      sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
    }
  }

/* --------------------------------------------------------
   5. Compute source terms (if any)
   -------------------------------------------------------- */

#if DIVB_CONTROL == EIGHT_WAVES
  HLL_DivBSource (sweep, Uhll, beg + 1, end, grid);
#endif

/* ----------------------------------------------------------
   6. Add CR flux contribution using simplified upwinding.
   ---------------------------------------------------------- */

#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux (stateL, beg, end);
  Particles_CR_Flux (stateR, beg, end);

  for (i = beg; i <= end; i++){
    if (sweep->flux[i][RHO] > 0.0) {
      for (nv = NFLX; nv--; ) sweep->flux[i][nv] += stateL->fluxCR[i][nv];
    }else if (sweep->flux[i][RHO] < 0.0){
      for (nv = NFLX; nv--; ) sweep->flux[i][nv] += stateR->fluxCR[i][nv];
    }else{
      for (nv = NFLX; nv--; ) {
        sweep->flux[i][nv] += 0.5*(stateL->fluxCR[i][nv] + stateR->fluxCR[i][nv]);
      }
    }
  }  
  #endif
#endif

}
