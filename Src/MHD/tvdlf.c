/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for MHD.

  Solve the Riemann problem for the adiabatic and isothermal MHD 
  equations using the Lax-Friedrichs Rusanov Riemann solver with local
  maximum characteristic speed:
  \f[
      \hat{F}_{i+\HALF} = \frac{F^L_{i+\HALF} + F^R_{i+\HALF}}{2}
   - c^{\rm max}_{i+\HALF}\frac{U^R_{i+\HALF} - U^L_{i-\HALF}}{2}
     \qquad{\rm where}\quad c^{\rm max} = \lambda_{\rm max}\left(
          \frac{V^R_{i+\HALF} + V^L_{i+\HALF}}{2}\right)
  \f]                          
  where \f$\lambda_{\rm \max}\f$ is a function of the arithmetic 
  average between the left and the right states.
  
  On input, this function takes left and right primitive state vectors 
  \c stateL->v and \c stateR->v at zone edge i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "Comparison of Some Flux Correced Transport and Total Variation 
       Diminishing Numerical Schemes for Hydrodynamics and Magnetohydrodynamics 
       Problems", Toth and Odstrcil, JCP (1996), 128,82
       
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 20, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const Sweep *sweep, int beg, int end, 
                double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic MHD equations using the 
 * Lax-Friedrichs (Rusanov) Riemann solver.
 * 
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int     nv, i;

  State        stateRL;  
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double  cRL, *uR, *uL, **bgf;
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;
  static double **vRL,  *cmin_RL, *cmax_RL;

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (vRL == NULL){
    vRL = ARRAY_2D(NMAX_POINT, NFLX, double);
    cmin_RL = ARRAY_1D(NMAX_POINT, double);
    cmax_RL = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
   1. Do some preliminary operations...
   -------------------------------------------------------- */
  
#if BACKGROUND_FIELD == YES
  GetBackgroundField (stateL, beg, end, FACE_CENTER, grid);
#endif

#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   2. Compute sound speed & fluxes at zone interfaces
   -------------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* --------------------------------------------------------
   3. Compute max eigenvalue and fluxes
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    NVAR_LOOP(nv) vRL[i][nv] = 0.5*(stateL->v[i][nv] + stateR->v[i][nv]);
    #if EOS == IDEAL
     g_maxMach = MAX(g_maxMach, 
                     fabs(vRL[i][VXn])/sqrt(g_gamma*vRL[i][PRS]/vRL[i][RHO])); 
    #elif EOS == BAROTROPIC
     cRL = sqrt(g_gamma*BAROTROPIC_PR(vRL[i][RHO])/vRL[i][RHO]);
     g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/cRL); 
    #elif EOS == ISOTHERMAL
     g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/g_isoSoundSpeed); 
    #endif
  }

  stateRL.v    = vRL;
  stateRL.a2   = a2R;   /* Reuse a2R for convenience */
  stateRL.Bbck = stateL->Bbck;
  SoundSpeed2    (&stateRL, beg, end, FACE_CENTER, grid);
  MaxSignalSpeed (&stateRL, cmin_RL, cmax_RL, beg, end);

  for (i = beg; i <= end; i++) {
    cRL = MAX(fabs(cmin_RL[i]), fabs(cmax_RL[i]));
    sweep->SL[i] = -cRL;
    sweep->SR[i] =  cRL;
    cmax[i] = cRL;
    uL = stateL->u[i];
    uR = stateR->u[i];
    for (nv = NFLX; nv--;    ) {
      sweep->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - cRL*(uR[nv] - uL[nv]));
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }

/* --------------------------------------------------------
   4. Compute source term
   -------------------------------------------------------- */

#if DIVB_CONTROL == EIGHT_WAVES
  Roe_DivBSource (sweep, beg + 1, end, grid);
#endif

/* ----------------------------------------------------------
   5. Add CR flux contribution using simplified upwinding.
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
