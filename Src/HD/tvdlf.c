/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for HD. 

  Solve the Riemann problem for the adiabatic and isothermal HD 
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
  \date    Feb 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const Sweep *sweep, int beg, int end, 
                double *cmax, Grid *grid)
/*!
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

  State stateRL;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double scrh;
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  static double **vRL;
  static double *cRL_min, *cRL_max;
  double *vR, *vL, *uR, *uL;
  
  if (vRL == NULL){
    vRL = ARRAY_2D(NMAX_POINT, NFLX, double);
    cRL_min = ARRAY_1D(NMAX_POINT, double);
    cRL_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* ---------------------------------------------------------------------
    Compute average state in order to get the local max characteristic
    velocity.
    In order to avoid underestimating this speed at strong reflecting
    boundaries (or when vxL = -vxR in general) , we average the
    absolute value of the normal velocity.
   ------------------------------------------------------------------- */

  for (i = beg; i <= end; i++){
    VAR_LOOP(nv) vRL[i][nv] = 0.5*(stateL->v[i][nv] + stateR->v[i][nv]);
    vRL[i][VXn] = 0.5*(fabs(stateL->v[i][VXn]) + fabs(stateR->v[i][VXn]));
  }

/* ---------------------------------------------------------------------
         Compute sound speed, max and min signal speeds
   --------------------------------------------------------------------- */

  stateRL.v  = vRL;
  stateRL.a2 = a2R;

  SoundSpeed2    (&stateRL, beg, end, FACE_CENTER, grid);
  MaxSignalSpeed (&stateRL, cRL_min, cRL_max, beg, end);

  for (i = beg; i <= end; i++) {
    uR   = stateR->u[i];
    uL   = stateL->u[i];
  
  /* -- compute max eigenvalue -- */

    cmax[i] = MAX(fabs(cRL_max[i]), fabs(cRL_min[i]));
    g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/sqrt(a2R[i]));

  /* -- compute fluxes -- */
  
    for (nv = NFLX; nv--;  ) {
      sweep->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - cmax[i]*(uR[nv] - uL[nv]));
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
