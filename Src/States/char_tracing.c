/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute time-centered interface states using characteristic 
         tracing.

  Advance 1-D left and right interface states, previously computed
  with any of the States functions, to the half time level
  (n+1/2) by extrapolating characteristic variables.
  This is done using an upwind selection rule that discards waves
  not reaching the interface in dt/2.
  
  \b References:
     - "The PLUTO Code for Adaptive Mesh Refinement in Astrophysical
        Fluid Dynamics"\n
        Mignone et al, ApJS (2012) 198, 7M
        
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 22, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if RECONSTRUCTION == PARABOLIC || RECONSTRUCTION == WENO3 || RECONSTRUCTION == LimO3

/* -------------------------------------------------------- */
/*! Flag to control the choice of the reference state in 
    the upwind selection rule.
    Set CHTR_REF_STATE to 1,2,3 to use 
    - cell centered value (1),
    - interpolated states (2),
    - fastest wave (3, the usual PPM rescription).          */
/*  ------------------------------------------------------- */

#if CHAR_LIMITING == YES   
 #ifndef CHTR_REF_STATE
  #define CHTR_REF_STATE  2
 #endif
#else
 #ifndef CHTR_REF_STATE
  #define CHTR_REF_STATE  3
 #endif
#endif

/* ********************************************************************* */
void CharTracingStep(const Sweep *sweep, int beg, int end, Grid *grid)
/*! 
 *  Compute interface states using characteristic tracing step.
 *
 * \param [in]  sweep  pointer to a Sweep structure
 * \param [in]    beg  initial index of computation
 * \param [in]    end  final   index of computation
 * \param [in]   grid  pointer to Grid structure

 * \return This function has no return value.
 *
 * Tasks are numbered below.
 *
 *********************************************************************** */
{
  int    i, nv, k;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double dx, dtdx;
  double dwh, d2w, dwp[NVAR], dwm[NVAR], dvp[NVAR], dvm[NVAR];
  double nup=1.0, num=-1.0, nu[NFLX];
  double Spp, Smm;
  double *vc, *vp, *vm, **L, **R, *lambda;
  static double **src;

/* -------------------------------------------------------
   0. Allocate memory
   ------------------------------------------------------- */

  if (src == NULL){
    src = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* -------------------------------------------------------
   1. Compute preliminary quantities such as sound speed,
      eigenvectors, etc...
      CR Flux is computed at time level n and will be
      added at the end of the function.
   ------------------------------------------------------- */

#if CHAR_LIMITING == NO
  SoundSpeed2 (stateC, beg, end, CELL_CENTER, grid);
  PrimEigenvectors (stateC, beg, end);
#endif
  PrimSource  (stateC, src, beg, end, grid);
#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux(stateL, beg, end);
  Particles_CR_Flux(stateR, beg-1, end-1);
  #endif
#endif

/* -------------------------------------------------------
   2. Sweep along 1D row of zones
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    dx   = grid->dx[g_dir][beg];
    dtdx = g_dt/dx;

    vc     = stateC->v[i]; 
    vp     = stateL->v[i];
    vm     = stateR->v[i-1];
    L      = stateC->Lp[i];
    R      = stateC->Rp[i];
    lambda = stateC->lambda[i];

  /* ------------------------------------------------------------ */
  /*! 2a) Define characteristic cfl coefficients.                 */
  /* ------------------------------------------------------------ */

    for (k = 0; k < NFLX; k++) nu[k] = dtdx*lambda[k];

  /* ------------------------------------------------------------- */
  /*! 2b) Obtain characteristic variable increments dwp and dwm.    */
  /* ------------------------------------------------------------- */

    for (nv = NVAR; nv--;  ) {
      dvp[nv] = vp[nv] - vc[nv];
      dvm[nv] = vm[nv] - vc[nv];
    }     
    PrimToChar(L, dvp, dwp);
    PrimToChar(L, dvm, dwm);

  /* ------------------------------------------------------------- */
  /*! 2c) Initialize vp and vm to the reference state.
          Since this is somewhat arbitrary we use the value of 
          ::CHTR_REF_STATE to select one of the following cases:
         
          - CHTR_REF_STATE==1: use cell-center value;
            
          - CHTR_REF_STATE==2: interpolated value at base 
                                    time level 
          - CHTR_REF_STATE==3: 
            traditional PPM reference state (fastest wave), minimize 
            the size of the term subject to characteristic limiting. 

         Passive scalars use always CHTR_REF_STATE == 2.      */
  /* ------------------------------------------------------------- */

    #if CHTR_REF_STATE == 1
     for (nv = NFLX; nv--;   ) vp[nv] = vm[nv] = vc[nv];
    #elif CHTR_REF_STATE == 3
     nup = MAX(nu[1], 0.0); num = MIN(nu[0], 0.0);
     for (nv = NFLX; nv--;   ){
       dwh = vp[nv] - vm[nv];
       d2w = vp[nv] + vm[nv] - 2.0*vc[nv];
       vp[nv] -= 0.5*nup*(dwh + d2w*(3.0 - 2.0*nup));
       vm[nv] -= 0.5*num*(dwh - d2w*(3.0 + 2.0*num));
     }
    #endif

  /* ------------------------------------------------------------- */
  /*! 2d) Compute left and right states in primitive variables. 
          This step also depends on the value of 
          ::CHTR_REF_STATE and include:

          - evolve characteristic variable increments by dt/2;
          - discard contributions from waves not reaching the 
            interface;
          - project characteristic differences dwp and dwm onto 
            right eigenvectors                                      */
  /* ------------------------------------------------------------- */

    for (k = 0; k < NFLX; k++){ 
      dwh = dwp[k] - dwm[k];
      d2w = dwp[k] + dwm[k]; 
      if (nu[k] >= 0.0) {
        #if CHTR_REF_STATE == 1
         Spp = dwp[k] - 0.5*nu[k]*(dwh + d2w*(3.0 - 2.0*nu[k]));
        #elif CHTR_REF_STATE == 2
         Spp =        - 0.5*nu[k]*(dwh + d2w*(3.0 - 2.0*nu[k]));
        #elif CHTR_REF_STATE == 3
         Spp = -0.5*(nu[k]-nup)*(dwh + 3.0*d2w) + (nu[k]*nu[k] - nup*nup)*d2w;
        #endif
        for (nv = NFLX; nv--;   ) vp[nv] += Spp*R[nv][k];
      } else {
        #if CHTR_REF_STATE == 1
         Smm = dwm[k] - 0.5*nu[k]*(dwh - d2w*(3.0 + 2.0*nu[k]));
        #elif CHTR_REF_STATE == 2
         Smm =        - 0.5*nu[k]*(dwh - d2w*(3.0 + 2.0*nu[k]));
        #elif CHTR_REF_STATE == 3
         Smm = -0.5*(nu[k]-num)*(dwh - 3.0*d2w) + (nu[k]*nu[k] - num*num)*d2w;
        #endif
        for (nv = NFLX; nv--;   ) vm[nv] += Smm*R[nv][k];
      }
    }

  /* ------------------------------------------------------- */
  /*! 2e) Add source term to L/R states                      */
  /* ------------------------------------------------------- */

    for (nv = NFLX; nv--;   ){   
      dwh = 0.5*g_dt*src[i][nv];
      vp[nv] += dwh;
      vm[nv] += dwh;
    }

  /* ------------------------------------------------------- */
  /*! 2f) Repeat construction for passive scalars            */
  /* ------------------------------------------------------- */

    #if NVAR != NFLX
     nu[0] = dtdx*vc[VXn];
     if (nu[0] >= 0.0) {   /* -- scalars all move at the flow speed -- */
       for (k = NFLX; k < NVAR; k++){
         dwh = dwp[k] - dwm[k];
         d2w = dwp[k] + dwm[k]; 
         vp[k] -= 0.5*nu[0]*(dwh + d2w*(3.0 - 2.0*nu[0]));
       }
     }else{
       for (k = NFLX; k < NVAR; k++){
         dwh = dwp[k] - dwm[k];
         d2w = dwp[k] + dwm[k]; 
         vm[k] -= 0.5*nu[0]*(dwh - d2w*(3.0 + 2.0*nu[0]));
       }
     }
    #endif
  } /* end loop on points */

/*  -------------------------------------------------------
    3. Assign face-centered magnetic field
    -------------------------------------------------------  */

#ifdef STAGGERED_MHD
  for (i = beg-1; i <= end; i++) {
    stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
  }
#endif

/* --------------------------------------------------------
   4. Add CR source term 
   -------------------------------------------------------- */

#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
  Particles_CR_StatesSource(sweep, 0.5*g_dt, beg, end, grid);
  #endif
#endif

  CheckPrimStates (stateL->v, stateR->v-1, stateC->v, beg, end);

/* --------------------------------------------------------
   5. Evolve cell-center values by dt/2
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    vp = stateL->v[i];
    vm = stateR->v[i-1];
    NVAR_LOOP(nv) stateC->v[i][nv] = 0.5*(vp[nv] + vm[nv]);
  }
#if RECONSTRUCT_4VEL
  ConvertTo3vel (stateC->v, beg, end);
#endif

}
#undef CHTR_REF_STATE

#else /* if RECONSTRUCTION == LINEAR */

/* ------------------------------------------------------------
    Set REF_STATE to 1,2,3 to use cell centered value (1), 
    interpolated states (2) or fastest wave (3, the standard 
    PPM prescription). Default is 3.
   ------------------------------------------------------------ */

#ifndef CHTR_REF_STATE
 #define CHTR_REF_STATE     3
#endif

/* ********************************************************************* */
void CharTracingStep(const Sweep *sweep, int beg, int end, Grid *grid)
/*
 * Same thing as before, but optimized for linear reconstruction.
 * 
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dx, dtdx, scrh;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double dv[NVAR], dw[NVAR];
  double nu_max=1.0, nu_min=-1.0, nu[NFLX];
  double *vc, *vp, *vm, **L, **R, *lambda;
  double dp, dm;
  double *a2 = stateC->a2;
  PLM_Coeffs plm_coeffs;
  #if GEOMETRY != CARTESIAN
   double betaL[NVAR], betaR[NVAR];
  #endif
  static double **src;

/* --------------------------------------------------------
   0. Allocate memory & obtain interpolation coefficients
   -------------------------------------------------------- */

  if (src == NULL) src = ARRAY_2D(NMAX_POINT, NVAR, double);

#if UNIFORM_CARTESIAN_GRID == NO
  PLM_CoefficientsGet (&plm_coeffs, g_dir);
#endif

/* --------------------------------------------------------
   1. Compute preliminary quantities such as sound speed,
      eigenvectors, etc...
      CR Flux is computed at time level n and will be
      added at the end of the function.
   -------------------------------------------------------- */


#if CHAR_LIMITING == NO
  SoundSpeed2 (stateC, beg, end, CELL_CENTER, grid);
  PrimEigenvectors(stateC, beg, end);
#endif
  PrimSource  (stateC, src, beg, end, grid);
#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux(stateL, beg, end);
  Particles_CR_Flux(stateR, beg-1, end-1);
  #endif
#endif

/* --------------------------------------------------------
   2. Sweep along 1D row of zones
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    #if UNIFORM_CARTESIAN_GRID == YES
    dp = dm = 0.5;
    #else
    dp = plm_coeffs.dp[i]; 
    dm = plm_coeffs.dm[i];
    #endif

    dx   = grid->dx[g_dir][i];
    dtdx = g_dt/dx;

    vc     = stateC->v[i];
    vp     = stateL->v[i];
    vm     = stateR->v[i-1];
    L      = stateC->Lp[i];
    R      = stateC->Rp[i];
    lambda = stateC->lambda[i];

  /* --------------------------------------------
     2a) Compute eigenvectors and eigenvalues if
         not yet defined
     -------------------------------------------- */

    for (k = 0; k < NFLX; k++) nu[k] = dtdx*lambda[k];
    nu_max = MAX(nu[1], 0.0); nu_min = MIN(nu[0], 0.0);

  /* --------------------------------------------
     2b) Geometry
     -------------------------------------------- */

    #if GEOMETRY == CYLINDRICAL
    if (g_dir == IDIR) {
      double *xR = grid->xr[IDIR];
      for (k = NFLX; k--;  ){
        betaR[k] = betaL[k] = 0.0;
        scrh = nu[k]*dx;
        betaR[k] = (nu[k] >  1.e-12 ?  scrh/(6.0*xR[i]   - 3.0*scrh):0.0);
        betaL[k] = (nu[k] < -1.e-12 ? -scrh/(6.0*xR[i-1] - 3.0*scrh):0.0);
      }
      nu_max *= (1.0 - betaR[1]);
      nu_min *= (1.0 + betaL[0]);
    }else{
      for (k = NFLX; k--;  )  betaR[k] = betaL[k] = 0.0;
    }
    #endif

  /* --------------------------------------------
     2c) Obtain characteristic variable
         increment dw
     -------------------------------------------- */

    for (nv = NVAR; nv--;  ) dv[nv] = vp[nv] - vm[nv];
    PrimToChar(L, dv, dw);

  /* ------------------------------------------------------
     2d)  Since this is somewhat arbitrary we use the value
          of ::CHTR_REF_STATE to select one of the
          following cases:
         
         - CHTR_REF_STATE==1: use cell-center value;
           
         - CHTR_REF_STATE==2: interpolated value at base  
                                    time level 
         - CHTR_REF_STATE==3: 
           traditional PPM reference state (fastest wave), minimize 
           the size of the term subject to characteristic limiting. 

        Passive scalars use always CHTR_REF_STATE == 2.      
     ------------------------------------------------------ */

    #if CHTR_REF_STATE == 1
    for (nv = NFLX; nv--; ) vp[nv] = vm[nv] = vc[nv];
    #elif CHTR_REF_STATE == 3
    for (nv = NFLX; nv--; ) {
      #if GEOMETRY == CARTESIAN
      vp[nv] = vc[nv] + 0.5*dv[nv]*(1.0 - nu_max);
      vm[nv] = vc[nv] - 0.5*dv[nv]*(1.0 + nu_min);
      #else
      vp[nv] = vc[nv] + dv[nv]*(dp - 0.5*nu_max);
      vm[nv] = vc[nv] - dv[nv]*(dm + 0.5*nu_min);
      #endif
    }
    #endif

  /* --------------------------------------------------------------------
     2e) Compute L/R states in primitive variables:

         Evolve interface states for dt/2 using characteristic tracing.
         Depending on the choice of the reference state (i.e. 
         REF_STATE = 1,2,3) we compute, e.g. Vp(t+dt/2), as

         Vp(t+dt/2) = \sum [w + dw*dfR - nu*dw*(1-beta)]*R =

    = V    + \sum [dfR - nu*(1-beta)]*dw*R                  (REF_STATE = 1)
    = Vp   + \sum [    - nu*(1-beta)]*dw*R                  (REF_STATE = 2)
    = Vmax + \sum [nu_max*(1-beta_max) - nu*(1-beta)]*dw*R  (REF_STATE = 3)

        where Vmax = Vp - nu_max*(1-beta_max) is the reference state
        of the original PPM algorithm.
     -------------------------------------------------------------------- */

    #ifdef GLM_MHD  /* -- hot fix for backward test compatibility 
                          must review this at some point !! -- */
    nu_max = nu[1]; nu_min = nu[0];
    #endif

    for (k = 0; k < NFLX; k++){
      if (nu[k] >= 0.0) {
        #if GEOMETRY != CARTESIAN
        nu[k] *= (1.0 - betaR[k]);
        #endif
        #if CHTR_REF_STATE == 1
        dw[k] *= 0.5*(1.0 - nu[k]);
        #elif CHTR_REF_STATE == 2
        dw[k] *= -0.5*nu[k];
        #elif CHTR_REF_STATE == 3
        dw[k] *= 0.5*(nu_max - nu[k]);
        #endif
        for (nv = 0; nv < NFLX; nv++) vp[nv] += dw[k]*R[nv][k];
      }else{     
        #if GEOMETRY != CARTESIAN
        nu[k] *= (1.0 + betaL[k]);
        #endif
        #if CHTR_REF_STATE == 1
        dw[k] *= -0.5*(1.0 + nu[k]);
        #elif CHTR_REF_STATE == 2
        dw[k] *= -0.5*nu[k];
        #elif CHTR_REF_STATE == 3
        dw[k] *= 0.5*(nu_min - nu[k]);
        #endif
        for (nv = 0; nv < NFLX; nv++) vm[nv] += dw[k]*R[nv][k];
      }
    }

  /* --------------------------------------------
     2f) Add source term to L/R states
     -------------------------------------------- */

    for (nv = NFLX; nv--;   ) {
      vp[nv] += 0.5*g_dt*src[i][nv];
      vm[nv] += 0.5*g_dt*src[i][nv];
    }

 /* ---------------------------------------------
     2g) Repeat construction for passive scalars
    --------------------------------------------- */

    #if NVAR != NFLX
    for (nv = NFLX; nv < NVAR; nv++) {
      #if GEOMETRY == CARTESIAN
      vp[nv] = vc[nv] + 0.5*dv[nv];
      vm[nv] = vc[nv] - 0.5*dv[nv];
      #else
      vp[nv] = vc[nv] + dv[nv]*dp;
      vm[nv] = vc[nv] - dv[nv]*dm;
      #endif
    }
     
    nu[0] = dtdx*vc[VXn];
    #if GEOMETRY == CYLINDRICAL
    if (g_dir == IDIR) {
      double *xR = grid->xr[IDIR];
      betaR[0] = betaL[0] = 0.0;
      scrh = nu[0]*dx;
       betaR[0] = (nu[0] >  1.e-12 ?  scrh/(6.0*xR[i]   - 3.0*scrh):0.0);
       betaL[0] = (nu[0] < -1.e-12 ? -scrh/(6.0*xR[i-1] - 3.0*scrh):0.0);
     }else{
       betaR[0] = betaL[0] = 0.0;
     }
    #endif

    if (nu[0] >= 0.0){  /* -- scalars all move at the flow speed -- */
      scrh = -0.5*nu[0];
      #if GEOMETRY != CARTESIAN
       scrh *= (1.0 - betaR[0]);
      #endif
      for (k = NFLX; k < NVAR; k++) vp[k] += dw[k]*scrh;
    }else {
      scrh = -0.5*nu[0];
      #if GEOMETRY != CARTESIAN
       scrh *= (1.0 + betaL[0]);
      #endif
      for (k = NFLX; k < NVAR; k++) vm[k] += dw[k]*scrh;
    }
    #endif

  }  /* -- end main loop on grid points -- */

/* --------------------------------------------------------
   3. Apply 1D shock flattening (only 1D)
   --------------------------------------------------------  */

#if SHOCK_FLATTENING == ONED
  Flatten (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   4. Assign face-centered magnetic field
   --------------------------------------------------------  */

#ifdef STAGGERED_MHD
  for (i = beg-1; i <= end; i++) {
   stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
  }
#endif

/* --------------------------------------------------------
   5. Add CR source term 
   -------------------------------------------------------- */

#ifdef PARTICLES
#if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
  Particles_CR_StatesSource(sweep, 0.5*g_dt, beg, end, grid);
#endif
#endif

/* --------------------------------------------------------
   6. Evolve cell-center values by dt/2
   -------------------------------------------------------- */

  CheckPrimStates (stateL->v, stateR->v-1, stateC->v, beg, end);

  for (i = beg; i <= end; i++){ 
    vp = stateL->v[i];
    vm = stateR->v[i-1];
    NVAR_LOOP(nv) stateC->v[i][nv] = 0.5*(vp[nv] + vm[nv]);
  }

#if RECONSTRUCT_4VEL
  ConvertTo3vel (stateC->v, beg, end);
#endif

}
#undef CHTR_REF_STATE 
#endif  /* RECONSTRUCTION == LINEAR */
