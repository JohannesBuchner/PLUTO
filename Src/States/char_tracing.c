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
  \date   Jan 28, 2014
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
void CharTracingStep(const State_1D *state, int beg, int end, Grid *grid)
/*! 
 *  Compute interface states using characteristic tracing step.
 *
 * \param [in]  state  pointer to a State_1D structure
 * \param [in]    beg  initial index of computation
 * \param [in]    end  final   index of computation
 * \param [in]   grid  pointer to an array of Grid structures

 * \return This function has no return value.
 *
 * Tasks are numbered below.
 *
 *********************************************************************** */
{
  int    i, nv, k;
  double dx, dtdx;
  double dwh, d2w, dwp[NVAR], dwm[NVAR], dvp[NVAR], dvm[NVAR];
  double nup=1.0, num=-1.0, nu[NFLX];
  double Spp, Smm;
  double *vc, *vp, *vm, **L, **R, *lambda;
  static double **src;

  if (src == NULL){
    src = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  #if CHAR_LIMITING == NO
   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
  #endif

  PrimSource (state, beg, end, state->a2, state->h, src, grid);
  for (i = beg; i <= end; i++){    

    dx   = grid[g_dir].dx[beg];
    dtdx = g_dt/dx;

    vc     = state->v[i]; 
    vp     = state->vp[i];
    vm     = state->vm[i];
    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

  /* ------------------------------------------------------------ */
  /*! 1) Compute eigenvectors and eigenvalues if not yet defined. */
  /* ------------------------------------------------------------ */

    #if CHAR_LIMITING == NO
     PrimEigenvectors (vc, state->a2[i], state->h[i], lambda, L, R);
    #endif
    for (k = 0; k < NFLX; k++) nu[k] = dtdx*lambda[k];

  /* ------------------------------------------------------------- */
  /*! 2) Obtain characteristic variable increments dwp and dwm.    */
  /* ------------------------------------------------------------- */

    for (nv = NVAR; nv--;  ) {
      dvp[nv] = vp[nv] - vc[nv];
      dvm[nv] = vm[nv] - vc[nv];
    }     
    PrimToChar(L, dvp, dwp);
    PrimToChar(L, dvm, dwm);

  /* ------------------------------------------------------------- */
  /*! 3) Initialize vp and vm to the reference state.
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
  /*! 4) Compute left and right states in primitive variables. 
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
  /*! 5) Add source term to L/R states                       */
  /* ------------------------------------------------------- */

    for (nv = NFLX; nv--;   ){   
      dwh = 0.5*g_dt*src[i][nv];
      vp[nv] += dwh;
      vm[nv] += dwh;
    }

  /* ------------------------------------------------------- */
  /*! 6) Repeat construction for passive scalars             */
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
  }

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

  CheckPrimStates (state->vm, state->vp, state->v, beg, end);

  for (i = beg; i <= end; i++) {
    vp = state->vp[i];
    vm = state->vm[i];
    NVAR_LOOP(nv) state->vh[i][nv] = 0.5*(vp[nv] + vm[nv]);
  }
#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->vh, beg, end);
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
void CharTracingStep(const State_1D *state, int beg, int end, Grid *grid)
/*
 * Same thing as before, but optimized for linear reconstruction.
 * 
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dx, dtdx, scrh;
  double dv[NVAR], dw[NVAR];
  double nu_max=1.0, nu_min=-1.0, nu[NFLX];
  double *vc, *vp, *vm, **L, **R, *lambda;
  double dp, dm;
  PLM_Coeffs plm_coeffs;
  #if GEOMETRY != CARTESIAN
   double betaL[NVAR], betaR[NVAR];
  #endif
  static double **src;

/* --------------------------------------------
    allocate memory and set pointer shortcuts
   -------------------------------------------- */

  if (src == NULL) src = ARRAY_2D(NMAX_POINT, NVAR, double);

  #if UNIFORM_CARTESIAN_GRID == NO
   PLM_CoefficientsGet (&plm_coeffs, g_dir);
  #endif

  #if CHAR_LIMITING == NO
   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
  #endif

  PrimSource  (state, beg, end, state->a2, state->h, src, grid);
  for (i = beg; i <= end; i++){    

    #if UNIFORM_CARTESIAN_GRID == YES
     dp = dm = 0.5;
    #else
     dp = plm_coeffs.dp[i]; 
     dm = plm_coeffs.dm[i];
    #endif

    dx   = grid[g_dir].dx[i];
    dtdx = g_dt/dx;

    vc     = state->v[i];
    vp     = state->vp[i];
    vm     = state->vm[i];
    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

  /* --------------------------------------------------------------
     1. Compute eigenvectors and eigenvalues if not yet defined
     -------------------------------------------------------------- */

    #if CHAR_LIMITING == NO
     PrimEigenvectors (vc, state->a2[i], state->h[i], lambda, L, R);
    #endif
    for (k = 0; k < NFLX; k++) nu[k] = dtdx*lambda[k];
    nu_max = MAX(nu[1], 0.0); nu_min = MIN(nu[0], 0.0);

  /* ------------------------------------------------------------
     1a. Geometry
     ------------------------------------------------------------ */

    #if GEOMETRY == CYLINDRICAL
     if (g_dir == IDIR) {
       double *xR = grid[IDIR].xr;
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

  /* --------------------------------------------------------------
     2. Obtain characteristic variable increment dw
     -------------------------------------------------------------- */

    for (nv = NVAR; nv--;  ) dv[nv] = vp[nv] - vm[nv];
    PrimToChar(L, dv, dw);

  /* ----------------------------------------------------------------
         Since this is somewhat arbitrary we use the value of 
         ::CHTR_REF_STATE to select one of the following cases:
         
         - CHTR_REF_STATE==1: use cell-center value;
           
         - CHTR_REF_STATE==2: interpolated value at base 
                                    time level 
         - CHTR_REF_STATE==3: 
           traditional PPM reference state (fastest wave), minimize 
           the size of the term subject to characteristic limiting. 

        Passive scalars use always CHTR_REF_STATE == 2.      
     ---------------------------------------------------------------- */

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
     4. Compute L/R states in primitive variables:

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

  /* -------------------------------------------------------
     5. Add source term to L/R states
     ------------------------------------------------------- */

    for (nv = NFLX; nv--;   ) {
      vp[nv] += 0.5*g_dt*src[i][nv];
      vm[nv] += 0.5*g_dt*src[i][nv];
    }

 /* -------------------------------------------------------
     6. Repeat construction for passive scalars
    ------------------------------------------------------- */

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
        double *xR = grid[IDIR].xr;
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

/*  -------------------------------------------
        Shock flattening (only 1D)
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
          evolve cell-center values by dt/2
   -------------------------------------------------------- */

  CheckPrimStates (state->vm, state->vp, state->v, beg, end);

  for (i = beg; i <= end; i++){ 
    vp = state->vp[i];
    vm = state->vm[i];
    NVAR_LOOP(nv) state->vh[i][nv] = 0.5*(vp[nv] + vm[nv]);
  }
#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->vh, beg, end);
#endif  
}
#undef CHTR_REF_STATE 
#endif  /* RECONSTRUCTION == LINEAR */
