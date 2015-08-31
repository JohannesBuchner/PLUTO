/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief MUSCL-Hancock predictor step.

  Use time-extrapolation to compute interface states and cell-centered 
  value at the half-time step, 
  \f$ V_{i,\pm}^{n+\HALF} = V_{i,\pm}^n + \partial_t V \Delta t^n/2 \f$

  This is done using the one-dimensional primitive form ot the 
  equations for the HD, RHD and MHD modules since we have at disposal
  \f$ \partial_t V = -A\partial V_x + S \f$.
  
  Conversely, for relativistic MHD, we adopt the one-dimensional 
  conservative form of the equations (conservative Hancock predictor) and 
  compute the primitive values using
  \f$ \partial_t V \approx (V^{n+\HALF}_i-V^n_i)/(\Delta t/2) \f$.
  This requires taking the following steps:
  
  - Convert to conservative: vp -> up, vm -> um
  - Evolve  u(n+1/2) = u(n) - dt/(2dx)*[F(vp) - F(vm)]
  - Convert to primitive  u(n+1/2) -> v(n+1/2)
  - Add time incerement v(n+1/2)-v(n) to left/right states

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS != RMHD
/* ********************************************************************* */
void HancockStep (const State_1D *state, int beg, int end, Grid *grid)
/*!
 * Use Taylor expansion to compute time-centered states in primitive
 * variables.
 *
 * \param [in,out] state  pointer to a State_1D structure
 * \param [in]     beg    initial index of computation
 * \param [in]     end    final index of computation 
 * \param [in      grid   pointer to array of Grid structures
 *
 * \b References
 *    - "Riemann Solvers and Numerical Methods for  Fluid Dynamics" \n
 *      E.F. Toro
 *
 *********************************************************************** */
{
  int    nv, i;
  double scrh, dt_2, ch2;
  double Adv[NVAR], dv[NVAR];
  double *vp, *vm, *vc;
  static double **src, *d_dl;

/* -----------------------------------------
         Check scheme compatibility
   ----------------------------------------- */

  #if RECONSTRUCTION != LINEAR
   print1 ("! MUSCL-Hancock scheme works with Linear reconstruction only\n");
   QUIT_PLUTO(1);
  #endif


  if (src == NULL){
    src  = ARRAY_2D(NMAX_POINT, NVAR, double);
    d_dl = ARRAY_1D(NMAX_POINT, double);
  }

/* -------------------------------------------------------
        compute geometric factor 1/(h*dx)
   ------------------------------------------------------- */

  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL

   for (i = beg; i <= end; i++) {
     d_dl[i] = 1.0/grid[g_dir].dx[i];
   }

  #elif GEOMETRY == SPHERICAL || GEOMETRY == POLAR

   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
       d_dl[i] = 1.0/grid[g_dir].dx[i];
     }
   }else if (g_dir == JDIR){
     for (i = beg; i <= end; i++) {
       d_dl[i] = 1.0/(grid[IDIR].x[g_i]*grid[JDIR].dx[i]);
     }
   }

  #endif

  #if CHAR_LIMITING == NO
   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
  #endif

  PrimSource  (state, beg, end, state->a2, state->h, src, grid);
  dt_2 = 0.5*g_dt;
  for (i = beg; i <= end; i++) {
    vc = state->v[i];
    vm = state->vm[i];
    vp = state->vp[i];

    for (nv = NVAR; nv--;  ) dv[nv] = vp[nv] - vm[nv];

    PrimRHS (vc, dv, state->a2[i], state->h[i], Adv);
    for (nv = NVAR; nv--;  ) {
      scrh = dt_2*(d_dl[i]*Adv[nv] - src[i][nv]);
      vp[nv] -= scrh;
      vm[nv] -= scrh;
    }
  }

  CheckPrimStates (state->vm, state->vp, state->v, beg, end);

/* -----------------------------------
     evolve center value by dt/2 
   ----------------------------------- */

  for (i = beg; i <= end; i++) {
    vc = state->vh[i];
    vp = state->vp[i];
    vm = state->vm[i];
    NVAR_LOOP(nv) vc[nv] = 0.5*(vp[nv] + vm[nv]);
  }
#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->vh, beg, end);
#endif  

}

#elif PHYSICS == RMHD

/* ********************************************************************* */
void HancockStep (const State_1D *state, int beg, int end, Grid *grid)
/*
 *
 *********************************************************************** */
#define CHAR_UPWIND  NO
 {
  int    nv, i, ifail;
  double dBx, dpsi, dtdx, dt, dtdV;
  double *A, *dx, *dV, r;
  double *vp, *vm, *vc, dvp[NVAR], dvm[NVAR];
  double P0, P1, P2; 
  double **rhs;
  static double **fp, **fm, **uh, **u;
  static double *pp, *pm, *hp, *hm;
  static double *lambda_max, *lambda_min;

#if GEOMETRY != CARTESIAN && GEOMETRY != CYLINDRICAL
  print1 ("! Hancock does not work in this geometry \n");
  QUIT_PLUTO(1);
#endif
#if BODY_FORCE != NO
  print ("! Conservative Hancock scheme does not support gravity\n");
  QUIT_PLUTO(1);
#endif
#if RECONSTRUCTION != LINEAR
  print1 ("! The MUSCL-Hancock scheme works with Linear reconstruction only\n");
  QUIT_PLUTO(1);
#endif
#if (PARABOLIC_FLUX & EXPLICIT)
  print1 ("! Conservative MUSCL-Hancock incompatible with Explicit Diffusion\n");
  QUIT_PLUTO(1);
#endif

  if (fp == NULL){
    fp   = ARRAY_2D(NMAX_POINT, NVAR, double);
    fm   = ARRAY_2D(NMAX_POINT, NVAR, double);
    pp   = ARRAY_1D(NMAX_POINT, double);
    pm   = ARRAY_1D(NMAX_POINT, double);
    u    = ARRAY_2D(NMAX_POINT, NVAR, double);
    uh   = ARRAY_2D(NMAX_POINT, NVAR, double);

    hp   = ARRAY_1D(NMAX_POINT, double);
    hm   = ARRAY_1D(NMAX_POINT, double);
    lambda_max = ARRAY_1D(NMAX_POINT, double);
    lambda_min = ARRAY_1D(NMAX_POINT, double);
  }
  rhs = state->rhs;

#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->v, beg-1, end+1);
  ConvertTo3vel (state->vp, beg, end);
  ConvertTo3vel (state->vm, beg, end);
#endif

/* --------------------------------
         make fluxes  
   -------------------------------- */

  PrimToCons (state->v, u, beg, end);
  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);

  Enthalpy (state->vp, hp, beg, end);
  Enthalpy (state->vm, hm, beg, end);

  Flux(state->up, state->vp, hp, fp, pp, beg, end);
  Flux(state->um, state->vm, hm, fm, pm, beg, end);

  dx = grid[g_dir].dx; A = grid[g_dir].A; dV = grid[g_dir].dV;

/* ---------------------------------------------------
          Compute passive scalar fluxes
   --------------------------------------------------- */

  #if NFLX != NVAR
   for (i = beg; i <= end; i++){
     vp  = state->vp[i]; vm  = state->vm[i];
     #if NSCL > 0
      for (nv = NFLX; nv < (NFLX + NSCL); nv++){
        fp[i][nv] = fp[i][RHO]*vp[nv];
        fm[i][nv] = fm[i][RHO]*vm[nv];
      }
     #endif
   }
  #endif

  #if CHAR_UPWIND == YES  
   MAX_CH_SPEED (state->v, lambda_min, lambda_max, beg, end);
  #endif

/* ------------------------------------------------------
      Initialize right-hand-side with flux difference  
   ------------------------------------------------------ */

  dt = 0.5*g_dt; 
  for (i = beg; i <= end; i++) { 

    dtdx = dt/dx[i]; 
    vp  = state->vp[i]; 
    vm  = state->vm[i];
    #if GEOMETRY == CARTESIAN
     for (nv = NVAR; nv--;   ) rhs[i][nv] = -dtdx*(fp[i][nv] - fm[i][nv]);
    #elif GEOMETRY == CYLINDRICAL
     dtdV = dt/dV[i];
     for (nv = NVAR; nv--;   ) {
       rhs[i][nv] = -dtdV*(A[i]*fp[i][nv] - A[i-1]*fm[i][nv]);
     }
     #ifdef GLM_MHD
      rhs[i][BXn] = -dtdx*(vp[PSI_GLM] - vm[PSI_GLM]);
     #endif
    #endif
    rhs[i][MXn] -= dtdx*(pp[i] - pm[i]);

  /* ---- add Powell source term ---- */

    #if (PHYSICS == MHD || PHYSICS == RMHD) && (DIVB_CONTROL == EIGHT_WAVES)
     dBx = (vp[BXn]*A[i] - vm[BXn]*A[i-1])/dV[i];
     EXPAND(rhs[i][BX1] -= dt*state->v[i][VX1]*dBx;  ,
            rhs[i][BX2] -= dt*state->v[i][VX2]*dBx;  ,
            rhs[i][BX3] -= dt*state->v[i][VX3]*dBx;)
    #endif

  /* ---- Extended GLM source term ---- */

    #ifdef GLM_MHD
     #if GLM_EXTENDED == YES
      dBx  = (vp[BXn]*A[i] - vm[BXn]*A[i-1])/dV[i];
      EXPAND(rhs[i][MX1] -= dt*state->v[i][BX1]*dBx;  ,
             rhs[i][MX2] -= dt*state->v[i][BX2]*dBx;  ,
             rhs[i][MX3] -= dt*state->v[i][BX3]*dBx;)
      #if HAVE_ENERGY
       dpsi = (vp[PSI_GLM] - vm[PSI_GLM])/dx[i];
       rhs[i][ENG] -= dt*state->v[i][BXn]*dpsi;
      #endif
     #endif
    #endif
  }
   
/* ------------------------------------------------
       Source terms: geometry
   ------------------------------------------------ */

  #if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
   if (g_dir == IDIR) {
     double vB, vel2, g_2, r_1, *uc;

     for (i = beg; i <= end; i++){
       r_1 = grid[IDIR].r_1[i];
       vc  = state->v[i];
       uc  = u[i];

       vB   = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
       vel2 = EXPAND(vc[VX1]*vc[VX1], + vc[VX2]*vc[VX2], + vc[VX3]*vc[VX3]);
       g_2  = 1.0 - vel2;

       rhs[i][iMR] += dt*((uc[MX3]*vc[VX3] - (vc[BX3]*g_2 + vB*vc[VX3])*vc[BX3]))*r_1;

       rhs[i][iMPHI] = -dtdV*(A[i]*A[i]*fp[i][iMPHI] - A[i-1]*A[i-1]*fm[i][iMPHI]);
       rhs[i][iMPHI] *= r_1;

       rhs[i][iBPHI] -= dt*(vc[VX3]*uc[BX1] - uc[BX3]*vc[VX1])*r_1;
     }
   }
  #endif

/* ------------------------------------------------
        Update solution vector to get u(n+1/2)
   ------------------------------------------------ */

  for (i = beg; i <= end; i++) {
    for (nv = NVAR; nv--;  ) uh[i][nv] = u[i][nv] + rhs[i][nv];
  }

/* -------------------------------------------------
     check if U(x_i, n+1/2) has physical meaning.
     Revert to 1st order otherwise.
   ------------------------------------------------- */

  if (ConsToPrim (uh, state->vh, beg, end, state->flag) != 0){
    for (i = beg; i <= end; i++){
      if (state->flag[i] != 0){
        for (nv = NVAR; nv--;   ){
          uh[i][nv]        = u[i][nv];
          state->vh[i][nv] = state->v[i][nv];
          state->vp[i][nv] = state->vm[i][nv] = state->v[i][nv];
        }
        #ifdef STAGGERED_MHD
         state->vp[i][BXn] = state->bn[i];
         state->vm[i][BXn] = state->bn[i-1];
        #endif
      }
    }
  }    

/* ----------------------------------------------------
        evolve interface values accordingly 
   ---------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    vp  = state->vp[i]; 
    vm  = state->vm[i];
    vc  = state->v[i];

    for (nv = NVAR; nv--;    ){
/*
      dvp[nv] = dvm[nv] = state->vh[i][nv] - state->v[i][nv]; 

      P0 = 1.5*vc[nv] - 0.25*(vp[nv] + vm[nv]);
      P1 = vp[nv] - vm[nv];
      P2 = 3.0*(vp[nv] - 2.0*vc[nv] + vm[nv]);

      dvp[nv] += P0 + 0.5*P1 + 0.25*P2 - vc[nv];
      dvm[nv] += P0 - 0.5*P1 + 0.25*P2 - vc[nv];
*/
      P1 = vp[nv] - vm[nv];
      dvp[nv] = dvm[nv] = state->vh[i][nv] - state->v[i][nv]; 
      dvp[nv] += 0.5*P1;
      dvm[nv] -= 0.5*P1;
    }

    #if CHAR_UPWIND == YES
  
     if (lambda_min[i] > 0.0) 
       for (nv = NVAR; nv--;    )  dvm[nv] = 0.0;
     else if (lambda_max[i] < 0.0)
       for (nv = NVAR; nv--;    )  dvp[nv] = 0.0;
        
    #endif

    for (nv = NVAR; nv--;    ){
      vp[nv] = vc[nv] + dvp[nv];
      vm[nv] = vc[nv] + dvm[nv];
    }
  }
  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++){
     state->vp[i][BXn] = state->vm[i+1][BXn] = state->bn[i];
   }
  #endif

  CheckPrimStates (state->vm, state->vp, state->v, beg, end);
}
#undef CHAR_UPWIND
#endif

