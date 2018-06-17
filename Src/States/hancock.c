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
  \date   April 16, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PRIMITIVE_HANCOCK == YES
/* ********************************************************************* */
void HancockStep (const Sweep *sweep, int beg, int end, Grid *grid)
/*!
 * Use Taylor expansion to compute time-centered states in primitive
 * variables.
 *
 * \param [in,out] sweep  pointer to a Sweep structure
 * \param [in]     beg    initial index of computation
 * \param [in]     end    final index of computation 
 * \param [in      grid   pointer to Grid structure
 *
 * \b References
 *    - "Riemann Solvers and Numerical Methods for  Fluid Dynamics" \n
 *      E.F. Toro
 *
 *********************************************************************** */
{
  int    nv, i;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double scrh, dt_2, ch2;
  double Adv[NVAR], dv[NVAR];
  double *vp, *vm, *vc;
  static double **src, *d_dl;
  double *a2 = stateC->a2;

/* ------------------------------------------------
   0. Check scheme compatibility, allocate memory
   ------------------------------------------------ */

#if RECONSTRUCTION != LINEAR
  print ("! MUSCL-Hancock scheme works with Linear reconstruction only\n");
  QUIT_PLUTO(1);
#endif
#if PHYSICS == RMHD
  print ("! Primitive MUSCL-Hancock scheme does not work with RMHD\n");
  QUIT_PLUTO(1);
#endif

  if (src == NULL){
    src  = ARRAY_2D(NMAX_POINT, NVAR, double);
    d_dl = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
   1. Compute geometric factor 1/(h*dx)
   -------------------------------------------------------- */

#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  for (i = beg; i <= end; i++) d_dl[i] = 1.0/grid->dx[g_dir][i];  
#elif GEOMETRY == SPHERICAL || GEOMETRY == POLAR
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) d_dl[i] = 1.0/grid->dx[g_dir][i];     
  }else if (g_dir == JDIR){
    for (i = beg; i <= end; i++) {
      d_dl[i] = 1.0/(grid->x[IDIR][g_i]*grid->dx[JDIR][i]);
    }
  }
#endif

/* --------------------------------------------------------
   2. Compute preliminary quantities such as sound speed,
      sources, etc...
      CR flux is computed at time level t^n and will be
      added at the end of the function.
   -------------------------------------------------------- */

#if CHAR_LIMITING == NO
  SoundSpeed2 (stateC, beg, end, CELL_CENTER, grid);
#endif
  PrimSource (stateC, src, beg, end, grid);
#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
  Particles_CR_Flux(stateL, beg, end);
  Particles_CR_Flux(stateR, beg-1, end-1);
  #endif
#endif

/* --------------------------------------------------------
   3. Loop over 1D row of zones
   -------------------------------------------------------- */

  dt_2 = 0.5*g_dt;
  for (i = beg; i <= end; i++) {
    vc = stateC->v[i];
    vp = stateL->v[i];
    vm = stateR->v[i-1];

    NVAR_LOOP(nv) dv[nv] = vp[nv] - vm[nv];

    PrimRHS (vc, dv, a2[i], stateC->h[i], Adv);
    for (nv = NVAR; nv--;  ) {
      scrh = dt_2*(d_dl[i]*Adv[nv] - src[i][nv]);
      vp[nv] -= scrh;
      vm[nv] -= scrh;
    }
  }

/* --------------------------------------------------------
   4. Add CR source term 
   -------------------------------------------------------- */

#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
  Particles_CR_StatesSource(sweep, 0.5*g_dt, beg, end, grid);
  #endif
#endif

/* --------------------------------------------------------
   5. Evolve center value by dt/2 
   -------------------------------------------------------- */

  CheckPrimStates (stateR->v-1, stateL->v, stateC->v, beg, end);
  for (i = beg; i <= end; i++) {
    vc = stateC->v[i];
    vp = stateL->v[i];
    vm = stateR->v[i-1];
    NVAR_LOOP(nv) vc[nv] = 0.5*(vp[nv] + vm[nv]);
  }
#if RECONSTRUCT_4VEL
  ConvertTo3vel (stateC->v, beg, end);
#endif  

}

#else

/* ********************************************************************* */
void HancockStep (const Sweep *sweep, int beg, int end, Grid *grid)
/*!
 * Conservative Hancock step.
 *
 *********************************************************************** */
 {
  int    nv, i, ifail;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double dBx, dpsi, dtdx, dt, dtdV;
  double *dx, r;
  double *vp, *vm, *vc, *uc, dvp[NVAR], dvm[NVAR];
  double P0, P1, P2; 
  double **rhs = sweep->rhs;
  double **fp = stateL->flux, **fm = stateR->flux - 1;
  double  *hp = stateL->h,     *hm = stateR->h    - 1;
  double  *pp = stateL->prs,   *pm = stateR->prs  - 1;

  static double **uh, **vn, *A, *dV;
  static double *lambda_max, *lambda_min;

/* --------------------------------------------------------
   0. Check scheme compatibility and allocate memory
   -------------------------------------------------------- */

#if GEOMETRY != CARTESIAN && GEOMETRY != CYLINDRICAL
  print ("! HancockStep(): conservative MUSCL-Hancock does not work\n");
  print ("                 in this geometry \n");
  QUIT_PLUTO(1);
#endif
#if BODY_FORCE != NO
  print ("! HancockStep(): conservative Hancock scheme does not support gravity\n");
  QUIT_PLUTO(1);
#endif
#if RECONSTRUCTION != LINEAR
  print ("! HancockStep(): conservative Hancock scheme works with \n");
  print ("                 linear reconstruction only\n");
  QUIT_PLUTO(1);
#endif
#if (PARABOLIC_FLUX & EXPLICIT)
  print ("! HancockStep(): conservative MUSCL-Hancock does not support \n");
  print ("                 explicit diffusion \n");
  QUIT_PLUTO(1);
#endif

  if (uh == NULL){
    uh         = ARRAY_2D(NMAX_POINT, NVAR, double); /* U^{n+1/2} */
    vn         = ARRAY_2D(NMAX_POINT, NVAR, double); /* V^n */
    lambda_max = ARRAY_1D(NMAX_POINT, double);
    lambda_min = ARRAY_1D(NMAX_POINT, double);
    A          = ARRAY_1D(NMAX_POINT, double);
    dV         = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
   1. Compute fluxes f(vp), f(vm)
   -------------------------------------------------------- */

  PrimToCons (stateC->v, stateC->u, beg, end);
  PrimToCons (stateL->v, stateL->u, beg, end);
  PrimToCons (stateR->v, stateR->u, beg-1, end-1);

  Enthalpy (stateL->v, stateL->h, beg, end);
  Enthalpy (stateR->v, stateR->h, beg-1, end-1);

  Flux(stateL, beg, end);
  Flux(stateR, beg-1, end-1);
#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux (stateL, beg, end);
  Particles_CR_Flux (stateR, beg-1, end-1);
  for (i = beg; i <= end; i++){
    NVAR_LOOP(nv){
      fp[i][nv] += stateL->fluxCR[i][nv];
      fm[i][nv] += stateR->fluxCR[i-1][nv];
    }
  }
  #endif
#endif

  dx = grid->dx[g_dir];

  if (g_dir == IDIR){
    #if GEOMETRY == CARTESIAN
    for (i = beg-1; i <= end; i++) { 
      A[i]  = 1.0;
      dV[i] = dx[i];
    }
    #elif GEOMETRY == CYLINDRICAL
    for (i = beg-1; i <= end; i++) { 
      A[i]  = fabs(grid->xr[IDIR][i]);
      dV[i] = fabs(grid->x[IDIR][i])*dx[i];
    }
    #endif
  }else{
    for (i = beg-1; i <= end; i++) { 
      A[i]  = 1.0;
      dV[i] = dx[i];
    }
  }

/* --------------------------------------------------------
   1a. Compute passive scalar fluxes
   -------------------------------------------------------- */

#if NFLX != NVAR
  for (i = beg; i <= end; i++){
    vp  = stateL->v[i];
    vm  = stateR->v[i-1];
    #if NSCL > 0
    for (nv = NFLX; nv < (NFLX + NSCL); nv++){
      fp[i][nv] = fp[i][RHO]*vp[nv];
      fm[i][nv] = fm[i][RHO]*vm[nv];
    }
   #endif
  }
#endif

/* --------------------------------------------------------
   2.  Initialize right-hand-side with flux difference 
       rhs = dt/2*[- (Ap*fp - Am*fm)/dV + S]
   -------------------------------------------------------- */

  dt = 0.5*g_dt; 
  for (i = beg; i <= end; i++) { 
    dtdx = dt/dx[i]; 
    vc  = stateC->v[i];
    vp  = stateL->v[i]; 
    vm  = stateR->v[i-1];

    #if GEOMETRY == CARTESIAN
    NVAR_LOOP(nv) rhs[i][nv] = -dtdx*(fp[i][nv] - fm[i][nv]);
    #elif GEOMETRY == CYLINDRICAL
    dtdV = dt/dV[i];
    NVAR_LOOP(nv) rhs[i][nv] = -dtdV*(A[i]*fp[i][nv] - A[i-1]*fm[i][nv]);
    #ifdef GLM_MHD
    rhs[i][BXn] = -dtdx*(vp[PSI_GLM] - vm[PSI_GLM]);
    #endif
    #endif
    rhs[i][MXn] -= dtdx*(pp[i] - pm[i]);

  /* -- 2a. Add Powell source term -- */

    #if (PHYSICS == MHD || PHYSICS == RMHD) && (DIVB_CONTROL == EIGHT_WAVES)
    dBx = (vp[BXn]*A[i] - vm[BXn]*A[i-1])/dV[i];
    EXPAND(rhs[i][BX1] -= dt*vc[VX1]*dBx;  ,
           rhs[i][BX2] -= dt*vc[VX2]*dBx;  ,
           rhs[i][BX3] -= dt*vc[VX3]*dBx;)
    #endif

  /* -- 2b- Add extended GLM source term -- */

    #if (defined GLM_MHD) && GLM_EXTENDED == YES
    dBx  = (vp[BXn]*A[i] - vm[BXn]*A[i-1])/dV[i];
    EXPAND(rhs[i][MX1] -= dt*vc[BX1]*dBx;  ,
           rhs[i][MX2] -= dt*vc[BX2]*dBx;  ,
           rhs[i][MX3] -= dt*vc[BX3]*dBx;)
    #if HAVE_ENERGY
    dpsi = (vp[PSI_GLM] - vm[PSI_GLM])/dx[i];
    rhs[i][ENG] -= dt*vc[BXn]*dpsi;
    #endif
    #endif
  }
   
/* --------------------------------------------------------
   3. Add source term in cylindrical geometry
   -------------------------------------------------------- */

#if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
  #if PHYSICS == RMHD
  if (g_dir == IDIR) {
    double vB, vel2, g_2, r_1, *uc;

    for (i = beg; i <= end; i++){
      r_1 = 1.0/grid->x[IDIR][i];
      vc  = stateC->v[i];
      uc  = stateC->u[i];

      vB   = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
      vel2 = EXPAND(vc[VX1]*vc[VX1], + vc[VX2]*vc[VX2], + vc[VX3]*vc[VX3]);
      g_2  = 1.0 - vel2;

      rhs[i][iMR] += dt*((uc[MX3]*vc[VX3] - (vc[BX3]*g_2 + vB*vc[VX3])*vc[BX3]))*r_1;

      rhs[i][iMPHI] = -dtdV*(A[i]*A[i]*fp[i][iMPHI] - A[i-1]*A[i-1]*fm[i][iMPHI]);
      rhs[i][iMPHI] *= r_1;

      rhs[i][iBPHI] -= dt*(vc[VX3]*uc[BX1] - uc[BX3]*vc[VX1])*r_1;
    }
  }
  #else
  print ("! Hancock(): conservative hancock step does not work in this geometry\n");
  QUIT_PLUTO(1);
  #endif 
#endif

/* --------------------------------------------------------
   4. Update solution array: U^{n+1/2} = U^n + rhs
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    NVAR_LOOP(nv) {
      uh[i][nv] = stateC->u[i][nv] + rhs[i][nv];
      vn[i][nv] = stateC->v[i][nv];  /* Save a copy before overwriting */
    }
    
  /* -- Add CR feedback (source term) to conservative cell-centered values -- */

    #ifdef PARTICLES
    #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES) 
    uh[i][MX1] -= dt*stateC->Fcr[i][IDIR];
    uh[i][MX2] -= dt*stateC->Fcr[i][JDIR];
    uh[i][MX3] -= dt*stateC->Fcr[i][KDIR];
    uh[i][ENG] -= dt*stateC->Fcr[i][3];
    #endif
    #endif
  }

/* --------------------------------------------------------
   5. Check if U^{n+1/2} has physical meaning.
      Revert to 1st order otherwise.
   -------------------------------------------------------- */

  if (ConsToPrim (uh, stateC->v, beg, end, sweep->flag) != 0){
    for (i = beg; i <= end; i++){
      if (sweep->flag[i] != 0){
        NVAR_LOOP(nv) {
          uh[i][nv] = stateC->u[i][nv];
          stateC->v[i][nv] = vn[i][nv];
          stateL->v[i][nv] = stateR->v[i-1][nv] = vn[i][nv];
        }
        #ifdef STAGGERED_MHD
        stateL->v[i][BXn]   = sweep->bn[i];
        stateR->v[i-1][BXn] = sweep->bn[i-1];
        #endif
      }
    }
  }    

/* --------------------------------------------------------
   6. Evolve interface values accordingly
      vp(n+1/2) = v(n) + dv/2 + (v(n+1/2) - v(n)) 
      vm(n+1/2) = v(n) - dv/2 + (v(n+1/2) - v(n)) 
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    vp  = stateL->v[i]; 
    vm  = stateR->v[i-1];
    
    NVAR_LOOP(nv){
      P1 = vp[nv] - vm[nv];
      dvp[nv] = dvm[nv] = stateC->v[i][nv] - vn[i][nv]; 
      dvp[nv] += 0.5*P1;
      dvm[nv] -= 0.5*P1;
    }

    NVAR_LOOP(nv) {
      vp[nv] = vn[i][nv] + dvp[nv];
      vm[nv] = vn[i][nv] + dvm[nv];
//      stateC->u[i][nv] = uh[i][nv];
    }
  }

#ifdef STAGGERED_MHD
  for (i = beg-1; i <= end; i++){
    stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
  }
#endif

  CheckPrimStates (stateL->v, stateR->v-1, vn, beg, end);
}
#endif

