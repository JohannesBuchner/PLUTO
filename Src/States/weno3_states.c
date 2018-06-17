#include "pluto.h"

static void WENO3_COEFF(double **, double **, double **, double **, Grid *);

/* ************************************************************* */
void States (const Sweep *sweep, int beg, int end, Grid *grid)
/*! 
 *
 * Provide a three-point stencil, third-order 
 * reconstruction algorithm based on the WENO3.
 *
 *
 * LAST MODIFIED
 *
 *   March 5, 2017
 *   by A. Mignone (mignone@ph.unito.it)
 *
 **************************************************************** */
{
  int    k, nv, i;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);
  double **v  = stateC->v;
  double **vp = stateL->v;
  double **vm = stateR->v-1;
  double **up = stateL->u;
  double **um = stateR->u-1;

  double dmm, dx2;
  double *dvp, *dvm, *dx;
  double **Lv, **Rv, *lambda;
  double *R, *L, *P, *M;
  double b0, S0, wp, tau;
  double b1, S1, wm;
  double dwp[NVAR], dwp_lim[NVAR];
  double dwm[NVAR], dwm_lim[NVAR];
  double dvpR, dvmR;
  static double **Rg, **Lg, **Pg, **Mg; /* -- interpolation coeffs -- */
  static double **dv;

/* -----------------------------------------------------
   0. Allocate memory and set pointer shortcuts
   ----------------------------------------------------- */
   
  if (dv == NULL) {
    dv = ARRAY_2D(NMAX_POINT, NVAR, double);
    Rg = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    Lg = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    Pg = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    Mg = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    WENO3_COEFF(Rg, Lg, Pg, Mg, grid);
  }

#if RECONSTRUCT_4VEL
  ConvertTo4vel (stateC->v, beg-1, end+1);
#endif

  R  = Rg[g_dir]; L = Lg[g_dir];
  P  = Pg[g_dir]; M = Mg[g_dir];
  dx = grid->dx[g_dir];

/* ----------------------------------------------------
   1. compute slopes and left and right interface values 
   ---------------------------------------------------- */

#if CHAR_LIMITING == NO  /* ----------------------------------------
                            a: Limiter on primitive variables
                            ----------------------------------------  */
  for (i = beg-1; i <= end; i++){
    NVAR_LOOP(nv)dv[i][nv] = v[i+1][nv] - v[i][nv];
  }

  for (i = beg; i <= end; i++){
    dvp = dv[i];  
    dvm = dv[i-1];

    #if SHOCK_FLATTENING == MULTID    
    if (sweep->flag[i] & FLAG_MINMOD){  
      for (nv = NVAR; nv--;    ) {
        dmm = MINMOD(dvp[nv], dvm[nv]);
        vp[i][nv] = v[i][nv] + 0.5*dmm;
        vm[i][nv] = v[i][nv] - 0.5*dmm;
      }
      continue;
    }
    #endif

    dx2 = dx[i]*dx[i];
    for (nv = 0; nv < NVAR; nv++){
      b0 = dvp[nv]*dvp[nv] + dx2;
      b1 = dvm[nv]*dvm[nv] + dx2;

      tau = dvp[nv] - dvm[nv];
      tau = tau*tau;

      S0 = 1.0 + tau/b0;
      S1 = 1.0 + tau/b1;
/*
S0 = 1.0/(b0 + 1.e-6)/(b0 + 1.e-6);
S1 = 1.0/(b1 + 1.e-6)/(b1 + 1.e-6);
       vp[nv] = v[i][nv] + (S0*dvp[nv] + 0.5*S1*dvm[nv])/(2.0*S0 + S1);
       vm[nv] = v[i][nv] - (S1*dvm[nv] + 0.5*S0*dvp[nv])/(2.0*S1 + S0);
*/

      vp[i][nv] = v[i][nv] + (S0*R[i]*dvp[nv] + P[i]*S1*R[i-1]*dvm[nv])
                            /(S0 + P[i]*S1);
      vm[i][nv] = v[i][nv] - (M[i]*S0*L[i]*dvp[nv] + S1*L[i-1]*dvm[nv])
                            /(M[i]*S0 + S1);
    }
  }       
  
#else       /* --------------------------------------------
               b: Limiter on characteristic variables
               --------------------------------------------  */

  SoundSpeed2 (stateC, beg, end, CELL_CENTER, grid);
  PrimEigenvectors (stateC, beg, end);
  i = beg-1;
  NVAR_LOOP(nv) dv[i][nv] = v[i+1][nv] - v[i][nv];

  for (i = beg; i <= end; i++){

    NVAR_LOOP(nv) dv[i][nv] = v[i+1][nv] - v[i][nv];
    Lv     = stateC->Lp[i];
    Rv     = stateC->Rp[i]; 
    lambda = stateC->lambda[i];
    dvp    = dv[i];   
    dvm    = dv[i-1];
    
  /* -------------------------------
      project undivided differences 
      onto characteristic space
     ------------------------------- */
     
    PrimToChar(Lv, dvp, dwp);
    PrimToChar(Lv, dvm, dwm);

    #if SHOCK_FLATTENING == MULTID    
    if (sweep->flag[i] & FLAG_MINMOD){  
      for (nv = NVAR; nv--;    ) {
        dmm = MINMOD(dvp[nv], dvm[nv]);
        vp[i][nv] = v[i][nv] + 0.5*dmm;
        vm[i][nv] = v[i][nv] - 0.5*dmm;
      }
      continue;
    }
    #endif

    dx2 = dx[i]*dx[i];
    for (k = NVAR; k--; ){
      b0 = dwp[k]*dwp[k] + dx2;
      b1 = dwm[k]*dwm[k] + dx2;

      tau = dwp[k] - dwm[k];
      tau = tau*tau;

      S0 = 1.0 + tau/b0;
      S1 = 1.0 + tau/b1;

      dwp_lim[k] = (S0*dwp[k] + 0.5*S1*dwm[k])/(2.0*S0 + S1);
      dwm_lim[k] = (S1*dwm[k] + 0.5*S0*dwp[k])/(2.0*S1 + S0);
/*
    dwp_lim[k] = (S0*R[i]*dwp[k] + P[i]*S1*R[i-1]*dwm[k])/(S0 + P[i]*S1);
    dwm_lim[k] = (M[i]*S0*L[i]*dwp[nv] + S1*L[i-1]*dwm[nv])/(M[i]*S0 + S1);
*/

    }
    for (nv = NFLX; nv--;   ){
      dvpR = dvmR = 0.0;
      #ifdef STAGGERED_MHD
      if (nv == BXn) continue;
      #endif
      for (k = NFLX; k--; ){
        dvpR += dwp_lim[k]*Rv[nv][k];
        dvmR += dwm_lim[k]*Rv[nv][k];
      }
      vp[i][nv] = v[i][nv] + dvpR;
      vm[i][nv] = v[i][nv] - dvmR;
    }

  /* -------------------------------------------------- 
      Passive scalars are treated in a separate loop
      since the characteristic structure is simpler,
      that is, L=R=diag(1) and the characteristic 
      variable coincides with the primitive one.
     -------------------------------------------------- */            

     #if NFLX != NVAR
     for (nv = NFLX; nv < NVAR; nv++ ){
       vp[i][nv] = v[i][nv] + dwp_lim[nv];
       vm[i][nv] = v[i][nv] - dwm_lim[nv];
     }
     #endif

  }
#endif /* CHAR_LIMITING == YES */
   
/* --------------------------------------------------
   2. Ensure positivity of density and pressure 
   -------------------------------------------------- */

  for (i = beg; i <= end; i++){
    dvp = dv[i];  
    dvm = dv[i-1]; 
    if (vp[i][RHO] < 0.0 || vm[i][RHO] < 0.0){
      for (nv = NFLX; nv--;  ){
        dmm = MINMOD(dvp[RHO], dvm[RHO]);
        vp[i][RHO] = v[i][RHO] + 0.5*dmm;
        vm[i][RHO] = v[i][RHO] - 0.5*dmm;  
      }
    }

    #if HAVE_ENERGY
    if (vp[i][PRS] < 0.0 || vm[i][PRS] < 0.0){
      for (nv = NFLX; nv--;  ){
        dmm = MINMOD(dvp[PRS], dvm[PRS]);
        vp[i][PRS] = v[i][PRS] + 0.5*dmm;
        vm[i][PRS] = v[i][PRS] - 0.5*dmm;  
      }
    }
    #endif
    #if ENTROPY_SWITCH
    if (vp[i][ENTR] < 0.0 || vm[i][ENTR] < 0.0){
      dmm = MINMOD(dvp[ENTR], dvm[ENTR]);
      vp[i][ENTR] = v[i][ENTR] + 0.5*dmm;
      vm[i][ENTR] = v[i][ENTR] - 0.5*dmm;
    }
    #endif      

  /* -- Relativistic Limiter -- */

    #if PHYSICS == RHD || PHYSICS == RMHD
    VelocityLimiter (v[i], vp[i], vm[i]);
    #endif
  }

/* -------------------------------------------
   3. Shock flattening 
   -------------------------------------------  */

#if SHOCK_FLATTENING == ONED 
  Flatten (sweep, beg, end, grid);
#endif

/* -------------------------------------------
   4. Assign face-centered magnetic field
   -------------------------------------------  */

#ifdef STAGGERED_MHD
  for (i = beg - 1; i <= end; i++) {
    stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
  }
#endif

/* --------------------------------------------------------
   5. Evolve L/R states and center value by dt/2
   -------------------------------------------------------- */

#if TIME_STEPPING == CHARACTERISTIC_TRACING
  CharTracingStep(sweep, beg, end, grid);
#elif TIME_STEPPING == HANCOCK
  HancockStep(sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   6. Convert back to 3-velocity
   -------------------------------------------------------- */

#if RECONSTRUCT_4VEL
  ConvertTo3vel (v, beg-1, end+1);
  ConvertTo3vel (vp, beg, end);
  ConvertTo3vel (vm, beg, end);  
#endif

/* -----------------------------------------------
   7. Obtain L/R states in conservative variables
   ----------------------------------------------- */

  PrimToCons (vp, up, beg, end);
  PrimToCons (vm, um, beg, end);
}

/* ************************************************************************* */
void WENO3_COEFF(double **R, double **L, double **P, double **M, Grid *grid)
/*
 *
 *
 *
 *
 **************************************************************************** */
{
  int  d, i, beg, end;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  double *x1p = grid->xr[IDIR];
  double *x2p = grid->xr[JDIR];
  double *x3p = grid->xr[KDIR];
  double *x1m = grid->xl[IDIR];
  double *x2m = grid->xl[JDIR];
  double *x3m = grid->xl[KDIR];

  static double *dV;

/* ---------------------------------------------
   0. Define 1D volume coordinate
   --------------------------------------------- */

  if (dV == NULL){
    dV = ARRAY_1D(NMAX_POINT, double);
  }
  
  for (d = 0; d < DIMENSIONS; d++){

    beg = 0; end = grid->np_tot[d] - 1;
    #if GEOMETRY == CARTESIAN
    for (i = beg; i <= end; i++){
      dV[i] = grid->dx[d][i];
    }
    #elif GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    if (d == IDIR) {
      for (i = beg; i <= end; i++){
        dV[i] = fabs(x1[i])*grid->dx[d][i];
      }
    }else if (d == JDIR || d == KDIR){
      for (i = beg; i <= end; i++){
        dV[i] = grid->dx[d][i];
      }
    }
    #elif GEOMETRY == SPHERICAL
    if (d == IDIR) {
      for (i = beg; i <= end; i++){
        dV[i] = (x1p[i]*x1p[i]*x1p[i] - x1m[i]*x1m[i]*x1m[i])/3.0;
      }
    }else if (d == JDIR){
      for (i = beg; i <= end; i++){
        dV[i] = fabs(cos(x2m[i] - x2p[i]));
      }
    }else if (d == KDIR){
      for (i = beg; i <= end; i++){
        dV[i] = grid->dx[d][i];
      }
    }
    #endif
    
    beg = 0; end = grid->np_tot[d] - 2;
    for (i = beg; i <= end; i++){
      R[d][i] = dV[i + 1]/(dV[i] + dV[i + 1]);
      L[d][i] = 1.0 - R[d][i];
    }
    beg = 1; end = grid->np_tot[d] - 2;
    for (i = beg; i <= end; i++){
      P[d][i] = dV[i + 1]/(dV[i] + dV[i - 1]);
      M[d][i] = dV[i - 1]/(dV[i] + dV[i + 1]);
    }
  }
}

