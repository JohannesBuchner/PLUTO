/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief LimO3 reconstruction.

  Provide a three-point stencil, third-order reconstruction algorithm 
  based on the limiter function of Cada & Torrilhon

  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 10, 2016

  \b References
     - "Compact third-order limiter functions for finite volume
        methods", Cada & Torrilhon, JCP (2009) 228, 4118.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static double LimO3Func (double, double, double);

/* ********************************************************************* */
void States (const Sweep *sweep, int beg, int end, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
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

  double dmm;
  double *dvp, *dvm, *dx;
  double **L, **R, *lambda;
  double dwp[NVAR], dwp_lim[NVAR];
  double dwm[NVAR], dwm_lim[NVAR];
  double dvpR, dvmR;
  static double **dv;

/* --------------------------------------------------------
   0. Allocate memory, set pointer shortcuts 
   -------------------------------------------------------- */

  if (dv == NULL) {
    dv = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  v  = stateC->v;
  vp = stateL->v;
  vm = stateR->v-1;

#if RECONSTRUCT_4VEL
  ConvertTo4vel (stateC->v, beg-1, end+1);
#endif
   
  dx = grid->dx[g_dir];

/* --------------------------------------------------------
   1. compute slopes and left and right interface values 
   -------------------------------------------------------- */

  #if CHAR_LIMITING == NO  /* ----------------------------------------
                                 Limiter on primitive variables
                              ----------------------------------------  */
   for (i = beg-1; i <= end; i++){
   for (nv = NVAR; nv--; ){
     dv[i][nv] = v[i+1][nv] - v[i][nv];
   }}

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

     for (nv = NVAR; nv--; ){
       vp[i][nv] = v[i][nv] + 0.5*dvp[nv]*LimO3Func(dvp[nv], dvm[nv], dx[i]);
       vm[i][nv] = v[i][nv] - 0.5*dvm[nv]*LimO3Func(dvm[nv], dvp[nv], dx[i]);
     }
   }       
  
  #else       /* --------------------------------------------
                    Limiter on characteristic variables
                 --------------------------------------------  */

   SoundSpeed2 (stateC, beg, end, CELL_CENTER, grid);
   PrimEigenvectors (stateC, beg, end);
   i = beg-1;
   NVAR_LOOP(nv) dv[i][nv] = v[i+1][nv] - v[i][nv];

   for (i = beg; i <= end; i++){

     NVAR_LOOP(nv) dv[i][nv] = v[i+1][nv] - v[i][nv];
     L      = stateC->Lp[i];
     R      = stateC->Rp[i];
     lambda = stateC->lambda[i];
     dvp = dv[i];  
     dvm = dv[i-1];
    
  /* -------------------------------
      project undivided differences 
      onto characteristic space
     ------------------------------- */
     
     PrimToChar(L, dvp, dwp);
     PrimToChar(L, dvm, dwm);

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

  /* -----------------------------
      limit undivided differences
     ----------------------------- */

     for (k = NFLX; k--; ){
       dwp_lim[k] = dwp[k]*LimO3Func(dwp[k], dwm[k], dx[i]);
       dwm_lim[k] = dwm[k]*LimO3Func(dwm[k], dwp[k], dx[i]);
     }
     for (nv = NFLX; nv--;   ){
       dvpR = dvmR = 0.0;
       #ifdef STAGGERED_MHD
        if (nv == BXn) continue;
       #endif
       for (k = NFLX; k--; ){
         dvpR += dwp_lim[k]*R[nv][k];
         dvmR += dwm_lim[k]*R[nv][k];
       }
       vp[i][nv] = v[i][nv] + 0.5*dvpR;
       vm[i][nv] = v[i][nv] - 0.5*dvmR;
     }

  /* -------------------------------------- 
      Compute limited slopes for tracers
      exploiting the simple characteristic 
      structure, L=R=diag(1).
     -------------------------------------- */            

    #if NFLX != NVAR
    for (nv = NFLX; nv < NVAR; nv++ ){
      dvpR = dvp[nv]*LimO3Func(dvp[nv], dvm[nv], dx[i]);
      dvmR = dvm[nv]*LimO3Func(dvm[nv], dvp[nv], dx[i]);
      vp[i][nv] = v[i][nv] + 0.5*dvpR;
      vm[i][nv] = v[i][nv] - 0.5*dvmR;
    }
    #endif
   }
   
  #endif /* CHAR_LIMITING == YES */
   
/* --------------------------------------------------------
   2. Since the third-order limiter is not TVD, we 
      need to ensure positivity of density and pressure 
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++){
    dvp = dv[i];   
    dvm = dv[i-1]; 
    if (vp[i][RHO] < 0.0 || vm[i][RHO] < 0.0){
      dmm = MINMOD(dvp[RHO], dvm[RHO]);
      vp[i][RHO] = v[i][RHO] + 0.5*dmm;
      vm[i][RHO] = v[i][RHO] - 0.5*dmm;
    }

    #if HAVE_ENERGY
    if (vp[i][PRS] < 0.0 || vm[i][PRS] < 0.0){
      dmm = MINMOD(dvp[PRS], dvm[PRS]);
      vp[i][PRS] = v[i][PRS] + 0.5*dmm;
      vm[i][PRS] = v[i][PRS] - 0.5*dmm;
    }
    #endif
    #if ENTROPY_SWITCH
    if (vp[i][ENTR] < 0.0 || vm[i][ENTR] < 0.0){
      dmm = MINMOD(dvp[ENTR], dvm[ENTR]);
      vp[i][ENTR] = v[i][ENTR] + 0.5*dmm;
      vm[i][ENTR] = v[i][ENTR] - 0.5*dmm;
    }
    #endif      

  /* -- relativistic limiter --*/

    #if PHYSICS == RHD || PHYSICS == RMHD
    VelocityLimiter (v[i], vp[i], vm[i]);
    #endif
  }

/* --------------------------------------------------------
    3. Shock flattening 
   -------------------------------------------------------- */

#if SHOCK_FLATTENING == ONED 
  Flatten (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   4. Assign face-centered magnetic field
   -------------------------------------------------------- */

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

/* --------------------------------------------------------
   7. Obtain L/R states in conservative variables
   -------------------------------------------------------- */

  PrimToCons (vp, up, beg, end);
  PrimToCons (vm, um, beg, end);
}

/* ********************************************************************* */
double LimO3Func (double dvp, double dvm, double dx)
/*!
 *  Implement the 3rd-order limiter function, 
 *  Eq. [3.22]
 *
 *
 *********************************************************************** */
{
  double r = 0.1;
  double a,b,c,q, th, lim;
  double eta, psi, eps = 1.e-12;

  th  = dvm/(dvp + 1.e-16);

  q = (2.0 + th)/3.0;

  a = MIN(1.5,2.0*th);
  a = MIN(q,a);
  b = MAX(-0.5*th,a);
  c = MIN(q,b);
  psi = MAX(0.0,c);

  eta = r*dx;
  eta = (dvm*dvm + dvp*dvp)/(eta*eta);
  if ( eta <= 1.0 - eps) {
    lim = q;
  }else if (eta >= 1.0 + eps){
    lim = psi;
  }else{
    psi =   (1.0 - (eta - 1.0)/eps)*q
          + (1.0 + (eta - 1.0)/eps)*psi;
    lim = 0.5*psi;
  }
  return (lim);
}
