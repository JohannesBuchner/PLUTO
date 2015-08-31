/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Two-shock Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the two-shock approximation.
  
      Reference:    
   - On input, it takes left and right primitive state
     vectors state->vL and state->vR at zone edge i+1/2;
     On output, return flux and pressure vectors at the
     same interface.
 
   - Also, compute maximum wave propagation speed (cmax) 
     for  explicit time step computation
   
 
  LAST_MODIFIED
 
    April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 
   
  \b Reference:
   -  "FLASH: an Adaptive Mesh Hydrodynamics Code for Modeling
      Astrophysical Thermonuclear Flashes"
      Fryxell et al, ApJS(2000) 131:273 (pages 292-294)
        
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 5, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER   5
#define small_p       1.e-9
#define small_rho     1.e-9

/* ***************************************************************************** */
void TwoShock_Solver (const State_1D *state, int beg, int end, 
                      double *cmax, Grid *grid)
/*!
 *
 ******************************************************************************* */
{
#if EOS == IDEAL
  int i, iter, nv;

  double   vxl, taul, cl, *ql;
  double   vxr, taur, cr, *qr;
  double   zs, taus, cs, *qs;
  double   pstar, ustar, rho_star, cstar;
  double   sigma, lambda_s, lambda_star, zeta, dp;
  double   g1_g, scrh1, scrh2, scrh3, scrh4;
  static double  **ws, **us;
  static double **fL, **fR, *pL, *pR, *a2L, *a2R;
  double *uL, *uR;

  if (ws == NULL){
    ws = ARRAY_2D(NMAX_POINT, NVAR, double);  
    us = ARRAY_2D(NMAX_POINT, NVAR, double);  

    fL = ARRAY_2D(NMAX_POINT, NFLX, double);  
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);  
    pL = ARRAY_1D(NMAX_POINT, double);  
    pR = ARRAY_1D(NMAX_POINT, double);  

    a2L = ARRAY_1D(NMAX_POINT, double);  
    a2R = ARRAY_1D(NMAX_POINT, double);  
  }

/*  ---------------------------------------------------------------
                       SOLVE RIEMANN PROBLEM
    ---------------------------------------------------------------   */

  g1_g = 0.5*(g_gamma + 1.0)/g_gamma;
  for (i = beg; i <= end; i++) {

    uL = state->uL[i];
    uR = state->uR[i];

#if SHOCK_FLATTENING == MULTID   
    if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){        
      a2L[i] = g_gamma*state->vL[i][PRS]/state->vL[i][RHO];
      a2R[i] = g_gamma*state->vR[i][PRS]/state->vR[i][RHO];
      HLL_Speed (state->vL, state->vR, a2L, a2R, &cl - i, &cr - i, i, i);
      Flux (state->uL, state->vL, a2L, fL, pL, i, i);
      Flux (state->uR, state->vR, a2R, fR, pR, i, i);

      cs    = MAX(fabs(cl), fabs(cr));
      cmax[i] = cs;
      cl    = MIN(0.0, cl);
      cr    = MAX(0.0, cr);
      scrh1  = 1.0/(cr - cl);
      for (nv = NFLX; nv--; ){
        state->flux[i][nv]  = cl*cr*(uR[nv] - uL[nv])
                           +  cr*fL[i][nv] - cl*fR[i][nv];
        state->flux[i][nv] *= scrh1;
      }
      state->press[i] = (cr*pL[i] - cl*pR[i])*scrh1;
      continue;
    }
#endif
  
    qr = state->vR[i];
    ql = state->vL[i];
      
    cl = sqrt(g_gamma*ql[PRS]*ql[RHO]);
    cr = sqrt(g_gamma*qr[PRS]*qr[RHO]);

    taul = 1.0/ql[RHO];
    taur = 1.0/qr[RHO];

  /* -- First guess -- */

    pstar = qr[PRS] - ql[PRS] - cr*(qr[VXn] - ql[VXn]);
    pstar = ql[PRS] + pstar*cl/(cl + cr);
    pstar = MAX(small_p , pstar);

  /* -- Begin to iterate --  */

    for (iter = 1; iter <= MAX_ITER; iter++) {
      vxl = cl*sqrt(1.0 + g1_g*(pstar - ql[PRS])/ql[PRS]);
      vxr = cr*sqrt(1.0 + g1_g*(pstar - qr[PRS])/qr[PRS]);

      scrh1 = vxl*vxl;
      scrh1 = 2.0*scrh1*vxl/(scrh1 + cl*cl);
      
      scrh2 = vxr*vxr;
      scrh2 = 2.0*scrh2*vxr/(scrh2 + cr*cr);

      scrh3 = ql[VXn] - (pstar - ql[PRS])/vxl;
      scrh4 = qr[VXn] + (pstar - qr[PRS])/vxr;

      dp = scrh1*scrh2/(scrh1 + scrh2)*(scrh4 - scrh3);

      pstar -= dp;
/*
      scrh1 = 4.0*taul*vxl*vxl;
      scrh1 = -scrh1*vxl/(scrh1 - (g_gamma + 1.0)*(pstar - ql[PRS]));
      scrh2 = 4.0*taur*vxr*vxr;
      scrh2 =  scrh2*vxr/(scrh2 - (g_gamma + 1.0)*(pstar - qr[PRS]));

      scrh3 = ql[VXn] - (pstar - ql[PRS])/vxl;
      scrh4 = qr[VXn] + (pstar - qr[PRS])/vxr;
      dp    = (scrh4 - scrh3)*(scrh1*scrh2)/(scrh2 - scrh1);
      pstar = pstar  + dp;
*/

      pstar = MAX (small_p, pstar);
      g_maxRiemannIter = MAX(g_maxRiemannIter,iter);
      if (fabs(dp/pstar) < 1.e-6) break;
/*
      if (iter == MAX_ITER) {
        print ("Rieman solver not converging  ps %12.6e  dp %12.6e %12.6e  %12.6e %12.6e %12.6e \n",
           pstar,dp, scrh1, scrh2, scrh3, scrh4);
      }
*/
    }

/*            End iterating           */

    scrh3 = ql[VXn] - (pstar - ql[PRS])/vxl;
    scrh4 = qr[VXn] + (pstar - qr[PRS])/vxr;
    ustar = 0.5*(scrh3 + scrh4);

    if (ustar > 0.0) {
      sigma = 1.0;
      taus  = taul;
      cs    = cl*taul;
      zs    = vxl;
      qs    = ql;
    }else{
      sigma = -1.0;
      taus  = taur;
      cs    = cr*taur;
      zs    = vxr;
      qs    = qr;
    }
        
    rho_star = taus - (pstar - qs[PRS])/(zs*zs);
    rho_star = MAX(small_rho,1.0/rho_star);
    cstar    = sqrt(g_gamma*pstar/rho_star);
     if (pstar < qs[PRS]){
      lambda_s    = cs    - sigma*qs[VXn];
      lambda_star = cstar - sigma*ustar;
    }else{
      lambda_s    = lambda_star = zs*taus  - sigma*qs[VXn];
    }      

  /* -- Now find solution -- */

    if (lambda_star > 0.0){
      
      ws[i][RHO] = rho_star;
      ws[i][VXn] = ustar;  
      ws[i][PRS] = pstar;
                
    } else if (lambda_s < 0.0){
      
      ws[i][RHO] = qs[RHO];
      ws[i][VXn] = qs[VXn];
      ws[i][PRS] = qs[PRS];

    } else {  /*   linearly interpolate rarefaction fan  */

      scrh1 = MAX(lambda_s - lambda_star, lambda_s + lambda_star);
      zeta  = 0.5*(1.0 + (lambda_s + lambda_star)/scrh1);

      ws[i][RHO] = zeta*rho_star + (1.0 - zeta)*qs[RHO];
      ws[i][VXn] = zeta*ustar    + (1.0 - zeta)*qs[VXn];
      ws[i][PRS] = zeta*pstar    + (1.0 - zeta)*qs[PRS];
    }  
        
  /* -- transverse velocities are advected --  */

    EXPAND(                    , 
           ws[i][VXt] = qs[VXt]; ,  
           ws[i][VXb] = qs[VXb];)

  /* -- compute flux -- */

    PrimToCons (ws, us, i, i);
    scrh2 = g_gamma*ws[i][PRS]/ws[i][RHO];
    Flux (us, ws, &scrh2 - i, state->flux, state->press, i, i);
    cstar = sqrt(scrh2);

  /* -- compute max speed -- */

    scrh1 = fabs(ws[i][VXn])/cstar;
    g_maxMach = MAX(scrh1, g_maxMach);
    scrh1 = fabs(ws[i][VXn]) + cstar;
    cmax[i] = scrh1;
  
  /* -- Add artificial viscosity -- */

#ifdef ARTIFICIAL_VISC
    scrh1 = ARTIFICIAL_VISC*(MAX(0.0, ql[VXn] - qr[VXn]));
    for (nv = 0; nv < NFLX; nv++) {
      state->flux[i][nv] += scrh1*(uL[nv] - uR[nv]);
    }
#endif    
    
  }

  
#else
  print1 ("! TwoShock_Solver: not defined for this EOS\n");
  QUIT_PLUTO(1);
#endif /* EOS == IDEAL */
}

#undef MAX_ITER   
#undef small_p     
#undef small_rho     
