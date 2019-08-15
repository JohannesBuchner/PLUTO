/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Two-shock Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the two-shock approximation.
  
  - On input, it takes left and right primitive state
    vectors stateL->v and stateR->v at zone edge i+1/2;
    On output, return flux and pressure vectors at the
    same interface.
 
  - Also, compute maximum wave propagation speed (cmax) 
    for  explicit time step computation   
 
  \b Reference:
   -  "FLASH: an Adaptive Mesh Hydrodynamics Code for Modeling
      Astrophysical Thermonuclear Flashes"
      Fryxell et al, ApJS(2000) 131:273 (pages 292-294)
        
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Jan 26, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER   5
#define small_p       1.e-9
#define small_rho     1.e-9

/* ***************************************************************************** */
void TwoShock_Solver (const Sweep *sweep, int beg, int end, 
                      double *cmax, Grid *grid)
/*!
 *
 ******************************************************************************* */
{
#if EOS == IDEAL
  int i, iter, nv;

  static State stateS;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double  vxl, taul, cl, *ql;
  double  vxr, taur, cr, *qr;
  double  zs, taus, cs, *qs;
  double  pstar, ustar, rho_star, cstar;
  double  sigma, lambda_s, lambda_star, zeta, dp;
  double  g1_g, scrh1, scrh2, scrh3, scrh4;
  double  *uL, *uR;

  if (stateS.v == NULL){
    StateStructAllocate (&stateS);
  }

/*  ---------------------------------------------------------------
                       SOLVE RIEMANN PROBLEM
    ---------------------------------------------------------------   */

  g1_g = 0.5*(g_gamma + 1.0)/g_gamma;
  for (i = beg; i <= end; i++) {

    uL = stateL->u[i];
    uR = stateR->u[i];

#if SHOCK_FLATTENING == MULTID   
    if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){        
      stateL->a2[i] = g_gamma*stateL->v[i][PRS]/stateL->v[i][RHO];
      stateR->a2[i] = g_gamma*stateR->v[i][PRS]/stateR->v[i][RHO];
      HLL_Speed (stateL, stateR, &cl - i, &cr - i, i, i);
      Flux (stateL, i, i);
      Flux (stateR, i, i);

      cs    = MAX(fabs(cl), fabs(cr));
      cmax[i] = cs;
      cl    = MIN(0.0, cl);
      cr    = MAX(0.0, cr);
      scrh1  = 1.0/(cr - cl);
      for (nv = NFLX; nv--; ){
        sweep->flux[i][nv]  = cl*cr*(uR[nv] - uL[nv])
                           +  cr*stateL->flux[i][nv] - cl*stateR->flux[i][nv];
        sweep->flux[i][nv] *= scrh1;
      }
      sweep->press[i] = (cr*stateL->prs[i] - cl*stateR->prs[i])*scrh1;
      continue;
    }
#endif
  
    qr = stateR->v[i];
    ql = stateL->v[i];
      
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
    } /* -- End iteration -- */

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

  /* -- Now find solution in primitive variables -- */

    if (lambda_star > 0.0){
      
      stateS.v[i][RHO] = rho_star;
      stateS.v[i][VXn] = ustar;  
      stateS.v[i][PRS] = pstar;
                
    } else if (lambda_s < 0.0){
      
      stateS.v[i][RHO] = qs[RHO];
      stateS.v[i][VXn] = qs[VXn];
      stateS.v[i][PRS] = qs[PRS];

    } else {  /*   linearly interpolate rarefaction fan  */

      scrh1 = MAX(lambda_s - lambda_star, lambda_s + lambda_star);
      scrh1 = MAX(1.e-12,scrh1); /* AM (26/01/2017): avoid division by 0.0 at
                                    sonic rarefactions when
                                    lambda_s = lambda_star = 0.0 */
      zeta  = 0.5*(1.0 + (lambda_s + lambda_star)/scrh1);

      stateS.v[i][RHO] = zeta*rho_star + (1.0 - zeta)*qs[RHO];
      stateS.v[i][VXn] = zeta*ustar    + (1.0 - zeta)*qs[VXn];
      stateS.v[i][PRS] = zeta*pstar    + (1.0 - zeta)*qs[PRS];
    }  
        
  /* -- transverse velocities are advected --  */

    EXPAND(                             , 
           stateS.v[i][VXt] = qs[VXt]; ,  
           stateS.v[i][VXb] = qs[VXb];)

  /* -- compute flux -- */

    PrimToCons (stateS.v, stateS.u, i, i);
    stateS.a2[i] = g_gamma*stateS.v[i][PRS]/stateS.v[i][RHO];
    Flux (&stateS, i,i);
    for (nv = 0; nv < NFLX; nv++) {
      sweep->flux[i][nv] = stateS.flux[i][nv];
    }
    sweep->press[i] = stateS.prs[i];
    
    cstar = sqrt(stateS.a2[i]);

  /* -- compute max speed -- */

    scrh1 = fabs(stateS.v[i][VXn])/cstar;
    g_maxMach = MAX(scrh1, g_maxMach);
    scrh1 = fabs(stateS.v[i][VXn]) + cstar;
    cmax[i] = scrh1;
  
  /* -- Add artificial viscosity -- */

#ifdef ARTIFICIAL_VISC
    scrh1 = ARTIFICIAL_VISC*(MAX(0.0, ql[VXn] - qr[VXn]));
    for (nv = 0; nv < NFLX; nv++) {
      sweep->flux[i][nv] += scrh1*(uL[nv] - uR[nv]);
    }
#endif

  }

  
#else
  print ("! TwoShock_Solver: not defined for this EOS\n");
  QUIT_PLUTO(1);
#endif /* EOS == IDEAL */
}

#undef MAX_ITER   
#undef small_p     
#undef small_rho     
