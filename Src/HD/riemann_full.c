#include "pluto.h"

#define MAX_ITER   20
#define small_p       1.e-9
#define small_rho     1.e-9

static void PFUN (real p, real *vL, real *vR, real *f, real *df);
static void FUN_LR (real p, real *v, real *fLR, real *dfLR);

/* ***************************************************************************** */
void TWO_SHOCK (const Sweep *sweep, real *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   RIEMANN
 *
 *
 * PURPOSE
 *
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 ******************************************************************************* */
{
  int  nv, k, i, beg, end;
  real pstar, ustar, dp, fp, dfp;
  real th, q, scrh, gmmp, gp1, gm1;
  real *vL, *vR;
  real fL, dfL, SL, STL, csL;
  real fR, dfR, SR, STR, csR;
  static real *s, **vs, **us, *cmax_loc;
  static int *shock;

  if (vs == NULL){
    vs       = array_2D(NMAX_POINT, NVAR);
    us       = array_2D(NMAX_POINT, NVAR);
    s        = array_1D(NMAX_POINT);
    cmax_loc = array_1D(NMAX_POINT);
    shock    = int_array_1D(NMAX_POINT);
  }
 
  beg = grid[g_dir].lbeg - 1;
  end = grid[g_dir].lend;

  gm1  = g_gamma - 1.0;
  gp1  = g_gamma + 1.0;
  gmmp = gm1/gp1;

  for (i = beg; i <= end; i++){

    vL = sweep->vL[i];
    vR = sweep->vR[i];
    s[i] = 0.0;

  /* -- guess here -- */

    pstar = 0.5*(vL[PRS] + vR[PRS]);
    
    for (k = 0; k < MAX_ITER; k++){

      PFUN (pstar, vL, vR, &fp, &dfp);
      dp     = fp/dfp;
      pstar -= dp;

      if (fabs(dp) < 1.e-7*pstar) break;
      if (k == (MAX_ITER-5)){
        print ("! Too many iterations in Rieman\n");
        Show(sweep->vL,i);
        Show(sweep->vR,i);

        QUIT_PLUTO(1);
      }
    }

    FUN_LR (pstar, vL, &fL, &dfL);
    FUN_LR (pstar, vR, &fR, &dfR);

    ustar = 0.5*(vL[VXn] + vR[VXn] + fR - fL);

  /* -- sample solution -- */

    if (s[i] <= ustar){  /* -- left of CD -- */
      q   = pstar/vL[PRS];
      csL = sqrt(g_gamma*vL[PRS]/vL[RHO]);

      if (q > 1.0) { /* -- left wave is a shock -- */

        scrh = gp1*q + gm1;          
        SL   = vL[VXn] - csL*sqrt(0.5/g_gamma*scrh);
        if (s[i] < SL){
          for (nv = NVAR; nv--; ) vs[i][nv] = vL[nv];
        }else { 
          vs[i][RHO] = vL[RHO]*(q + gmmp)/(gmmp*q + 1.0);
          vs[i][VXn] = ustar;
          vs[i][PRS] = pstar;
        }

      }else{  /* -- left wave is a rarefaction -- */

        SL = vL[VXn] - csL;
 
        if (s[i] < SL) {
          for (nv = NVAR; nv--; ) vs[i][nv] = vL[nv];
        }else { 
          vs[i][RHO] = vL[RHO]*pow(q, 1.0/g_gamma);
          vs[i][VXn] = ustar;
          vs[i][PRS] = pstar;
          STL = ustar - sqrt(g_gamma*pstar/vs[i][RHO]);
          if (s[i] < STL){ /* -- sol inside rarefaction -- */
            scrh = 2.0 + gm1/csL*(vL[VXn] - s[i]);
            vs[i][RHO] = vL[RHO]*pow(scrh/gp1, 2.0/gm1);
            vs[i][PRS] = vL[PRS]*pow(scrh/gp1, 2.0*g_gamma/gm1);
            vs[i][VXn] = 2.0/gp1*(csL + 0.5*gm1*vL[VXn] + s[i]);
          }
        }
      } 

    }else{  /* -- right of CD -- */

      q   = pstar/vR[PRS];
      csR = sqrt(g_gamma*vR[PRS]/vR[RHO]);

      if (q > 1.0) { /* -- right wave is a shock -- */

        scrh = gp1*q + gm1;          
        SR   = vR[VXn] + csR*sqrt(0.5/g_gamma*scrh);
        if (s[i] > SR){
          for (nv = NVAR; nv--; ) vs[i][nv] = vR[nv];
        }else { 
          vs[i][RHO] = vR[RHO]*(q + gmmp)/(gmmp*q + 1.0);
          vs[i][VXn] = ustar;
          vs[i][PRS] = pstar;
        }

      }else{  /* -- right wave is a rarefaction -- */

        SR = vR[VXn] + csR;
 
        if (s[i] > SR) {
          for (nv = NVAR; nv--; ) vs[i][nv] = vR[nv];
        }else { 
          vs[i][RHO] = vR[RHO]*pow(q, 1.0/g_gamma);
          vs[i][VXn] = ustar;
          vs[i][PRS] = pstar;
          STR = ustar + sqrt(g_gamma*pstar/vs[i][RHO]);
          if (s[i] > STR){ /* -- sol inside rarefaction -- */
            scrh = 2.0 - gm1/csR*(vR[VXn] - s[i]);
            vs[i][RHO] = vR[RHO]*pow(scrh/gp1, 2.0/gm1);
            vs[i][PRS] = vR[PRS]*pow(scrh/gp1, 2.0*g_gamma/gm1);
            vs[i][VXn] = 2.0/gp1*(-csR + 0.5*gm1*vR[VXn] + s[i]);
          }
        }
      } 
    }   

    if (ustar > 0.0) {
      EXPAND(                     ,
             vs[i][VXt] = vL[VXt];  ,
             vs[i][VXb] = vL[VXb];)
    }else{
      EXPAND(                     ,
             vs[i][VXt] = vR[VXt];  ,
             vs[i][VXb] = vR[VXb];)
    }
      
  }

/* -- compute fluxes -- */

  MAX_CH_SPEED (vs, cmax_loc, grid, beg, end);
  for (i = beg; i <= end; i++) {
    *cmax = MAX(cmax_loc[i], *cmax);
  }
  FLUX(us, vs, sweep->flux, sweep->press, beg, end);

}

/* ************************************************ */
void PFUN (real p, real *vL, real *vR, real *f, real *df)
/*
 *
 *
 *
 *
 ************************************************** */
{
  real fL , fR; 
  real dfL, dfR;
  
  FUN_LR (p, vL, &fL, &dfL);
  FUN_LR (p, vR, &fR, &dfR);

  *f  = fL  + fR + vR[VXn] - vL[VXn];
  *df = dfL + dfR; 
}

/* ************************************************ */
void FUN_LR (real p, real *v, real *fLR, real *dfLR)
/*
 *
 *
 *
 *
 ************************************************** */
{
  real A, B, scrh, cs;
  real q;

  if (p > v[PRS]) {  /* -- (shock) -- */

    A = 2.0/(g_gamma + 1.0)/v[RHO];
    B = (g_gamma - 1.0)/(g_gamma + 1.0)*v[PRS];
    scrh  = A/(p + B);
    *fLR  = (p - v[PRS])*sqrt(scrh);
    *dfLR = sqrt(scrh)*(1.0 - 0.5*(p - v[PRS])/(B + p));

  }else{   /* -- (rarefaction) -- */

    cs = sqrt(g_gamma*v[PRS]/v[RHO]);
    q  = p/v[PRS];
    scrh = pow(q, 0.5*(g_gamma - 1.0)/g_gamma) - 1.0;
    *fLR  = 2.0*cs/(g_gamma - 1.0)*scrh;
    scrh  = pow(q, -0.5*(g_gamma + 1.0)/g_gamma);
    *dfLR = 1.0/(v[RHO]*cs)*scrh;
  }

}


