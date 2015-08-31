/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implementation of the two-shock Riemann solver for the 
         relativistic hydro equations.

  Solve the Riemann problem for the relativistic hydrodynamics 
  equations using the two-shock approach described by Mignone, Plewa 
  and Bodo (2005).
  The formulation works with IDEAL or TAUB equation of state.

  On input, this function takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge \c i+1/2.
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "The Piecewise Parabolic Method for  Multidimensional Relativistic 
       Fluid Dynamics", Mignone, Plewa and Bodo, ApJS (2005) 160,199.
       
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 5, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER   20
#define small_p    1.e-12
#define small_rho  1.e-12
#define INTERPOLATE_RAREFACTION   YES

#define accuracy   1.e-6

static double TwoShock_Lorentz (double *U, int n);
static void   TwoShock_Shock (double, double, double, double, double, double,
                              double, double *, double *, double *, int);
static double TwoShock_RarefactionSpeed (double *u, int side);

static double qglob_r[NFLX], qglob_l[NFLX], gmmr;

/* *********************************************************************  */
void TwoShock_Solver (const State_1D *state, int beg, int end, 
                      double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the relativistic HD equations using the 
 * two-shock Riemann solver of Mignone et al. (2005).
 * 
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int     iter, i, nv;
  int     k, nfail = 0, izone_fail[1024];
  double  Ustar[NFLX];
  double  gL, gR, gL1, gR1;
  double  uxR, uxR1, pR, j2R, wR, VR, VR1;
  double  uxL, uxL1, pL, j2L, wL, VL, VL1;
  double  duR, duL, p1, u1, dp, a0, a1;
  double  tauR, tauL, am, ap;
  double  *ql, *qr, *qs;
  static double  **fl, **us, **ws;
  static double *hR, *hL, *a2L, *a2R, *cmax_loc, *cmin_loc;

  static double **fL, **fR, *prL, *prR;
  double bmin, bmax;
  double SL, SR;

  #if EOS == IDEAL
   gmmr = g_gamma/(g_gamma - 1.0);
  #endif

  if (fl == NULL){
    fl       = ARRAY_2D(NMAX_POINT, NFLX, double); 
    us       = ARRAY_2D(NMAX_POINT, NVAR, double);  
    ws       = ARRAY_2D(NMAX_POINT, NVAR, double);  
    cmax_loc = ARRAY_1D(NMAX_POINT, double);
    cmin_loc = ARRAY_1D(NMAX_POINT, double);

    a2R      = ARRAY_1D(NMAX_POINT, double);
    a2L      = ARRAY_1D(NMAX_POINT, double);
    hR       = ARRAY_1D(NMAX_POINT, double);
    hL       = ARRAY_1D(NMAX_POINT, double);

    fL  = ARRAY_2D(NMAX_POINT, NVAR, double);
    fR  = ARRAY_2D(NMAX_POINT, NVAR, double);
    prL = ARRAY_1D(NMAX_POINT, double);
    prR = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
        compute sound speed at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, hR, beg, end, FACE_CENTER, grid);

/*  ---------------------------------------------------------------
                       SOLVE RIEMANN PROBLEM
    ---------------------------------------------------------------   */

  for (i = beg; i <= end; i++) {

#if SHOCK_FLATTENING == MULTID
    if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){        
      HLL_Speed (state->vL, state->vR, a2L, a2R, &gL - i, &gR - i, i, i);
      Flux  (state->uL, state->vL, a2L, fL, prL, i, i);
      Flux  (state->uR, state->vR, a2R, fR, prR, i, i);

      a0 = MAX(fabs(gR), fabs(gL));
      cmax[i] = a0;

      gL = MIN(0.0, gL);
      gR = MAX(0.0, gR);
      a0 = 1.0/(gR - gL);
      for (nv = NFLX; nv--; ){
        state->flux[i][nv]  = gL*gR*(state->uR[i][nv] - state->uL[i][nv])
                             + gR*fL[i][nv] - gL*fR[i][nv];
        state->flux[i][nv] *= a0;
      }
      state->press[i] = (gR*prL[i] - gL*prR[i])*a0;
      continue;
    }
#endif

    ql = state->vL[i];
    qr = state->vR[i];

    gR = TwoShock_Lorentz (qr,0);
    gL = TwoShock_Lorentz (ql,1);

    tauR = 1.0/qr[RHO];
    tauL = 1.0/ql[RHO];

    uxR = qr[VXn];
    uxL = ql[VXn];

    pR = qr[PRS];
    pL = ql[PRS];

    wR = pR*tauR;
    wL = pL*tauL;

    VR = tauR/gR;
    VL = tauL/gL;

/*  ---------------------------------------------
       First iteration is done out of the loop
    ---------------------------------------------  */

    #if EOS == IDEAL
     j2R = tauR*(hR[i]*(gmmr - 2.0) + 1.0) / (gmmr*pR);
     j2L = tauL*(hL[i]*(gmmr - 2.0) + 1.0) / (gmmr*pL);
    #endif
    #if EOS == TAUB
     a0  = (5.0*hR[i] - 8.0*wR)/(2.0*hR[i] - 5.0*wR);
     j2R = -tauR*(a0*wR + hR[i]*(1.0 - a0))/(a0*pR); 

     a0  = (5.0*hL[i] - 8.0*wL)/(2.0*hL[i] - 5.0*wL);
     j2L = -tauL*(a0*wL + hL[i]*(1.0 - a0))/(a0*pL);
    #endif

    gR1 = sqrt(VR*VR + (1.0 - uxR*uxR)*j2R);   /*    RIGHT    */
    wR  = VR*uxR + gR1;

    gL1 = -sqrt(VL*VL + (1.0 - uxL*uxL)*j2L);  /*    LEFT    */
    wL  = VL*uxL + gL1;

    duR = gR1/(hR[i]*gR);
    duL = gL1/(hL[i]*gL);

    p1  = ql[VXn] - qr[VXn] + duR*pR - duL*pL;
    p1 /= duR - duL;

    if (p1 < 0.0) {
      p1 = MIN(pR,pL);
    }

/*  -----------------------------------------------
             BEGIN iteration  loop   
    -----------------------------------------------  */

    for (iter = 1; iter < MAX_ITER; iter++)  {

      TwoShock_Shock (tauL, uxL, pL, gL, VL, hL[i], p1, &uxL1, &duL, &wL, -1);
      TwoShock_Shock (tauR, uxR, pR, gR, VR, hR[i], p1, &uxR1, &duR, &wR, 1);

  /*  ----  Find  next  approximate solution  ----  */

      dp  = (uxR1 - uxL1)/(duL - duR);
      if (-dp > p1) dp = -0.5*p1;
      p1 += dp;

      if (p1 < 0.0){
        p1 -= dp;
        p1 *= 0.5;
      }
      g_maxRiemannIter = MAX(g_maxRiemannIter,iter);
      if ( (fabs(dp) < accuracy*p1) ) break;
    }

/*  -----------------------------------------------
             END iteration  loop   
    -----------------------------------------------  */

    u1 = 0.5*(uxR1 + uxL1); 
    
  /*  ---- Check possible failures, and flag zones  ----  */

    if (u1 != u1 || iter == MAX_ITER) {

      u1 = 0.5*(uxL + uxR);
      p1 = 0.5*(pL  + pR);
      nfail++;
      izone_fail[nfail] = i;
      for (nv = NFLX; nv--; ){
        ws[i][nv] = ql[nv];
      }
      continue;
    } 

/*  This switch works to prevent shear smearing, by
    filtering noise in the velocity  */

/*
    u1 = (fabs(u1)<1.e-13 ? 0.0:u1);
*/

    Ustar[VXn] = u1;
    Ustar[PRS] = p1;

/*  ------------------------------------------------------------------ 
             Sample solution on  x/t = 0 axis       
    ------------------------------------------------------------------ */

    if (u1 >= 0.0) {  /* --  Solution is sampled to the LEFT of contact -- */

      EXPAND(                        ;   ,
             gL1   = gL*hL[i]/(gL*hL[i] + (p1 - pL)*(VL + wL*uxL));
             Ustar[VXt] = ql[VXt]*gL1;  ,
             Ustar[VXb] = ql[VXb]*gL1;)

      VL1        = VL - (u1 - uxL)*wL;
      gL1        = TwoShock_Lorentz (Ustar,2);
      Ustar[RHO] = MAX(small_rho,1.0/(VL1*gL1));

/*     Shock or rarefaction ?   */

      #if INTERPOLATE_RAREFACTION == YES
       if (p1 < ql[PRS]){
         am  = TwoShock_RarefactionSpeed (ql,    -1);  /*  rarefaction tail  */
         ap  = TwoShock_RarefactionSpeed (Ustar, -1);  /*  rarefaction head  */
         am  = MIN(am, ap);
       }else{        
         ap = am = VL/wL + uxL;                /* -- Shock speed --  */
       }
      #else
       am = ap = VL/wL + uxL;                  /*  -- Shock speed --  */
      #endif

      if (am >= 0.0) {                      /* --  region L  -- */
        for (nv = NFLX; nv--;) ws[i][nv] = ql[nv];

      }else if (ap <= 0.0) {                /* -- region L1 -- */

        for (nv = NFLX; nv--;) ws[i][nv] = Ustar[nv];

      }else{    /*  Solution is inside rarefaction fan, --> interpolate  */

        for (nv = NFLX; nv--;) {
          ws[i][nv] = (am*Ustar[nv] - ap*ql[nv])/(am - ap);
        }
      }
          
    } else {   /* -- Solution is sampled to the RIGHT of contact -- */

      EXPAND(                        ;   ,
             gR1 = gR*hR[i]/(gR*hR[i] + (p1 - pR)*(VR + wR*uxR));
             Ustar[VXt] = qr[VXt]*gR1;  ,
             Ustar[VXb] = qr[VXb]*gR1;)

      gR1        = TwoShock_Lorentz (Ustar,3);
      VR1        = VR - (u1 - uxR)*wR;
      Ustar[RHO] = MAX(small_rho,1.0/(VR1*gR1));

      #if INTERPOLATE_RAREFACTION == YES
       if (p1 < qr[PRS]){
         am  = TwoShock_RarefactionSpeed (Ustar, 1);  /* rarefaction tail  */
         ap  = TwoShock_RarefactionSpeed (qr,    1);  /* rarefaction head  */
         ap  = MAX(am, ap);
       }else{        
         am = ap = VR/wR + uxR;               /* -- Shock speed -- */
       }
      #else
       am = ap = VR/wR + uxR;                 /* -- Shock speed -- */
      #endif

      if (ap <= 0.0) {                        /* -- Shock speed -- */

        for (nv = NFLX; nv--;) ws[i][nv] = qr[nv];

      }else if (am >= 0.0){              /*   region R1   */

        for (nv = NFLX; nv--;) ws[i][nv] = Ustar[nv];

      }else{      /*  Solution is inside rarefaction fan, --> interpolate  */

        for (nv = NFLX; nv--;) {
          ws[i][nv] = (ap*Ustar[nv] - am*qr[nv])/(ap - am);
        }
      }
    }
 
/* ----------------------------------------------------------
                      Compute Fluxes               
   ---------------------------------------------------------- */

    PrimToCons (ws, us, i, i);
    SoundSpeed2 (ws, a2R, hR, i, i, FACE_CENTER, grid);
    Flux (us, ws, a2R, state->flux, state->press, i, i);
    MaxSignalSpeed (ws, a2R, cmin_loc, cmax_loc, i, i);
    a0 = MAX(fabs(cmax_loc[i]), fabs(cmin_loc[i]));
    cmax[i] = a0;

  }

    /* -- Add artificial viscosity -- */

#ifdef ARTIFICIAL_VISC
  for (i = beg; i <= end; i++){
    a0 = ARTIFICIAL_VISC*(MAX(0.0, state->vL[i][VXn] - state->vR[i][VXn]));
    for (nv = 0; nv < NFLX; nv++) {
      state->flux[i][nv] += a0*(state->uL[i][nv] - state->uR[i][nv]);
    }
  }  
#endif

/* ---------------------------------------------------------
    Compute HLL flux in those zones where the Riemann 
    Solver has failed (izone_fail[k] = 1, k > 0)
   --------------------------------------------------------- */

  for (k = 1; k <= nfail; k++){ 
    print1 ("! Failure in Riemann - substituting HLL flux: ");
    Where (izone_fail[k],NULL);
    HLL_Solver (state, izone_fail[k]-2, izone_fail[k]+3, cmax, grid);
  }

}
#undef  MAX_ITER

/* ********************************************************************* */
double TwoShock_Lorentz (double *U, int n)
/*!
 * Compute Lorentz gamma factor.
 *
 *********************************************************************** */
{

   double wl2, beta_fix=0.9999, scrh;

   scrh = EXPAND(U[VX1]*U[VX1], + U[VX2]*U[VX2],  + U[VX3]*U[VX3]);

   if (scrh >= 1.0){

     print ("! u2 > 1 (%12.6e) in TwoShock_Lorentz\n", scrh);
     scrh = beta_fix/sqrt(scrh);
     EXPAND(U[VX1] *= scrh;  ,
            U[VX2] *= scrh;  ,
            U[VX3] *= scrh;)
     scrh = beta_fix*beta_fix;
     exit(1);
   }

   wl2 = 1.0/(1.0 - scrh);

   return(sqrt(wl2));
}

/* ********************************************************************* */
void TwoShock_Shock (double tau0, double u0, double p0, double g0, 
            double V0, double h0, double p1, double *u1,
            double *dudp, double *zeta, int istate)
/*!
 * Compute post shock quantities  u1, dudp, zeta, for a given
 * value of the post-shock pressure p1
 *
 *********************************************************************** */
{
  double a,b,c, da,db,dc, dp;
  double tau1,h1,d_htau1,g1;
  double j2, dx;

/* ***********************************************************
        Use Taub Adiabat to find post-shock Enthalpy
        and mass flux; here j2 --> 1/(j)^2
   *********************************************************** */

  dp = p1 - p0;

  #if EOS == IDEAL
   a = 1.0 - dp/(gmmr*p1);
   b = 1.0 - a;
   c = -h0*(h0 + tau0*dp);
   
   h1 = 0.5/a*(-b + sqrt(b*b - 4.0*a*c));
   tau1 = (h1 - 1.0)/(gmmr*p1);
   g1   = 2.0*h1*gmmr/(2.0*h1 - 1.0);

   j2  = h0*gmmr*tau0 + (h1*tau1 + h0*tau0)*(1.0/(h0 + h1) - 1.0);
   j2 /= gmmr*p1;

   d_htau1  = (h1*tau1 + h0*tau0 - g1*h1*tau1);
   d_htau1 /= (g1*p1 - dp);
  #endif
  #if EOS == TAUB
   a = p0*(3.0*p1 + p0);
   b = -h0*h0*(3.0*p1 + 2.0*p0)
       -h0*tau0*(3.0*p1*p1 - 7.0*p1*p0 - 4.0*p0*p0) - dp;
   c = -h0*tau0*(h0*h0 + 2.0*h0*tau0*p1 + 2.0);
   
   dx = 2.0*c*dp/(-b + sqrt(b*b - 4.0*a*c*dp));/* = [h\tau] */

   g1 = 2.0*c/(-b + sqrt(b*b - 4.0*a*c*dp));/* = [h\tau]/[dp] */
   h1 = sqrt(h0*h0 + (dx + 2.0*h0*tau0)*dp);

   j2 = -g1;

   da = 3.0*p0;
   db = -h0*h0*3.0 - h0*tau0*(6.0*p1 - 7.0*p0) - 1.0;
   dc = c - h0*tau0*dp*2.0*h0*tau0;

   d_htau1 = -(da*dx*dx + db*dx + dc)/(2.0*a*dx + b);
  #endif

  g1    = ((double)istate)*sqrt(V0*V0 + (1.0 - u0*u0)*j2);
  *zeta = (V0*u0 + g1) / (1.0 - u0*u0);

/* ***********************************************************
                  Get post-shock velocity
   *********************************************************** */

  b   = 1.0/(h0*g0 + (p1 - p0)*((*zeta)*u0 + V0));
  *u1 = (h0*g0*u0 + (*zeta)*(p1 - p0)) * b;

/* ***********************************************************
                  Get du/dp for next iteration
   *********************************************************** */

  a     = -0.5*(d_htau1 + j2) / g1;
  *dudp = ((*zeta) + a - (*u1)*((*zeta)*u0 + V0 + a*u0)) * b;
}
/* ********************************************************************* */
double TwoShock_RarefactionSpeed (double w[], int iside)
/*!
 *  Compute head or tail characteristic speeds enclosing
 *  the rarefaction fan:
 *
 * \param[in]  w          a vector of primitive quantities containing 
 *                        the three velocities.
 * \param[in]  iside(IN)  an integer specifying a left (-1) or right (+1)
 *                        rarefaction wave.
 *
 *********************************************************************** */
{
  int    nv;
  double   vx, vt2, vel2;
  double   sroot, delta2, cs2[1], h[1];
  static double **q;
  
  if (q == NULL) q = ARRAY_2D(10, NFLX, double);
  
  for (nv = 0; nv < NFLX; nv++) q[0][nv] = w[nv];
   
  SoundSpeed2(q, cs2, h, 0, 0, FACE_CENTER, NULL);
  
  vx   = w[VXn];
  vt2  = EXPAND(0.0, + w[VXt]*w[VXt], + w[VXb]*w[VXb]);
  vel2 = vx*vx + vt2;

  sroot = cs2[0]*(1.0 - vx*vx - vt2*cs2[0])*(1.0 - vel2);   /* this is eta/gamma */
    
  if (sroot < 0.0){
    print ("! sroot < 0 in RIemann \n");
    for (nv = 0; nv < NFLX; nv++){
      print ("%d  %12.6e  %12.6e\n",nv, qglob_l[nv], qglob_r[nv]);
    }
    QUIT_PLUTO(1);
  }
  sroot = sqrt(sroot);   /*this is eta/gamma */
  
  delta2   = 1.0 - vel2*cs2[0];
  return( (vx*(1.0 - cs2[0]) + (double)iside*sroot)/delta2);
}

#undef  MAX_ITER
#undef  small_p
#undef  small_rho
#undef  accuracy 
#undef  INTERPOLATE_RAREFACTION 
