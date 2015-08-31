#include"pluto.h"

#define  sqrt_1_2  (0.70710678118654752440)

/* **************************************************************************** */
void Roe_Solver (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
/*
 *
 *
 *
 **************************************************************************** */
 {
  int    nv, i, j, k;
  double gmm1;
  double scrh0, scrh1, scrh2, scrh3, scrh4;
  double rho, u, v, w, vel2, bx, by, bz, pr;
  double b1x, b1y, b1z;
  double a2, a, ca2, cf2, cs2;
  double cs, ca, cf, b2, A2, At2;
  double S;
  double alpha_f, alpha_s, beta_y, beta_z;
  double dw[NVAR], *vL, *vR, *uL, *uR;
  double *SL, *SR;
  double eta[NVAR], lambda[NVAR];

  double tau, sqrt_rho;
  double delta, delta_inv, gmm1_inv;

  static double *a2L, *a2R, *pR, *pL;
  static double **fL, **fR, **Rp, **Rc;
  static double **VL, **VR, **UL, **UR;

  double **bgf;
  
  gmm1      = g_gamma - 1.0;
  delta     = 1.e-7;
  delta_inv = 1.0 / delta;
  gmm1_inv  = 1.0 / gmm1;

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NVAR, double);
    fR = ARRAY_2D(NMAX_POINT, NVAR, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif

    Rc = ARRAY_2D(NMAX_POINT, NVAR, double);
    Rp = ARRAY_2D(NMAX_POINT, NVAR, double);
    
  }

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg, end, FACE_CENTER, grid);
  #endif

  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, pL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, pR, beg, end);

  SL = state->SL; SR = state->SR;

/*  Some eigenvector components will always be zero so set
      R = 0 initially   */
      
  for (k = NFLX; k--;  ) {
  for (j = NFLX; j--;  ) {
    Rp[k][j] = Rc[k][j] = 0.0;
  }}


  for (i = beg; i <= end; i++) {

    vL = VL[i]; vR = VR[i];
    uL = UL[i]; uR = UR[i];

    #if SHOCK_FLATTENING == MULTID

  /* -- revert to HLL in proximity of strong shock -- */

     if (CheckZone(i, FLAG_HLL) || CheckZone(i+1,FLAG_HLL)){
       HLL_SPEED (VL, VR, NULL, SL, SR, i, i);

       scrh0 = MAX(fabs(SL[i]), fabs(SR[i]));
       cmax[i] = scrh0; 

       if (SL[i] > 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
         state->press[i] = pL[i];
       } else if (SR[i] < 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
         state->press[i] = pR[i];
       }else{
         scrh0 = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh0;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh0;
       }
       continue;
     }
    #endif

    rho = 0.5*(vL[RHO] + vR[RHO]);
    EXPAND (u = 0.5*(vL[VXn] + vR[VXn]); ,
            v = 0.5*(vL[VXt] + vR[VXt]); ,
            w = 0.5*(vL[VXb] + vR[VXb]); )
    EXPAND (bx = b1x = 0.5*(vL[BXn] + vR[BXn]);  ,
            by = b1y = 0.5*(vL[BXt] + vR[BXt]);  ,
            bz = b1z = 0.5*(vL[BXb] + vR[BXb]);)

    #if BACKGROUND_FIELD == YES
     EXPAND (bx += bgf[i][BXn];   ,
             by += bgf[i][BXt];   ,
             bz += bgf[i][BXb];)
    #endif

    pr = 0.5*(vL[PRS] + vR[PRS]);
  
    for (nv = 0; nv < NFLX; nv++) dw[nv] = vR[nv] - vL[nv];

    tau = 1.0/rho;
    sqrt_rho = sqrt(rho);

    vel2 = EXPAND(u*u, +v*v, +w*w);
    a2   = g_gamma*pr*tau;

    scrh2 = bx*bx;                               /* > 0 */
    scrh3 = SELECT(0.0, by*by, by*by + bz*bz);   /* > 0 */

    b2  = scrh2 + scrh3;   /*  Magnetic field module,  >0  */
    ca2 = scrh2*tau;       /*                          >0  if tau >0 */
    A2  = b2*tau;          /*                          >0  ''   ''  */
    At2 = scrh3*tau;       /*                          >0  */

    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);      /*   >0   */

/*   Now get fast and slow speeds   */
    
    cf2 = 0.5*(a2 + A2 + scrh0);
    cs2 = a2*ca2/cf2;

    cf = sqrt(cf2);  /*  > 0  */
    cs = sqrt(cs2);  /*  > 0  */
    ca = sqrt(ca2);  /*  > 0  */
    a  = sqrt(a2);   /*  > 0  */

    scrh0   = 1.0/scrh0;
    alpha_f = (a2 - cs2)*scrh0;
    alpha_s = (cf2 - a2)*scrh0;

    alpha_f = MAX(0.0, alpha_f);
    alpha_s = MAX(0.0, alpha_s);

    alpha_f = sqrt(alpha_f);
    alpha_s = sqrt(alpha_s);

    scrh0 = sqrt(scrh3);

    if (scrh0 > 1.e-8) {
      SELECT (, beta_y = DSIGN(by);  ,
                beta_y = by / scrh0; beta_z = bz / scrh0;)
    } else {
      SELECT (, beta_y = 1.0;  ,
                beta_z = beta_y = sqrt_1_2;)
    }

    S = (bx >= 0.0 ? 1.0 : -1.0);

  /* --------------------------------------------------------
      Define non-zero entries of primitive right 
      eigenvectors (Rp), wave strength l*dw (=eta) for all 
      8 (or 7) waves.
      left eigenvectors for fast & slow waves can be defined
      in terms of right eigenvectors (see page 296)    
     -------------------------------------------------------- */

  /*  ----  FAST WAVE  (u - c_f)  ----  */

    k         = KFASTM;  
    lambda[k] = u - cf;

    scrh0 = alpha_s*cs*S;
    scrh1 = alpha_s*sqrt_rho*a;
    scrh2 = 0.5 / a2;
    scrh3 = scrh2*tau;
    scrh4 = alpha_f*g_gamma*pr;

    Rp[RHO][k] = rho*alpha_f;         /* right eigenvectors */
    EXPAND(Rp[VXn][k] = -cf*alpha_f;   ,
           Rp[VXt][k] = scrh0*beta_y;  ,
           Rp[VXb][k] = scrh0*beta_z;)
    EXPAND(                            ,                            
           Rp[BXt][k] = scrh1*beta_y;  ,
           Rp[BXb][k] = scrh1*beta_z;)
    Rp[PRS][k] = scrh4;

    eta[k] = (EXPAND(  Rp[VXn][k]*dw[VXn],
                     + Rp[VXt][k]*dw[VXt],
                     + Rp[VXb][k]*dw[VXb]))*scrh2;
    eta[k] += (EXPAND(   0.0               ,
                       + Rp[BXt][k]*dw[BXt] ,
                       + Rp[BXb][k]*dw[BXb])
                + alpha_f*dw[PRS])*scrh3;

  /*  ----  FAST WAVE (u + c_f)  ----  */

    k         = KFASTP;
    lambda[k] = u + cf;

    Rp[RHO][k] = rho*alpha_f;
    EXPAND(Rp[VXn][k] = cf*alpha_f;    ,
           Rp[VXt][k] = -scrh0*beta_y; ,
           Rp[VXb][k] = -scrh0*beta_z;)
    EXPAND(, Rp[BXt][k] = scrh1*beta_y;  ,
             Rp[BXb][k] = scrh1*beta_z;)
    Rp[PRS][k] = scrh4;

    eta[k] = (EXPAND(  Rp[VXn][k]*dw[VXn],
                     + Rp[VXt][k]*dw[VXt],
                     + Rp[VXb][k]*dw[VXb]))*scrh2;
    eta[k] += (EXPAND(   0.0,
                       + Rp[BXt][k]*dw[BXt],
                       + Rp[BXb][k]*dw[BXb])
                + alpha_f*dw[PRS])*scrh3;

  /*  ----  ENTROPY WAVE  (u)  ----  */

    k         = KENTRP;
    lambda[k] = u;

    Rp[RHO][k] = 1.0;

    eta[k] = dw[RHO] - dw[PRS]/a2;

  /*  ----  MAGNETIC FLUX, for 8-wave formulation only (u)  ----  */

    #ifdef GLM_MHD
     lambda[KPSI_GLMP] =  glm_ch;
     lambda[KPSI_GLMM] = -glm_ch;
     eta[KPSI_GLMP] = eta[KPSI_GLMM] = 0.0;
    #else
     k = KDIVB;
     lambda[k] = u;
     #if DIVB_CONTROL == EIGHT_WAVES
      Rc[BXn][k] = 1.0;
      eta[k]    = dU[BXn];
     #else
      Rc[BXn][k] = eta[k] = 0.0;
     #endif
    #endif

    #if COMPONENTS > 1       /*  ----  SLOW WAVE (u - c_s)    ----  */
     k         = KSLOWM; 
     lambda[k] = u - cs;

     scrh0 = alpha_f*cf*S;
     scrh1 = alpha_f*sqrt_rho*a;
     scrh4 = alpha_s*g_gamma*pr;

     Rp[RHO][k] = rho*alpha_s;
     EXPAND(Rp[VXn][k] = -cs*alpha_s;   ,
            Rp[VXt][k] = -scrh0*beta_y; ,
            Rp[VXb][k] = -scrh0*beta_z;)
     EXPAND(, Rp[BXt][k] = -scrh1*beta_y; ,
              Rp[BXb][k] = -scrh1*beta_z;)
     Rp[PRS][k] = scrh4;

     eta[k] = (EXPAND(  Rp[VXn][k]*dw[VXn],
                      + Rp[VXt][k]*dw[VXt],
                      + Rp[VXb][k]*dw[VXb]))*scrh2;
     eta[k] +=(EXPAND(   0.0 ,
                       + Rp[BXt][k]*dw[BXt],
                       + Rp[BXb][k]*dw[BXb])
              +alpha_s*dw[PRS])*scrh3;

    /*  ----   SLOW WAVE (u + c_s)  ----  */

     k         = KSLOWP;
     lambda[k] = u + cs;

     Rp[RHO][k] = rho*alpha_s;
     EXPAND(Rp[VXn][k] = cs*alpha_s;   ,
            Rp[VXt][k] = scrh0*beta_y; ,
            Rp[VXb][k] = scrh0*beta_z;)

     EXPAND (  , Rp[BXt][k] = -scrh1*beta_y;  ,
                 Rp[BXb][k] = -scrh1*beta_z;)
     Rp[PRS][k] = scrh4; 

     eta[k] = (EXPAND(  Rp[VXn][k]*dw[VXn],
                      + Rp[VXt][k]*dw[VXt],
                      + Rp[VXb][k]*dw[VXb]))*scrh2;
     eta[k] +=(EXPAND(   0.0 ,
                       + Rp[BXt][k]*dw[BXt],
                       + Rp[BXb][k]*dw[BXb])
              +alpha_s*dw[PRS])*scrh3;
    #endif

    #if COMPONENTS == 3

    /*  ----  ALFVEN WAVE (u-c_a)  ----  */

     k         = KALFVM;
     lambda[k] = u - ca;

     scrh2 = beta_y*sqrt_1_2;
     scrh3 = beta_z*sqrt_1_2;
     Rp[VXt][k] = -scrh3;  
     Rp[VXb][k] =  scrh2;
     Rp[BXt][k] = -scrh3*sqrt_rho*S;   
     Rp[BXb][k] =  scrh2*sqrt_rho*S;

     eta[k] =   Rp[VXt][k]*dw[VXt] + Rp[VXb][k]*dw[VXb]
              +(Rp[BXt][k]*dw[BXt] + Rp[BXb][k]*dw[BXb])*tau;

    /*  ----  ALFVEN WAVE (u+c_a)  ----  */

     k = KALFVP; 
     lambda[k] = u + ca;
     Rp[VXt][k] = -scrh3; 
     Rp[VXb][k] =  scrh2;
     Rp[BXt][k] =  scrh3*sqrt_rho*S;
     Rp[BXb][k] = -scrh2*sqrt_rho*S; 

     eta[k] =   Rp[VXt][k]*dw[VXt] + Rp[VXb][k]*dw[VXb]
              +(Rp[BXt][k]*dw[BXt] + Rp[BXb][k]*dw[BXb])*tau;
    #endif

  /*  -------------------------------------------------------------------
       Find conservative eigenvectors; this is done by vector 
       multiplication, as shown on the reference paper ("A solution 
       adaptive upwind scheme for MHD", JCP 154, 284 (1999)), on page
       297:

                         Rc = dU/dW * Rp
 
       Primitive left eigenvectors are not necessary, since the jump
       is invariant, i.e.:  

                          Lp*dW = Lc*dU
      ------------------------------------------------------------------- */

    for (k = 0; k < NFLX; k++) {
      Rc[RHO][k] = Rp[RHO][k];
      EXPAND (Rc[MXn][k] = u*Rp[RHO][k] + rho*Rp[VXn][k]; ,
              Rc[MXt][k] = v*Rp[RHO][k] + rho*Rp[VXt][k]; ,
              Rc[MXb][k] = w*Rp[RHO][k] + rho*Rp[VXb][k];)
      EXPAND (Rc[BXn][k] = Rp[BXn][k];  ,
              Rc[BXt][k] = Rp[BXt][k];  ,
              Rc[BXb][k] = Rp[BXb][k];)
      
      scrh0 = EXPAND(u*Rp[VXn][k], + v*Rp[VXt][k], + w*Rp[VXb][k]);
      scrh1 = EXPAND(b1x*Rp[BXn][k], + b1y*Rp[BXt][k], + b1z*Rp[BXb][k]);
      Rc[ENG][k] =  0.5*vel2*Rp[RHO][k] + rho*scrh0 
                  + scrh1 + Rp[PRS][k]*gmm1_inv;                  
    }

  /*  ----  COMPUTE MAX EIGENVALUE  ----  */

    cmax[i] = fabs (u) + cf;  
    g_maxMach = MAX (fabs (u / a), g_maxMach);

    SL[i] = lambda[KFASTM];
    SR[i] = lambda[KFASTP];

  /*  ----  ADD 'VISCOUS' FLUX FIRST  ---- */

    for (nv = 0; nv < NFLX; nv++) {
      scrh0 = 0.0;
      for (k = 0; k < NFLX; k++) {
        scrh1 = fabs(lambda[k]);
        if ((k == KFASTM || k == KFASTP || k == KSLOWM || k == KSLOWP) 
             && scrh1 < 0.5*delta){        /* entropy fix     */
           scrh1 = scrh1*scrh1/delta + 0.25*delta;
        }
        scrh0 += scrh1*eta[k]*Rc[nv][k];
      }
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - scrh0);
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }

  #if DIVB_CONTROL == EIGHT_WAVES
   ROE_DivBSource (state, grid, beg, end);
  #endif
}
#undef sqrt_1_2
