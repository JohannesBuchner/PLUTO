/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implementation of the Roe Riemann solver for the MHD equations.

  Solve the Riemann problem for the adiabatic and isothermal MHD 
  equations using the linearized Riemann solver of Roe.
  The implementation follows the approach outlined by 
  Cargo & Gallice (1997). 
  
  The isothermal version is recovered by taking the limit of 
  \f$ \bar{a}^2 \f$ for \f$ \gamma \to 1\f$ which gives 
  (page 451) \f$\bar{a}^2 \to {\rm g\_isoSoundSpeed2} + X\f$
  where \c X is defined as in the adiabatic case.
  This follows by imposing zero jump across the entropy wave 
  (first of Eq. 4.20, page 452) giving 
  \f$ \Delta p = (\bar{a}^2 - X)\Delta\rho\f$.
  Then all the terms like 
  \f$ [X\Delta\rho + \Delta p] \to \bar{a}^2\Delta\rho \f$.
 
  When the macro HLL_HYBRIDIZATION is set to YES, we revert to the HLL 
  solver whenever an unphysical state appear in the solution.
  
  The macro CHECK_ROE_MATRIX can be used to verify that the 
  characteristic decomposition reproduces the Roe matrix.
  This can be done also when BACKGROUND_FIELD is set to YES.

  On input, this function takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "Roe Matrices for Ideal MHD and Systematic Construction 
       of Roe Matrices for Systems of Conservation Laws",
       P. Cargo, G. Gallice, JCP (1997), 136, 446.
       
  \authors A. Mignone (mignone@ph.unito.it)
           C. Zanni   (zanni@oato.inaf.it)
  \date    April 11, 2014 
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define sqrt_1_2  (0.70710678118654752440)
#define HLL_HYBRIDIZATION    NO
#define CHECK_ROE_MATRIX     NO

/* ********************************************************************* */
void Roe_Solver (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic MHD equations using the 
 * Roe Riemann solver of Cargo & Gallice (1997).
 * 
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int  nv, i, j, k;
  int  ifail;
  double rho, u, v, w, vel2, bx, by, bz, pr;
  double a2, a, ca2, cf2, cs2;
  double cs, ca, cf, b2;
  double alpha_f, alpha_s, beta_y, beta_z, beta_v, scrh, sBx;
  double dV[NFLX], dU[NFLX], *vL, *vR, *uL, *uR;
  double *SL, *SR;
  double Rc[NFLX][NFLX], eta[NFLX], lambda[NFLX];
  double alambda[NFLX], Uv[NFLX];

  double sqrt_rho;
  double delta, delta_inv;
  
  double g1, sl, sr, H, Hgas, HL, HR, Bx, By, Bz, X;
  double vdm, BdB, beta_dv, beta_dB;
  double bt2, Btmag, sqr_rho_L, sqr_rho_R;

  static double *pR, *pL, *a2L, *a2R;
  static double **fL, **fR;
  static double **VL, **VR, **UL, **UR;
  double **bgf;
  #if BACKGROUND_FIELD == YES
   double B0x, B0y, B0z, B1x, B1y, B1z;
  #endif      
  double Us[NFLX];
  delta    = 1.e-6;

/* -----------------------------------------------------------
   1. Allocate static memory areas
   ----------------------------------------------------------- */
   
  if (fL == NULL){

    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

/* -----------------------------------------------------------
   2. Background field
   ----------------------------------------------------------- */

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg, end, FACE_CENTER, grid);
   #if (DIVB_CONTROL == EIGHT_WAVES)
    print1 ("! Roe_Solver: background field and 8wave not tested\n");
    QUIT_PLUTO(1);
   #endif
  #endif

/* -----------------------------------------------------------
   3. GLM pre-Rieman solver
   ----------------------------------------------------------- */
   
  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
   4. Compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, pL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, pR, beg, end);

  SL = state->SL; SR = state->SR;

  #if EOS == IDEAL
   g1 = g_gamma - 1.0;
  #endif

/* -------------------------------------------------
   5. Some eigenvectors components will always be 
      zero so set Rc = 0 initially  
   -------------------------------------------------- */
     
  for (k = NFLX; k--;  ) {
  for (j = NFLX; j--;  ) {
    Rc[k][j] = 0.0;
  }}

/* ---------------------------------------------------------------------
   6. Begin main loop
   --------------------------------------------------------------------- */
   
  for (i = beg; i <= end; i++) {

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

  /* -----------------------------------------------------------
     6a.  switch to HLL in proximity of strong shock 
     ----------------------------------------------------------- */

    #if SHOCK_FLATTENING == MULTID
     if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){
       HLL_Speed (VL, VR, a2L, a2R, NULL, SL, SR, i, i);

       scrh = MAX(fabs(SL[i]), fabs(SR[i]));
       cmax[i] = scrh;

       if (SL[i] > 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
         state->press[i] = pL[i];
       } else if (SR[i] < 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
         state->press[i] = pR[i];
       }else{
         scrh = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
       }
       continue;
     }
    #endif

  /* ----------------------------------------------------------------
     6b. Compute jumps in primitive and conservative variables.
      
         Note on the velocity jump:
         Formally the jump in velocity should be written using
         conservative variables:
      
         Delta(v) = Delta(m/rho) = (Delta(m) - v*Delta(rho))/rho
      
         with rho and v defined by the second and first of [4.6]
         Still, since we reconstruct v and not m, the following MAPLE 
         script shows that the previous expression simplifies to
         Delta(v) = vR - vL.
      
         restart;
         sR   := sqrt(rho[R])/(sqrt(rho[R]) + sqrt(rho[L])):
         sL   := sqrt(rho[L])/(sqrt(rho[R]) + sqrt(rho[L])):
         rhoa := sL*rho[R] + sR*rho[L]:
         va   := sL*v[L]   + sR*v[R]:
         dV := (rho[R]*v[R] - rho[L]*v[L] - va*(rho[R]-rho[L]))/rhoa:
         simplify(dV);
     ---------------------------------------------------------------- */

    for (nv = 0; nv < NFLX; nv++) { 
      dV[nv] = vR[nv] - vL[nv];
      dU[nv] = uR[nv] - uL[nv];
    }

  /* ---------------------------------
     6c. Compute Roe averages 
     --------------------------------- */

    sqr_rho_L = sqrt(vL[RHO]);
    sqr_rho_R = sqrt(vR[RHO]);

    sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
    sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);

/*      sl = sr = 0.5;    */
    
    rho = sr*vL[RHO] + sl*vR[RHO];

    sqrt_rho = sqrt(rho);

    EXPAND (u = sl*vL[VXn] + sr*vR[VXn];  ,
            v = sl*vL[VXt] + sr*vR[VXt];  ,
            w = sl*vL[VXb] + sr*vR[VXb];)

    EXPAND (Bx = sr*vL[BXn] + sl*vR[BXn];  ,
            By = sr*vL[BXt] + sl*vR[BXt];  ,
            Bz = sr*vL[BXb] + sl*vR[BXb];)

   #if BACKGROUND_FIELD == YES
   /* -- Define field B0 and total B. B1 is the deviation -- */   
    EXPAND (B0x = bgf[i][BXn]; B1x = sr*vL[BXn] + sl*vR[BXn]; Bx = B0x + B1x;  ,
            B0y = bgf[i][BXt]; B1y = sr*vL[BXt] + sl*vR[BXt]; By = B0y + B1y;  ,
            B0z = bgf[i][BXb]; B1z = sr*vL[BXb] + sl*vR[BXb]; Bz = B0z + B1z;)
   #else
    EXPAND (Bx = sr*vL[BXn] + sl*vR[BXn];  ,
            By = sr*vL[BXt] + sl*vR[BXt];  ,
            Bz = sr*vL[BXb] + sl*vR[BXb];)
   #endif

    sBx = (Bx >= 0.0 ? 1.0 : -1.0);

    EXPAND(bx = Bx/sqrt_rho;  ,
           by = By/sqrt_rho;  ,
           bz = Bz/sqrt_rho; )
    
    bt2   = EXPAND(0.0  , + by*by, + bz*bz);
    b2    = bx*bx + bt2;
    Btmag = sqrt(bt2*rho);

    X  = EXPAND(dV[BXn]*dV[BXn], + dV[BXt]*dV[BXt], + dV[BXb]*dV[BXb]);
    X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;   

    vdm = EXPAND(u*dU[MXn],  + v*dU[MXt],  + w*dU[MXb]);
    #if BACKGROUND_FIELD == YES /* BdB = B1.dB1 (deviation only) */
     BdB = EXPAND(B1x*dU[BXn], + B1y*dU[BXt], + B1z*dU[BXb]); 
    #else
     BdB = EXPAND(Bx*dU[BXn], + By*dU[BXt], + Bz*dU[BXb]);
    #endif
    
  /* ---------------------------------------
     6d. Compute enthalpy and sound speed.
     --------------------------------------- */

    #if EOS == ISOTHERMAL 
     a2 = 0.5*(a2L[i] + a2R[i]) + X;  /* in most cases a2L = a2R
                                         for isothermal MHD */
    #elif EOS == BAROTROPIC
     print ("! Roe_Solver: not implemented for barotropic EOS\n");
     QUIT_PLUTO(1);
    #elif EOS == IDEAL
     vel2    = EXPAND(u*u, + v*v, + w*w);
     dV[PRS] = g1*((0.5*vel2 - X)*dV[RHO] - vdm + dU[ENG] - BdB); 
     
     HL   = (uL[ENG] + pL[i])/vL[RHO];
     HR   = (uR[ENG] + pR[i])/vR[RHO];
     H    = sl*HL + sr*HR;   /* total enthalpy */

     #if BACKGROUND_FIELD == YES
      scrh = EXPAND(B1x*Bx, + B1y*By, + B1z*Bz);
      Hgas = H - scrh/rho;   /* gas enthalpy */
     #else
      Hgas = H - b2;         /* gas enthalpy */
     #endif

     a2 = (2.0 - g_gamma)*X + g1*(Hgas - 0.5*vel2);
     if (a2 < 0.0) {
      printf ("! Roe: a2 = %12.6e < 0.0 !! \n",a2);
      Show(VL,i);
      Show(VR,i);
      QUIT_PLUTO(1);
     }      
    #endif /* EOS == IDEAL */
    
  /* ------------------------------------------------------------
     6e. Compute fast and slow magnetosonic speeds.

      The following expression appearing in the definitions
      of the fast magnetosonic speed 
    
       (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

      is always positive and avoids round-off errors.
      
      Note that we always use the total field to compute the 
      characteristic speeds.
     ------------------------------------------------------------ */
        
    scrh = a2 - b2;
    ca2  = bx*bx;
    scrh = scrh*scrh + 4.0*bt2*a2;    
    scrh = sqrt(scrh);    

    cf2 = 0.5*(a2 + b2 + scrh); 
    cs2 = a2*ca2/cf2;   /* -- same as 0.5*(a2 + b2 - scrh) -- */
    
    cf = sqrt(cf2);
    cs = sqrt(cs2);
    ca = sqrt(ca2);
    a  = sqrt(a2); 
    
    if (cf == cs) {
      alpha_f = 1.0;
      alpha_s = 0.0;
    }else if (a <= cs) {
      alpha_f = 0.0;
      alpha_s = 1.0;
    }else if (cf <= a){
      alpha_f = 1.0;
      alpha_s = 0.0;
    }else{ 
      scrh    = 1.0/(cf2 - cs2);
      alpha_f = (a2  - cs2)*scrh;
      alpha_s = (cf2 -  a2)*scrh;
      alpha_f = MAX(0.0, alpha_f);
      alpha_s = MAX(0.0, alpha_s);
      alpha_f = sqrt(alpha_f);
      alpha_s = sqrt(alpha_s);
    }

    if (Btmag > 1.e-9) {
      SELECT(                     , 
             beta_y = DSIGN(By);  ,
             beta_y = By/Btmag; 
             beta_z = Bz/Btmag;)
    } else {
      SELECT(                       , 
             beta_y = 1.0;          ,
             beta_z = beta_y = sqrt_1_2;)
    }

  /* -------------------------------------------------------------------
     6f. Compute non-zero entries of conservative eigenvectors (Rc), 
         wave strength L*dU (=eta) for all 8 (or 7) waves using the
         expressions given by Eq. [4.18]--[4.21]. 
         Fast and slow eigenvectors are multiplied by a^2 while
         jumps are divided by a^2.
      
         Notes:
         - the expression on the paper has a typo in the very last term 
           of the energy component: it should be + and not - !
         - with background field splitting: additional terms must be 
           added to the energy component for fast, slow and Alfven waves.
           To obtain energy element, conservative eigenvector (with 
           total field) must be multiplied by | 0 0 0 0 -B0y -B0z 1 |.
           Also, H - b2 does not give gas enthalpy. A term b0*btot must 
           be added and eta (wave strength) should contain total field 
           and deviation's delta.
     ------------------------------------------------------------------- */

  /* -----------------------
      Fast wave:  u - c_f
     ----------------------- */

    k = KFASTM;
    lambda[k] = u - cf;

    scrh    = alpha_s*cs*sBx;
    beta_dv = EXPAND(0.0, + beta_y*dV[VXt], + beta_z*dV[VXb]);
    beta_dB = EXPAND(0.0, + beta_y*dV[BXt], + beta_z*dV[BXb]);
    beta_v  = EXPAND(0.0, + beta_y*v,       + beta_z*w);

    Rc[RHO][k] = alpha_f;
    EXPAND(Rc[MXn][k] = alpha_f*lambda[k];          ,
           Rc[MXt][k] = alpha_f*v + scrh*beta_y;    ,
           Rc[MXb][k] = alpha_f*w + scrh*beta_z;) 
    EXPAND(                                         ,
           Rc[BXt][k] = alpha_s*a*beta_y/sqrt_rho;  ,
           Rc[BXb][k] = alpha_s*a*beta_z/sqrt_rho;)

    #if EOS == IDEAL
     Rc[ENG][k] =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;
     #if BACKGROUND_FIELD == YES
      Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
     #endif

     eta[k] =   alpha_f*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
              - rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
    #elif EOS == ISOTHERMAL
     eta[k] =   alpha_f*(0.0*X + a2)*dV[RHO] + rho*scrh*beta_dv
              - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*beta_dB;
    #endif
    
    eta[k] *= 0.5/a2;

  /* -----------------------
      Fast wave:  u + c_f
     ----------------------- */

    k = KFASTP;
    lambda[k] = u + cf;

    Rc[RHO][k] = alpha_f;
    EXPAND( Rc[MXn][k] = alpha_f*lambda[k];         ,
            Rc[MXt][k] = alpha_f*v - scrh*beta_y;  ,
            Rc[MXb][k] = alpha_f*w - scrh*beta_z; ) 
    EXPAND(                              ,                                
            Rc[BXt][k] = Rc[BXt][KFASTM];  ,
            Rc[BXb][k] = Rc[BXb][KFASTM]; )

    #if EOS == IDEAL
     Rc[ENG][k] =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;

     #if BACKGROUND_FIELD == YES
      Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
     #endif

     eta[k] =   alpha_f*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
    #elif EOS == ISOTHERMAL
     eta[k] =   alpha_f*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[VXn]      + sqrt_rho*alpha_s*a*beta_dB;
    #endif

    eta[k] *= 0.5/a2;

  /* -----------------------
      Entropy wave:  u
     ----------------------- */

    #if EOS == IDEAL
     k = KENTRP;
     lambda[k] = u;

     Rc[RHO][k] = 1.0;
     EXPAND( Rc[MXn][k] = u; ,
             Rc[MXt][k] = v; ,
             Rc[MXb][k] = w; )
     Rc[ENG][k] = 0.5*vel2 + (g_gamma - 2.0)/g1*X;

     eta[k] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
    #endif

  /* -----------------------------------------------------------------
      div.B wave (u): this wave exists when: 

       1) 8 wave formulation
       2) CT, since we always have 8 components, but it 
          carries zero jump.

      With GLM, KDIVB is replaced by KPSI_GLMM, KPSI_GLMP and these
      two waves should not enter in the Riemann solver (eta = 0.0) 
      since the 2x2 linear system formed by (B,psi) has already 
      been solved.
     ----------------------------------------------------------------- */

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
    
    #if COMPONENTS > 1    

   /* -----------------------
       Slow wave:  u - c_s
      ----------------------- */

     scrh = alpha_f*cf*sBx;
     
     k = KSLOWM;
     lambda[k] = u - cs;

     Rc[RHO][k] = alpha_s;
     EXPAND( Rc[MXn][k] = alpha_s*lambda[k];        ,
             Rc[MXt][k] = alpha_s*v - scrh*beta_y;  ,
             Rc[MXb][k] = alpha_s*w - scrh*beta_z;) 
     EXPAND(                                            ,                                
             Rc[BXt][k] = - alpha_f*a*beta_y/sqrt_rho;  ,
             Rc[BXb][k] = - alpha_f*a*beta_z/sqrt_rho; )

     #if EOS == IDEAL
      Rc[ENG][k] =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                   - alpha_f*a*Btmag/sqrt_rho; 
      #if BACKGROUND_FIELD == YES
       Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
      #endif

      eta[k] =   alpha_s*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
               - rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB;
     #elif EOS == ISOTHERMAL
      eta[k] =   alpha_s*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
               - rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB;
     #endif

     eta[k] *= 0.5/a2;

   /* -----------------------
       Slow wave:  u + c_s
      ----------------------- */

     k = KSLOWP;
     lambda[k] = u + cs; 

     Rc[RHO][k] = alpha_s;
     EXPAND(Rc[MXn][k] = alpha_s*lambda[k];         ,
            Rc[MXt][k] = alpha_s*v + scrh*beta_y;   ,
            Rc[MXb][k] = alpha_s*w + scrh*beta_z; ) 
     EXPAND(                                , 
            Rc[BXt][k] = Rc[BXt][KSLOWM];   ,
            Rc[BXb][k] = Rc[BXb][KSLOWM];)

     #if EOS == IDEAL
      Rc[ENG][k] =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                   - alpha_f*a*Btmag/sqrt_rho;
      #if BACKGROUND_FIELD == YES
       Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
      #endif

      eta[k] =   alpha_s*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
               + rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB; 
     #elif EOS == ISOTHERMAL
      eta[k] =   alpha_s*(0.*X + a2)*dV[RHO] + rho*scrh*beta_dv
               + rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB; 
     #endif

     eta[k] *= 0.5/a2;

    #endif

    #if COMPONENTS == 3

   /* ------------------------
       Alfven wave:  u - c_a
      ------------------------ */

     k = KALFVM;
     lambda[k] = u - ca;

     Rc[MXt][k] = - rho*beta_z;  
     Rc[MXb][k] = + rho*beta_y;
     Rc[BXt][k] = - sBx*sqrt_rho*beta_z;   
     Rc[BXb][k] =   sBx*sqrt_rho*beta_y;
     #if EOS == IDEAL
      Rc[ENG][k] = - rho*(v*beta_z - w*beta_y);
      #if BACKGROUND_FIELD == YES
       Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
      #endif
     #endif

     eta[k] = + beta_y*dV[VXb]               - beta_z*dV[VXt] 
              + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

     eta[k] *= 0.5;

   /* -----------------------
       Alfven wave:  u + c_a 
      ----------------------- */

     k = KALFVP;
     lambda[k] = u + ca;

     Rc[MXt][k] = - Rc[MXt][KALFVM];  
     Rc[MXb][k] = - Rc[MXb][KALFVM];
     Rc[BXt][k] =   Rc[BXt][KALFVM];   
     Rc[BXb][k] =   Rc[BXb][KALFVM];
     #if EOS == IDEAL
      Rc[ENG][k] = - Rc[ENG][KALFVM];
      #if BACKGROUND_FIELD == YES
       Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
      #endif
     #endif

     eta[k] = - beta_y*dV[VXb]               + beta_z*dV[VXt] 
              + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

     eta[k] *= 0.5;
    #endif

   /* -----------------------------------------
      6g. Compute maximum signal velocity
      ----------------------------------------- */

    cmax[i] = fabs (u) + cf;
    g_maxMach = MAX (fabs (u / a), g_maxMach);
    for (k = 0; k < NFLX; k++) alambda[k] = fabs(lambda[k]);

   /* --------------------------------
      6h. Entropy Fix 
      -------------------------------- */
      
    if (alambda[KFASTM] < 0.5*delta) {
      alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
    }
    if (alambda[KFASTP] < 0.5*delta) {
      alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
    }
    #if COMPONENTS > 1
     if (alambda[KSLOWM] < 0.5*delta) {
       alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
     }
     if (alambda[KSLOWP] < 0.5*delta) {
       alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta; 
     }
    #endif
   
  /*  ---------------------------------
      6i. Compute Roe numerical flux 
      --------------------------------- */

    for (nv = 0; nv < NFLX; nv++) {
      scrh = 0.0;
      for (k = 0; k < NFLX; k++) {
        scrh += alambda[k]*eta[k]*Rc[nv][k];
      }
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - scrh);
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);

  /* --------------------------------------------------------
     6j. Check the Roe matrix condition, FR - FL = A*(UR - UL)
         where A*(UR - UL) = R*lambda*eta.
    -------------------------------------------------------- */

    #if CHECK_ROE_MATRIX == YES
     for (nv = 0; nv < NFLX; nv++){
       dV[nv] = fR[i][nv] - fL[i][nv]; 
       if (nv == MXn) dV[MXn] += pR[i] - pL[i];
       for (k = 0; k < NFLX; k++){
         dV[nv] -= Rc[nv][k]*eta[k]*lambda[k];
       }
       if (fabs(dV[nv]) > 1.e-4){
         printf (" ! Roe matrix condition not satisfied, var = %d\n", nv);
         printf (" ! Err = %12.6e\n",dV[nv]); 
         Show(VL, i);
         Show(VR, i);
         exit(1);
       }
     } 
    #endif

  /* -------------------------------------------------------------
     6k. Save max and min Riemann fan speeds for EMF computation.
     ------------------------------------------------------------- */

    SL[i] = lambda[KFASTM];
    SR[i] = lambda[KFASTP];

  /* -----------------------------------------------------------------
     6l. Hybridize with HLL solver: replace occurences of unphysical 
        states (p < 0, rho < 0) with HLL Flux. Reference:
      
        "A Positive Conservative Method for MHD based based on HLL 
         and Roe methods", P. Janhunen, JCP (2000), 160, 649
     ----------------------------------------------------------------- */

    #if HLL_HYBRIDIZATION == YES

     if (SL[i] < 0.0 && SR[i] > 0.0){
       ifail = 0;    

      /* -----------------------
           check left state 
         ----------------------- */

       #if EOS == ISOTHERMAL
        Uv[RHO] = uL[RHO] + (state->flux[i][RHO] - fL[i][RHO])/SL[i];        
        ifail  = (Uv[RHO] < 0.0);
       #else
        for (nv = NFLX; nv--; ){
          Uv[nv] = uL[nv] + (state->flux[i][nv] - fL[i][nv])/SL[i];        
        }
        Uv[MXn] += (state->press[i] - pL[i])/SL[i];    
 
        vel2  = EXPAND(Uv[MX1]*Uv[MX1], + Uv[MX2]*Uv[MX2], + Uv[MX3]*Uv[MX3]);
        b2    = EXPAND(Uv[BX1]*Uv[BX1], + Uv[BX2]*Uv[BX2], + Uv[BX3]*Uv[BX3]);    
        pr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
        ifail = (pr < 0.0) || (Uv[RHO] < 0.0);
       #endif

      /* -----------------------
           check right state 
         ----------------------- */

       #if EOS == ISOTHERMAL
        Uv[RHO] = uR[RHO] + (state->flux[i][RHO] - fR[i][RHO])/SR[i];
        ifail  = (Uv[RHO] < 0.0);
       #else
        for (nv = NFLX; nv--;  ){
          Uv[nv] = uR[nv] + (state->flux[i][nv] - fR[i][nv])/SR[i];
        }
        Uv[MXn] += (state->press[i] - pR[i])/SR[i];

        vel2  = EXPAND(Uv[MX1]*Uv[MX1], + Uv[MX2]*Uv[MX2], + Uv[MX3]*Uv[MX3]);
        b2    = EXPAND(Uv[BX1]*Uv[BX1], + Uv[BX2]*Uv[BX2], + Uv[BX3]*Uv[BX3]);    
        pr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
        ifail = (pr < 0.0) || (Uv[RHO] < 0.0);
       #endif

    /* -------------------------------------------------------
        Use the HLL flux function if the interface lies 
        within a strong shock. The effect of this switch is 
        visible in the Mach reflection test.
       ------------------------------------------------------- */

       #if DIMENSIONS > 1   
        #if EOS == ISOTHERMAL
         scrh  = fabs(vL[RHO] - vR[RHO]);
         scrh /= MIN(vL[RHO], vR[RHO]);
        #else       
         scrh  = fabs(vL[PRS] - vR[PRS]);
         scrh /= MIN(vL[PRS], vR[PRS]);
        #endif
        if (scrh > 1.0 && (vR[VXn] < vL[VXn])) ifail = 1;
       #endif
      
       if (ifail){
         scrh = 1.0/(SR[i] - SL[i]);
         for (nv = 0; nv < NFLX; nv++) {
           state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                                SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
       }
     }
    #endif
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
  
  #if DIVB_CONTROL == EIGHT_WAVES
   Roe_DivBSource (state, beg + 1, end, grid);
  #endif

}
#undef sqrt_1_2
#undef HLL_HYBRIDIZATION
#undef CHECK_ROE_MATRIX
