/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLLD Riemann solver for the classical MHD equations.

  Solve the Riemann problem for the classical MHD equations with 
  an ideal or isothermal equation of state.
  For the ideal case, we implement the four-state (of five-wave) 
  HLLD solver of   Miyoshi & Kusano (2005).
  For the isothermal case, we use the three-state approximation
  described by Mignone (2007).
  
  The macro VERIFY_CONSISTENCY_CONDITION can be turned to YES to verify
  the correctness of the implementation.
  
  On input, this function takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "A muLti-state HLL approximate Riemann Solver for ideal MHD", 
       Miyoshi, T., Kusano, K., JCP 2005 
    - "A simple and accurate Riemann solver for isothermal MHD",
       Mignone, JCP (2007) 225, 1427.
       
  \authors A. Mignone (mignone@ph.unito.it)
           C. Zanni   (zanni@oato.inaf.it)
  \date    April 14, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define VERIFY_CONSISTENCY_CONDITION NO

#if EOS == IDEAL || EOS == PVTE_LAW
/* ********************************************************************* */
void HLLD_Solver (const State_1D *state, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic MHD equations using the 
 * four-state HLLD Riemann solver of Miyoshi & Kusano (2005).
 * 
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;
  int    revert_to_hllc;
  double   scrh, Uhll[NFLX];

  double usL[NFLX], ussl[NFLX];
  double usR[NFLX], ussr[NFLX];
  
  double vsL, wsL, scrhL, S1L, sqrL, duL;
  double vsR, wsR, scrhR, S1R, sqrR, duR;
  double Bx, Bx1, SM, sBx, pts;
  double vss, wss;
  double *vL, *vR, *uL, *uR, *SL, *SR;
  static double *ptL, *ptR;
  static double **fL, **fR;
  static double **VL, **VR, **UL, **UR;
  static double *a2L, *a2R;
  #if BACKGROUND_FIELD == YES
   double B0n, B0t, B0b;
  #endif
  double **bgf;

  if (fL == NULL){
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);

    ptR = ARRAY_1D(NMAX_POINT, double);
    ptL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg, end, FACE_CENTER, grid);
  #endif

  #if DIVB_CONTROL == EIGHT_WAVES
   print ("! hlld Riemann solver does not work with Powell's 8-wave\n");
   QUIT_PLUTO(1);
  #endif

  #ifdef GLM_MHD
/*
double **uL0, **uR0;
static double *BxL, *BxR, *psiL, *psiR;
if (BxL == NULL){
 BxL = ARRAY_1D(NMAX_POINT, double);
 BxR = ARRAY_1D(NMAX_POINT, double);
 psiL = ARRAY_1D(NMAX_POINT, double);
 psiR = ARRAY_1D(NMAX_POINT, double);
}
uL0 = state->uL;  
uR0 = state->uR;
for (i = beg; i <= end; i++){
 BxL[i]  = state->vL[i][BXn];
 BxR[i]  = state->vR[i][BXn];
 psiL[i] = state->vL[i][PSI_GLM];
 psiR[i] = state->vL[i][PSI_GLM];
}
state->uL = UL;
state->uR = UR;

   GLM_Solve (state, beg, end, grid);
*/

   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);

  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     Compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, ptL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, ptR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, bgf, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {
    
    #if BACKGROUND_FIELD == YES
     EXPAND (B0n = bgf[i][BXn];  ,
             B0t = bgf[i][BXt];  ,
             B0b = bgf[i][BXb];)
    #endif
    
  /* ----------------------------------------
      get max propagation speed for dt comp.
     ---------------------------------------- */             

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

/* ---------------------------------------------------------- 
                COMPUTE FLUXES and STATES
   ---------------------------------------------------------- */

    if (SL[i] >= 0.0){                     /*  ----  Region L  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = ptL[i];

    }else if (SR[i] <= 0.0) {              /*  ----  Region R  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = ptR[i];
 
    } else {

#if SHOCK_FLATTENING == MULTID
    if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){
      scrh = 1.0/(SR[i] - SL[i]);
      for (nv = NFLX; nv--; ){
        state->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                           +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
      continue;
    }
#endif

    /* ---------------------------
              Compute U*  
       --------------------------- */

      scrh = 1.0/(SR[i] - SL[i]);
      Bx1  = Bx = (SR[i]*vR[BXn] - SL[i]*vL[BXn])*scrh; 
      #if BACKGROUND_FIELD == YES
       Bx += B0n;   /* Bx will be now the (normal) total field */
      #endif
      sBx  = (Bx > 0.0 ? 1.0 : -1.0);

      duL  = SL[i] - vL[VXn];
      duR  = SR[i] - vR[VXn];

      scrh = 1.0/(duR*uR[RHO] - duL*uL[RHO]);
      SM   = (duR*uR[MXn] - duL*uL[MXn] - ptR[i] + ptL[i])*scrh;

      pts  = duR*uR[RHO]*ptL[i] - duL*uL[RHO]*ptR[i] + 
             vL[RHO]*vR[RHO]*duR*duL*(vR[VXn]- vL[VXn]);
      pts *= scrh;

      usL[RHO] = uL[RHO]*duL/(SL[i] - SM);
      usR[RHO] = uR[RHO]*duR/(SR[i] - SM);

      sqrL = sqrt(usL[RHO]);
      sqrR = sqrt(usR[RHO]);

      S1L = SM - fabs(Bx)/sqrL;
      S1R = SM + fabs(Bx)/sqrR;

    /* -------------------------------------------------------------
        When S1L -> SL or S1R -> SR a degeneracy occurs. 
        Although Miyoshi & Kusano say that no jump exists, we don't 
        think this is actually true. 
        Indeed, vy*, vz*, By*, Bz* cannot be solved independently. 
        In this case we revert to the HLLC solver of Li (2005), 
        except for the term v.B in the region, which we compute in 
        our own way.
        Note, that by comparing the expressions of Li (2005) and 
        Miyoshi & Kusano (2005), the only change involves a 
        re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).
       ------------------------------------------------------------- */

      revert_to_hllc = 0;

      if ( (S1L - SL[i]) <  1.e-4*(SM - SL[i]) ) revert_to_hllc = 1;
      if ( (S1R - SR[i]) > -1.e-4*(SR[i] - SM) ) revert_to_hllc = 1;

      if (revert_to_hllc){

        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){  
          Uhll[nv]  = SR[i]*uR[nv] - SL[i]*uL[nv] + fL[i][nv] - fR[i][nv];
          Uhll[nv] *= scrh;
        }

        EXPAND(usL[BXn] = usR[BXn] = Uhll[BXn];   ,
               usL[BXt] = usR[BXt] = Uhll[BXt];   ,
               usL[BXb] = usR[BXb] = Uhll[BXb];)
 
        S1L = S1R = SM; /* region ** should never be computed since */ 
                        /* fluxes are given in terms of UL* and UR* */

      }else{

        scrhL = (uL[RHO]*duL*duL - Bx*Bx)/(uL[RHO]*duL*(SL[i] - SM) - Bx*Bx);
        scrhR = (uR[RHO]*duR*duR - Bx*Bx)/(uR[RHO]*duR*(SR[i] - SM) - Bx*Bx);
 
        EXPAND(usL[BXn]  = Bx1;            ,
               usL[BXt]  = uL[BXt]*scrhL;  ,
               usL[BXb]  = uL[BXb]*scrhL;)           

        #if BACKGROUND_FIELD == YES
         EXPAND(                              ;  ,
                usL[BXt]  += B0t*(scrhL - 1.0);  ,  /* Eq. [40] of         */
                usL[BXb]  += B0b*(scrhL - 1.0);)    /* Miyoshi etal (2010) */
        #endif

        EXPAND(usR[BXn] = Bx1;            ,
               usR[BXt] = uR[BXt]*scrhR;  ,   
               usR[BXb] = uR[BXb]*scrhR;)     

        #if BACKGROUND_FIELD == YES
         EXPAND(                              ;  ,
                usR[BXt]  += B0t*(scrhR - 1.0);  , /* Eq. [40] of         */
                usR[BXb]  += B0b*(scrhR - 1.0);)   /* Miyoshi etal (2010) */
        #endif
        
      }

      scrhL = Bx/(uL[RHO]*duL);
      scrhR = Bx/(uR[RHO]*duR);

      EXPAND(                                          ;  ,
             vsL = vL[VXt] - scrhL*(usL[BXt] - uL[BXt]);
             vsR = vR[VXt] - scrhR*(usR[BXt] - uR[BXt]);  ,

             wsL = vL[VXb] - scrhL*(usL[BXb] - uL[BXb]);
             wsR = vR[VXb] - scrhR*(usR[BXb] - uR[BXb]); )
         
      EXPAND(usL[MXn] = usL[RHO]*SM; 
             usR[MXn] = usR[RHO]*SM;   ,
    
             usL[MXt] = usL[RHO]*vsL;
             usR[MXt] = usR[RHO]*vsR;  ,

             usL[MXb] = usL[RHO]*wsL;
             usR[MXb] = usR[RHO]*wsR;)

      scrhL  = EXPAND(vL[VXn]*Bx1, + vL[VXt]*uL[BXt], + vL[VXb]*uL[BXb]);
      scrhL -= EXPAND(     SM*Bx1, +    vsL*usL[BXt], +    wsL*usL[BXb]);
     
      usL[ENG]  = duL*uL[ENG] - ptL[i]*vL[VXn] + pts*SM + Bx*scrhL;
      usL[ENG] /= SL[i] - SM;

      scrhR  = EXPAND(vR[VXn]*Bx1, + vR[VXt]*uR[BXt], + vR[VXb]*uR[BXb]);
      scrhR -= EXPAND(     SM*Bx1, +    vsR*usR[BXt], +    wsR*usR[BXb]);
     
      usR[ENG]  = duR*uR[ENG] - ptR[i]*vR[VXn] + pts*SM + Bx*scrhR;
      usR[ENG] /= SR[i] - SM;

      #ifdef GLM_MHD
       usL[PSI_GLM] = usR[PSI_GLM] = vL[PSI_GLM];
      #endif

  /* ------------------------------
         compute HLLD flux 
     ------------------------------ */

      if (S1L >= 0.0){       /*  ----  Region L*  ---- */

        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = fL[i][nv] + SL[i]*(usL[nv] - uL[nv]);
        }
        state->press[i] = ptL[i];

      }else if (S1R <= 0.0) {    /*  ----  Region R*  ---- */
    
        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = fR[i][nv] + SR[i]*(usR[nv] - uR[nv]);
        }
        state->press[i] = ptR[i];
         
      } else {   /* -- This state exists only if B_x != 0 -- */
      
  /* ---------------------------
           Compute U**
     --------------------------- */

        ussl[RHO] = usL[RHO];
        ussr[RHO] = usR[RHO];
 
        EXPAND(                           ,
       
               vss  = sqrL*vsL + sqrR*vsR + (usR[BXt] - usL[BXt])*sBx;       
               vss /= sqrL + sqrR;        ,
            
               wss  = sqrL*wsL + sqrR*wsR + (usR[BXb] - usL[BXb])*sBx;
               wss /= sqrL + sqrR;)
           
        EXPAND(ussl[MXn] = ussl[RHO]*SM;
               ussr[MXn] = ussr[RHO]*SM;    ,
     
               ussl[MXt] = ussl[RHO]*vss;
               ussr[MXt] = ussr[RHO]*vss;  ,
           
               ussl[MXb] = ussl[RHO]*wss;
               ussr[MXb] = ussr[RHO]*wss;)           
    
        EXPAND(ussl[BXn] = ussr[BXn] = Bx1;   ,

               ussl[BXt]  = sqrL*usR[BXt] + sqrR*usL[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
               ussl[BXt] /= sqrL + sqrR;        
               ussr[BXt]  = ussl[BXt];        ,
           
               ussl[BXb]  = sqrL*usR[BXb] + sqrR*usL[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
               ussl[BXb] /= sqrL + sqrR;        
               ussr[BXb]  = ussl[BXb];)
          
        scrhL  = EXPAND(SM*Bx1, +  vsL*usL [BXt], +  wsL*usL [BXb]);
        scrhL -= EXPAND(SM*Bx1, +  vss*ussl[BXt], +  wss*ussl[BXb]);

        scrhR  = EXPAND(SM*Bx1, +  vsR*usR [BXt], +  wsR*usR [BXb]);
        scrhR -= EXPAND(SM*Bx1, +  vss*ussr[BXt], +  wss*ussr[BXb]);

        ussl[ENG] = usL[ENG] - sqrL*scrhL*sBx;
        ussr[ENG] = usR[ENG] + sqrR*scrhR*sBx;

        #ifdef GLM_MHD
         ussl[PSI_GLM] = ussr[PSI_GLM] = vL[PSI_GLM];
        #endif
    
  /* --------------------------------------
      verify consistency condition 
     -------------------------------------- */

/*
      for (nv = 0; nv < NFLX; nv++){
        scrh = (SR[i] - S1R)*usR[nv]  + (S1R - SM)*ussr[nv] +
               (SM - S1L)*ussl[nv] + (S1L - SL[i])*usL[nv] -
               SR[i]*UR[nv] + SL[i]*UL[nv] + fr[i][nv] - fl[i][nv];

        if (fabs(scrh) > 1.e-2){
          printf (" ! Consistency condition violated, pt %d, nv %d, %12.6e \n", 
                   i,nv,scrh);
          printf ("%f %f %f %f %f %f\n",UL[BXn], usL[BXn], ussl[BXn],
                                   ussr[BXn], usR[BXn], UR[BXn]);
          printf ("%f %f %f %f %f\n",SL[i], S1L, SM, S1R, SR[i]);
          printf ("%f %f  %d\n",fr[i][nv],fl[i][nv],BXn);
          exit(1);
        }
      }
*/

        if (SM >= 0.0){           /*  ----  Region L**  ---- */

          for (nv = NFLX; nv--; ){
            state->flux[i][nv] = fL[i][nv] + S1L*(ussl[nv]  - usL[nv])
                                           + SL[i]*(usL[nv] - uL[nv]);
          }
          state->press[i] = ptL[i];
        }else{                   /*  ----  Region R**  ---- */

          for (nv = NFLX; nv--; ){
            state->flux[i][nv] = fR[i][nv] + S1R*(ussr[nv]  - usR[nv])
                                           + SR[i]*(usR[nv] - uR[nv]);
          }
          state->press[i] = ptR[i];
        }
      }
    } 
  }
}
#endif

#if EOS == ISOTHERMAL
/* ********************************************************************* */
void HLLD_Solver (const State_1D *state, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the isothermal MHD equations using the 
 * three-state HLLD Riemann solver of Mignone (2007).
 * 
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int  nv, i;
  int  revert_to_hll;
  double scrh;
  double usL[NFLX], *SL;
  double usR[NFLX], *SR;
  double usc[NFLX];
  
  double scrhL, S1L, duL;
  double scrhR, S1R, duR;
  double Bx, Bx1, SM, sBx, rho, sqrho;

  double *vL, *vR, *uL, *uR;
  static double *ptL, *ptR, *a2L, *a2R;
  static double **fL, **fR;
  static double **VL, **VR, **UL, **UR;
  #if BACKGROUND_FIELD == YES
   double B0n, B0t, B0b;
  #endif
  double **bgf;

  if (fL == NULL){
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);

    ptR = ARRAY_1D(NMAX_POINT, double);
    ptL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg, end, FACE_CENTER, grid);
  #endif

  #if DIVB_CONTROL == EIGHT_WAVES
   print ("! hlld Riemann solver does not work with Powell\n");
   QUIT_PLUTO(1);
  #endif

  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, ptL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, ptR, beg, end);

/* ----------------------------------------
      get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, bgf, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {
    
    #if BACKGROUND_FIELD == YES
     EXPAND (B0n = bgf[i][BXn];  ,
             B0t = bgf[i][BXt];  ,
             B0b = bgf[i][BXb];)
    #endif

  /* ----------------------------------------
      get max propagation speed for dt comp.
     ---------------------------------------- */             

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

/* ---------------------------------------------------------- 
                COMPUTE FLUXES and STATES
   ---------------------------------------------------------- */

    if (SL[i] >= 0.0){                     /*  ----  Region L  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = ptL[i];

    }else if (SR[i] <= 0.0) {              /*  ----  Region R   ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = ptR[i];
 
    } else {

      scrh = 1.0/(SR[i] - SL[i]);
      duL = SL[i] - vL[VXn];
      duR = SR[i] - vR[VXn];

      Bx1 = Bx = (SR[i]*vR[BXn] - SL[i]*vL[BXn])*scrh; 
      #if BACKGROUND_FIELD == YES
       Bx += B0n;   /* total field */
      #endif

      rho                = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
      state->flux[i][RHO] = (SL[i]*uR[RHO]*duR - SR[i]*uL[RHO]*duL)*scrh;
           
  /* ---------------------------
          compute S*
     --------------------------- */

      sqrho = sqrt(rho);

      SM  = state->flux[i][RHO]/rho;
      S1L = SM - fabs(Bx)/sqrho;
      S1R = SM + fabs(Bx)/sqrho;

    /* ---------------------------------------------
        Prevent degeneracies when S1L -> SL or 
        S1R -> SR. Revert to HLL if necessary.
       --------------------------------------------- */

      revert_to_hll = 0;

      if ( (S1L - SL[i]) <  1.e-4*(SR[i] - SL[i]) ) revert_to_hll = 1;
      if ( (S1R - SR[i]) > -1.e-4*(SR[i] - SL[i]) ) revert_to_hll = 1;

      if (revert_to_hll){
        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                               SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          state->flux[i][nv] *= scrh;
        }
        state->press[i] = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
        continue;
      }

      state->flux[i][MXn] = (SR[i]*fL[i][MXn] - SL[i]*fR[i][MXn] 
                            + SR[i]*SL[i]*(uR[MXn] - uL[MXn]))*scrh;

      state->press[i] = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
      #ifdef GLM_MHD
       state->flux[i][BXn]     = fL[i][BXn];
       state->flux[i][PSI_GLM] = fL[i][PSI_GLM];
      #else
       state->flux[i][BXn] = SR[i]*SL[i]*(uR[BXn] - uL[BXn])*scrh;
      #endif

  /* ---------------------------
             Compute U*  
     --------------------------- */
       
      scrhL = 1.0/((SL[i] - S1L)*(SL[i] - S1R));
      scrhR = 1.0/((SR[i] - S1L)*(SR[i] - S1R));

      #if BACKGROUND_FIELD == YES
       EXPAND(                                                              ;  ,
              usL[MXt] = rho*vL[VXt] - Bx*(uL[BXt]+B0t)*(SM - vL[VXn])*scrhL;
              usR[MXt] = rho*vR[VXt] - Bx*(uR[BXt]+B0t)*(SM - vR[VXn])*scrhR;  ,  

              usL[MXb] = rho*vL[VXb] - Bx*(uL[BXb]+B0b)*(SM - vL[VXn])*scrhL;
              usR[MXb] = rho*vR[VXb] - Bx*(uR[BXb]+B0b)*(SM - vR[VXn])*scrhR;)
      #else
       EXPAND(                                                        ;  ,
              usL[MXt] = rho*vL[VXt] - Bx*uL[BXt]*(SM - vL[VXn])*scrhL;
              usR[MXt] = rho*vR[VXt] - Bx*uR[BXt]*(SM - vR[VXn])*scrhR;  ,

              usL[MXb] = rho*vL[VXb] - Bx*uL[BXb]*(SM - vL[VXn])*scrhL;
              usR[MXb] = rho*vR[VXb] - Bx*uR[BXb]*(SM - vR[VXn])*scrhR;)
      #endif

      #if BACKGROUND_FIELD == YES
       EXPAND(                                                                ; ,
              usL[BXt] = (uL[BXt]+B0t)/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL-B0t;
              usR[BXt] = (uR[BXt]+B0t)/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR-B0t; ,
 
              usL[BXb] = (uL[BXb]+B0b)/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL-B0b;
              usR[BXb] = (uR[BXb]+B0b)/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR-B0b;)
      #else
       EXPAND(                                                      ;  ,
              usL[BXt] = uL[BXt]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL; 
              usR[BXt] = uR[BXt]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;  ,

              usL[BXb] = uL[BXb]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL;           
              usR[BXb] = uR[BXb]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;)           
      #endif

      if (S1L >= 0.0){       /*  ----  Region L*  ---- */

        EXPAND(                                                    ;  ,
          state->flux[i][MXt] = fL[i][MXt] + SL[i]*(usL[MXt] - uL[MXt]);  ,
          state->flux[i][MXb] = fL[i][MXb] + SL[i]*(usL[MXb] - uL[MXb]);  
        ) 
        EXPAND(                                                    ;  ,
          state->flux[i][BXt] = fL[i][BXt] + SL[i]*(usL[BXt] - uL[BXt]);  ,
          state->flux[i][BXb] = fL[i][BXb] + SL[i]*(usL[BXb] - uL[BXb]);  
        ) 

      }else if (S1R <= 0.0) {    /*  ----  Region R*  ---- */
    
        EXPAND(                                                    ;  ,
          state->flux[i][MXt] = fR[i][MXt] + SR[i]*(usR[MXt] - uR[MXt]);  ,
          state->flux[i][MXb] = fR[i][MXb] + SR[i]*(usR[MXb] - uR[MXb]);  
        ) 
        EXPAND(                                                    ;  ,
          state->flux[i][BXt] = fR[i][BXt] + SR[i]*(usR[BXt] - uR[BXt]);  ,
          state->flux[i][BXb] = fR[i][BXb] + SR[i]*(usR[BXb] - uR[BXb]);  
        ) 
         
      } else {
                      
       /* ---------------------------
               Compute U** = Uc
          --------------------------- */

        sBx = (Bx > 0.0 ? 1.0 : -1.0);

        EXPAND(                                                  ;  ,
               usc[MXt] = 0.5*(usR[MXt] + usL[MXt] 
                               + (usR[BXt] - usL[BXt])*sBx*sqrho);  ,     
               usc[MXb] = 0.5*(   usR[MXb] + usL[MXb] 
                               + (usR[BXb] - usL[BXb])*sBx*sqrho);)
           
        EXPAND(                                                  ;  ,
               usc[BXt] = 0.5*(   usR[BXt] + usL[BXt]  
                               + (usR[MXt] - usL[MXt])*sBx/sqrho);  ,
               usc[BXb] = 0.5*(   usR[BXb] + usL[BXb] 
                               + (usR[MXb] - usL[MXb])*sBx/sqrho);)

        EXPAND(                                               ;  ,
               state->flux[i][MXt] = usc[MXt]*SM - Bx*usc[BXt];  ,
               state->flux[i][MXb] = usc[MXb]*SM - Bx*usc[BXb]; )
        #if BACKGROUND_FIELD == YES
         EXPAND(                              ;  ,
                state->flux[i][MXt] -= Bx1*B0t;  ,
                state->flux[i][MXb] -= Bx1*B0b; )
        #endif
               
        EXPAND(                                                   ;  ,
               state->flux[i][BXt] = usc[BXt]*SM - Bx*usc[MXt]/rho;  ,
               state->flux[i][BXb] = usc[BXb]*SM - Bx*usc[MXb]/rho;)
        #if BACKGROUND_FIELD == YES
         EXPAND(                             ;  ,
                state->flux[i][BXt] += B0t*SM;  ,
                state->flux[i][BXb] += B0b*SM;)
        #endif

    /* --------------------------------------
          verify consistency condition 
       -------------------------------------- */

        #if VERIFY_CONSISTENCY_CONDITION == YES
         for (nv = NFLX; nv--; ){
           if (nv == RHO || nv == MXn || nv == BXn) continue;
           scrh = (S1L - SL[i])*usL[nv]  + (S1R - S1L)*usc[nv] +
                  (SR[i] - S1R)*usR[nv] -
                  SR[i]*uR[nv] + SL[i]*uL[nv] + fR[i][nv] - fL[i][nv];

           if (fabs(scrh) > 1.e-6){
             printf (" ! Consistency condition violated, pt %d, nv %d, %12.6e \n", 
                     i,nv,scrh);
             printf (" scrhL = %12.6e   scrhR = %12.6e\n",scrhL, scrhR);
             printf (" SL = %12.6e, S1L = %12.6e, S1R = %12.6e, SR = %12.6e\n",
                     SL[i],S1L,S1R, SR[i]);
             Show(state->vL,i);
             Show(state->vR,i);

             exit(1);
           }
         }
        #endif	  
                 
      }
    }
  }
}
#endif /* end #if on EOS  */

#undef VERIFY_CONSISTENCY_CONDITION 
