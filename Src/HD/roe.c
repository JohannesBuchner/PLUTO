/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Roe Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the standard Roe solver with local characteristic decomposition.
  Eigenvectors are identical to the ones given in the book by Toro 
  and can also be derived from the maple script "eigenv.maple"
  in Src/HD/.
  The solver can be used for adiabatic and isothermal hydrodynamics.
  
  The macro ROE_AVERAGE specifies how the averaging process is done:
    - ROE_AVERAGE == YES use Roe average (default);
    - ROE_AVERAGE == NO  use arithmetic average.

  The macro CHECK_ROE_MATRIX can be used to verify that the 
  characteristic decomposition reproduces the Roe matrix.

  On input, it takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
   -   "Riemann Solver and Numerical Methods for Fluid Dynamics"
        by E.F. Toro (Chapter 11)
        
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Dec 10, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define ROE_AVERAGE       YES
#define CHECK_ROE_MATRIX  NO

/* ********************************************************************* */
void Roe_Solver (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve the Riemann problem using the Roe solver.
 *
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int   nv, i, j, k, nn;
  double  scrh;
  double  um[NFLX], vel2;
  double  a2, a, h;
  double  dv[NFLX], eta[NFLX];
  double  Rc[NFLX][NFLX], alambda[NFLX], lambda[NFLX];
  double  delta, delta_inv, gmm1, gmm1_inv;
#if ROE_AVERAGE == YES
  double s, c, hl, hr;
#endif
  double *ql, *qr, *uL, *uR;
  static double  **fL, **fR, *pL, *pR, *a2L, *a2R;

  double bmin, bmax, scrh1;
  double Us[NFLX];

  delta     = 1.e-7;
  delta_inv = 1.0/delta;
  #if EOS == IDEAL
   gmm1      = g_gamma - 1.0;
   gmm1_inv  = 1.0/gmm1;
  #endif

  if (fL == NULL){
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    pL  = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
  }

  for (i = NFLX; i--;  ) {
  for (j = NFLX; j--;  ) {
    Rc[i][j] = 0.0;
  }}

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

  for (i = beg; i <= end; i++)  {

    uR = state->uR[i];
    uL = state->uL[i];

#if SHOCK_FLATTENING == MULTID   

    /* ---------------------------------------------
       HLL switching function as in Quirk (1994).
       Since the problem is related to multidimensional 
       pathologies, it works in more than 1-D only.	 
       Use the HLL flux function if the interface 
       lies within a strong shock.
       The effect of this switch is visible
       in the Mach reflection test.
      --------------------------------------------- */

    if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){        
      HLL_Speed (state->vL, state->vR, a2L, a2R, 
                 &bmin - i, &bmax - i, i, i);
      a     = MAX(fabs(bmin), fabs(bmax));
      cmax[i] = a;
      bmin  = MIN(0.0, bmin);
      bmax  = MAX(0.0, bmax);
      scrh  = 1.0/(bmax - bmin);
      for (nv = NFLX; nv--; ){
        state->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                           +  bmax*fL[i][nv] - bmin*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh;
      continue;
    }
#endif

    ql = state->vL[i];
    qr = state->vR[i];

  /*  ----  Define Wave Jumps  ----  */

    for (nv = NFLX; nv--;   ) dv[nv] = qr[nv] - ql[nv];

    #if ROE_AVERAGE == YES    
     s       = sqrt(qr[RHO]/ql[RHO]);
     um[RHO]  = ql[RHO]*s;
     s       = 1.0/(1.0 + s); 
     c       = 1.0 - s;
  
     EXPAND(um[VX1] = s*ql[VX1] + c*qr[VX1];  ,
            um[VX2] = s*ql[VX2] + c*qr[VX2];  ,
            um[VX3] = s*ql[VX3] + c*qr[VX3];)

     #if EOS == IDEAL
      vel2 = EXPAND(um[VX1]*um[VX1], + um[VX2]*um[VX2], + um[VX3]*um[VX3]);

      hl  = 0.5*(EXPAND(ql[VX1]*ql[VX1], + ql[VX2]*ql[VX2], + ql[VX3]*ql[VX3]));    
      hl += a2L[i]*gmm1_inv;
     
      hr = 0.5*(EXPAND(qr[VX1]*qr[VX1], + qr[VX2]*qr[VX2], + qr[VX3]*qr[VX3]));    
      hr += a2R[i]*gmm1_inv;

      h = s*hl + c*hr;

  /* -------------------------------------------------
       the following should be  equivalent to 
    
       scrh = EXPAND(   dv[VX1]*dv[VX1],
                      + dv[VX2]*dv[VX2],
                      + dv[VX3]*dv[VX3]);

       a2 = s*a2L + c*a2R + 0.5*gmm1*s*c*scrh;

       and therefore always positive.
       just work out the coefficiendnts...
     -------------------------------------------------- */
     
      a2 = gmm1*(h - 0.5*vel2);
      a  = sqrt(a2);
     #endif
    #else
     for (nv = NFLX; nv--;   ) um[nv] = 0.5*(ql[nv] + qr[nv]);  
     #if EOS == IDEAL
      a2   = g_gamma*um[PRS]/um[RHO];
      a    = sqrt(a2);
     
      vel2 = EXPAND(um[VX1]*um[VX1], + um[VX2]*um[VX2], + um[VX3]*um[VX3]);
      h    = 0.5*vel2 + a2/gmm1;
     #endif /* EOS == IDEAL */
    #endif /* ROE_AVERAGE == YES/NO */
  
    #if EOS == ISOTHERMAL
     a2 = 0.5*(a2L[i] + a2R[i]);
     a  = sqrt(a2);
    #endif

  /* ----------------------------------------------------------------
      define non-zero components of conservative eigenvectors Rc, 
      eigenvalues (lambda) and wave strenght eta = L.du     
     ----------------------------------------------------------------  */

  /*  ---- (u - c_s)  ----  */ 

    nn         = 0;
    lambda[nn] = um[VXn] - a;
    #if EOS == IDEAL
     eta[nn] = 0.5/a2*(dv[PRS] - dv[VXn]*um[RHO]*a);
    #elif EOS == ISOTHERMAL
     eta[nn] = 0.5*(dv[RHO] - um[RHO]*dv[VXn]/a);
    #endif

    Rc[RHO][nn]        = 1.0;
    EXPAND(Rc[MXn][nn] = um[VXn] - a;   ,
           Rc[MXt][nn] = um[VXt];       ,
           Rc[MXb][nn] = um[VXb];)
    #if EOS == IDEAL
     Rc[ENG][nn] = h - um[VXn]*a;
    #endif

  /*  ---- (u + c_s)  ----  */ 

    nn         = 1;
    lambda[nn] = um[VXn] + a;
    #if EOS == IDEAL
     eta[nn]    = 0.5/a2*(dv[PRS] + dv[VXn]*um[RHO]*a);
    #elif EOS == ISOTHERMAL
     eta[nn] = 0.5*(dv[RHO] + um[RHO]*dv[VXn]/a);
    #endif

    Rc[RHO][nn]        = 1.0;
    EXPAND(Rc[MXn][nn] = um[VXn] + a;   ,
           Rc[MXt][nn] = um[VXt];       ,
           Rc[MXb][nn] = um[VXb];)
    #if EOS == IDEAL
     Rc[ENG][nn] = h + um[VXn]*a;
    #endif

  /*  ----  (u)  ----  */ 
     
    #if EOS == IDEAL
     nn         = 2;
     lambda[nn] = um[VXn];
     eta[nn]    = dv[RHO] - dv[PRS]/a2;
     Rc[RHO][nn]        = 1.0;
     EXPAND(Rc[MX1][nn] = um[VX1];   ,
            Rc[MX2][nn] = um[VX2];   ,
            Rc[MX3][nn] = um[VX3];)
     Rc[ENG][nn]        = 0.5*vel2;
    #endif
    
    #if COMPONENTS > 1

  /*  ----  (u)  ----  */ 

     nn++;
     lambda[nn] = um[VXn];
     eta[nn]    = um[RHO]*dv[VXt];
     Rc[MXt][nn] = 1.0;
     #if EOS == IDEAL
      Rc[ENG][nn] = um[VXt];  
     #endif
    #endif

    #if COMPONENTS > 2

  /*  ----  (u)  ----  */ 

     nn++;
     lambda[nn] = um[VXn];
     eta[nn]    = um[RHO]*dv[VXb];
     Rc[MXb][nn] = 1.0;
     #if EOS == IDEAL
      Rc[ENG][nn] = um[VXb];  
     #endif
    #endif

  /*  ----  get max eigenvalue  ----  */

    cmax[i] = fabs(um[VXn]) + a;
    g_maxMach = MAX(fabs(um[VXn]/a), g_maxMach);

    #if DIMENSIONS > 1
   
    /* ---------------------------------------------
         use the HLL flux function if the interface 
         lies within a strong shock.
         The effect of this switch is visible
         in the Mach reflection test.
      --------------------------------------------- */

     #if EOS == IDEAL
      scrh  = fabs(ql[PRS] - qr[PRS]);
      scrh /= MIN(ql[PRS],qr[PRS]);
     #elif EOS == ISOTHERMAL
      scrh  = fabs(ql[RHO] - qr[RHO]);
      scrh /= MIN(ql[RHO],qr[RHO]);
      scrh *= a*a;
     #endif
     if (scrh > 0.5 && (qr[VXn] < ql[VXn])){   /* -- tunable parameter -- */
       bmin = MIN(0.0, lambda[0]);
       bmax = MAX(0.0, lambda[1]);
       scrh1 = 1.0/(bmax - bmin);
       for (nv = NFLX; nv--;   ){
        state->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                          +   bmax*fL[i][nv] - bmin*fR[i][nv];
        state->flux[i][nv] *= scrh1;
       }
       state->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh1;
       continue;
     } 
    #endif

    #if CHECK_ROE_MATRIX == YES
     for (nv = 0; nv < NFLX; nv++){
       um[nv] = 0.0;
       for (k = 0; k < NFLX; k++){
       for (j = 0; j < NFLX; j++){
         um[nv] += Rc[nv][k]*(k==j)*lambda[k]*eta[j];
       }}
     }
     for (nv = 0; nv < NFLX; nv++){
       scrh = fR[i][nv] - fL[i][nv] - um[nv];
       if (nv == MXn) scrh += pR[i] - pL[i];
       if (fabs(scrh) > 1.e-6){
         print ("! Matrix condition not satisfied %d, %12.6e\n", nv, scrh);
         Show(state->vL, i);
         Show(state->vR, i);
         exit(1);
       }
     }
    #endif

  /* -----------------------------------------------------------
                      compute Roe flux 
     ----------------------------------------------------------- */
      
    for (nv = NFLX; nv--;   ) alambda[nv]  = fabs(lambda[nv]);

  /*  ----  entropy fix  ----  */

    if (alambda[0] <= delta) {
      alambda[0] = 0.5*lambda[0]*lambda[0]/delta + 0.5*delta;
    }
    if (alambda[1] <= delta) {
      alambda[1] = 0.5*lambda[1]*lambda[1]/delta + 0.5*delta;
    }

    for (nv = NFLX; nv--;   ) {
      state->flux[i][nv] = fL[i][nv] + fR[i][nv];
      for (k  = NFLX; k-- ;   ) {
        state->flux[i][nv] -= alambda[k]*eta[k]*Rc[nv][k];
      }
      state->flux[i][nv] *= 0.5;
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
#undef ROE_AVERAGE        
#undef CHECK_ROE_MATRIX
