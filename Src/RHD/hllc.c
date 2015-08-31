/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implement the HLLC Riemann solver for relativistic HD.

  Solve the Riemann problem for the relativistic hydro (RHD) equations 
  using the HLLC solver of Mignone & Bodo (2005).
   
  On input, it takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics",
       Mignone and Bodo, MNRAS (2005) 364,1126.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    Dec 10, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLLC_Solver (const State_1D *state, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve the RHD Riemann problem using the HLLC Riemann solver.
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
  double scrh;

  double usl[NFLX], usr[NFLX], vm[NFLX], us, ps;
  double AL, BL, AR, BR, a,b,c;
  double vxl, vxr;
  
  double *vL, *vR, *uL, *uR;
  static double *SL, *SR;
  static double **Uhll, **Fhll;
  static double **fL, **fR;
  static double *pR, *pL;
  static double *a2L, *a2R, *hL, *hR;

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    Uhll = ARRAY_2D(NMAX_POINT, NFLX, double);
    Fhll = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR   = ARRAY_1D(NMAX_POINT, double);
    pL   = ARRAY_1D(NMAX_POINT, double);
    SR   = ARRAY_1D(NMAX_POINT, double);
    SL   = ARRAY_1D(NMAX_POINT, double);

    a2L  = ARRAY_1D(NMAX_POINT, double);
    a2R  = ARRAY_1D(NMAX_POINT, double);

    hL  = ARRAY_1D(NMAX_POINT, double);
    hR  = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

  HLL_Speed (state->vL, state->vR, a2L, a2R, SL, SR, beg, end);
  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

/* --------------------------------------------------
      compute HLL state and flux
   -------------------------------------------------- */
/*   
    scrh = 1.0/(Sr - Sl);
    for (nv = NFLX; nv--; ){  
      Uhll[i][nv]  =   Sr*ur[i][nv] - Sl*ul[i][nv] 
                     + fl[i][nv] - fr[i][nv];
      Uhll[i][nv] *= scrh;
      
      Fhll[i][nv]  =   Sl*Sr*(ur[i][nv] - ul[i][nv])
                     + Sr*fl[i][nv] - Sl*fr[i][nv];
      Fhll[i][nv] *= scrh;
    }
    Uhll[i][MXn] += (pl[i] - pr[i])*scrh;
    Fhll[i][MXn] += (Sr*pl[i] - Sl*pr[i])*scrh;
*/
/* --------------------------------------------------
        compute HLLC  flux
   -------------------------------------------------- */

    if (SL[i] >= 0.0){
    
      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = pL[i];

    }else if (SR[i] <= 0.0){

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = pR[i];

    }else{

      vL = state->vL[i]; uL = state->uL[i];
      vR = state->vR[i]; uR = state->uR[i];

#if SHOCK_FLATTENING == MULTID   
      if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){        
        scrh  = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){
          state->flux[i][nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                             +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          state->flux[i][nv] *= scrh;
        }
        state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
        continue;
      }
#endif

      vxl = vL[VXn];
      vxr = vR[VXn];
    
/* ---------------------------------------
                   get u* 
   --------------------------------------- */    

      AL = SL[i]*uL[ENG] - fL[i][ENG];
      AR = SR[i]*uR[ENG] - fR[i][ENG];
    
      BL = SL[i]*uL[MXn] - fL[i][MXn] - pL[i];
      BR = SR[i]*uR[MXn] - fR[i][MXn] - pR[i];
    
      a = AR*SL[i] - AL*SR[i];
      b = AL + BL*SR[i] - AR - BR*SL[i];
      c = BR - BL;
/*
      if (fabs(a) > 1.e-9){       
        usp = 0.5*(- b + sqrt(b*b - 4.0*a*c))/a; 
        usm = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a; 
      }else{
        usp = usm = -c/b;
      }
*/
      scrh = -0.5*(b + DSIGN(b)*sqrt(b*b - 4.0*a*c));
      us   = c/scrh;

      ps = (AL*us - BL)/(1.0 - us*SL[i]);
    
      usl[RHO] = uL[RHO]*(SL[i] - vxl)/(SL[i] - us);
      usr[RHO] = uR[RHO]*(SR[i] - vxr)/(SR[i] - us);
      EXPAND(usl[MXn] = (SL[i]*(uL[ENG] + ps) - uL[MXn])*us/(SL[i] - us); 
             usr[MXn] = (SR[i]*(uR[ENG] + ps) - uR[MXn])*us/(SR[i] - us);  ,
             usl[MXt] =  uL[MXt]*(SL[i] - vxl)/(SL[i] - us); 
             usr[MXt] =  uR[MXt]*(SR[i] - vxr)/(SR[i] - us);              ,
             usl[MXb] =  uL[MXb]*(SL[i] - vxl)/(SL[i] - us); 
             usr[MXb] =  uR[MXb]*(SR[i] - vxr)/(SR[i] - us);)
           
      usl[ENG] = uL[ENG] + (usl[MXn] - uL[MXn])/SL[i];
      usr[ENG] = uR[ENG] + (usr[MXn] - uR[MXn])/SR[i];

     /*  ----  Compute HLLC flux  ----  */

      if (us >= 0.0) {
        for (nv = NFLX; nv--;  ) {
          state->flux[i][nv] = fL[i][nv] + SL[i]*(usl[nv] - uL[nv]);
        }
        state->press[i] = pL[i];
      }else {
        for (nv = NFLX; nv--; ) {
          state->flux[i][nv] = fR[i][nv] + SR[i]*(usr[nv] - uR[nv]);
        }
        state->press[i] = pR[i];
      }
    }   /* -- end loop on speed signs  -- */
  }   /* -- end loop on points -- */
}
