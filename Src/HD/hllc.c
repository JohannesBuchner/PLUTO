/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLLC Riemann solver for MHD.

  Solve the Riemann problem for the HD equations using the 
  two-state HLLC solver by Toro.
  
  On input, this function takes left and right primitive state vectors 
  \c stateL->v and \c stateR->v at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
   -   "Riemann Solver and Numerical Methods for Fluid Dynamics"
        by E.F. Toro (Chapter 10)
       
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Oct 12, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLLC_Solver (const Sweep *sweep, int beg, int end, 
                  real *cmax, Grid *grid)
/*!
 * Solve Riemann problem using the HLLC Riemann solver.
 * 
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double **fL = stateL->flux, **fR = stateR->flux;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  double scrh, vxr, vxl;
  double usL[NFLX], usR[NFLX], vs;
  double qL, qR, wL, wR;
  double *vL, *vR, *uL, *uR, *SL, *SR;
  #if EOS == ISOTHERMAL
   double rho, mx;
  #endif

/* ----------------------------------------------------
    Compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i]  = scrh;

    if (SL[i] > 0.0){
    
      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fL[i][nv];
      sweep->press[i] = pL[i];

    }else if (SR[i] < 0.0){

      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fR[i][nv];
      sweep->press[i] = pR[i];

    }else{

      vR = stateR->v[i]; uR = stateR->u[i];
      vL = stateL->v[i]; uL = stateL->u[i];

      vxr = vR[VXn];
      vxl = vL[VXn];
 
#if SHOCK_FLATTENING == MULTID   
      if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){        
         scrh  = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           sweep->flux[i][nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           sweep->flux[i][nv] *= scrh;
         }
         sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
         continue;
       }
#endif

  /* ---------------------------------------
                   get u* 
     --------------------------------------- */    

      #if HAVE_ENERGY
       qL = vL[PRS] + uL[MXn]*(vL[VXn] - SL[i]);
       qR = vR[PRS] + uR[MXn]*(vR[VXn] - SR[i]);

       wL = vL[RHO]*(vL[VXn] - SL[i]);
       wR = vR[RHO]*(vR[VXn] - SR[i]);

       vs = (qR - qL)/(wR - wL); /* wR - wL > 0 since SL < 0, SR > 0 */
/*
      vs = vR[PRS] - vL[PRS] + uL[MXn]*(SL[i] - vxl) 
                           - uR[MXn]*(SR[i] - vxr);
      vs /= vL[RHO]*(SL[i] - vxl) - vR[RHO]*(SR[i] - vxr);
*/

       usL[RHO] = uL[RHO]*(SL[i] - vxl)/(SL[i] - vs);
       usR[RHO] = uR[RHO]*(SR[i] - vxr)/(SR[i] - vs);
       EXPAND(usL[MXn] = usL[RHO]*vs;     usR[MXn] = usR[RHO]*vs;      ,
              usL[MXt] = usL[RHO]*vL[VXt]; usR[MXt] = usR[RHO]*vR[VXt];  ,
              usL[MXb] = usL[RHO]*vL[VXb]; usR[MXb] = usR[RHO]*vR[VXb];)
           
       usL[ENG] =    uL[ENG]/vL[RHO] 
                  + (vs - vxl)*(vs + vL[PRS]/(vL[RHO]*(SL[i] - vxl)));
       usR[ENG] =    uR[ENG]/vR[RHO] 
                  + (vs - vxr)*(vs + vR[PRS]/(vR[RHO]*(SR[i] - vxr)));

       usL[ENG] *= usL[RHO];
       usR[ENG] *= usR[RHO];
      #elif EOS == ISOTHERMAL
       scrh = 1.0/(SR[i] - SL[i]);
       rho  = (SR[i]*uR[RHO] - SL[i]*uL[RHO] - fR[i][RHO] + fL[i][RHO])*scrh;
       mx   = (SR[i]*uR[MXn] - SL[i]*uL[MXn] - fR[i][MXn] + fL[i][MXn])*scrh;
       
       usL[RHO] = usR[RHO] = rho;
       usL[MXn] = usR[MXn] = mx;
       vs  = (  SR[i]*fL[i][RHO] - SL[i]*fR[i][RHO] 
              + SR[i]*SL[i]*(uR[RHO] - uL[RHO]));
       vs *= scrh;
       vs /= rho;
       EXPAND(                                            ,
              usL[MXt] = rho*vL[VXt]; usR[MXt] = rho*vR[VXt]; ,
              usL[MXb] = rho*vL[VXb]; usR[MXb] = rho*vR[VXb];)
      #endif

        
 /* ---- verify consistency condition ------ */
 /*
      for (nv = 0; nv < NFLX; nv++) {
        scrh  = (vs - SL[i])*usL[nv] + (SR[i] - vs)*usR[nv];
        scrh -= SR[i]*uR[i][nv] - SL[i]*uL[i][nv] + 
                fl[i][nv] - fr[i][nv];
          if (fabs(scrh) > 1.e-9){
          printf (" Consistency condition violated\n");
          printf ("%d %d  %12.6e\n",i,nv, scrh);
        }
      }
 */          
/*  ----  Compute HLLC flux  ----  */

      if (vs >= 0.0){
        for (nv = NFLX; nv--;   ) {
          sweep->flux[i][nv] = fL[i][nv] + SL[i]*(usL[nv] - uL[nv]);
        }
        sweep->press[i] = pL[i];
      } else {
        for (nv = NFLX; nv--;   ) {
          sweep->flux[i][nv] = fR[i][nv] + SR[i]*(usR[nv] - uR[nv]);
        }
        sweep->press[i] = pR[i];
      }
    }
  } /* end loops on points */
}
