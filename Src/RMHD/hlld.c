/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLLD Riemann solver for the relativistic MHD equations.

  Solve the Riemann problem for the relativistic MHD equations 
  using the HLLD solver of Mignone, Ugliano & Bodo (2009).
  The solver requires the solution of a nonlinear equation (Eq. [48]) 
  and the maximum number of iterations before convergence can be 
  set using the macro MAX_ITER (default 20).
  A physical relevant solution is accepted if it satisfies
  the constraints by Eq. [54], implemented in the function HLLD_Fstar().
  If this is not the case or if other failure conditions occur during 
  the iteration loop, we switch to the HLL Riemann solver.
  
  The macro ::COUNT_FAILURES can be set to YES if one wishes to count
  how many times the solver fails.
  
  On input, it takes left and right primitive sweep vectors 
  \c sweep->vL and \c sweep->vR at zone edge i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "A five wave Harte-Lax-van Leer Riemann solver for relativistic
       magnetohydrodynamics" Mignone et al, MNRAS (2009) 393,1141.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    Dec 10, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"
#define MAX_ITER  20
#define COUNT_FAILURES  NO /**< When set to YES, count number of failures and 
        write the count to "hlld_fails.dat" (works in parallel as well).
        The number of failures is normalized to the total number of zones: 
        \c nf/nz where \c nz is incremented during the main loop and \c nf
        is incremented when a failure occur. */
 
typedef struct RIEMANN_STATE{
  int    fail;
  double vx, vy, vz;
  double Bx, By, Bz;
  double Kx, Ky, Kz, K2;
  double w, sw, p, rho;
  double u[NFLX], R[NVAR];
  double S, Sa;
/*
  double fun1, fun2, fun3; 
  double denL, denR;
*/
} Riemann_State;


static double Sc, Bx;
static double HLLD_Fstar (Riemann_State *, Riemann_State *, double);

static int  HLLD_GetRiemannState (Riemann_State *, double, int);
static void HLLD_GetAState (Riemann_State *,  double p);
static void HLLD_GetCState (Riemann_State *, Riemann_State *, double, double *);

static double HLLD_TotalPressure  (double *);

static void HLLD_PrintStates(double *, double *);
static void HLLD_PrintWhatsWrong(Riemann_State *, Riemann_State *, 
                                 double, double *, double *);

#define HLLD_DEBUG  NO
/* ********************************************************************* */
void HLLD_Solver (const Sweep *sweep, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve the Riemann problem using the HLLD Riemann solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i, k;
  int    switch_to_hll;
  double   scrh;
  static double **Uhll, **Fhll, **Vhll;
  double *vL, *vR, *fL, *fR, *uL, *uR, *SL, *SR;
 
  double p0, f0, p, f, dp, dS_1;
  double Uc[NVAR];
  Riemann_State PaL, PaR;
  double pguess;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double **fluxL = stateL->flux, **fluxR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *hL = stateL->h,     *hR = stateR->h;
  double  *pL = stateL->prs,   *pR = stateR->prs;

/* -------------------------------------------------------
     for information purposes, show how many failures
     have been encountered in every step
   ------------------------------------------------------- */

  #if COUNT_FAILURES == YES
   static double totfail=0.0,totzones=1.0; 
   static int last_step = -1;
   FILE *fp;   

   if (g_stepNumber == 0 && prank == 0) {
     fp = fopen("hlld_fails.dat","w");
     fprintf (fp, "#  Step    Failures (x 100) \n");
     fprintf (fp, "# --------------------------\n");
     fclose(fp);
   }
   if (last_step != g_stepNumber){ /* -- at the beginning of new step -- */
     #ifdef PARALLEL
      MPI_Allreduce (&totfail, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      totfail = scrh;
      MPI_Allreduce (&totzones, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      totzones = scrh;
     #endif
     if (prank == 0){
       fp = fopen("hlld_fails.dat","a+");
       fprintf (fp,"    %d      %f\n",g_stepNumber, totfail/totzones*100.);
       fclose(fp);
     }
     totfail = totzones = 0.0; 
     last_step = g_stepNumber;
   } 
  #endif
 
/* ----------------------------------------------------------- */

  if (Uhll == NULL){
    Uhll = ARRAY_2D(NMAX_POINT, NVAR, double);
    Fhll = ARRAY_2D(NMAX_POINT, NVAR, double);
    Vhll = ARRAY_2D(NMAX_POINT, NVAR, double);
  }
/*
  #if DIVB_CONTROL == EIGHT_WAVES
   print ("! hlld Riemann solver does not work with Powell's 8-wave\n");
   QUIT_PLUTO(1);
  #endif
*/

#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* -- compute speeds -- */

  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

/* -------------------------------------------------------
               Begin main loop 
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++) {

#if HLLD_DEBUG == YES
    if (!(grid[IDIR].x[i] < 0.5 && grid[IDIR].x[i+1] > 0.5)) continue;
#endif

    #if COUNT_FAILURES == YES
     totzones += 1.0;
    #endif

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    vL = stateL->v[i]; uL = stateL->u[i]; fL = stateL->flux[i];
    vR = stateR->v[i]; uR = stateR->u[i]; fR = stateR->flux[i];

    if (SL[i] >= 0.0){

      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fL[nv];
      sweep->press[i] = pL[i];

    }else if (SR[i] <= 0.0){

      for (nv = NFLX; nv--; ) sweep->flux[i][nv] = fR[nv];
      sweep->press[i] = pR[i];

    }else{

  /* ---- build the HLL average sweep ---- */

      dS_1 = 1.0/(SR[i] - SL[i]);
      for (nv = NFLX; nv--;  ){  /* -- we use NVAR and not NFLX  since 
                                       ConsToPrim may need entropy  -- */
        Uhll[i][nv]  = SR[i]*uR[nv] - SL[i]*uL[nv]  + fL[nv] - fR[nv];
        Uhll[i][nv] *= dS_1;

        Fhll[i][nv]  = SR[i]*fL[nv] - SL[i]*fR[nv] + SL[i]*SR[i]*(uR[nv] - uL[nv]);
        Fhll[i][nv] *= dS_1;
      }
      Uhll[i][MXn] += (pL[i] - pR[i])*dS_1;
#if NSCL > 0 
      NSCL_LOOP(nv) {
        double vxR = vR[VXn];
        double vxL = vL[VXn];
        Uhll[i][nv]  = (SR[i] - vxR)*uR[nv] - (SL[i] - vxL)*uL[nv];
        Uhll[i][nv] *= dS_1;

        Fhll[i][nv]  =   SR[i]*uL[nv]*vL[VXn] - SL[i]*uR[nv]*vR[VXn] 
                       + SL[i]*SR[i]*(uR[nv] - uL[nv]);
        Fhll[i][nv] *= dS_1;
      }
#endif

  /* ---- revert to HLL in proximity of strong shocks ---- */

#if SHOCK_FLATTENING == MULTID
      if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){        
        for (nv = NFLX; nv--; ) sweep->flux[i][nv] = Fhll[i][nv];
        sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*dS_1;
        continue;
      }
#endif

  /* ---- proceed normally otherwise ---- */

      PaL.S = SL[i];
      PaR.S = SR[i];
      Bx    = Uhll[i][BXn];
      for (nv = 0; nv < NFLX; nv++){
        PaL.R[nv] = SL[i]*uL[nv] - fL[nv];
        PaR.R[nv] = SR[i]*uR[nv] - fR[nv];
      }
      PaL.R[MXn] -= pL[i];
      PaR.R[MXn] -= pR[i];
      #if RMHD_REDUCED_ENERGY == YES
       PaL.R[ENG] += PaL.R[RHO];
       PaR.R[ENG] += PaR.R[RHO];
      #endif

  /* ---- provide an initial guess ---- */

      scrh = MAX(pL[i], pR[i]);
      if (Bx*Bx/scrh < 0.01) { /* -- try the B->0 limit -- */

        double a,b,c;
        a = SR[i] - SL[i];
        b = PaR.R[ENG] - PaL.R[ENG] + SR[i]*PaL.R[MXn] - SL[i]*PaR.R[MXn];
        c = PaL.R[MXn]*PaR.R[ENG] - PaR.R[MXn]*PaL.R[ENG];
        scrh = b*b - 4.0*a*c;
        scrh = MAX(scrh,0.0);
        p0 = 0.5*(- b + sqrt(scrh))*dS_1;

      }else{  /* ----  use HLL average ---- */
                         
        ConsToPrim(Uhll, Vhll, i, i, sweep->flag);
        p0    = HLLD_TotalPressure(Vhll[i]);
      }
      
  /* ---- check if guess makes sense ---- */

      pguess = p0;
      switch_to_hll = 0;
      f0 = HLLD_Fstar(&PaL, &PaR, p0);
      if (f0 != f0 || PaL.fail) switch_to_hll = 1;

  /* ---- Root finder ---- */

      k = 0;
      if (fabs(f0) > 1.e-12 && !switch_to_hll){
        p  = 1.025*p0; f  = f0;
        for (k = 1; k < MAX_ITER; k++){

#if HLLD_DEBUG == YES
 printf ("k = %d, p = %12.6e,  f = %12.6e\n",k,p,f);
#endif

          f  = HLLD_Fstar(&PaL, &PaR, p);
          if (f != f  || PaL.fail || (k > 7) || (fabs(f) > fabs(f0) && k > 4)) {
            switch_to_hll = 1;
            break;
          }

          dp = (p - p0)/(f - f0)*f;

          p0 = p; f0 = f;
          p -= dp;
          if (p < 0.0) p = 1.e-6;
          if (fabs(dp) < 1.e-5*p || fabs(f) < 1.e-6) break;
        }
      }else p = p0;

#if HLLD_DEBUG == YES
HLLD_PrintWhatsWrong(&PaL, &PaR, phll, phll, p0Bx, p, vL, vR);
#endif

    /* ----  too many iter ? --> use HLL ---- */

      if (PaL.fail) switch_to_hll = 1;
      if (switch_to_hll){

        #if COUNT_FAILURES == YES
         totfail += 1.0;
        #endif

        for (nv = NFLX; nv--; ) sweep->flux[i][nv] = Fhll[i][nv];
        sweep->press[i]  = (SR[i]*pL[i] - SL[i]*pR[i])*dS_1;
        continue;
      }

  /* -- ok, solution should be reliable -- */

      g_maxRiemannIter = MAX(g_maxRiemannIter, k);

      #ifdef GLM_MHD
       PaL.u[PSI_GLM] = PaR.u[PSI_GLM] = Uc[PSI_GLM] = uL[PSI_GLM];
      #endif

      if (PaL.Sa >= -1.e-6){      

        HLLD_GetAState (&PaL, p);
        #if RMHD_REDUCED_ENERGY == YES
         PaL.u[ENG] -= PaL.u[RHO];
        #endif
        for (nv = NFLX; nv--;   ) {
          sweep->flux[i][nv] = fL[nv] + SL[i]*(PaL.u[nv] - uL[nv]);
        }
        sweep->press[i] = pL[i];

      }else if (PaR.Sa <= 1.e-6){

        HLLD_GetAState (&PaR, p);
        #if RMHD_REDUCED_ENERGY == YES
         PaR.u[ENG] -= PaR.u[RHO];
        #endif
        for (nv = NFLX; nv--;   ) {
          sweep->flux[i][nv] = fR[nv] + SR[i]*(PaR.u[nv] - uR[nv]);
        }
        sweep->press[i] = pR[i];

      }else{

        HLLD_GetCState (&PaL, &PaR, p, Uc);
        if (Sc > 0.0){
          #if RMHD_REDUCED_ENERGY == YES
           PaL.u[ENG] -= PaL.u[RHO];
           Uc[ENG]    -= Uc[RHO];
          #endif
          for (nv = NFLX; nv--;   ) {
            sweep->flux[i][nv] = fL[nv] + SL[i]*(PaL.u[nv] - uL[nv]) 
                                        + PaL.Sa*(Uc[nv] - PaL.u[nv]);
          }
          sweep->press[i] = pL[i];

        }else{
          #if RMHD_REDUCED_ENERGY == YES
           PaR.u[ENG] -= PaR.u[RHO];
           Uc[ENG]    -= Uc[RHO];
          #endif
          for (nv = NFLX; nv--;   ) {
            sweep->flux[i][nv] = fR[nv] + SR[i]*(PaR.u[nv] - uR[nv]) 
                                        + PaR.Sa*(Uc[nv] - PaR.u[nv]);
          }
          sweep->press[i] = pR[i];
        }  
      }
    } /* --- end if on SL, SR -- */
  } /* --- end loop on points -- */

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
 
  #if DIVB_CONTROL == EIGHT_WAVES
   HLL_DIVB_SOURCE (sweep, Uhll, beg + 1, end, grid);
  #endif

}
#undef MAX_ITER
#undef COUNT_FAILURES
/* ********************************************************************* */
double HLLD_Fstar (Riemann_State *PaL, Riemann_State *PaR, double p)
/*!
 * Return the velocity difference across the contact mode as a 
 * function of the total pressure p.
 *
 *
 *********************************************************************** */
{
  int    success = 1;
  double dK, Bxc, Byc, Bzc;
  double sBx, fun;
  double vxcL, KLBc;
  double vxcR, KRBc;

  sBx = Bx > 0.0 ? 1.0:-1.0;

  success *= HLLD_GetRiemannState (PaL, p, -1);
  success *= HLLD_GetRiemannState (PaR, p,  1);

/* -- compute B from average sweep -- */

  dK  = PaR->Kx - PaL->Kx + 1.e-12;

  EXPAND(Bxc = Bx*dK;                            ,

         Byc =   PaR->By*(PaR->Kx - PaR->vx) 
               - PaL->By*(PaL->Kx - PaL->vx)
               + Bx*(PaR->vy - PaL->vy);         ,

         Bzc =   PaR->Bz*(PaR->Kx - PaR->vx) 
               - PaL->Bz*(PaL->Kx - PaL->vx)
               + Bx*(PaR->vz - PaL->vz);)
  
  KLBc = EXPAND(PaL->Kx*Bxc, + PaL->Ky*Byc, + PaL->Kz*Bzc);
  KRBc = EXPAND(PaR->Kx*Bxc, + PaR->Ky*Byc, + PaR->Kz*Bzc);

  vxcL = PaL->Kx - dK*Bx*(1.0 - PaL->K2)/(PaL->sw*dK - KLBc);
  vxcR = PaR->Kx - dK*Bx*(1.0 - PaR->K2)/(PaR->sw*dK - KRBc);

  PaL->Sa = PaL->Kx;
  PaR->Sa = PaR->Kx;
  Sc      = 0.5*(vxcL + vxcR);
  fun     = vxcL - vxcR;   /* -- essentially, Eq. [48] -- */

/*
fun = dK*(1.0 - Bx*(  (1.0 - PaR->K2)/(PaR->sw*dK - KRBc)
                    - (1.0 - PaL->K2)/(PaL->sw*dK - KLBc)) );
*/

  /* -- check if sweep makes physically sense -- */

  success  = (vxcL - PaL->Kx)  > -1.e-6;
  success *= (PaR->Kx  - vxcR) > -1.e-6;

  success *= (PaL->S - PaL->vx) < 0.0;
  success *= (PaR->S - PaR->vx) > 0.0;

  success *= (PaR->w - p) > 0.0;
  success *= (PaL->w - p) > 0.0;
  success *= (PaL->Sa - PaL->S)  > -1.e-6;
  success *= (PaR->S  - PaR->Sa) > -1.e-6;

  PaL->fail = !success;

/*
scrh  = (1.0 - PaR->K2)*(PaL->sw*dK - KLBc);
scrh -= (1.0 - PaL->K2)*(PaR->sw*dK - KRBc);

PaL->fun1 = (PaR->sw*dK - KRBc)*(PaL->sw*dK - KLBc) - Bx*scrh; 
PaL->fun2 = scrh;
PaL->denL = (PaL->sw*dK - KLBc);
PaL->denR = (PaR->sw*dK - KRBc);
*/
  return (fun);
}

/* ********************************************************************* */
int HLLD_GetRiemannState (Riemann_State *Pv, double p, int side)
/*!
 * Express the sweep behind a wave as function of the total
 * pressure p and the right hand side on the other side of the wave.
 *  
 * On output, return 1 if succesful, 0 if w < 0 is encountered.
 * side = -1 : means left
 * side =  1 : means right
 * 
 *********************************************************************** */
{
  double A, C, G, X, s;
  double vx, vy, vz, scrh, S, *R;

  S = Pv->S;
  R = Pv->R;

  A = R[MXn] + p*(1.0 - S*S) - S*R[ENG];             /* Eq. [26] */      
  G = EXPAND(0.0, + R[BXt]*R[BXt], + R[BXb]*R[BXb]); /* Eq. [27] */
  C = EXPAND(0.0, + R[BXt]*R[MXt], + R[BXb]*R[MXb]); /* Eq. [28] */
  X = Bx*(A*S*Bx + C) - (A + G)*(S*p + R[ENG]);      /* Eq. [30] */

/* -- compute the numerators of Eqs. [23,23,45] -- */

  EXPAND(vx = ( Bx*(A*Bx + C*S) - (R[MXn] + p)*(G + A) );   , 

         vy = ( - (A + G - Bx*Bx*(1.0 - S*S))*R[MXt]     
                  + R[BXt]*(C + Bx*(S*R[MXn] - R[ENG])) );  , 
 
         vz = ( - (A + G - Bx*Bx*(1.0 - S*S))*R[MXb]          
                 + R[BXb]*(C + Bx*(S*R[MXn] - R[ENG])) );)

  scrh = EXPAND(vx*R[MXn], + vy*R[MXt], + vz*R[MXb]);
  scrh = X*R[ENG] - scrh;
  Pv->w = p + scrh/(X*S - vx);  /* Eq [31] */

  if (Pv->w < 0.0) return(0);  /* -- failure -- */

  EXPAND(Pv->vx = vx/X;  , 
         Pv->vy = vy/X;  ,
         Pv->vz = vz/X;)
         
/*  -- original Eq. [21] -- */
/*
  EXPAND(Pv->Bx = Bx;                            , 
         Pv->By = (R[BXt]*X - Bx*vy)/(X*S - vx);  ,
         Pv->Bz = (R[BXb]*X - Bx*vz)/(X*S - vx);)
*/
/* --------------------------------------------------------------------- */
/*! When computing Bx, By and Bz, we use Eq. [21] with 
    vx, vy and vz replaced by by Eq. [23,24,25]. 
    The following short MAPLE script is used to verify the 
    correctness.
 \code
   restart;
   A    := R[mx] - lambda*R[E] + p*(1 - lambda^2):
   G    := R[By]*R[By] + R[Bz]*R[Bz]:
   C    := R[my]*R[By] + R[mz]*R[Bz]:
   Q    := - A - G + B[x]*B[x]*(1-lambda^2):
   X    := B[x]*(A*lambda*B[x]+ C) - (A + G)*(lambda*p + R[E]):
   v[x] := (B[x]*(A*B[x] + lambda*C) - (A + G)*(p + R[mx]))/X:
   v[y] := (Q*R[my] + R[By]*(C + B[x]*(lambda*R[mx] - R[E])))/X:
   B[y] := (R[By] - B[x]*v[y])/(lambda - v[x]):
   simplify(B[y]);
 \endcode
*/
/* --------------------------------------------------------------------- */

  EXPAND(Pv->Bx = Bx;                                      , 
         Pv->By = -(R[BXt]*(S*p + R[ENG]) - Bx*R[MXt])/A;  ,
         Pv->Bz = -(R[BXb]*(S*p + R[ENG]) - Bx*R[MXb])/A;)
  
  s  = Bx > 0.0 ? 1.0:-1.0;
  if (side < 0) s *= -1.0;

  Pv->sw = s*sqrt(Pv->w);

  scrh = 1.0/(S*p +  R[ENG] + Bx*Pv->sw);
  EXPAND(Pv->Kx = scrh*(R[MXn] + p + R[BXn]*Pv->sw);  ,
         Pv->Ky = scrh*(R[MXt]     + R[BXt]*Pv->sw);  ,
         Pv->Kz = scrh*(R[MXb]     + R[BXb]*Pv->sw);)

  Pv->K2 = EXPAND(Pv->Kx*Pv->Kx, + Pv->Ky*Pv->Ky, + Pv->Kz*Pv->Kz);
  return(1); /* -- success -- */
}
/* ********************************************************************* */
void HLLD_GetAState (Riemann_State *Pa,  double p)
/*!
 *  Compute sweeps aL and aR behind fast waves.
 *
 *********************************************************************** */
{
  double vB, *ua, *R, S;
  double scrh;
  
  ua = Pa->u;
  S  = Pa->S;
  R  = Pa->R;

  scrh = 1.0/(S - Pa->vx);

  ua[RHO] =  R[RHO]*scrh;
  EXPAND(ua[BXn] =  Bx;                       ,
         ua[BXt] = (R[BXt] - Bx*Pa->vy)*scrh;  ,
         ua[BXb] = (R[BXb] - Bx*Pa->vz)*scrh;)

  vB     = EXPAND(Pa->vx*ua[BXn], + Pa->vy*ua[BXt], + Pa->vz*ua[BXb]);
  ua[ENG] = (R[ENG] + p*Pa->vx - vB*Bx)*scrh;

  EXPAND(ua[MXn] = (ua[ENG] + p)*Pa->vx - vB*ua[BXn];  ,
         ua[MXt] = (ua[ENG] + p)*Pa->vy - vB*ua[BXt];  ,
         ua[MXb] = (ua[ENG] + p)*Pa->vz - vB*ua[BXb];)
}

/* ********************************************************************* */
void HLLD_GetCState (Riemann_State *PaL, Riemann_State *PaR, double p,
                     double *Uc)
/*!
 *  Compute sweeps cL and cR across contact mode.
 *
 *************************************************************** */
{
  double dK, *ua;
  double vxcL, vycL, vzcL, KLBc;
  double vxcR, vycR, vzcR, KRBc;
  double vxc, vyc, vzc, vBc;
  double Bxc, Byc, Bzc, Sa, vxa;
  double scrhL, scrhR;

  HLLD_GetAState (PaL, p);
  HLLD_GetAState (PaR, p);
  dK = (PaR->Kx - PaL->Kx) + 1.e-12;

  EXPAND(Bxc = Bx*dK;                          ,

         Byc =   PaR->By*(PaR->Kx - PaR->vx)    /* Eq. [45] */
               - PaL->By*(PaL->Kx - PaL->vx)
                    + Bx*(PaR->vy - PaL->vy);  ,

         Bzc =   PaR->Bz*(PaR->Kx - PaR->vx)   /* Eq.[45] */
               - PaL->Bz*(PaL->Kx - PaL->vx)
                    + Bx*(PaR->vz - PaL->vz);)
   
  EXPAND(Bxc  = Bx;  ,
         Byc /= dK;  ,
         Bzc /= dK;)

  EXPAND(Uc[BXn] = Bxc;  ,
         Uc[BXt] = Byc;  ,
         Uc[BXb] = Bzc;)

  KLBc = EXPAND(PaL->Kx*Bxc, + PaL->Ky*Byc, + PaL->Kz*Bzc);
  KRBc = EXPAND(PaR->Kx*Bxc, + PaR->Ky*Byc, + PaR->Kz*Bzc);

  scrhL = (1.0 - PaL->K2)/(PaL->sw - KLBc);
  scrhR = (1.0 - PaR->K2)/(PaR->sw - KRBc);

  EXPAND(vxcL = PaL->Kx - Uc[BXn]*scrhL;     /* Eq [47] */
         vxcR = PaR->Kx - Uc[BXn]*scrhR;  ,

         vycL = PaL->Ky - Uc[BXt]*scrhL;  
         vycR = PaR->Ky - Uc[BXt]*scrhR;  ,

         vzcL = PaL->Kz - Uc[BXb]*scrhL;
         vzcR = PaR->Kz - Uc[BXb]*scrhR;)

  EXPAND(vxc = 0.5*(vxcL + vxcR);  ,  
         vyc = 0.5*(vycL + vycR);  ,
         vzc = 0.5*(vzcL + vzcR);)

  if (vxc > 0.0) {
    HLLD_GetAState (PaL, p);
    ua  = PaL->u;
    Sa  = PaL->Sa;
    vxa = PaL->vx;
  } else {
    HLLD_GetAState (PaR, p);
    ua  = PaR->u;
    Sa  = PaR->Sa;
    vxa = PaR->vx;
  }
  
  vBc = EXPAND(vxc*Uc[BXn], + vyc*Uc[BXt], + vzc*Uc[BXb]);

  Uc[RHO] = ua[RHO]*(Sa - vxa)/(Sa - vxc);                      /* Eq [32] */
  Uc[ENG] = (Sa*ua[ENG] - ua[MXn] + p*vxc - vBc*Bx)/(Sa - vxc); /* Eq [33] */

  EXPAND(Uc[MXn] = (Uc[ENG] + p)*vxc - vBc*Bx;       ,          /* Eq [34] */
         Uc[MXt] = (Uc[ENG] + p)*vyc - vBc*Uc[BXt];  ,
         Uc[MXb] = (Uc[ENG] + p)*vzc - vBc*Uc[BXb];)
}

/* ********************************************************************* */
double HLLD_TotalPressure (double *v)
/*!
 * Compute total pressure
 *
 *********************************************************************** */
{
  double vel2, Bmag2, vB, lor2;
  double pt;

  Bmag2 = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
  vel2  = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
  vB    = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
  pt   = v[PRS] + 0.5*(Bmag2*(1.0 - vel2) + vB*vB);
  return(pt);  
}



void HLLD_PrintStates(double *vL, double *vR)
{
  printf ("RHO_LEFT  %18.9e\n",vL[RHO]);
  printf ("VX_LEFT   %18.9e\n",vL[VXn]);
  printf ("VY_LEFT   %18.9e\n",vL[VXt]);
  #if COMPONENTS == 3
   printf ("VZ_LEFT   %18.9e \n" ,vL[VXb]);
  #else
   printf ("VZ_LEFT   %18.9e \n" ,0.0);
  #endif
  printf ("BY_LEFT   %18.9e\n",vL[BXt]);
  #if COMPONENTS == 3
   printf ("BZ_LEFT   %18.9e \n" ,vL[BXb]);
  #else
   printf ("BZ_LEFT   %18.9e \n" ,0.0);
  #endif
  printf ("PR_LEFT   %18.9e \n" ,vL[PRS]);

  printf ("RHO_RIGHT  %18.9e\n",vR[RHO]);
  printf ("VX_RIGHT   %18.9e\n",vR[VXn]);
  printf ("VY_RIGHT   %18.9e\n",vR[VXt]);
  #if COMPONENTS == 3
   printf ("VZ_RIGHT   %18.9e \n", vR[VXb]);
  #else
   printf ("VZ_RIGHT   %18.9e \n", 0.0);
  #endif
  printf ("BY_RIGHT   %18.9e\n", vR[BXt]);
  #if COMPONENTS == 3
   printf ("BZ_RIGHT   %18.9e\n", vR[BXb]);
  #else
   printf ("BZ_RIGHT   %18.9e\n",0.0);
  #endif
  printf ("PR_RIGHT   %18.9e \n",vR[PRS]);
  printf ("BX_CONST   %18.9e \n",vR[BXn]);
}


void HLLD_PrintWhatsWrong(Riemann_State *PaL, Riemann_State *PaR, 
                      double pguess, double *vL, double *vR)
{
  double f,p;
  FILE *fp;

  HLLD_PrintStates(vL, vR);

  HLLD_GetAState (PaL, pguess);
  HLLD_GetAState (PaR, pguess);

  printf ("pguess  = %f, F = %f\n",pguess, HLLD_Fstar(PaL, PaR, pguess));
  printf ("S = %f  %f  %f  %f  %f\n", PaL->S, PaL->Sa, Sc, PaR->Sa, PaR->S);

  fp = fopen("new.dat","w");
  for (p = 0.01*pguess; p <= 10.0*pguess; p += 1.e-4*pguess){
    f = HLLD_Fstar(PaL, PaR, p);
    fprintf (fp,"%f  %f %f  %f\n", p, f, Sc - PaL->Sa, PaR->Sa - Sc);
  }
  fclose(fp);
  exit(1);
}


