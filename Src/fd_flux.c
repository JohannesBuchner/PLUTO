/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute finite difference fluxes.
  
  The function FD_Flux() computes a high-order approximation to 
  the flux divergence using WENO or Monotonicity-preserving 
  reconstruction in characteristic variables.
  
  It is called as a regular Riemann solver function except 
  that the left and right states do not need to be defined 
  since high-order reconstruction is done on the characteristic
  projection of the cell-centered fluxes.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Fen 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void FD_Flux (const Sweep *sweep, int beg, int end, 
              double *cmax, Grid *grid)
/*!
 * Compute interface flux by suitable high-order finite
 * difference non-oscillatory interpolants.
 *
 * \param [in]  sweep   pointer to Sweep structure
 * \param [in]    beg   initial index of computation
 * \param [in]    end   final   index of computation
 * \param [out]  cmax   array of maximum characteristic speeds
 * \param [in]    grid  pointer to Grid structure
 *
 * \b Reference
 *    - "High-order conservative finite difference GLM-MHD schemes
 *     for cell-centered MHD"\n
 *      Mignone et al., JCP (2010), 229, 5896
 *
 * \return This function has no return value.    
 * 
 ****************************************************************** */
{
  int   nv, i, j, k, np_tot;
  int    Kmax, S;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);
  static State *stateLR;

  double **v = stateC->v;
  double **u = stateC->u;
  double **F = stateC->flux;
  double  *p = stateC->prs;
  double **Uave, **Vave;

  double **flux, *press;
  double fs, fp, fm, flx[NFLX], dx;

  static double *Fp, *Fm;
  static double *psim, *Bm, **lp;
  double **L, **R;
  double (*REC)(double *, double, int);

/* -----------------------------------------------
   0. Check compatibility with other modules 
   ----------------------------------------------- */

#if TIME_STEPPING != RK3 && TIME_STEPPING != SSP_RK4
  print ("! Finite Difference schemes work with RK3 only \n");
  QUIT_PLUTO(1);
#endif
#if GEOMETRY != CARTESIAN
  print ("! Finite Difference schemes work in Cartesian coordinates only\n");
  QUIT_PLUTO(1);
#endif
#if PHYSICS == RMHD || PHYSICS == RHD
  print ("! Finite difference schemes work only for HD od MHD modules\n");
  QUIT_PLUTO(1);
#endif   
#ifdef STAGGERED_MHD
  print ("! Finite difference schemes work only with cell-centered schemes\n");
  QUIT_PLUTO(1);
#endif   
#if PHYSICS == MHD && BACKGROUND_FIELD == YES
  print ("! Background field splitting not supported with FD schemes\n");
  QUIT_PLUTO(1);
#endif  

/* --------------------------------------------------------------
   1. Define pointer to reconstruction function and 
      compute the stencil for interpolation.
      For given S, the stencil is: i-S <= i <= i+S+1
   -------------------------------------------------------------- */

  dx = grid->dx[g_dir][beg];
#if RECONSTRUCTION == WENO3_FD
  REC = WENO3_Reconstruct;
  S   = 1;
#elif RECONSTRUCTION == LIMO3_FD
  REC = LIMO3_Reconstruct;
  S   = 1;
#elif RECONSTRUCTION == WENOZ_FD
  REC = WENOZ_Reconstruct;
  S   = 2;
#elif RECONSTRUCTION == MP5_FD
  REC = MP5_Reconstruct;
  S   = 2;
#endif
  
/* ----------------------------------------------------
   2. Allocate memory / set pointer shortcuts
   ---------------------------------------------------- */

  if (Fp == NULL){
    stateLR = malloc(sizeof(struct State_));
    StateStructAllocate (stateLR);

    Fp   = ARRAY_1D(NMAX_POINT, double);
    Fm   = ARRAY_1D(NMAX_POINT, double);
    lp   = ARRAY_2D(NMAX_POINT, NVAR, double);
    psim = ARRAY_1D(NMAX_POINT, double);
    Bm   = ARRAY_1D(NMAX_POINT, double);
  }

  np_tot = grid->np_tot[g_dir];
  flux   = sweep->flux;
  press  = sweep->press;
  Uave   = stateLR->u;
  Vave   = stateLR->v;
  memset ((void *)stateLR->Rp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)stateLR->Lp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));

/* ---------------------------------------------------------------------
   3. Compute cell-centered fluxes. Pressure is added to the normal
      component of momentum.
      To save computational time, we also solve the linear GLM-MHD
      subsystem here and use the solution in the nonlinear fluxes later.
   --------------------------------------------------------------------- */

  PrimToCons (v, u, 0, np_tot - 1);
  SoundSpeed2 (stateC, 0, np_tot-1, CELL_CENTER, grid);
  Flux (stateC, 0, np_tot-1);

  for (i = 0; i < np_tot; i++){
    F[i][MXn] += p[i];
    p[i]       = press[i] = 0.0;
    #if EOS != ISOTHERMAL
    fs = g_gamma*v[i][PRS]/v[i][RHO];
    #else
    fs = g_isoSoundSpeed*g_isoSoundSpeed;
    #endif 
    fs = fabs(v[i][VXn])/sqrt(fs);
    g_maxMach = MAX(g_maxMach, fs);
    #ifdef GLM_MHD
    j = np_tot - 1 - i;
    Fp[i] = 0.5*(u[i][PSI_GLM] + glm_ch*u[i][BXn]);
    Fm[j] = 0.5*(u[i][PSI_GLM] - glm_ch*u[i][BXn]);
    #endif
  }

#ifdef GLM_MHD
  for (i = beg; i <= end; i++){
    fp = REC(Fp, dx, i); 
    fm = REC(Fm, dx, np_tot - i - 2); 
    #if SHOCK_FLATTENING == MULTID
    if (sweep->flag[i] & FLAG_MINMOD){
      fp = WENO3_Reconstruct(Fp, dx, i);
      fm = WENO3_Reconstruct(Fm, dx, np_tot - i - 2); 
    }
    #endif
    Bm[i]   = (fp - fm)/glm_ch;
    psim[i] =  fp + fm;
    #if GLM_COMPUTE_DIVB == YES
    sweep->bn[i] = Bm[i];
    #endif
  }
  #if GLM_COMPUTE_DIVB == YES
  if (g_intStage == 1) GLM_ComputeDivB(sweep, grid);
  #endif

  Kmax = KPSI_GLMM;
#else
  Kmax = NFLX;
#endif 

/* ----------------------------------------------------
   4. Compute average state
   ---------------------------------------------------- */

  for (i = beg; i <= end; i++){
    for (nv = NFLX; nv--;    ){
      Uave[i][nv] = 0.5*(u[i][nv] + u[i+1][nv]);
      Vave[i][nv] = 0.5*(v[i][nv] + v[i+1][nv]);
    }
    #ifdef GLM_MHD
    Uave[i][BXn]     = Vave[i][BXn]     = Bm[i];
    Uave[i][PSI_GLM] = Vave[i][PSI_GLM] = psim[i];
    #endif 
  }

/* ---------------------------------------
   5. Compute sound speed
   --------------------------------------- */

  SoundSpeed2 (stateLR, beg, end, FACE_CENTER, grid);

/* ----------------------------------------------------
   6. Compute eigenvectors and eigenvalues
   ---------------------------------------------------- */

  for (i = beg; i <= end; i++){
    L = stateLR->Lp[i];
    R = stateLR->Rp[i];

    ConsEigenvectors(Uave[i], Vave[i], stateLR->a2[i], L, R, lp[i]);

    fs = 0.0;
    for (k = Kmax; k--;     ){
      fs = MAX(fs, fabs(sweep->lmax[k]));
    }  
    cmax[0] = MAX(cmax[0], fs);
    cmax[i] = fs;
  }

/* ---------------------------------------------------------
   6. Main spatial loop.
      At each interface construct, for each characteristic 
      field k, the positive and negative part of the flux:
   
        F_{+,j} = 1/2 L.[F(j)  + a*u(j)]
        F_{-,j} = 1/2 L.[F(j') - a*u(j')]

      a is the global maximum eigenvalue for the k-th 
      characteristic.
      Reconstruct F_{+,j} and F_{-,j}.
   --------------------------------------------------------- */

  for (i = beg; i <= end; i++){
    L = stateLR->Lp[i];
    R = stateLR->Rp[i];
    for (k = 0; k < Kmax; k++){
      fs = sweep->lmax[k];
      for (j = i-S; j <= i+S; j++) {
        Fp[j] = Fm[j] = 0.0;
        for (nv = NFLX; nv--;   ){ /* -- LF flux splitting -- */
          Fp[j] += 0.5*L[k][nv]*(F[j][nv]       + fs*u[j][nv]);
          Fm[j] += 0.5*L[k][nv]*(F[2*i-j+1][nv] - fs*u[2*i-j+1][nv]);
        }
      }
      #if SHOCK_FLATTENING == MULTID
      if ( (sweep->flag[i] & FLAG_MINMOD) || (sweep->flag[i+1] & FLAG_MINMOD)){
        fs = LIN_Reconstruct(Fp, dx, i) + LIN_Reconstruct(Fm, dx, i); 
      }else
      #endif
      flx[k] = REC(Fp, dx, i) + REC(Fm, dx,i);
    }

    for (nv = 0; nv < NFLX; nv++){
      fs = 0.0;
      for (k=0; k < Kmax; k++) fs += flx[k]*R[nv][k];
      flux[i][nv] = fs;
    }

    #ifdef GLM_MHD
    flux[i][BXn]     = psim[i];
    flux[i][PSI_GLM] = glm_ch*glm_ch*Bm[i];
    #endif

  /* -------------------------------------------
      Passive scalars are treated separately in 
      a way similar to the finite volume case:
      
        F(T) = F(\rho)*TL  if F(rho) >= 0
        F(T) = F(\rho)*TR  if F(rho) <  0
 
       (inside adv_flux).
       This seems to preserve better the 
       original tracer bounds (e.g., 0<T<1).
      ------------------------------------------ */

    #if NSCL > 0 
    if (flux[i][RHO] >= 0.0){
      NSCL_LOOP(k){
        for (j = i-S; j <= i+S; j++) Fp[j] = v[j][k];
        stateL->v[i][k] = REC(Fp, dx, i);
        stateR->v[i][k] = v[i+1][k];
      } 
    }else{
      NSCL_LOOP(k){
        for (j = i-S; j <= i+S; j++) Fm[j] = v[2*i-j+1][k];
        stateL->v[i][k] = v[i][k];
        stateR->v[i][k] = REC(Fm, dx, i);
      } 
    }
    #endif
  }  /* End loop on i=beg,end */
}

/* ********************************************************************* */
void FD_GetMaxEigenvalues (const Data *d, Sweep *sweep, Grid *grid)
/*
 * PURPOSE:
 *
 *   Compute, for each characteristic wave, its maximum over
 *   the entire computational domain
 *
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *x1, *x2, *x3;
  static double **v, **lambda, *a2;
  State state;
  
  if (v == NULL){
    v      = ARRAY_2D(NMAX_POINT, NFLX, double);
    lambda = ARRAY_2D(NMAX_POINT, NFLX, double);
    a2     = ARRAY_1D(NMAX_POINT, double);
  }

  NVAR_LOOP(nv) sweep->lmax[nv] = 0.0;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  state.v = v;
  state.lambda = lambda;
  state.a2 = a2;
  KDOM_LOOP(k) JDOM_LOOP(j){
    IDOM_LOOP(i) NVAR_LOOP(nv) v[i][nv] = d->Vc[nv][k][j][i];
    SoundSpeed2  (&state, IBEG, IEND, CELL_CENTER, grid);
    Eigenvalues  (v, a2, lambda, IBEG, IEND);
    IDOM_LOOP(i){
      for (nv = NFLX; nv--;  ) 
        sweep->lmax[nv] = MAX(fabs(lambda[i][nv]), sweep->lmax[nv]);
    }
  }

  #ifdef PARALLEL
   MPI_Allreduce (sweep->lmax, lambda[0], NFLX, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   for (nv = 0; nv < NFLX; nv++) sweep->lmax[nv] = lambda[0][nv];
  #endif
}
