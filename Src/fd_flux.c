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
  \date   April 02, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void FD_Flux (const State_1D *state, int beg, int end, 
               double *cmax, Grid *grid)
/*!
 * Compute interface flux by suitable high-order finite
 * difference non-oscillatory interpolants.
 *
 * \param [in]  state   pointer to State_1D structure
 * \param [in]    beg   initial index of computation
 * \param [in]    end   final   index of computation
 * \param [out]  cmax   array of maximum characteristic speeds
 * \param [in]    grid  pointer to an array of Grid structures
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
  int    nv, i, j, k, np_tot;
  int    Kmax, S;
  double **v, **flux, *press;
  double fs, fp, fm, flx[NFLX], dx;

  static double *Fp, *Fm, **F, *a2, **u;
  static double **Uave, **Vave, **lp;
  double **L, **R;
  static double *psim, *Bm; 
  double (*REC)(double *, double, int);

/* ------------------------------------------------------------------
           Check compatibility with other modules 
   ------------------------------------------------------------------ */

  #if TIME_STEPPING != RK3 && TIME_STEPPING != SSP_RK4
   print1 ("! Finite Difference schemes work with RK3 only \n");
   QUIT_PLUTO(1);
  #endif
  #if GEOMETRY != CARTESIAN
   print1 ("! Finite Difference schemes work in Cartesian coordinates only\n");
   QUIT_PLUTO(1);
  #endif
  #if PHYSICS == RMHD || PHYSICS == RHD
   print1 ("! Finite difference schemes work only for HD od MHD modules\n");
   QUIT_PLUTO(1);
  #endif   
  #ifdef STAGGERED_MHD
   print1 ("! Finite difference schemes work only with cell-centered schemes\n");
   QUIT_PLUTO(1);
  #endif   
  #if PHYSICS == MHD && BACKGROUND_FIELD == YES
   print1 ("! Background field splitting not supported with FD schemes\n");
   QUIT_PLUTO(1);
  #endif  

/* --------------------------------------------------------------
     - Define pointer to reconstruction function
     - Determine the stencil for interpolation
       For a given S, the stencil is: i-S <= i <= i+S+1
   -------------------------------------------------------------- */

  dx = grid[g_dir].dx[beg];
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
  
  if (F == NULL){
    Fp   = ARRAY_1D(NMAX_POINT, double);
    Fm   = ARRAY_1D(NMAX_POINT, double);
    F    = ARRAY_2D(NMAX_POINT, NVAR, double);
    lp   = ARRAY_2D(NMAX_POINT, NVAR, double);
    psim = ARRAY_1D(NMAX_POINT, double);
    Bm   = ARRAY_1D(NMAX_POINT, double);
    a2   = ARRAY_1D(NMAX_POINT, double);
    Uave = ARRAY_2D(NMAX_POINT, NVAR, double);
    Vave = ARRAY_2D(NMAX_POINT, NVAR, double);
    u    = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  np_tot = grid[g_dir].np_tot;
  flux   = state->flux;
  press  = state->press;
  v      = state->v;

/* ------------------------------------------------------
    Compute cell-centered fluxes. Pressure is added to 
    the normal component of momentum.
    To save computational time, we also solve the linear 
    GLM-MHD subsystem here and use the solution in the 
    nonlinear fluxes later.
   ------------------------------------------------------ */

  PrimToCons (v, u, 0, np_tot - 1);
  SoundSpeed2 (v, a2, NULL, 0, np_tot-1, CELL_CENTER, grid);
  #if PHYSICS == MHD
   Flux (u, v, a2, NULL, F, press, 0, np_tot-1);
  #else
   Flux (u, v, a2, F, press, 0, np_tot-1);
  #endif

  for (i = 0; i < np_tot; i++){
    F[i][MXn] += press[i];
    press[i]  = 0.0;
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
      if (state->flag[i] & FLAG_MINMOD){
        fp = WENO3_Reconstruct(Fp, dx, i);
        fm = WENO3_Reconstruct(Fm, dx, np_tot - i - 2); 
      }
     #endif
     Bm[i]   = (fp - fm)/glm_ch;
     psim[i] =  fp + fm;
     #if COMPUTE_DIVB == YES
      state->bn[i] = Bm[i];
     #endif
   }
   #if COMPUTE_DIVB == YES
    if (g_intStage == 1) GLM_ComputeDivB(state, grid);
   #endif

   Kmax = KPSI_GLMM;
  #else
   Kmax = NFLX;
  #endif

/* ----------------------------------------------------
              Compute average state
   ---------------------------------------------------- */

  for (i = beg; i <= end; i++){
    for (nv = NFLX; nv--;    ){
      Uave[i][nv] = 0.5*(u[i][nv] + u[i + 1][nv]);
      Vave[i][nv] = 0.5*(v[i][nv] + v[i + 1][nv]);
    }

    #ifdef GLM_MHD
     Uave[i][BXn]      = Vave[i][BXn]      = Bm[i];
     Uave[i][PSI_GLM] = Vave[i][PSI_GLM] = psim[i];
    #endif 
  }

/* ---------------------------------------
     compute sound speed
   --------------------------------------- */

  SoundSpeed2 (Vave, a2, NULL, beg, end, FACE_CENTER, grid);

/* ----------------------------------------------------
         Compute eigenvectors and eigenvalues
   ---------------------------------------------------- */

  for (i = beg; i <= end; i++){
    L = state->Lp[i];
    R = state->Rp[i];

    ConsEigenvectors(Uave[i], Vave[i], a2[i], L, R, lp[i]);

    fs = 0.0;
    for (k = Kmax; k--;     ){
      fs = MAX(fs, fabs(state->lmax[k]));
    }  
    cmax[0] = MAX(cmax[0], fs);
    cmax[i] = fs;
  }

/* ---------------------------------------------------------
    Main spatial loop.
    At each interface construct, for each characteristic 
    field k, the positive and negative part of the flux:
   
      F_{+,j} = 1/2 L.[F(j)  + a*u(j)]
      F_{-,j} = 1/2 L.[F(j') - a*u(j')]

    a is the global maximum eigenvalue for the k-th 
    characteristic.
    Reconstruct F_{+,j} and F_{-,j}.
   --------------------------------------------------------- */

  for (i = beg; i <= end; i++){
    L = state->Lp[i];
    R = state->Rp[i];

    for (k = 0; k < Kmax; k++){
      fs = state->lmax[k];
      for (j = i-S; j <= i+S; j++) {
        Fp[j] = Fm[j] = 0.0;
        for (nv = NFLX; nv--;   ){ /* -- LF flux splitting -- */
          Fp[j] += 0.5*L[k][nv]*(F[j][nv]       + fs*u[j][nv]);
          Fm[j] += 0.5*L[k][nv]*(F[2*i-j+1][nv] - fs*u[2*i-j+1][nv]);
        }
      }
      #if SHOCK_FLATTENING == MULTID
       if ( (state->flag[i] & FLAG_MINMOD) || (state->flag[i+1] & FLAG_MINMOD)){
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
     flux[i][BXn]      = psim[i];
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
        state->vL[i][k] = REC(Fp, dx, i);
        state->vR[i][k] = v[i+1][k];
      } 
    }else{
      NSCL_LOOP(k){
        for (j = i-S; j <= i+S; j++) Fm[j] = v[2*i-j+1][k];
        state->vL[i][k] = v[i][k];
        state->vR[i][k] = REC(Fm, dx, i);
      } 
    }
#endif
  }
}

/* ********************************************************************* */
void FD_GetMaxEigenvalues (const Data *d, State_1D *state, Grid *grid)
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
  
  if (v == NULL){
    v      = ARRAY_2D(NMAX_POINT, NFLX, double);
    lambda = ARRAY_2D(NMAX_POINT, NFLX, double);
    a2     = ARRAY_1D(NMAX_POINT, double);
  }

  for (nv = NVAR; nv--;  ) state->lmax[nv] = 0.0;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  KDOM_LOOP(k) JDOM_LOOP(j){
    IDOM_LOOP(i) for (nv = NVAR; nv--;  ) v[i][nv] = d->Vc[nv][k][j][i];
    SoundSpeed2  (v, a2, NULL, IBEG, IEND, CELL_CENTER, grid);
    Eigenvalues  (v, a2, lambda, IBEG, IEND);
    IDOM_LOOP(i){
      for (nv = NFLX; nv--;  ) 
        state->lmax[nv] = MAX(fabs(lambda[i][nv]), state->lmax[nv]);
    }
  }

  #ifdef PARALLEL
   MPI_Allreduce (state->lmax, lambda[0], NFLX, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   for (nv = 0; nv < NFLX; nv++) state->lmax[nv] = lambda[0][nv];
  #endif
}
