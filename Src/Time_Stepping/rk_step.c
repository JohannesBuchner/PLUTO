/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK split/unsplit integrations and finite difference
  methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Nov 21, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* Weight factor for 2nd stage of RK integrators */

#if TIME_STEPPING == RK2
 #define w0  0.5
 #define wc  0.5
#elif TIME_STEPPING == RK3
 #define w0 0.75
 #define wc 0.25
#endif

/* ********************************************************************* */
int AdvanceStep (Data *d, Riemann_Solver *Riemann, 
                 timeStep *Dts, Grid *grid)
/*!
 * Advance the equations by a single time step using unsplit 
 * integrators based on the method of lines.
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in]    Riemann  pointer to a Riemann solver function
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *    
 *********************************************************************** */
{
  int  i, j, k, nv;
  static double  one_third = 1.0/3.0;
  static Data_Arr U0, Bs0;
  RBox   box;

static Data_Arr Uhalf;
  
  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

#if (defined PARTICLES) && DIMENSIONAL_SPLITTING == YES
  print ("! AdvanceStep(): particles require DIMENSIONAL_SPLITTING == NO\n");
  QUIT_PLUTO(1);
#endif

/* --------------------------------------------------------
   0. Allocate memory 
   -------------------------------------------------------- */

  if (U0 == NULL){
    U0 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
Uhalf = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

    #ifdef STAGGERED_MHD
     Bs0 = ARRAY_4D(DIMENSIONS, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

/* --------------------------------------------------------
   1. Predictor step (EULER, RK2, RK3, MH (split), 
                      ChTr (split))

      After baoundaries have been set we flag zones lying 
      in a shock. 
      This is useful for shock flattening or 
      entropy/energy selective update.

      Note: when using FARGO, boundary condition must be 
      set on the *total* velocity while the update step is 
      performed on the *residual* velocity.
      The addition and subtraction operations are 
      automatically performed in Boundary() function.
   -------------------------------------------------------- */

/* -- 1a. Set boundary conditions -- */

  g_intStage = 1;  
  Boundary (d, ALL_DIR, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (d, grid);
  #endif

/* -- 1b. Convert primitive to conservative, save initial stage  -- */

  PrimToCons3D(d->Vc, d->Uc, &box);
  KDOM_LOOP(k) JDOM_LOOP(j){
    memcpy ((void *)U0[k][j][IBEG], d->Uc[k][j][IBEG], NX1*NVAR*sizeof(double));
  }

  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) Bs0[nv][k][j][i] = d->Vs[nv][k][j][i];
  #endif

/* -- 1c. Compute Particles feedback / Advance with pred. step -- */

#ifdef PARTICLES
  #if PARTICLES_TYPE == COSMIC_RAYS
  Particles_CR_ComputeForce (d->Vc, d, grid);
  #elif PARTICLES_TYPE == DUST
  Particles_Dust_ComputeForce (d->Vc, d, grid);
  #elif PARTICLES_TYPE == LAGRANGIAN
  Particles_LP_Predictor(d, Dts, g_dt, grid);
  #endif
#endif

/* -- 1d. Advance conservative variables array -- */

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
#ifdef STAGGERED_MHD
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
#endif

/* -- 1e. Advance particles by a full step -- */

#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) || (PARTICLES_TYPE == DUST)
  TOT_LOOP(k,j,i) NVAR_LOOP(nv) {
    Uhalf[k][j][i][nv] = 0.5*(U0[k][j][i][nv] + d->Uc[k][j][i][nv]);
  }
  ConsToPrim3D (Uhalf, d->Vc, d->flag,&box);
  #endif
  #if PARTICLES_TYPE == COSMIC_RAYS
  Particles_CR_Update(d, Dts, g_dt, grid);
  #elif PARTICLES_TYPE == DUST
  Particles_Dust_Update(d, Dts, g_dt, grid);
  #endif
#endif
 
/* -- 1f. Convert to primitive vars -- */

  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);

/* --------------------------------------------------------
   2. Corrector step (RK2, RK3)
   -------------------------------------------------------- */

#if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)

/* -- 2a. Set boundary conditions -- */

  g_intStage = 2;
  Boundary (d, ALL_DIR, grid);

/* -- 2b. Advance Particles -- */

  #if (defined PARTICLES) && (PARTICLES_TYPE == LAGRANGIAN)
  Particles_LP_Corrector(d, Dts, g_dt, grid);
  #endif

/* -- 2c. need an extra conversion if INTERNAL_BOUNDARY is enabled 
         [note: done only with dimensional splitting for backward compat.] -- */

  #if (INTERNAL_BOUNDARY == YES) && (DIMENSIONAL_SPLITTING == YES)
  PrimToCons3D (d->Vc, d->Uc, &box);
  #endif   

/* -- 2d. Advance solution array -- */

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
  DOM_LOOP(k, j, i) NVAR_LOOP(nv){
    d->Uc[k][j][i][nv] = w0*U0[k][j][i][nv] + wc*d->Uc[k][j][i][nv];
  }
  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) {
    d->Vs[nv][k][j][i] = w0*Bs0[nv][k][j][i] + wc*d->Vs[nv][k][j][i];
  }
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
  #endif

/* -- 2e. Apply FARGO orbital shift -- */

  #if (defined FARGO) && (TIME_STEPPING == RK2)
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif

/* -- 2f. Convert to Primitive -- */

  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);

/* -- 2g. Inject particles or update spectra-- */
  
  #ifdef PARTICLES
    #if (PARTICLES_TYPE == COSMIC_RAYS)
    Particles_Inject(d,grid);
    #elif (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == YES)
    Particles_LP_UpdateSpectra (d, g_dt, grid);
    #endif
  #endif

#endif  /* TIME_STEPPING == RK2/RK3 */

/* --------------------------------------------------------
   3. Last corrector step (RK3 only) 
   -------------------------------------------------------- */

#if TIME_STEPPING == RK3
  #ifdef PARTICLES 
  #if PARTICLES_TYPE == COSMIC_RAYS || PARTICLES_TYPE == DUST
  print ("! AdvanceStep(): RK3 algorithm not permitted\n");
  QUIT_PLUTO(1);
  #endif
  #endif

/* -- 3a. Set Boundary conditions -- */

  g_intStage = 3;
  Boundary (d, ALL_DIR, grid);

/* -- 3b. need an extra conversion if INTERNAL_BOUNDARY is enabled 
         [note: done only with dimensional splitting for backward compat.] -- */

  #if (INTERNAL_BOUNDARY == YES) && (DIMENSIONAL_SPLITTING == YES)
  PrimToCons3D (d->Vc, d->Uc, &box);
  #endif

/* -- 3c. Update solution array -- */

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
  DOM_LOOP(k,j,i) NVAR_LOOP(nv){
    d->Uc[k][j][i][nv] = one_third*(U0[k][j][i][nv] + 2.0*d->Uc[k][j][i][nv]);
  }
  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i){
    d->Vs[nv][k][j][i] = (Bs0[nv][k][j][i] + 2.0*d->Vs[nv][k][j][i])/3.0;
  }
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
  #endif

/* -- 3d. Apply FARGO orbital shift -- */

  #ifdef FARGO
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
#endif /* TIME_STEPPING == RK3 */

/* --------------------------------------------------------
   4. Add velocity if FARGO is used
   -------------------------------------------------------- */

#ifdef FARGO
  FARGO_AddVelocity (d,grid);
#endif

  return 0; /* -- step has been achieved, return success -- */
}

