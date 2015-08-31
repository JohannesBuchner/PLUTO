/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK split/unsplit integrations and finite difference
  methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Dec 18, 2014
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
int AdvanceStep (const Data *d, Riemann_Solver *Riemann, 
                 Time_Step *Dts, Grid *grid)
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
  RBox *box = GetRBox (DOM, CENTER);
  
/* ----------------------------------------------------
   0. Allocate memory 
   ---------------------------------------------------- */

  if (U0 == NULL){
    U0 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    #ifdef STAGGERED_MHD
     Bs0 = ARRAY_4D(DIMENSIONS, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

#ifdef FARGO
  FARGO_SubtractVelocity (d,grid);
#endif

/* ---------------------------------------------------------------
   1. Predictor step (EULER, RK2, RK3)

      After baoundaries have been set we flag zones lying in a
      shock. This is useful for shock flattening or entropy/energy
      selective update.

      Note: when using FARGO, boundary condition must be set on
      the *total* velocity while the update step is performed on
      the *residual* velocity.
      For this reason we add/subtract the mean velocity field
      immediately before/after setting the boundary conditions.
   --------------------------------------------------------------- */

  g_intStage = 1;  
  Boundary (d, ALL_DIR, grid);
#if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (d, grid);
#endif

/* -- Convert primitive to conservative, save initial stage  -- */

  PrimToCons3D(d->Vc, d->Uc, box);
  KDOM_LOOP(k) JDOM_LOOP(j){
    memcpy ((void *)U0[k][j][IBEG], d->Uc[k][j][IBEG], NX1*NVAR*sizeof(double));
  }
#ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) Bs0[nv][k][j][i] = d->Vs[nv][k][j][i];
#endif

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
#ifdef STAGGERED_MHD
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
#endif
  ConsToPrim3D (d->Uc, d->Vc, d->flag, box);

/* ----------------------------------------------------
   2. Corrector step (RK2, RK3)
   ---------------------------------------------------- */

#if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)

   g_intStage = 2;
   Boundary (d, ALL_DIR, grid);

/* -- need an extra conversion if INTERNAL_BOUNDARY is enabled 
      [note: done only with dimensional splitting for backward compat.] -- */

  #if (INTERNAL_BOUNDARY == YES) && (DIMENSIONAL_SPLITTING == YES)
  PrimToCons3D (d->Vc, d->Uc, box);
  #endif   

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
  DOM_LOOP(k, j, i) VAR_LOOP(nv){
    d->Uc[k][j][i][nv] = w0*U0[k][j][i][nv] + wc*d->Uc[k][j][i][nv];
  }
  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) {
    d->Vs[nv][k][j][i] = w0*Bs0[nv][k][j][i] + wc*d->Vs[nv][k][j][i];
  }
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
  #endif

  #if (defined FARGO) && (TIME_STEPPING == RK2)
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif 
  ConsToPrim3D (d->Uc, d->Vc, d->flag, box);

#endif  /* TIME_STEPPING == RK2/RK3 */

/* ----------------------------------------------------
   3. Last corrector step (RK3 only) 
   ---------------------------------------------------- */

#if TIME_STEPPING == RK3
  g_intStage = 3;
  Boundary (d, ALL_DIR, grid);

/* -- need an extra conversion if INTERNAL_BOUNDARY is enabled -- */

  #if (INTERNAL_BOUNDARY == YES) && (DIMENSIONAL_SPLITTING == YES)
  PrimToCons3D (d->Vc, d->Uc, box);
  #endif

  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
  DOM_LOOP(k,j,i) VAR_LOOP(nv){
    d->Uc[k][j][i][nv] = one_third*(U0[k][j][i][nv] + 2.0*d->Uc[k][j][i][nv]);
  }
  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i){
    d->Vs[nv][k][j][i] = (Bs0[nv][k][j][i] + 2.0*d->Vs[nv][k][j][i])/3.0;
  }
  CT_AverageMagneticField (d->Vs, d->Uc, grid);
  #endif

  #ifdef FARGO
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif
  ConsToPrim3D (d->Uc, d->Vc, d->flag, box);
#endif /* TIME_STEPPING == RK3 */


#ifdef FARGO
  FARGO_AddVelocity (d,grid);
#endif

  return 0; /* -- step has been achieved, return success -- */
}