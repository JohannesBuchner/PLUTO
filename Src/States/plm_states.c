/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Piecewise linear reconstruction.

  Compute interface states using piecewise linear reconstruction 
  inside each zone.
  Reconstruction is performed in primitive variable when 
   <tt> CHAR_LIMITING == NO </tt>) or characteristic variables 
  when <tt> CHAR_LIMITING == YES </tt>.
  
  The convention used throughout is that \b vL and \b vR are left and 
  right states with respect to the \e interface, while \b vm and \b vp 
  refer to the cell \e center:
  \verbatim
                      vL(i)-> <-vR(i)
       |----------*----------|----------*----------|
        <-vm(i)  (i)  vp(i)->         (i+1)
  \endverbatim    
 
  The default setting (LIMITER == DEFAULT) applies a different limiter
  to each variable: 
  - MC for density
  - VanLeer for velcoity and magnetic field
  - MinMod for pressure.
  
  Otherwise the same limiter can be imposed to all variables from 
  definitions.h.

  The approach followed here is taken from Mignone (JCP, 2014) "High-order 
  conservative reconstruction schemes for finite volume methods in cylindrical
  and spherical geometries" where interface states are constructed as
  \f[
    V_i^{\pm} = V_i \;\pm\; \overline{\Delta V}_i \; d^\pm_i\,,
    \qquad
   \overline{\Delta V}_i = {\rm Lim}\left(\Delta V_i^+,\,\Delta V_i^-,\,
                                          c^+_i, c^-_i\right)\,,
   \qquad
     \Delta V^\pm_i = \pm\left(V_{i\pm1}-V_i\right)w^\pm_i
  \f]
  where
  \f[
     w^\pm_i = \left|\frac{\Delta\xi_i}{\bar{\xi}_{i\pm1} - \bar{\xi}_i}\right| \,,
    \qquad
     d^\pm_i = \left|\frac{\xi_{i\pm\HALF} - \bar{\xi}_i}{\Delta\xi_i}\right|\,,
    \qquad
     c^\pm_i = \left|\frac{\bar{\xi}_{i\pm1} - \bar{\xi}_i}
                          {\xi_{i\pm\HALF}- \bar{\xi}_i}\right|
  \f]
  are interpolation coefficients returned by the PLM_Coeffs structure
  (see the ::PLM_CoefficientsGet() function).
  The slope limiter \f$\overline{\Delta V}_i=\f$ <tt> dv\_lim[i][nv] </tt> 
  is a function of the forward  and backward (undivided) derivatives 
  and geometrical coefficients in case of non-uniform or non-Cartesian grids.
  Different limiters are implemented through small macros defined in 
  plm_coeffs.h.
  Note that \f$ \bar{\xi} \f$ defines the centroid of volume which differs 
  from the cell center in cylindrical and spherical geometires.

  A stencil of 3 zones is required for all limiters except for
  the FOURTH_ORDER_LIM which requires  5 zones. 

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 16, 2015

  \b References
     - "High-order conservative reconstruction schemes for finite
        volume methods in cylindrical and spherical coordinates"
        A. Mignone, JCP (2014), 270, 784.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void MonotonicityTest(double **, double **, double **, int, int);

#if CHAR_LIMITING == NO
static void FourthOrderLinear(const State_1D *, int, int, Grid *);

/* ********************************************************************* */
void States (const State_1D *state, int beg, int end, Grid *grid)
/*! 
 * Compute states using piecewise linear interpolation.
 *
 * \param [in] state pointer to a State_1D structure
 * \param [in] beg   starting point where vp and vm must be computed
 * \param [in] end   final    point where vp and vm must be computed
 * \param [in] grid  pointer to array of Grid structures
 *
 * \return This function has no return value.
 *
 ************************************************************************ */
{
  int    nv, i;
  double **v, **vp, **vm;
  double dv_lim[NVAR], dvp[NVAR], dvm[NVAR];
  double cp, cm, wp, wm, dp, dm;
  PLM_Coeffs plm_coeffs;
  static double **dv;

  #if LIMITER == FOURTH_ORDER_LIM
   FourthOrderLinear(state, beg, end, grid);
   return;
  #endif
  
/* -----------------------------------------------------------
   0. Memory allocation and pointer shortcuts, geometrical
      coefficients and conversion to 4vel (if required)
   ----------------------------------------------------------- */

  if (dv == NULL) {
    dv = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  v  = state->v;
  vp = state->vp;
  vm = state->vm;

#if UNIFORM_CARTESIAN_GRID == NO
  PLM_CoefficientsGet (&plm_coeffs, g_dir);
#endif

#if RECONSTRUCT_4VEL
  ConvertTo4vel (state->v, beg-1, end+1);
#endif

/* -------------------------------------------
        compute undivided differences
   ------------------------------------------- */

  for (i = beg-1; i <= end; i++){
    NVAR_LOOP(nv)  dv[i][nv] = v[i+1][nv] - v[i][nv];
  }

/* -------------------------------------------
    2. Main spatial loop
   ------------------------------------------- */

  for (i = beg; i <= end; i++){

  /* ---------------------------------------------------------
     2a. compute forward (dvp) and backward (dvm) derivatives
     --------------------------------------------------------- */

    #if UNIFORM_CARTESIAN_GRID == YES
     cp = cm = 2.0;
     wp = wm = 1.0;
     dp = dm = 0.5;
     NVAR_LOOP(nv) {
       dvp[nv] = dv[i][nv];
       dvm[nv] = dv[i-1][nv];
     }
    #else
     cp = plm_coeffs.cp[i]; cm = plm_coeffs.cm[i];
     wp = plm_coeffs.wp[i]; wm = plm_coeffs.wm[i];
     dp = plm_coeffs.dp[i]; dm = plm_coeffs.dm[i];
     NVAR_LOOP(nv) {
       dvp[nv] = dv[i][nv]*wp;
       dvm[nv] = dv[i-1][nv]*wm;
     }
    #endif

  /* ---------------------------------------------------------------
     2b. if shock-flattening is enabled, revert to minmod limiter
     --------------------------------------------------------------- */
     
#if SHOCK_FLATTENING == MULTID
    if (state->flag[i] & FLAG_FLAT) {
      NVAR_LOOP(nv)vp[i][nv] = vm[i][nv] = v[i][nv];
      continue;
    }else if (state->flag[i] & FLAG_MINMOD) {
      NVAR_LOOP(nv){
        SET_MM_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
        vp[i][nv] = v[i][nv] + dv_lim[nv]*dp;
        vm[i][nv] = v[i][nv] - dv_lim[nv]*dm;
      }
  #if (PHYSICS == RHD || PHYSICS == RMHD)
      VelocityLimiter (v[i], vp[i], vm[i]);
  #endif
      continue;
    }
#endif    

  /* ---------------------------------------------------------
     2c. DEFAULT limiting: combination of different limiters
        (this has some hystorical reasons)
     --------------------------------------------------------- */

    #if LIMITER == DEFAULT
     SET_MC_LIMITER(dv_lim[RHO], dvp[RHO], dvm[RHO], cp, cm);
     #if PHYSICS != ADVECTION
      EXPAND(SET_VL_LIMITER(dv_lim[VX1], dvp[VX1], dvm[VX1], cp, cm);  ,
             SET_VL_LIMITER(dv_lim[VX2], dvp[VX2], dvm[VX2], cp, cm);  ,
             SET_VL_LIMITER(dv_lim[VX3], dvp[VX3], dvm[VX3], cp, cm);)
     #endif

     #if PHYSICS == MHD || PHYSICS == RMHD
      EXPAND(SET_VL_LIMITER(dv_lim[BX1], dvp[BX1], dvm[BX1], cp, cm);  ,
             SET_VL_LIMITER(dv_lim[BX2], dvp[BX2], dvm[BX2], cp, cm);  ,
             SET_VL_LIMITER(dv_lim[BX3], dvp[BX3], dvm[BX3], cp, cm);)
      #ifdef GLM_MHD
       SET_MC_LIMITER(dv_lim[PSI_GLM], dvp[PSI_GLM], dvm[PSI_GLM], cp, cm);
      #endif
     #endif

     #if HAVE_ENERGY
      SET_MM_LIMITER(dv_lim[PRS], dvp[PRS], dvm[PRS], cp, cm);
     #endif
 
     #if NFLX != NVAR /* -- scalars: MC lim  -- */
      for (nv = NFLX; nv < NVAR; nv++){
        SET_MC_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
      }
     #endif
    #endif /* LIMITER == DEFAULT */

  /* ----------------------------------------
     2d. construct (+) and (-) states
     ---------------------------------------- */

    for (nv = 0; nv < NVAR; nv++){
      #if LIMITER != DEFAULT  /* -- same limiter for all variables -- */
       SET_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
      #endif

      vp[i][nv] = v[i][nv] + dv_lim[nv]*dp;
      vm[i][nv] = v[i][nv] - dv_lim[nv]*dm;
    }

  /* --------------------------------------
      2e.  Relativistic Limiter
     -------------------------------------- */

    #if (PHYSICS == RHD || PHYSICS == RMHD) 
     VelocityLimiter (v[i], vp[i], vm[i]);
    #endif
  } /* -- end loop on zones -- */

/* ----------------------------------------------------
   3a. Check monotonicity
   ---------------------------------------------------- */ 

#if CHECK_MONOTONICITY == YES
  MonotonicityTest(state->v, state->vp, state->vm, beg, end);
#endif

/* -------------------------------------------
   3b. Shock flattening 
   -------------------------------------------  */

#if SHOCK_FLATTENING == ONED
  Flatten (state, beg, end, grid);
#endif

/* -------------------------------------------
   4.  Assign face-centered magnetic field
   -------------------------------------------  */

#ifdef STAGGERED_MHD
  for (i = beg - 1; i <= end; i++) {
    state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
  }
#endif

/* --------------------------------------------------------
   5. Evolve L/R states and center value by dt/2
   -------------------------------------------------------- */

#if TIME_STEPPING == CHARACTERISTIC_TRACING
  CharTracingStep(state, beg, end, grid);
#elif TIME_STEPPING == HANCOCK
  HancockStep(state, beg, end, grid);
  /* The conservative hancock scheme is used only for RMHD with linear
     interpolation in primitive variables.
     In this case, the conversion to 3-vel is already carried out in
     the predictor step and we need to skip it here. */
  #if PRIMITIVE_HANCOCK == NO
    PrimToCons (state->vp, state->up, beg, end);
    PrimToCons (state->vm, state->um, beg, end);
    return;
  #endif
#endif

/* --------------------------------------------------------
   6. Convert back to 3-velocity
   -------------------------------------------------------- */

#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->v, beg-1, end+1);
  ConvertTo3vel (state->vp, beg, end);
  ConvertTo3vel (state->vm, beg, end);  
#endif

/* ----------------------------------------------
   7. Obtain L/R states in conservative variables
   ---------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
/* ********************************************************************** */
void FourthOrderLinear(const State_1D *state, int beg, int end, Grid *grid)
/*
 * PURPOSE
 *
 *   Compute interface states using Colella's fourth-order slope limiter
 * 
 *  Ref:  Miller, G.H and P. COlella, 
 *        "A high-order Eulerian Godunov Method for 
 *         Elastic-Plastic Flow in Solids", JCP 167,131-176 (2001)
 *    
 *                             +
 *
 *        Saltzman, J, " An Unsplit 3D Upwind Method for 
 *                       Hyperbolic Conservation Laws", 
 *                       JCP 115, 153-168 (1994)
 *
 *********************************************************************** */
{
  int    i, nv;
  static double **s;
  static double **dv, **dvf, **dvc, **dvlim; 
  double scrh, dvp, dvm, dvl;
  double **v, **vp, **vm;

  v  = state->v;
  vp = state->vp;
  vm = state->vm;

  if (s == NULL){
    s     = ARRAY_2D(NMAX_POINT, NVAR, double);
    dv    = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvf   = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvc   = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvlim = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  #if TIME_STEPPING == HANCOCK && PHYSICS != RMHD
   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
  #endif
/*
  #if GEOMETRY != CARTESIAN
   print1 ("! FourthOrderLinear: only Cartesian geometry supported\n");
   QUIT_PLUTO(1);  
  #endif
*/
/* -----------------------------------------------------------
               compute undivided differences
   ----------------------------------------------------------- */

  for (i = beg-2; i <= end+1; i++){
    for (nv = 0; nv < NVAR; nv++) dv[i][nv] = v[i+1][nv] - v[i][nv];
  }

  for (i = beg - 1; i <= end + 1; i++){
  for (nv = 0; nv < NVAR; nv++){
    dvp = dv[i][nv]; dvm = dv[i-1][nv];
    dvc[i][nv] = 0.5*(dvp + dvm);
      s[i][nv] =  (dvp > 0.0 ? 0.5:-0.5) 
                + (dvm > 0.0 ? 0.5:-0.5);
    dvlim[i][nv] = 2.0*MIN(fabs(dvp), fabs(dvm));
    dvf[i][nv]   = MIN(fabs(dvc[i][nv]), dvlim[i][nv])*s[i][nv];
  }}

  for (i = beg; i <= end; i++){
    for (nv = 0; nv < NVAR; nv++){
      if (dv[i][nv]*dv[i-1][nv] > 0.0) {
        scrh = 4.0/3.0*dvc[i][nv] - (dvf[i+1][nv] + dvf[i-1][nv])/6.0; 
        dvlim[i][nv]  = MIN(fabs(scrh), dvlim[i][nv])*s[i][nv];
      }else{
        dvlim[i][nv] = 0.0;
      }
    }

    for (nv = NVAR; nv--;  ) {
      vp[i][nv] = v[i][nv] + 0.5*dvlim[i][nv];
      vm[i][nv] = v[i][nv] - 0.5*dvlim[i][nv];
    }
    #if (PHYSICS == RHD || PHYSICS == RMHD)
     VelocityLimiter(v[i], vp[i], vm[i]);
    #endif
  }

/*  -------------------------------------------
               Shock flattening
    -------------------------------------------  */

  #if SHOCK_FLATTENING == MULTID && CHAR_LIMITING == NO
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
        Shock flattening
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg - 1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
      evolve center values by dt/2
   -------------------------------------------------------- */


  #if TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}

#endif

#if CHAR_LIMITING == YES
/* *********************************************************************** */
void States (const State_1D *state, int beg, int end,  Grid *grid)
/*! 
 *   Compute 1D left and right interface states using piecewise
 *   linear reconstruction and the characteristic decomposition of the
 *   quasi-linear form of the equations.
 *
 *   This is done by first extrapolating the cell center value to the 
 *   interface using piecewise limited linear reconstruction
 *   on the characteristic variables.
 *
 *   Left and right states are then evolved for the half time step 
 *   using characteristic tracing if necessary.
 *
 ************************************************************************* */
{
  int    i, j, k, nv;
  double dvp[NVAR], dvm[NVAR], dv_lim[NVAR], dvc[NVAR], d2v;
  double dw_lim[NVAR], dwp[NVAR], dwm[NVAR];
  double dp, dm, dc;
  double **vp, **vm, **v;
  double **L, **R, *lambda;
  double cp, cm, wp, wm, cpk[NVAR], cmk[NVAR];
  double kstp[NVAR];
  PLM_Coeffs plm_coeffs;
  static double **dv;

/* ---------------------------------------------
   0. Allocate memory and set pointer shortcuts
   --------------------------------------------- */

  if (dv == NULL) {
    dv = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  v  = state->v; 
  vp = state->vp;
  vm = state->vm;

  #if UNIFORM_CARTESIAN_GRID == NO
   PLM_CoefficientsGet(&plm_coeffs, g_dir);
  #endif

/* ---------------------------------------------
    Compute sound speed and then convert 
    to 4-vel if necessary. 
   --------------------------------------------- */
 
  SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);

#if RECONSTRUCT_4VEL
   ConvertTo4vel (state->v, beg-1, end+1); /* From now up to the end of the
                                            * function, v contains the 4-vel */
#endif

  for (i = beg-1; i <= end; i++){
    NVAR_LOOP(nv) dv[i][nv] = v[i+1][nv] - v[i][nv];
  }

/* --------------------------------------------------------------
    Set the amount of steepening for each characteristic family.
    Default is 2, but nonlinear fields may be safely set to 1
    for strongly nonlinear problems.
   -------------------------------------------------------------- */

  for (k = NVAR; k--;  ) kstp[k] = 2.0;
  #if PHYSICS == HD || PHYSICS == RHD
   kstp[0] = kstp[1] = 1.0;
  #elif PHYSICS == MHD
   kstp[KFASTP] = kstp[KFASTM] = 1.0;
   #if COMPONENTS > 1
    kstp[KSLOWP] = kstp[KSLOWM] = 1.0; 
   #endif
  #endif

/* --------------------------------------------------------------
   2. Start main spatial loop
   -------------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

    PrimEigenvectors(v[i], state->a2[i], state->h[i], lambda, L, R);

  /* ---------------------------------------------------------------
     2a. Project forward, backward (and centered) undivided 
         differences of primitive variables along characteristics, 
         dw(k) = L(k).dv
     --------------------------------------------------------------- */

    #if UNIFORM_CARTESIAN_GRID == YES
     cp = cm = 2.0;
     wp = wm = 1.0;
     dp = dm = 0.5;
     NVAR_LOOP(nv) {
       dvp[nv] = dv[i][nv];
       dvm[nv] = dv[i-1][nv];

       k = nv;
       cpk[k] = cmk[k] = kstp[k];
     }
    #else
     cp = plm_coeffs.cp[i]; cm = plm_coeffs.cm[i];
     wp = plm_coeffs.wp[i]; wm = plm_coeffs.wm[i];
     dp = plm_coeffs.dp[i]; dm = plm_coeffs.dm[i];
     NVAR_LOOP(nv) {
       dvp[nv] = dv[i][nv]*wp;
       dvm[nv] = dv[i-1][nv]*wm;

     /* -- Map 1 < kstp < 2  ==>  1 < ck < c -- */

       k = nv;
       cpk[k] = (2.0 - cp) + (cp - 1.0)*kstp[k]; /* used only for general */
       cmk[k] = (2.0 - cm) + (cm - 1.0)*kstp[k]; /* minmod limiter        */
     }
    #endif

    PrimToChar(L, dvm, dwm);
    PrimToChar(L, dvp, dwp);

  /* ----------------------------------------------------------
     2b. Apply slope limiter to characteristic differences for
         nv < NFLX, dw = Limiter(dwp, dwm). 
     ------------------------------------------------------- */

    #if SHOCK_FLATTENING == MULTID
     if (state->flag[i] & FLAG_FLAT) {
       for (k = NFLX; k--;   )   dw_lim[k] = 0.0;
     } else if (state->flag[i] & FLAG_MINMOD) {
       for (k = NFLX; k--;   ){
         SET_MM_LIMITER(dw_lim[k], dwp[k], dwm[k], cp, cm);
       }
     } else
    #endif
    for (k = NFLX; k--;    ){
      #if LIMITER == DEFAULT
       SET_GM_LIMITER(dw_lim[k], dwp[k], dwm[k], cpk[k], cmk[k]);
      #else 
       SET_LIMITER(dw_lim[k], dwp[k], dwm[k], cp, cm);
      #endif
    }

  /* ------------------------------------------------------------------
     2c. Project limited slopes in characteristic variables on right 
         eigenvectors to obtain primitive slopes: dv = \sum dw.R

         Also, enforce monotonicity in primitive variables as well.
     ------------------------------------------------------------------ */

    for (nv = NFLX; nv--;   ){
      dc = 0.0;
      for (k = 0; k < NFLX; k++) dc += dw_lim[k]*R[nv][k];

      if (dvp[nv]*dvm[nv] > 0.0){
        d2v        = ABS_MIN(cp*dvp[nv], cm*dvm[nv]);
        dv_lim[nv] = MINMOD(d2v, dc);
      }else dv_lim[nv] = 0.0;
    }

  /* -----------------------------------------------------------------
     2d. Repeat construction for passive scalars (tracers).
         For a passive scalar, the primitive variable is the same as
         the characteristic one. We use the MC limiter always.
     ----------------------------------------------------------------- */
     
    #if NFLX != NVAR
     for (nv = NFLX; nv < NVAR; nv++ ){
       #if LIMITER == DEFAULT
        SET_MC_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
       #else
        SET_LIMITER(dv_lim[nv], dvp[nv], dvm[nv], cp, cm);
       #endif
     }
    #endif

  /* --------------------------------------------------------------------
     2e. Build L/R states at time level t^n
     -------------------------------------------------------------------- */

    for (nv = NVAR; nv--;   ) {
      vp[i][nv] = v[i][nv] + dv_lim[nv]*dp;
      vm[i][nv] = v[i][nv] - dv_lim[nv]*dm;
    }

    #if (PHYSICS == RHD || PHYSICS == RMHD)
     VelocityLimiter (v[i], vp[i], vm[i]);
    #endif

  }  /* -- end main loop on grid points -- */

/* ----------------------------------------------------
   3a. Check monotonicity
   ---------------------------------------------------- */ 

  #if CHECK_MONOTONICITY == YES
   MonotonicityTest(state->v, state->vp, state->vm, beg, end);
  #endif

/*  -------------------------------------------
    3b. Shock flattening (only 1D)
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
    4. Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
   5. Evolve L/R states and center value by dt/2
   -------------------------------------------------------- */

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #elif TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* --------------------------------------------------------
   6. Convert back to 3-velocity
   -------------------------------------------------------- */

#if RECONSTRUCT_4VEL
  ConvertTo3vel (state->v, beg-1, end+1);
  ConvertTo3vel (state->vp, beg, end);
  ConvertTo3vel (state->vm, beg, end);  
#endif

/* ----------------------------------------------
   7. Obtain L/R states in conservative variables
   ---------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
#endif  /* CHAR_LIMITING == YES */


/* ********************************************************************* */
void MonotonicityTest(double **v, double **vp, double **vm, int beg, int end)
/*
 *
 *********************************************************************** */
{
  int    i, nv;
  double vp_max, vp_min, vm_max, vm_min;
  
  for (i = beg; i <= end; i++){ VAR_LOOP(nv) {
    vp_max = MAX(v[i][nv], v[i+1][nv]) + 1.e-9;
    vp_min = MIN(v[i][nv], v[i+1][nv]) - 1.e-9;

    vm_max = MAX(v[i][nv], v[i-1][nv]) + 1.e-9;
    vm_min = MIN(v[i][nv], v[i-1][nv]) - 1.e-9;
  
    if (vp[i][nv] > vp_max || vp[i][nv] < vp_min){
       print ("! States: monotonicity violated at + interface, i = %d\n",i);
       print ("vp = %12.6e; vmax = %12.6e, vmin = %12.6e\n", 
               vp[i][nv], vp_max, vp_min);
       QUIT_PLUTO(1);
    }

    if (vm[i][nv] > vm_max || vm[i][nv] < vm_min){
      print ("! States: monotonicity violated at - interface, i = %d\n",i);
      print ("vm = %12.6e; vmax = %12.6e, vmin = %12.6e\n",
              vm[i][nv], vm_max, vm_min);
      QUIT_PLUTO(1);
    }
  }}
}
      
