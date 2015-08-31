/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Enforce conservation at the X1 boundaries in the shearing-box
         module.

  This file provides functions to store and then modify the upwind 
  fluxes computed during the Riemann solver at the leftmost 
  and righmost physical boundaries.
  These tasks are performed only when either SB_SYMMETRIZE_HYDRO, 
  SB_SYMMETRIZE_EY, SB_SYMMETRIZE_EZ flags are set to YES in 
  Src/MHD/shearingbox.h and are useful to avoid loss of conservation 
  in the hydrodynamical variables (density, momentum, energy) 
  and/or magnetic fields. 
  
  This is first achieved by calling the SB_SaveFluxes() function during
  the time stepping scheme, with the purpose of storing the leftmost and 
  rightmost conservative fluxes in the x-direction into the static
  arrays FluxL[] and FluxR[]. 

  These fluxes are then subsequently used by SB_CorrectFluxes() which 
  interpolates the fluxes and properly correct leftmost and rightmost 
  cell-centered flow quantities to ensure conservation.
  
  The treatment of staggered magnetic field is done similarly by 
  SB_CorrectEMF().
  
 \b References
   - "??" \n
     Mignone et al, in preparation

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.musicanisi@cineca.it)\n
           G. Bodo (bodo@oato.inaf.it)
           
  \date   Feb 14, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double ***FluxL; /**< Array of fluxes at the left x-boundary. */ 
static double ***FluxR; /**< Array of fluxes at the right x-boundary. */ 

#define NVLAST  (NVAR-1)

/*! Sets the weight coefficient in the range [0,1] used during flux
    symmetrization at the left and right boundaries, 
    
        FL -> swL*FL    + swR*I(FR)
        FR -> swL*I(FL) + swR*FR 
            
    where swR=1-swL and I() stands for conservative interpolation. */       
#define swL  0.0
#define swR  (1.0 - swL)

/* ********************************************************************* */
void SB_SaveFluxes (State_1D *state, Grid *grid)
/*!
 * Store leftmost and rightmost conservative fluxes in the x direction 
 * into FluxL and FluxR. 
 *
 * \param [in] state     pointer to a State_1D structure
 * \param [in] grid      pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int nv;

  #if SB_SYMMETRIZE_HYDRO == NO
   return;
  #endif

  if (FluxL == NULL){
    FluxL = ARRAY_3D(NVAR, NX3_TOT, NX2_TOT, double);
    FluxR = ARRAY_3D(NVAR, NX3_TOT, NX2_TOT, double);
  }

/* ---------------------------------------------------- 
    Save Fluxes on the LEFT physical boundary
   ---------------------------------------------------- */
  
  if (g_dir == IDIR && grid[IDIR].lbound != 0){
    state->flux[IBEG - 1][MX1] += state->press[IBEG - 1];
    state->press[IBEG - 1]      = 0.0;  
    for (nv = 0; nv <= NVLAST; nv++){
      FluxL[nv][g_k][g_j] = state->flux[IBEG - 1][nv];
    } 
  }

/* ---------------------------------------------------- 
    Save Fluxes on the RIGHT pysical boundary
   ---------------------------------------------------- */
  
  if (g_dir == IDIR && grid[IDIR].rbound != 0){
    state->flux[IEND][MX1] += state->press[IEND];    
    state->press[IEND]      = 0.0;  
    for (nv = 0; nv <= NVLAST; nv++){
      FluxR[nv][g_k][g_j] = state->flux[IEND][nv];
    }
  }
}

/* ********************************************************************* */
void SB_CorrectFluxes (Data_Arr U, double t, double dt, Grid *grid)
/*!
 * Interpolate x-fluxes FLuxL and FluxR and properly correct leftmost
 * and rightmost cells to ensure conservation.
 *
 *
 * \param [in,out] U   data array containing cell-centered quantities
 * \param [in] t       the time step at which fluxes have been computed
 * \param [in] dt      the time step being used in the integrator stage 
 * \param [in] grid    pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dtdx;
  static double **fL, **fR;
  static double ****Ftmp;
  
  #if SB_SYMMETRIZE_HYDRO == NO
   return;
  #endif

  if (fL == NULL){
    fL   = ARRAY_2D(NVAR, NMAX_POINT, double);
    fR   = ARRAY_2D(NVAR, NMAX_POINT, double);
    Ftmp = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, 1, double);
  }

  dtdx = dt/grid[IDIR].dx[IBEG];

{
t = g_time;
#if TIME_STEPPING == RK2
 if (g_intStage == 2) t += g_dt;
#elif TIME_STEPPING == RK3
 if (g_intStage == 2) t += 0.5*g_dt;
 if (g_intStage == 3) t += g_dt;
#endif
#ifdef CTU
 if (g_intStage == 2) t += 0.5*g_dt;
#endif
}

/* -----------------------------------------------------------------------
    Exchange left and right fluxes if they were initially computed on
    different processors.
    After this step, both the left- and right- processor- will share
    the same fluxes FluxL[] and FluxR[].
   ----------------------------------------------------------------------- */

  #ifdef PARALLEL
   ExchangeX (FluxL[0][0], FluxR[0][0], NVAR*NX2_TOT*NX3_TOT, grid);
  #endif

/* --------------------------------------------------------------------- */
/*! \note 
    Modify the values of FluxR (on the left) and FluxL (on the right)
    to properly account for the shift induced by the moving boundaries.
    For safety reason (avoid overwriting arrays when a processor owns 
    both the left and the right side) we first copy, then interpolate and
    finally average:

    - On the left boundary:

      FluxR -> FluxLR -> I(FluxLR) -> fL = wL*FluxL + wR*FluxLR

    - On the right boundary:

      FluxL -> FluxLR -> I(FluxLR) -> fR = wL*FluxLR + wR*FluxR
      
   --------------------------------------------------------------------- */

  if (grid[IDIR].lbound != 0){  /* ---- Left x-boundary ---- */

    RBox box;

  /* -- set grid ranges of FluxR and exchange b.c. -- */

    box.ib = 0; box.ie = 0; 
    box.jb = 0; box.je = NX2_TOT-1; 
    box.kb = 0; box.ke = NX3_TOT-1;

  /* -- copy FluxR on temporary storage and shift -- */

    for (nv = 0; nv <= NVLAST; nv++) {
      BOX_LOOP((&box), k, j, i) Ftmp[nv][k][j][i] = FluxR[nv][k][j];
      SB_SetBoundaryVar(Ftmp[nv], &box, X1_BEG, t, grid);
    }

  /* -- symmetrize fluxes -- */
  
    for (k = KBEG; k <= KEND; k++){
      for (nv = 0; nv <= NVLAST; nv++){
      for (j = JBEG; j <= JEND; j++){
        fL[nv][j] = swL*FluxL[nv][k][j] + swR*Ftmp[nv][k][j][0];
      }}
      for (j = JBEG; j <= JEND; j++){
        #if HAVE_ENERGY
         fL[ENG][j] += swR*sb_vy*(fL[MX2][j] + 0.5*sb_vy*fL[RHO][j]);
        #endif
        #ifdef GLM_MHD
         fL[BX2][j] -= swR*sb_vy*fL[PSI_GLM][j]/glm_ch/glm_ch;
/*         fL[BX2][j] -= 0.5*sb_vy*FluxL[PSI_GLM][k][j]/glm_ch/glm_ch; */
        #endif
        fL[MX2][j] += swR*sb_vy*fL[RHO][j];
        
      /* -- update solution vector -- */

        for (nv = 0; nv <= NVLAST; nv++){
          U[k][j][IBEG][nv] += dtdx*(fL[nv][j] - FluxL[nv][k][j]); 
        } 
      }
    }
  }   

  if (grid[IDIR].rbound != 0){  /* ---- Right x-boundary ---- */

  /* -- set grid ranges of FluxL and exchange b.c. -- */

    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;

  /* -- copy FluxL on temporary storage and shift -- */

    for (nv = 0; nv <= NVLAST; nv++) {
      BOX_LOOP((&box), k, j, i) Ftmp[nv][k][j][i] = FluxL[nv][k][j];
      SB_SetBoundaryVar(Ftmp[nv], &box, X1_END, t, grid);
    }

  /* -- symmetrize fluxes -- */
  
    for (k = KBEG; k <= KEND; k++){
      for (nv = 0; nv <= NVLAST; nv++){
      for (j = JBEG; j <= JEND; j++){
        fR[nv][j] = swL*Ftmp[nv][k][j][0] + swR*FluxR[nv][k][j];
      }}
      for (j = JBEG; j <= JEND; j++){
        #if HAVE_ENERGY
         fR[ENG][j] += swL*sb_vy*(-fR[MX2][j] + 0.5*sb_vy*fR[RHO][j]);
        #endif
        #ifdef GLM_MHD
         fR[BX2][j] += swL*sb_vy*fR[PSI_GLM][j]/glm_ch/glm_ch;
/*         fR[BX2][j] += 0.5*sb_vy*FluxR[PSI_GLM][k][j]/glm_ch/glm_ch; */
        #endif
        fR[MX2][j] -= swL*sb_vy*fR[RHO][j];

        for (nv = 0; nv <= NVLAST; nv++){
          U[k][j][IEND][nv] -= dtdx*(fR[nv][j] - FluxR[nv][k][j]); 
        }
      }
    }  
  }
}

#ifdef STAGGERED_MHD
/* ********************************************************************* */
void SB_CorrectEMF (EMF *emf, Data_Arr V0, Grid *grid)
/*! 
 * Correct the electromotive force by enforcing symmetrization in order 
 * to guarantee that Bx on the left is mapped from Bx on the right.
 * This is done by re-defining the y- and z-components of the electric 
 * field as follows:
 *
 * \f[ 
 *    E_{z,L} \to \frac{1}{2}\left[E_{z,L} + \Aop^{t_E}_{>}(E_{z,R})
 *          + wF_{+}\right]
 *  \,,\qquad
 *    E_{z,R} \to \frac{1}{2}\left[\Aop^{t_E}_<(E_{z,L}) + E_{z,R}
 *          - wF_-\right]
 * \f]
 * where \f$t_E =\f$ is time level of the electric field, $\Aop$ are 
 * the forward (>) or backward (<) shift operators.
 * The corrective terms \f$F_{\pm}\f$ are the upwind numerical flux 
 * computed in the advection of the x-component of magnetic field:
 * \f[
 *    F_{+} = B_{x,L} + \frac{\widetilde{\Delta B}_{xL}}{2}
 *            \left(1 - \frac{\delta t}{\Delta y}w\right)
 *  \,,\qquad
 *    F_{-} = B_{x,R} - \frac{\widetilde{\Delta B}_{xR}}{2}
 *            \left(1 - \frac{\delta t}{\Delta y}w\right)
 * \f]
 * where \f$ \delta t\ne \Delta t\f$ is the time increment of the shift
 * operator applied to the magnetic field, \f$B_{xL}\f$ and \f$B_{xR}\f$ 
 * are the x component of magnetic field computed at some intermediate
 * time \c tB while 
 * \f[
 *   \widetilde{\Delta B}_{xL} = 
 *   \frac{                \overline{\Delta B}_{xL} 
 *           + \Aop^{tB}_{>}(\overline{\Delta B}_{xR})}{2}
 *    \;,\qquad
 *   \widetilde{\Delta B}_{xR} = 
 *   \frac{ \Aop^{tB}_{<}(\overline{\Delta B}_{xL}) 
 *                      + \overline{\Delta B}_{xR}}{2}         
 * \f]
 * are the arithmetic average of the limited slopes of \c Bx computed
 * at the left and right boundaries.
 * 
 * \param [in,out] emf    pointer to EMF structure
 * \param [in]     V0     4D array of staggered fields
 * \param [in]   grid     pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nghost;
  int    dimx[3] = {1, 0, 0};
  int    dimy[3] = {0, 1, 0};
  int    dimz[3] = {0, 0, 1};
  double   fE, esym;
  double   tB, tE, w, dt, dtdy;
  static double **BxL0, **BxL, ***dBxL_lim, *dBxL, ***eyL, ***ezL;
  static double **BxR0, **BxR, ***dBxR_lim, *dBxR, ***eyR, ***ezR;
  static double ***etmp, ***dBtmp;

  #if    (SB_SYMMETRIZE_EY == NO) && (SB_SYMMETRIZE_EZ == NO) \
      && (SB_FORCE_EMF_PERIODS == NO)
   return;
  #endif

  if (ezL == NULL){
    eyL  = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    eyR  = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);

    ezL  = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    ezR  = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);

    BxL  = ARRAY_2D(NX3_TOT, NX2_TOT, double);
    BxR  = ARRAY_2D(NX3_TOT, NX2_TOT, double);

    BxL0 = ARRAY_2D(NX3_TOT, NX2_TOT, double);
    BxR0 = ARRAY_2D(NX3_TOT, NX2_TOT, double);

    dBxL     = ARRAY_1D(NMAX_POINT, double);
    dBxR     = ARRAY_1D(NMAX_POINT, double);
    dBxL_lim = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    dBxR_lim = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);

    etmp  = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double); /* temporary storage */
    dBtmp = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double); /* temporary storage */
  }

  w  = sb_vy;         /* shear velocity */
  dt = g_dt;          /* default increment of the shift operator A(B) */
  tB = g_time;        /* time level of the magnetic field V0 */
  tE = g_time + g_dt; /* increment for shift operator A(E) */

/* --------------------------------------------------------------
       Store initial Bx on the left and on the right
   -------------------------------------------------------------- */

  if (g_intStage == 1){
    if (grid[IDIR].lbound != 0){
      KTOT_LOOP(k) JTOT_LOOP(j) BxL[k][j] = BxL0[k][j] = V0[BX1s][k][j][IBEG - 1]; 
    }
    if (grid[IDIR].rbound != 0){
      KTOT_LOOP(k) JTOT_LOOP(j) BxR[k][j] = BxR0[k][j] = V0[BX1s][k][j][IEND];
    } 
  }

  #ifdef CTU
   if (grid[IDIR].lbound != 0){
     KTOT_LOOP(k) JTOT_LOOP(j) BxL[k][j] = BxL0[k][j] = V0[BX1s][k][j][IBEG - 1]; 
   }
   if (grid[IDIR].rbound != 0){
     KTOT_LOOP(k) JTOT_LOOP(j) BxR[k][j] = BxR0[k][j] = V0[BX1s][k][j][IEND];
   } 
  #elif TIME_STEPPING == RK2
   if (g_intStage == 2) {
     tB = g_time + 0.5*g_dt;
     dt = 0.5*g_dt;
     if (grid[IDIR].lbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxL[k][j] = 0.5*(BxL0[k][j] + V0[BX1s][k][j][IBEG - 1]);
     }
     if (grid[IDIR].rbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxR[k][j] = 0.5*(BxR0[k][j] + V0[BX1s][k][j][IEND]);
     }
   } 
  #elif TIME_STEPPING == RK3
   if (g_intStage == 2) {
     tB = g_time + 0.25*g_dt;
     tE = g_time + 0.5*g_dt;
     dt = 0.25*g_dt;
     if (grid[IDIR].lbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxL[k][j] = 0.75*BxL0[k][j] + 0.25*V0[BX1s][k][j][IBEG - 1];
     }
     if (grid[IDIR].rbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxR[k][j] = 0.75*BxR0[k][j] + 0.25*V0[BX1s][k][j][IEND];
     }
   }else if (g_intStage == 3) {
     tB = g_time + g_dt/3.0;
     tE = g_time + g_dt;
     dt = g_dt/1.5;
     if (grid[IDIR].lbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxL[k][j] = (BxL0[k][j] + 2.0*V0[BX1s][k][j][IBEG - 1])/3.0;
     }
     if (grid[IDIR].rbound != 0) KTOT_LOOP(k) JTOT_LOOP(j) {
       BxR[k][j] = (BxR0[k][j] + 2.0*V0[BX1s][k][j][IEND])/3.0;
     }
   } 
  #endif

  dtdy   = dt/grid[JDIR].dx[JBEG];
  nghost = grid[IDIR].nghost;

/* -- compute Bx slopes on left side -- */

  if (grid[IDIR].lbound != 0){

  /* -- Store Ey, Ez, Bx on the left and compute slopes -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
      D_EXPAND(                                         ,
               ezL[k][j][0] = emf->ez[k][j][IBEG - 1];  ,
               eyL[k][j][0] = emf->ey[k][j][IBEG - 1];)
    }

    for (k = KBEG; k <= KEND; k++){
      for (j = 1; j <= JEND + nghost; j++) dBxL[j] = BxL[k][j] - BxL[k][j-1]; 
      for (j = JBEG-1; j <= JEND+1; j++){
        dBxL_lim[k][j][0] = VAN_LEER(dBxL[j+1], dBxL[j]);
      }
    }
  }

/* -- compute Bx slopes on right side -- */

  if (grid[IDIR].rbound != 0){

  /* -- Store Ey, Ez, Bx on the right and compute slopes -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
      D_EXPAND(                                      ,
               ezR[k][j][0] = emf->ez[k][j][IEND];   ,
               eyR[k][j][0] = emf->ey[k][j][IEND]; )
    }

    for (k = KBEG; k <= KEND; k++){
      for (j = 1; j <= JEND + nghost; j++) dBxR[j] = BxR[k][j] - BxR[k][j-1];
      for (j = JBEG-1; j <= JEND+1; j++){
        dBxR_lim[k][j][0] = VAN_LEER(dBxR[j+1], dBxR[j]);
      }
    }
  }

/* ---------------------------------------------------------------------
            exchange data between processors
   --------------------------------------------------------------------- */

  #ifdef PARALLEL
   D_EXPAND(                                                                   ,
            ExchangeX (ezL[0][0], ezR[0][0], NX3_TOT*NX2_TOT, grid); 
            ExchangeX (dBxL_lim[0][0], dBxR_lim[0][0], NX3_TOT*NX2_TOT, grid); ,
            ExchangeX (eyL[0][0], eyR[0][0], NX3_TOT*NX2_TOT, grid); )
  #endif

/* ----------------------------------------------------------------
    Symmetrize Ey and Ez on the left.
    We copy one of the two arrays before doing interpolation since
    its original value is lost after b.c. have been assigned.
   ---------------------------------------------------------------- */

  if (grid[IDIR].lbound != 0){
    RBox box;
   
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;
    #if SB_SYMMETRIZE_EY == YES

   /* -- Interpolate eyR --> eyLR -- */

     BOX_LOOP((&box), k, j, i) etmp[k][j][i] = eyR[k][j][i];
     SB_SetBoundaryVar(etmp, &box, X1_BEG, tE, grid);
     for (k = KBEG - 1; k <= KEND; k++){
       for (j = JBEG; j <= JEND; j++) {
         emf->ey[k][j][IBEG - 1] = swL*eyL[k][j][0] + swR*etmp[k][j][0];
       }
     }
    #endif

/* --------------------------------------------------------------------- */
/*! \note 
    Symmetrization of Ez can only be done using weight coefficients of
    1/2 (and not swL, swR) since limited slopes (dbxL_lim and dbxR_lim) 
    do not satisfy the property that their sum is zero (as expected for
    a periodic function).                                                */
/* --------------------------------------------------------------------- */

    #if SB_SYMMETRIZE_EZ == YES

    /* -- interpolate dbxR_lim, ezR --> dbxLR_lim, ezLR -- */

     BOX_LOOP((&box), k, j, i) {
       dBtmp[k][j][i] = dBxR_lim[k][j][i];
        etmp[k][j][i] =      ezR[k][j][i];
     }
     SB_SetBoundaryVar(dBtmp, &box, X1_BEG, tB, grid);
     SB_SetBoundaryVar( etmp, &box, X1_BEG, tE, grid);
     for (k = KBEG; k <= KEND; k++){ 
       for (j = JBEG - 1; j <= JEND + 1; j++){
         dBxL[j] = 0.5*(dBxL_lim[k][j][0] + dBtmp[k][j][0]);
       }
       for (j = JBEG - 1; j <= JEND; j++){
         fE = w*(BxL[k][j] + 0.5*dBxL[j]*(1.0 - w*dtdy));
         emf->ez[k][j][IBEG-1] = 0.5*(ezL[k][j][0] + etmp[k][j][0] + fE);
       }
     }
    #endif
  }

/* ------------------------------------------------------
         Symmetrize Ey, and Ez on the right
   ------------------------------------------------------ */

  if (grid[IDIR].rbound != 0){
    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;

    #if SB_SYMMETRIZE_EY == YES
     SB_SetBoundaryVar(eyL, &box, X1_END, tE, grid);
     for (k = KBEG - 1; k <= KEND; k++){
       for (j = JBEG; j <= JEND; j++) {
         emf->ey[k][j][IEND] = swL*eyL[k][j][0] + swR*eyR[k][j][0];
       }
     } 
    #endif

    #if SB_SYMMETRIZE_EZ == YES
     SB_SetBoundaryVar(dBxL_lim, &box, X1_END, tB, grid);
     SB_SetBoundaryVar(     ezL, &box, X1_END, tE, grid);

     for (k = KBEG; k <= KEND; k++){ 
       for (j = JBEG - 1; j <= JEND + 1; j++){
         dBxR[j] = 0.5*(dBxL_lim[k][j][0] + dBxR_lim[k][j][0]);
       }
       for (j = JBEG - 1; j <= JEND; j++){
         fE = w*(BxR[k][j + 1] - 0.5*dBxR[j + 1]*(1.0 - w*dtdy));
         emf->ez[k][j][IEND] = 0.5*(ezR[k][j][0] + ezL[k][j][0] - fE);
       }
     }
    #endif
  }

/* --------------------------------------------------
                Force Periodicity 
   -------------------------------------------------- */

  #if DIMENSIONS == 3 && SB_FORCE_EMF_PERIODS == YES
               
  /* -- Ex at Z faces: force periodicty  -- */
                                                                     
   #ifdef PARALLEL
    MPI_Barrier (MPI_COMM_WORLD);
    AL_Exchange_dim (emf->ex[0][0], dimz, SZ);
    MPI_Barrier (MPI_COMM_WORLD);
   #else
    for (j = JBEG - 1; j <= JEND; j++){
    for (i = IBEG    ; i <= IEND; i++){
      esym = 0.5*(emf->ex[KBEG - 1][j][i] + emf->ex[KEND][j][i]);
      emf->ex[KBEG - 1][j][i] = esym;
      emf->ex[KEND    ][j][i] = esym;
    }}
   #endif

   /*  Ex at Y faces: force periodicity  */

   for (k = KBEG - 1; k <= KEND; k++){
   for (i = IBEG    ; i <= IEND; i++){
     esym = 0.5*(emf->ex[k][JBEG - 1][i] + emf->ex[k][JEND][i]);
     emf->ex[k][JBEG - 1][i] = esym;
     emf->ex[k][JEND    ][i] = esym;
   }}
     
   /*  Ey at Z faces: force periodicity   */

   #ifdef PARALLEL
    MPI_Barrier (MPI_COMM_WORLD);
    AL_Exchange_dim (emf->ey[0][0], dimz, SZ);
    MPI_Barrier (MPI_COMM_WORLD);
   #else
    for (j = JBEG    ; j <= JEND; j++){
    for (i = IBEG - 1; i <= IEND; i++){
      esym = 0.5*(emf->ey[KBEG - 1][j][i] + emf->ey[KEND][j][i]);
      emf->ey[KBEG - 1][j][i] = esym;
      emf->ey[KEND    ][j][i] = esym;
    }}
   #endif

   /*  Ez at Y faces: force periodicity   */
 
   for (k = KBEG    ; k <= KEND; k++){
   for (i = IBEG - 1; i <= IEND; i++){
     esym = 0.5*(ez[k][JBEG - 1][i] + ez[k][JEND][i]);
     ez[k][JBEG - 1][i] = esym;
     ez[k][JEND    ][i] = esym;
   }}
 
  #endif 
}
#endif

#ifdef PARALLEL
/* ********************************************************************* */
void ExchangeX (double *bufL, double *bufR, int nel, Grid *grid)
/*!
 * Send bufL owned by processor at X1_BEG to the processor at X1_END;
 * Send bufR owned by processor at X1_END to the processor at X1_BEG;
 *
 *********************************************************************** */
{
  static int dest = -1;
  int stag = 1, rtag = 1;
  MPI_Comm cartcomm;
  MPI_Status istat;
  MPI_Request req;
  static  int nprocs[3], periods[3], coords[3];

  AL_Get_cart_comm(SZ, &cartcomm);
  MPI_Barrier (MPI_COMM_WORLD);

/* --------------------------------------
     get rank of the processor lying 
     on the opposite side of x-domain
   -------------------------------------- */

  if (dest == -1){
    if (grid[IDIR].lbound != 0){
      MPI_Cart_get(cartcomm, 3, nprocs, periods, coords);
      coords[0] += nprocs[0] - 1;
      MPI_Cart_rank (cartcomm, coords, &dest);
    }

    if (grid[IDIR].rbound != 0){
      MPI_Cart_get(cartcomm, 3, nprocs, periods, coords);
      coords[0] += - nprocs[0] + 1;
      MPI_Cart_rank (cartcomm, coords, &dest); 
    }
  }

  if (grid[IDIR].lbound != 0){
    if (prank != dest){
      MPI_Sendrecv (bufL, nel, MPI_DOUBLE, dest, stag,
                    bufR, nel, MPI_DOUBLE, dest, rtag,
                    MPI_COMM_WORLD, &istat);
    }
  }

  if (grid[IDIR].rbound != 0){
    if (prank != dest){
      MPI_Sendrecv (bufR, nel, MPI_DOUBLE, dest, stag,
                    bufL, nel, MPI_DOUBLE, dest, rtag,
                    MPI_COMM_WORLD, &istat);
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);
}

#endif
