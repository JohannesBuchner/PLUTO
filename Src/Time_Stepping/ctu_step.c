/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations using Corner Transport Upwind (CTU).

  Implement the dimensionally unsplit, Corner Transport Upwind method 
  (CTU) of Colella. It consists of
 
  <b> A predictor step </b>: 
  -# start from cell center values of \c V and compute time-centered 
     interface states using either HANCOCK or CHARACTERISTIC_TRACING:
     \f[
       V^{n+\HALF}_{i,\pm} = V^n_{i,\pm} + \pd{V}{t}\frac{\Delta t}{2}
     \f] 
     Store the resulting states by converting \f$ V_\pm \f$ in \f$ U_\pm \f$
     (= \c UP, \c UM).
     Normal (n) and trasnverse (t) indexes for this step should span
     \f$ n \in [n_{\rm beg} - 1, n_{\rm end} + 1] \f$, 
     \f$ t \in [t_{\rm beg} - 1, t_{\rm end} + 1] \f$ 
     for cell-centered schemes and
     \f$ n \in [n_{\rm beg} - 2, n_{\rm end} + 2] \f$,  
     \f$ t \in [t_{\rm beg} - 2, t_{\rm end} + 2] \f$ for staggered MHD. 
  -# Solve Riemann problems between \f$ V_{i,+}\f$ and \f$ V_{i+1,-}\f$  
     and store the flux difference in RHS.
  -# Compute cell and time-centered values UH = U + dt/2*RHS with index
     range \f$ n \in [n_{\rm beg}, n_{\rm end}]\f$, 
           \f$ t \in [t_{\rm beg}, t_{\rm end}]\f$ for cell-centered schemes
     and   \f$ n \in [n_{\rm beg}-1, n_{\rm end}+1]\f$, 
           \f$ t \in [t_{\rm beg}-1, t_{\rm end}+1]\f$ for staggered MHD.
   
  <b> A corrector step </b>: 
  -# correct left and right states \c UP and \c UM computed in the previous 
     step by adding the transverse RHS contribution, e.g.
     \f$ U^{n+\HALF}_{i\pm} += ({\cal RHS}^n_y + {\cal RHS}^n_z)\Delta t/2 \f$.
     This step should cover \f$ [n_{\rm beg}-1, n_{\rm end} + 1]\f$, 
     \f$ [t_{\rm beg}, t_{\rm end}]\f$ for cell-centered MHD and
     \f$ [n_{\rm beg}-2, n_{\rm end} + 2]\f$, \f$ 
     [t_{\rm beg}-1, tend+1]\f$ for staggered MHD.
 
  -# Solve a normal Riemann problem with left and right states
     \c UP and \c UM, get RHS and evolve to final stage,
     \f$ U^{n+1} = U^n + \Delta t^n {\cal RHS}^{n+\HALF} \f$
   
 
  This integrator performs an integration in the ghost boundary zones, 
  in order to recover appropriate information to build the transverse 
  predictors.
 
  \note If explicit parabolic flux have to be included, the predictor
        step is modified by computing the RHS from 1st order states at t^n.
        This allows to obtain space and time centered states in one 
        row of boundary zones required during the following corrector step.
  
  \b References
     - "The Piecewise Parabolic Method for Multidimensional 
        Relativistic Fluid Dynamics" \n
       Mignone, A.; Plewa, T.; Bodo, G. ApJS (2005), 160..199M 
     - "An unsplit Godunov method for ideal MHD via constrained transport" \n
        Gardiner & Stone, JCP (2005), 205, 509
     - "A second-order unsplit Godunov scheme for cell-centered MHD:
        the CTU-GLM scheme" \n
        Mignone & Tzeferacos  JCP (2010), 229, 2117
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    June 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#ifdef STAGGERED_MHD
 #if TIME_STEPPING == CHARACTERISTIC_TRACING
  #define CTU_MHD_SOURCE YES
 #elif TIME_STEPPING == HANCOCK && PRIMITIVE_HANCOCK == YES
  #define CTU_MHD_SOURCE YES
 #else
  #define CTU_MHD_SOURCE NO
 #endif
#else
 #define CTU_MHD_SOURCE NO
#endif

#if CTU_MHD_SOURCE == YES
 static void CTU_CT_Source (double **, double **, double **,
                            double *, int, int, Grid *);
#endif
/* ********************************************************************* */
int AdvanceStep (const Data *d, Riemann_Solver *Riemann, 
                 Time_Step *Dts, Grid *grid)
/*!
 * Advance equations using the corner transport upwind method
 * (CTU)
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in]    Riemann  pointer to a Riemann solver function
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *         
 *********************************************************************** */
{
  int i, j,k, nv, *in;
  int errp, errm;
  double dt2, inv_dtp, *inv_dl;
  Index indx;
  static unsigned char *flagm, *flagp;
  static Data_Arr UM[DIMENSIONS], UP[DIMENSIONS];  

  static Data_Arr dU, UH, Bs0;
  static State_1D state;
  static double **dtdV, **dcoeff, ***T;
  double *dtdV2, **rhs;
  RBox *box = GetRBox(DOM, CENTER);
  
/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL)
   print1 ("! AdvanceStep(): ");
   print1 ("CTU only works in Cartesian or cylindrical coordinates\n");
   QUIT_PLUTO(1);
  #endif

/*  --------------------------------------------------------
                   Allocate static memory areas   
   --------------------------------------------------------- */

  if (dU == NULL){

    dtdV = ARRAY_2D(DIMENSIONS,NMAX_POINT, double);

    MakeState (&state);

    flagp = ARRAY_1D(NMAX_POINT, unsigned char);
    flagm = ARRAY_1D(NMAX_POINT, unsigned char);

    UH  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    dU  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    #ifdef STAGGERED_MHD
     Bs0 = ARRAY_4D(DIMENSIONS, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
    
  /* ---------------------------------------------------------
      corner-coupled multidimensional arrays are stored into 
      memory following the same conventions adopted when 
      sweeping along the coordinate directions, i.e., 

         (z,y,x)->(z,x,y)->(y,x,z).

      This allows 1-D arrays to conveniently point at the 
      fastest running indexes of the respective multi-D ones.
     --------------------------------------------------------- */  
     
    UM[IDIR] = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    UP[IDIR] = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

    UM[JDIR] = ARRAY_4D(NX3_TOT, NX1_TOT, NX2_TOT, NVAR, double);
    UP[JDIR] = ARRAY_4D(NX3_TOT, NX1_TOT, NX2_TOT, NVAR, double);

    #if DIMENSIONS == 3
     UM[KDIR] = ARRAY_4D(NX2_TOT, NX1_TOT, NX3_TOT, NVAR, double);
     UP[KDIR] = ARRAY_4D(NX2_TOT, NX1_TOT, NX3_TOT, NVAR, double);
    #endif
 
    #if (PARABOLIC_FLUX & EXPLICIT)
     dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
    #if THERMAL_CONDUCTION == EXPLICIT
     T = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

/* -- save pointers -- */

  rhs = state.rhs;
/*  memset (dU[0][0][0], 0.0, NX3_TOT*NX2_TOT*NX1_TOT*NVAR*sizeof(double));   */

  #ifdef FARGO
   FARGO_SubtractVelocity (d,grid); 
  #endif

/* -------------------------------------------------------
    Set boundary conditions and flag shocked regions for
    shock flattening or energy/entropy selective update.
   ------------------------------------------------------- */

  g_intStage = 1;
  Boundary (d, ALL_DIR, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
   FlagShock (d, grid);
  #endif

  dt2 = 0.5*g_dt;
  D_EXPAND(
    ITOT_LOOP(i) dtdV[IDIR][i] = dt2/grid[IDIR].dV[i]; ,
    JTOT_LOOP(j) dtdV[JDIR][j] = dt2/grid[JDIR].dV[j]; ,
    KTOT_LOOP(k) dtdV[KDIR][k] = dt2/grid[KDIR].dV[k];
  )

  #if THERMAL_CONDUCTION == EXPLICIT
   TOT_LOOP(k,j,i) T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
  #endif
  
/* ------------------------------------------------
    Compute current arrays
   ------------------------------------------------ */

  #if (RESISTIVITY == EXPLICIT) && (defined STAGGERED_MHD)
   GetCurrent (d, -1, grid);
  #endif
 
/* ----------------------------------------------------
    Convert primitive to conservative and reset arrays
   ---------------------------------------------------- */

  KTOT_LOOP(k) JTOT_LOOP(j){
    ITOT_LOOP(i){
      VAR_LOOP(nv) state.v[i][nv] = d->Vc[nv][k][j][i];
    }
    PrimToCons(state.v, d->Uc[k][j], 0, NX1_TOT-1);
  }

  #ifdef STAGGERED_MHD
   nv = NX1_TOT*sizeof(double);
   KTOT_LOOP(k) JTOT_LOOP(j){
     D_EXPAND(memcpy(Bs0[BX1s][k][j], d->Vs[BX1s][k][j], nv);  ,
              memcpy(Bs0[BX2s][k][j], d->Vs[BX2s][k][j], nv);  ,
              memcpy(Bs0[BX3s][k][j], d->Vs[BX3s][k][j], nv);)
   }
  #endif
  
/* ----------------------------------------------------
     1. Compute Normal predictors and
        solve normal Riemann problems. 
        Store computations in UP, UM, RHS (X,Y,Z)
   ---------------------------------------------------- */

  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);
    ResetState (d, &state, grid);
    #if (RESISTIVITY == EXPLICIT) && !(defined STAGGERED_MHD)
     GetCurrent(d, g_dir, grid);
    #endif
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      g_i = i; g_j = j; g_k = k;

    /* ---------------------------------------------
        save computational time by having state.up 
        and state.um pointing at the fastest 
        running indexes of UP and UM. 
        Also, during the x-sweep we initialize Uc 
        and dU by changing the memory address of 
        state.u and state.rhs. 
       --------------------------------------------- */
        
      state.up = UP[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uL = state.up;
      state.um = UM[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uR = state.um + 1;

    /* ---- get a 1-D array of primitive quantities ---- */

      for (*in = 0; (*in) < indx.ntot; (*in)++) {
        VAR_LOOP(nv) state.v[*in][nv] = d->Vc[nv][k][j][i];
        state.flag[*in] = d->flag[k][j][i];
        #ifdef STAGGERED_MHD
         state.bn[*in] = d->Vs[g_dir][k][j][i];
        #endif
      }

      CheckNaN (state.v, 0, indx.ntot - 1, 0);

#if !(PARABOLIC_FLUX & EXPLICIT)   /* adopt this formulation when there're no
                                      explicit diffusion  flux terms */
      States  (&state, indx.beg - 1, indx.end + 1, grid);
      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
       CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, dt2, grid);

      #if CTU_MHD_SOURCE == YES
       CTU_CT_Source (state.v, state.up, state.um,
                           dtdV[g_dir], indx.beg - 1, indx.end + 1, grid);
      #endif

   /* ------------------------------------------------------------------
       At this point we have at disposal the normal predictors U^_\pm.
       To save memory, we compute corner coupled states by first
       subtracting the normal contribution and then by adding the total
       time increment. For example, in the x-direction, we do

           U^{n+1/2}_\pm = U^_\pm - RX + (RX + RY + RZ)

       where the first subtraction is done here while dU = (RX + RY + RZ) is
       the total time increment which will be added later (step 2 below).
       In a 2-D domain dU  contains the following (X,Y,Z) contributions:

        0    X     0        0    X     0
         +--------+          +--------+
         |        |          |        |
        0|   X    |0   -->  Y|   XY   |Y  
         |        |          |        |
         +--------+          +--------+
        0    X     0        0    X     0

         (X sweep)           (Y sweep)

       Also, evolve cell center values by dt/2. 
       For staggered mhd, this step should be carried also in one rows of
       boundary zones in order to provide the cell  and time centered
       e.m.f.
      ------------------------------------------------------------------ */

      if (g_dir == IDIR){
          
        for (nv = NVAR; nv--; ) {
          state.rhs[indx.beg - 1][nv] = 0.0;
          state.rhs[indx.end + 1][nv] = 0.0;
        }

        for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv]  = state.rhs[*in][nv];
          UH[k][j][i][nv]  = d->Uc[k][j][i][nv] + state.rhs[*in][nv];
          state.up[*in][nv]   -= state.rhs[*in][nv];
          state.um[*in][nv]   -= state.rhs[*in][nv];
        }}
      }else{
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] += state.rhs[*in][nv];
          UH[k][j][i][nv] += state.rhs[*in][nv];
          state.up[*in][nv]   -= state.rhs[*in][nv];
          state.um[*in][nv]   -= state.rhs[*in][nv];
        }}
      }

#else 

   /* -------------------------------------------------------------
       When parabolic terms have to be included explicitly, we use 
       a slightly different formulation where the transverse 
       predictors are computed using 1st order states. 
       This still yields a second-order accurate scheme but allow 
       to compute the solution at the half time step in one rows 
       of ghost zones, which is essential to update the parabolic 
       terms using midpoint rule:

            U^{n+1} = U^n - RHS(hyp, n+1/2) - RHS(par, n+1/2) 
      
           
       Since RHS(par) has a stencil 3-point wide.
       Corner coupled states are modified to account for parabolic 
       terms in the following way:

            U^{n+1/2}_\pm = U^_\pm - RX,h + (RX,h + RY,h + RZ,h) 
                                          + (RX,p + RY,p + RZ,p) 
      -------------------------------------------------------------- */

      for ((*in) = 0; (*in) < indx.ntot; (*in)++) {
      for (nv = NVAR; nv--;  ) {
        state.vp[*in][nv] = state.vm[*in][nv] = state.vh[*in][nv] = state.v[*in][nv];
      }}
      #ifdef STAGGERED_MHD
       for (*in = 0; (*in) < indx.ntot-1; (*in)++) {
         state.vR[*in][BXn] = state.vL[*in][BXn] = state.bn[*in];
       }
      #endif
      PrimToCons(state.vm, state.um, 0, indx.ntot-1);
      PrimToCons(state.vp, state.up, 0, indx.ntot-1);
      
      Riemann (&state, indx.beg-1, indx.end, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
       CT_StoreEMF (&state, indx.beg-1, indx.end, grid);
      #endif

  /* -----------------------------------------------------------
          compute rhs using the hyperbolic fluxes only
     ----------------------------------------------------------- */

      #if (VISCOSITY == EXPLICIT)
       for ((*in) = 0; (*in) < indx.ntot; (*in)++) for (nv = NVAR; nv--;  )
         state.par_src[*in][nv] = 0.0;
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, dt2, grid);
      ParabolicFlux (d->Vc, d->J, T, &state, dcoeff, indx.beg-1, indx.end, grid);

  /* ----------------------------------------------------------------
      compute LR states and subtract normal (hyperbolic) 
      rhs contribution.
      NOTE: states are computed from IBEG - 1 (= indx.beg) up to
            IEND + 1 (= indx.end) since EMF has already been evaluated
            and stored using 1st order states above.
      NOTE: States should be called after ParabolicFlux since some
            terms (e.g. thermal conduction) may depend on left and
            right normal states.
     ---------------------------------------------------------------- */

      States  (&state, indx.beg, indx.end, grid);
      #if CTU_MHD_SOURCE == YES
       CTU_CT_Source (state.v, state.up, state.um,
                           dtdV[g_dir], indx.beg, indx.end, grid);
      #endif

      for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
      for (nv = NVAR; nv--; ){
        state.up[*in][nv] -= state.rhs[*in][nv];
        state.um[*in][nv] -= state.rhs[*in][nv];
      }}

  /* -----------------------------------------------------------
       re-compute the full rhs using the total (hyp+par) rhs
     ----------------------------------------------------------- */

      RightHandSide (&state, Dts, indx.beg, indx.end, dt2, grid);

      if (g_dir == IDIR){
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] = state.rhs[*in][nv];
          UH[k][j][i][nv] = d->Uc[k][j][i][nv] + state.rhs[*in][nv];
        }}
      }else{
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] += state.rhs[*in][nv];
          UH[k][j][i][nv] += state.rhs[*in][nv];
        }}
      }
#endif
    } /* -- end loop on transverse directions -- */

  } /* -- end loop on dimensions -- */

/* -------------------------------------------------------
     2a. Advance staggered magnetic fields by dt/2 
   ------------------------------------------------------- */

  #ifdef STAGGERED_MHD     
   CT_Update (d, d->Vs, 0.5*g_dt, grid);
   CT_AverageMagneticField (d->Vs, UH, grid);
  #endif

/* ----------------------------------------------------------
     2b. Convert cell and time centered values to primitive. 
         UH is transformed from conservative to primitive 
         variables for efficiency purposes.
   ---------------------------------------------------------- */

  g_dir = IDIR;
  SetIndexes (&indx, grid);
  TRANSVERSE_LOOP(indx,in,i,j,k){
    g_i = i; g_j = j; g_k = k;
    errp = ConsToPrim(UH[k][j], state.v, indx.beg, indx.end, d->flag[k][j]);
    for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
      VAR_LOOP(nv) d->Vc[nv][k][j][i] = state.v[*in][nv];
    }
    #if THERMAL_CONDUCTION == EXPLICIT
     for ((*in) = indx.beg; (*in) <= indx.end; (*in)++){
       T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
     }
    #endif
  }
  #if (RESISTIVITY == EXPLICIT) && (defined STAGGERED_MHD)
   GetCurrent (d, -1, grid);
  #endif
  
/* ---------------------------------------------------- 
     3. Final conservative update
   ---------------------------------------------------- */

  g_intStage = 2;
  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);
    ResetState (d, &state, grid);
    #if (RESISTIVITY == EXPLICIT) && !(defined STAGGERED_MHD)
     GetCurrent(d, g_dir, grid);
    #endif
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      g_i = i;  g_j = j;  g_k = k;

    /* ---------------------------------------------------------
        Convert conservative corner-coupled states to primitive 
       --------------------------------------------------------- */
       
      state.up = UP[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uL = state.up;
      state.um = UM[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uR = state.um + 1;
      
      for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++){
        for (nv = NVAR; nv--; ){
          state.up[*in][nv] += dU[k][j][i][nv];
          state.um[*in][nv] += dU[k][j][i][nv];
        }
        flagp[*in] = flagm[*in] = d->flag[k][j][i];
      }

    /* -------------------------------------------------------
        compute time centered, cell centered state. 
        Useful for source terms like gravity, 
        curvilinear terms and Powell's 8wave 
       ------------------------------------------------------- */

      for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++) {
        NVAR_LOOP(nv) state.vh[*in][nv] = d->Vc[nv][k][j][i];
        state.flag[*in] = d->flag[k][j][i];
      }

      #ifdef STAGGERED_MHD
      for ((*in) = indx.beg - 2; (*in) <= indx.end + 1; (*in)++){
        state.uL[*in][BXn] = state.uR[*in][BXn] = d->Vs[BXs + g_dir][k][j][i];
        state.bn[*in] = state.uL[*in][BXn];
      }
      #endif

      errm = ConsToPrim (state.um, state.vm, indx.beg - 1, indx.end + 1, flagm);
      errp = ConsToPrim (state.up, state.vp, indx.beg - 1, indx.end + 1, flagp);
/*
      if (errm || errp){
        WARNING(print ("! Corner coupled states not physical: reverting to 1st order\n");)
        for (in = indx.beg - 1; in <= indx.end + 1; in++){
          if (   (flagm[*in] & FLAG_CONS2PRIM_FAIL)
              || (flagm[in]  & FLAG_CONS2PRIM_FAIL) ) {
            for (nv = 0; nv < NVAR; nv++) {
              state.v[in][nv] = d->Vc[nv][*k][*j][*i];
            }
            PrimToCons(state.v, state.u, in, in);
            for (nv = 0; nv < NVAR; nv++) {
              state.vm[in][nv] = state.vp[in][nv] = state.vh[in][nv] = state.v[in][nv];
              state.um[in][nv] = state.up[in][nv] = state.uh[in][nv] = state.u[in][nv];
            }
          }
        }
      }
*/

   /* ------------------------------------------
           compute flux & righ-hand-side
      ------------------------------------------ */

      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
       CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
      #endif
      #if (PARABOLIC_FLUX & EXPLICIT)

     /* -------------------------------------------------------- 
         since integration in the transverse directions extends 
         one point below and one point above the computational 
         box (when using CT), we avoid computing parabolic 
         operators in one row of boundary zones below and above 
         since the 3-point wide stencil may not be available.
        -------------------------------------------------------- */

       errp = 1;
       #ifdef STAGGERED_MHD
       D_EXPAND(                                                        ,
         errp = (*(indx.pt1) < indx.t1_end) && (*(indx.pt1) > indx.t1_beg);     ,
         errp = errp && (*(indx.pt2) < indx.t2_end) && (*(indx.pt2) > indx.t2_beg);)
       #endif
       if (errp) {
         ParabolicFlux (d->Vc, d->J, T, &state, dcoeff, indx.beg - 1, indx.end, grid);
         inv_dl  = GetInverse_dl(grid);
         for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
           inv_dtp = 0.0;
           #if VISCOSITY == EXPLICIT
            inv_dtp = MAX(inv_dtp, dcoeff[*in][MX1]);
           #endif
           #if RESISTIVITY == EXPLICIT
            EXPAND(inv_dtp = MAX(inv_dtp, dcoeff[*in][BX1]);  ,
                   inv_dtp = MAX(inv_dtp, dcoeff[*in][BX2]);  ,
                   inv_dtp = MAX(inv_dtp, dcoeff[*in][BX3]);)
           #endif
           #if THERMAL_CONDUCTION == EXPLICIT
            inv_dtp = MAX(inv_dtp, dcoeff[*in][ENG]);
           #endif
           inv_dtp *= inv_dl[*in]*inv_dl[*in];

           Dts->inv_dtp = MAX(Dts->inv_dtp, inv_dtp);
         }
       }
      #endif
      #if UPDATE_VECTOR_POTENTIAL == YES
       VectorPotentialUpdate (d, NULL, &state, grid);
      #endif
      #ifdef SHEARINGBOX
       SB_SaveFluxes(&state, grid);
      #endif

      RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);

      for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        NVAR_LOOP(nv) d->Uc[k][j][i][nv] += state.rhs[*in][nv];
      }           
    }
  }

  #ifdef SHEARINGBOX
   SB_CorrectFluxes(d->Uc, g_time+0.5*g_dt, g_dt, grid);
  #endif

  #ifdef STAGGERED_MHD
   CT_Update (d, Bs0, g_dt, grid);
   CT_AverageMagneticField (d->Vs, d->Uc, grid);
  #endif

/* ----------------------------------------------
           convert to primitive 
   ---------------------------------------------- */

  #ifdef FARGO
   FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif 

  ConsToPrim3D(d->Uc, d->Vc, d->flag, box);

  #ifdef FARGO
   FARGO_AddVelocity (d,grid); 
  #endif

  return(0); /* -- step has been achieved, return success -- */
}

#if CTU_MHD_SOURCE == YES
/* ********************************************************************* */
void CTU_CT_Source (double **v, double **up, double **um,
                    double *dtdV, int beg, int end, Grid *grid)
/*!
 * Add source terms to conservative left and right states obtained from 
 * the primitive form  of the equations. The source terms are:
 *
 * - m  += dt/2 *  B  * dbx/dx
 * - Bt += dt/2 * vt  * dbx/dx   (t = transverse component)
 * - E  += dt/2 * v*B * dbx/dx
 *
 * These terms are NOT accounted for when the primitive form of the 
 * equations is used (see Gardiner & Stone JCP (2005), Crockett et al. 
 * JCP(2005)). This is true for both the Charactheristic Tracing AND the 
 * primitive Hancock scheme when the constrained transport is used, since 
 * the  resulting system is 7x7. To better understand this, you can 
 * consider the stationary solution rho = p = 1, v = 0  and Bx = x, 
 * By = -y. If these terms were not included the code would generate 
 * spurious velocities.
 *
 *********************************************************************** */
{
  int    i;
  double   scrh, *dx, *A;
  static double *db;

  if (db == NULL) db = ARRAY_1D(NMAX_POINT, double);

/* ----------------------------------------
              comput db/dx
   ---------------------------------------- */

  #if GEOMETRY == CARTESIAN
   for (i = beg; i <= end; i++){
     db[i] = dtdV[i]*(up[i][BXn] - um[i][BXn]); 
   }
  #elif GEOMETRY == CYLINDRICAL
   if (g_dir == IDIR){
     A = grid[IDIR].A;
     for (i = beg; i <= end; i++){
       db[i] = dtdV[i]*(up[i][BXn]*A[i] - um[i][BXn]*A[i - 1]);
     }
   }else{
     for (i = beg; i <= end; i++){
       db[i] = dtdV[i]*(up[i][BXn] - um[i][BXn]); 
     }
   }
  #else
   print1 (" ! CTU-MHD does not work in this geometry\n");
   QUIT_PLUTO(1);
  #endif

/* --------------------------------------------
         Add source terms
   -------------------------------------------- */

  for (i = beg; i <= end; i++){
    EXPAND( up[i][MX1] += v[i][BX1]*db[i];
            um[i][MX1] += v[i][BX1]*db[i];   ,
            up[i][MX2] += v[i][BX2]*db[i];
            um[i][MX2] += v[i][BX2]*db[i];   ,
            up[i][MX3] += v[i][BX3]*db[i];
            um[i][MX3] += v[i][BX3]*db[i]; ) 

    EXPAND(                            ;   ,
            up[i][BXt] += v[i][VXt]*db[i]; 
            um[i][BXt] += v[i][VXt]*db[i];   ,
            up[i][BXb] += v[i][VXb]*db[i]; 
            um[i][BXb] += v[i][VXb]*db[i];)

    #if HAVE_ENERGY 
     scrh = EXPAND(   v[i][VX1]*v[i][BX1]  , 
                    + v[i][VX2]*v[i][BX2]  , 
                    + v[i][VX3]*v[i][BX3]);
     up[i][ENG] += scrh*db[i];
     um[i][ENG] += scrh*db[i];
    #endif
  }

}
#endif
