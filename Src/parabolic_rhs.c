/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute right hand side for diffusion terms.

  The ParabolicRHS function computes the right hand side, in 
  conservative or divergence form, of the parabolic (diffusion)
  operators only:
  \f[
      \pd{q}{t} = - \nabla\cdot\Big(D\nabla q\Big)
  \f]  
  Here \c q is a generic cell-centered quantity (like momentum, magnetic
  field or total energy) and \f$\vec{F} = D\nabla q\f$ is the corresponding
  flux computed elsewhere. 
  Contributions may simultaneously come from viscosity, magnetic resistivity 
  and thermal conduction for which the flux functions are computed in the 
  function ViscousFlux(), ResistiveFlux() or TC_Flux().
  
  This function is intended for operator split algorithms (STS, RKC) 
  in order to solve for the parabolic part of the equations only.
  
  This function also computes the inverse diffusion time step by adding, 
  for each diffusion equation, contributions coming from different directions:
  \f[ \Delta t_{p}^{-1} = 
      \max\left[\frac{\eta_{i+\HALF} + \eta_{i-\HALF}}{2\Delta x^2} +
                \frac{\eta_{j+\HALF} + \eta_{j-\HALF}}{2\Delta y^2} +
                \frac{\eta_{k+\HALF} + \eta_{k-\HALF}}{2\Delta z^2}\right]
  \f]
  where \f$\eta_{x,y,z}\f$ are the diffusion coefficients available at 
  cell interfaces in the three directions, respectively, and the maximum
  is taken over the local processor grid.
  
  \date    Sept 1, 2014
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define ADD_RESISTIVITY (RESISTIVITY == SUPER_TIME_STEPPING || \
                         RESISTIVITY == RK_CHEBYSHEV)

#define ADD_VISCOSITY (VISCOSITY == SUPER_TIME_STEPPING || \
                       VISCOSITY == RK_CHEBYSHEV)

#define ADD_TC (THERMAL_CONDUCTION == SUPER_TIME_STEPPING || \
                THERMAL_CONDUCTION == RK_CHEBYSHEV)

/* ********************************************************************* */
double ParabolicRHS (const Data *d, Data_Arr dU, double dt, Grid *grid)
/*! 
 * \param [in]  V    3D array containing primitive variables
 * \param [out] dU   3D array containing the conservative right hand sides 
 * \param [in]  dt   the time step
 * \param [in]  grid a pointer to an array of Grid structure
 *
 * \return On output it returns the maximum diffusion coefficients 
 *         among all dissipative term over the local processor grid.
 *********************************************************************** */
{
  int i, j, k, nv;
  double r, th, r_1, s_1;
  double dtdl, dtdV, dt_rdr;
  double *inv_dl, inv_dl2;
  double *rp, *A, *du, **vh, *C, **tc_flx;
  double dvp[NVAR], dvm[NVAR], dvl;
  double ****V;
  static double **res_flx; 
  static double **dcoeff, ***T;
  static double **ViF, **ViS;
  
  static Data_Arr C_dtp;
  static State_1D state;
  double max_inv_dtp[3],inv_dtp;

  i = j = k = 0;
  
  if (C_dtp == NULL) {
    C_dtp = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    #if ADD_TC 
     MakeState (&state);
     #ifdef CH_SPACEDIM 
      i = j = k = 1;
      D_EXPAND(i = NMAX_POINT;, j = NMAX_POINT;, k = NMAX_POINT;)
      T = ARRAY_3D(k, j, i, double); 
     #else
      T = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     #endif
    #endif

    ViF = ARRAY_2D(NMAX_POINT, NVAR, double);
    ViS = ARRAY_2D(NMAX_POINT, NVAR, double); 
    res_flx = ARRAY_2D(NMAX_POINT, NVAR, double);
    dcoeff  = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  V = d->Vc;  /* -- set a pointer to the primitive vars array -- */
  
  #if ADD_RESISTIVITY && (defined STAGGERED_MHD)
   GetCurrent (d, -1, grid);
  #endif

/* ------------------------------------------------
    compute the temperature array if TC is needed
   ------------------------------------------------ */

  #if ADD_TC
   TOT_LOOP(k,j,i) T[k][j][i] = V[PRS][k][j][i]/V[RHO][k][j][i];
   vh     = state.vh;
   tc_flx = state.tc_flux;
  #endif

max_inv_dtp[0] = max_inv_dtp[1] = max_inv_dtp[2] = 0.0;

/* ----------------------------------------------------------------------   
         L O O P     O N    X1   D I R E C T I O N  ( I D I R )
   ---------------------------------------------------------------------- */

  g_dir = IDIR;
  A   = grid[g_dir].A;
  rp  = grid[g_dir].xr;
  #if ADD_RESISTIVITY  && !(defined STAGGERED_MHD)
   GetCurrent(d, g_dir, grid);
  #endif
  KDOM_LOOP (k){ 
  JDOM_LOOP (j){ 

    g_j = j; g_k = k;

  /* -------------------------------------------------------------------
      Compute viscous, resistive and thermal cond. fluxes in the x1 dir
     ------------------------------------------------------------------- */
     
    #if ADD_VISCOSITY
     ViscousFlux (V, ViF, ViS, dcoeff, IBEG-1, IEND, grid);
    #endif

    #if ADD_RESISTIVITY
     ResistiveFlux (V, d->J, res_flx, dcoeff, IBEG - 1, IEND, grid);
    #endif

    #if ADD_TC
     ITOT_LOOP (i) for (nv = NVAR; nv--; ) vh[i][nv] = V[nv][k][j][i];  

     for (nv = NVAR; nv--; ) dvm[nv] = vh[1][nv] - vh[0][nv];
     for(i=1; i<NX1_TOT-1; i++){
     for (nv = NVAR; nv--; ) {
       dvp[nv] = vh[i+1][nv] - vh[i][nv];
       dvl     = VAN_LEER(dvp[nv], dvm[nv]); 
       state.vp[i][nv] = vh[i][nv] + 0.5*dvl;
       state.vm[i][nv] = vh[i][nv] - 0.5*dvl;
       dvm[nv] = dvp[nv];
     }}
     TC_Flux (T, &state, dcoeff, IBEG - 1, IEND, grid);
    #endif

  /* ---------------------------
      compute inverse time-step
     --------------------------- */
     
    if (g_intStage == 1){ 
      inv_dl = GetInverse_dl(grid);
      IDOM_LOOP(i){
        C       = C_dtp[k][j][i];
        inv_dl2 = 0.5*inv_dl[i]*inv_dl[i];
        #if ADD_VISCOSITY
         C[MX1] = (dcoeff[i-1][MX1] + dcoeff[i][MX1])*inv_dl2;
        #endif
        #if ADD_RESISTIVITY
         EXPAND(C[BX1] = (dcoeff[i-1][BX1] + dcoeff[i][BX1])*inv_dl2;  ,
                C[BX2] = (dcoeff[i-1][BX2] + dcoeff[i][BX2])*inv_dl2;  ,
                C[BX3] = (dcoeff[i-1][BX3] + dcoeff[i][BX3])*inv_dl2;)
        #endif
        #if ADD_TC
         C[ENG] = (dcoeff[i-1][ENG] + dcoeff[i][ENG])*inv_dl2;
         inv_dtp = dcoeff[i][ENG]*inv_dl[i]*inv_dl[i];
         max_inv_dtp[g_dir] = MAX(max_inv_dtp[g_dir], inv_dtp);
        #endif
      }
    }

  /* ---------------------------
          Main X1-sweep
     --------------------------- */

    IDOM_LOOP (i){  

      du = dU[k][j][i];
      r  = grid[IDIR].x[i];
      dtdV = dt/grid[IDIR].dV[i];
      dtdl = dt/grid[IDIR].dx[i];

      #if HAVE_ENERGY
      du[ENG] = 0.0; /* energy contribution may come from any of the
                        dissipative terms. It's better to initialize it
                        to zero and then add contribution separately. */
      #endif

    /* -------------------------------------------
        VISCOSITY: build rhs in the x1 direction
       ------------------------------------------- */

      #if ADD_VISCOSITY
       r_1 = 1.0/grid[IDIR].x[i];
       #if GEOMETRY == CARTESIAN
        EXPAND(du[MX1] = (ViF[i][MX1] - ViF[i-1][MX1])*dtdl + dt*ViS[i][MX1];,
               du[MX2] = (ViF[i][MX2] - ViF[i-1][MX2])*dtdl + dt*ViS[i][MX2];,   
               du[MX3] = (ViF[i][MX3] - ViF[i-1][MX3])*dtdl + dt*ViS[i][MX3];)
        #if HAVE_ENERGY
         du[ENG] += (ViF[i][ENG] - ViF[i - 1][ENG])*dtdl;
        #endif
       #elif GEOMETRY == CYLINDRICAL
        EXPAND(du[MX1] = (A[i]*ViF[i][MX1] - A[i-1]*ViF[i-1][MX1])*dtdV 
                        + dt*ViS[i][MX1]; ,
               du[MX2] = (A[i]*ViF[i][MX2] - A[i-1]*ViF[i-1][MX2])*dtdV 
                        + dt*ViS[i][MX2]; ,   
               du[MX3] = (A[i]*A[i]*ViF[i][MX3] - A[i-1]*A[i-1]*ViF[i-1][MX3])*r_1*dtdV;)
        #if HAVE_ENERGY
         du[ENG] += (A[i]*ViF[i][ENG] - A[i-1]*ViF[i-1][ENG])*dtdV;
        #endif
       #elif GEOMETRY == POLAR
        EXPAND(du[MX1] = (A[i]*ViF[i][MX1] - A[i-1]*ViF[i-1][MX1])*dtdV 
                         + dt*ViS[i][MX1]; ,
               du[MX2] = (A[i]*A[i]*ViF[i][MX2] - A[i-1]*A[i-1]*ViF[i-1][MX2])*r_1*dtdV; ,   
               du[MX3] = (A[i]*ViF[i][MX3] - A[i-1]*ViF[i-1][MX3])*dtdV 
                         + dt*ViS[i][MX3];)
        #if HAVE_ENERGY
         du[ENG] += (A[i]*ViF[i][ENG] - A[i-1]*ViF[i-1][ENG])*dtdV;
        #endif
       #elif GEOMETRY == SPHERICAL
        EXPAND(du[MX1] = (A[i]*ViF[i][MX1] - A[i-1]*ViF[i-1][MX1])*dtdV 
                         + dt*ViS[i][MX1];,
               du[MX2] = (A[i]*ViF[i][MX2] - A[i-1]*ViF[i-1][MX2])*dtdV 
                         + dt*ViS[i][MX2];,   
               du[MX3] = (rp[i]*A[i]*ViF[i][MX3] - rp[i-1]*A[i-1]*ViF[i-1][MX3])*r_1*dtdV;)
        #if HAVE_ENERGY
         du[ENG] += (A[i]*ViF[i][ENG] - A[i-1]*ViF[i-1][ENG])*dtdV;
        #endif
       #endif
      #endif /* -- VISCOSITY -- */

    /* ----------------------------------------------
        RESISTIVITY: build rhs in the x1 direction
       ---------------------------------------------- */

      #if ADD_RESISTIVITY
       #if GEOMETRY == CARTESIAN    /* -- x coordinate -- */
        EXPAND(
          du[BX1] = 0.0;                                        ,
          du[BX2] = -(res_flx[i][BX2] - res_flx[i-1][BX2])*dtdl;  ,
          du[BX3] = -(res_flx[i][BX3] - res_flx[i-1][BX3])*dtdl; )
        #if EOS != ISOTHERMAL
         du[ENG] += -(res_flx[i][ENG] - res_flx[i-1][ENG])*dtdl; 
        #endif 

       #elif GEOMETRY == CYLINDRICAL  /* -- r coordinate -- */
        EXPAND(
          du[BX1] = 0.0;                                                   ,
          du[BX2] = -(A[i]*res_flx[i][BX2] - A[i-1]*res_flx[i-1][BX2])*dtdV; ,
          du[BX3] = -(res_flx[i][BX3] - res_flx[i-1][BX3])*dtdl;)
        #if EOS != ISOTHERMAL
         du[ENG] += -(A[i]*res_flx[i][ENG] - A[i-1]*res_flx[i-1][ENG])*dtdV; 
        #endif 

       #elif GEOMETRY == POLAR    /* -- r coordinate -- */
        EXPAND(
          du[BX1] = 0.0;                                                      ,
          du[BX2] = -(res_flx[i][BX2] - res_flx[i-1][BX2])*dtdl;                ,
          du[BX3] = -(A[i]*res_flx[i][BX3] - A[i-1]*res_flx[i-1][BX3])*dtdV; )
        #if EOS != ISOTHERMAL
         du[ENG] += -(A[i]*res_flx[i][ENG] - A[i-1]*res_flx[i-1][ENG])*dtdV; 
        #endif 

       #elif GEOMETRY == SPHERICAL  /* -- r coordinate -- */
        dt_rdr = dtdl/r;
        EXPAND(
          du[BX1] = 0.0;                                                      ,
          du[BX2] = -(rp[i]*res_flx[i][BX2] - rp[i-1]*res_flx[i-1][BX2])*dt_rdr; ,
          du[BX3] = -(rp[i]*res_flx[i][BX3] - rp[i-1]*res_flx[i-1][BX3])*dt_rdr; )
        #if EOS != ISOTHERMAL
         du[ENG] += -(A[i]*res_flx[i][ENG] - A[i-1]*res_flx[i-1][ENG])*dtdV; 
        #endif 
       #endif
      #endif /* -- RESISTIVITY -- */

    /* ---------------------------------------------------
        THERMAL_CONDUCTION: build rhs in the x1 direction
       --------------------------------------------------- */

      #if ADD_TC
       #if GEOMETRY == CARTESIAN 
        du[ENG] += (tc_flx[i][ENG] - tc_flx[i-1][ENG])*dtdl;
       #else
        du[ENG] += (A[i]*tc_flx[i][ENG] - A[i-1]*tc_flx[i-1][ENG])*dtdV;
       #endif
      #endif
    }
  }}

/* ----------------------------------------------------------------------   
         L O O P     O N    X2   D I R E C T I O N  ( J D I R )
   ---------------------------------------------------------------------- */

#if DIMENSIONS > 1
  g_dir = JDIR;
  A     = grid[g_dir].A;
  #if ADD_RESISTIVITY  && !(defined STAGGERED_MHD)
   GetCurrent(d, g_dir, grid);
  #endif
  KDOM_LOOP (k){
  IDOM_LOOP (i){

    g_i = i; g_k = k;

  /* -------------------------------------------------------------------
      Compute viscous, resistive and thermal cond. fluxes in the x2 dir
     ------------------------------------------------------------------- */

    #if ADD_VISCOSITY
     ViscousFlux (V, ViF, ViS, dcoeff, JBEG-1, JEND, grid);
    #endif
    #if ADD_RESISTIVITY
     ResistiveFlux (V, d->J, res_flx, dcoeff, JBEG - 1, JEND, grid);
    #endif

    #if ADD_TC
     JTOT_LOOP (j) for (nv = NVAR; nv--; ) vh[j][nv] = V[nv][k][j][i];  

     for (nv = NVAR; nv--; ) dvm[nv] = vh[1][nv] - vh[0][nv];
     for(j = 1; j < NX2_TOT-1; j++){
     for (nv = NVAR; nv--; ) {
       dvp[nv] = vh[j+1][nv] - vh[j][nv];
       dvl     = VAN_LEER(dvp[nv], dvm[nv]);
       state.vp[j][nv] = vh[j][nv] + 0.5*dvl;
       state.vm[j][nv] = vh[j][nv] - 0.5*dvl;
       dvm[nv] = dvp[nv];
     }}
     TC_Flux (T, &state, dcoeff, JBEG - 1, JEND, grid);
    #endif

  /* ---------------------------
      compute inverse time-step
     --------------------------- */

    if (g_intStage == 1){
      inv_dl = GetInverse_dl(grid);
      JDOM_LOOP(j){
        C       = C_dtp[k][j][i];
        inv_dl2 = 0.5*inv_dl[j]*inv_dl[j];
        #if ADD_VISCOSITY
         C[MX1] += (dcoeff[j-1][MX1] + dcoeff[j][MX1])*inv_dl2;
        #endif
        #if ADD_RESISTIVITY
         EXPAND(C[BX1] += (dcoeff[j-1][BX1] + dcoeff[j][BX1])*inv_dl2;  ,
                C[BX2] += (dcoeff[j-1][BX2] + dcoeff[j][BX2])*inv_dl2;  ,
                C[BX3] += (dcoeff[j-1][BX3] + dcoeff[j][BX3])*inv_dl2;)
        #endif
        #if ADD_TC
         C[ENG] += (dcoeff[j-1][ENG] + dcoeff[j][ENG])*inv_dl2;

         inv_dtp = dcoeff[j][ENG]*inv_dl[j]*inv_dl[j];
         max_inv_dtp[g_dir] = MAX(max_inv_dtp[g_dir], inv_dtp);
        #endif
      }
    }

  /* ---------------------------
          Main X2-sweep
     --------------------------- */

    r = grid[IDIR].x[i];
    JDOM_LOOP (j){   

      du   = dU[k][j][i];
      dtdV = dt/grid[JDIR].dV[j]; 
      dtdl = dt/grid[JDIR].dx[j];

      #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
       dtdl /= r;
       dtdV /= r;
      #endif

    /* ------------------------------------------
        VISCOSITY: build rhs in the x2 direction
       ------------------------------------------ */

      #if ADD_VISCOSITY
       s_1 = 1.0/sin(grid[JDIR].x[j]);
       #if GEOMETRY != SPHERICAL
        EXPAND(du[MX1] += (ViF[j][MX1] - ViF[j-1][MX1])*dtdl + dt*ViS[j][MX1];,
               du[MX2] += (ViF[j][MX2] - ViF[j-1][MX2])*dtdl + dt*ViS[j][MX2];,   
               du[MX3] += (ViF[j][MX3] - ViF[j-1][MX3])*dtdl + dt*ViS[j][MX3];)
        #if HAVE_ENERGY
         du[ENG] += (ViF[j][ENG] - ViF[j-1][ENG])*dtdl;
        #endif
       #elif GEOMETRY == SPHERICAL
        EXPAND(du[MX1] += (A[j]*ViF[j][MX1] - A[j-1]*ViF[j-1][MX1])*dtdV 
                          + dt*ViS[j][MX1];,
               du[MX2] += (A[j]*ViF[j][MX2] - A[j-1]*ViF[j-1][MX2])*dtdV 
                          + dt*ViS[j][MX2];,   
               du[MX3] += (A[j]*A[j]*ViF[j][MX3] - A[j-1]*A[j-1]*ViF[j-1][MX3])*s_1*dtdV;)
        #if HAVE_ENERGY
         du[ENG] += (A[j]*ViF[j][ENG] - A[j-1]*ViF[j-1][ENG])*dtdV;
        #endif
       #endif 
      #endif/* -- VISCOSITY -- */

    /* ----------------------------------------------
        RESISTIVITY: build rhs in the x2 direction
       ---------------------------------------------- */

      #if ADD_RESISTIVITY
       #if GEOMETRY != SPHERICAL  /* -- y coordinate -- */
        EXPAND(
          du[BX1] -= (res_flx[j][BX1] - res_flx[j-1][BX1])*dtdl;  ,
                                                              ,
          du[BX3] -= (res_flx[j][BX3] - res_flx[j-1][BX3])*dtdl; )
        #if HAVE_ENERGY
         du[ENG] -= (res_flx[j][ENG] - res_flx[j-1][ENG])*dtdl; 
        #endif 
       #elif GEOMETRY == SPHERICAL  /* -- theta coordinate -- */
        EXPAND(
          du[BX1] -= (A[j]*res_flx[j][BX1] - A[j-1]*res_flx[j-1][BX1])*dtdV;  ,
                                                                          ,
          du[BX3] -= (res_flx[j][BX3] - res_flx[j-1][BX3])*dtdl;)
        #if HAVE_ENERGY
         du[ENG] -= (A[j]*res_flx[j][ENG] - A[j-1]*res_flx[j-1][ENG])*dtdV; 
        #endif 
       #endif
      #endif /* -- RESISTIVE MHD -- */

    /* ---------------------------------------------------
        THERMAL_CONDUCTION: build rhs in the x2 direction
       --------------------------------------------------- */

      #if ADD_TC
       #if GEOMETRY != SPHERICAL
        du[ENG] += (tc_flx[j][ENG] - tc_flx[j-1][ENG])*dtdl;
       #else
        du[ENG] += (A[j]*tc_flx[j][ENG] - A[j-1]*tc_flx[j-1][ENG])*dtdV;
       #endif
      #endif
    }

  }}
#endif

/* ----------------------------------------------------------------------   
         L O O P     O N    X3   D I R E C T I O N  ( K D I R )
   ---------------------------------------------------------------------- */

#if DIMENSIONS == 3
  g_dir = KDIR;
  A   = grid[g_dir].A;
  #if ADD_RESISTIVITY && !(defined STAGGERED_MHD)
   GetCurrent(d, g_dir, grid);
  #endif
  JDOM_LOOP (j){
  IDOM_LOOP (i){

    g_i = i; g_j = j;

  /* -------------------------------------------------------------------
      Compute viscous, resistive and thermal cond. fluxes in the x3 dir
     ------------------------------------------------------------------- */

    #if ADD_VISCOSITY
     ViscousFlux (V, ViF, ViS, dcoeff, KBEG-1, KEND, grid);
    #endif
    #if ADD_RESISTIVITY
     ResistiveFlux (V, d->J, res_flx, dcoeff, KBEG - 1, KEND, grid);
    #endif

    #if ADD_TC
     KTOT_LOOP (k) for (nv = NVAR; nv--; ) vh[k][nv] = V[nv][k][j][i];  

     for (nv = NVAR; nv--; ) dvm[nv] = vh[1][nv] - vh[0][nv];
     for(k=1; k<NX3_TOT-1; k++){
     for (nv = NVAR; nv--; ) {
       dvp[nv] = vh[k+1][nv]- vh[k][nv];
       dvl     = VAN_LEER(dvp[nv], dvm[nv]);
       state.vp[k][nv] = vh[k][nv] + 0.5*dvl;
       state.vm[k][nv] = vh[k][nv] - 0.5*dvl;
       dvm[nv] = dvp[nv];
     }}
     TC_Flux (T, &state, dcoeff, KBEG - 1, KEND, grid);
    #endif

  /* ---------------------------
      compute inverse time-step
     --------------------------- */

    if (g_intStage == 1){
      inv_dl = GetInverse_dl(grid);
      KDOM_LOOP(k){
        C       = C_dtp[k][j][i];
        inv_dl2 = 0.5*inv_dl[k]*inv_dl[k];
        #if ADD_VISCOSITY
         C[MX1] += (dcoeff[k-1][MX1] + dcoeff[k][MX1])*inv_dl2;
        #endif
        #if ADD_RESISTIVITY
         EXPAND(C[BX1] += (dcoeff[k-1][BX1] + dcoeff[k][BX1])*inv_dl2;  ,
                C[BX2] += (dcoeff[k-1][BX2] + dcoeff[k][BX2])*inv_dl2;  ,
                C[BX3] += (dcoeff[k-1][BX3] + dcoeff[k][BX3])*inv_dl2;)
        #endif
        #if ADD_TC
         C[ENG] += (dcoeff[k-1][ENG] + dcoeff[k][ENG])*inv_dl2;

         inv_dtp = dcoeff[k][ENG]*inv_dl[k]*inv_dl[k];
         max_inv_dtp[g_dir] = MAX(max_inv_dtp[g_dir], inv_dtp);

        #endif
      }
    }

  /* ---------------------------
          Main X3-sweep
     --------------------------- */

    r  = grid[IDIR].x[i];
    th = grid[JDIR].x[j];
    KDOM_LOOP (k){  

      du   = dU[k][j][i];
      dtdV = dt/grid[KDIR].dV[k];
      dtdl = dt/grid[KDIR].dx[k];
      #if GEOMETRY == SPHERICAL
       dtdl /= r*sin(th);
      #endif

    /* ------------------------------------------
        VISCOSITY: build rhs in the x3 direction
       ------------------------------------------ */

      #if ADD_VISCOSITY
       du[MX1] += (ViF[k][MX1] - ViF[k-1][MX1])*dtdl + dt*ViS[k][MX1];
       du[MX2] += (ViF[k][MX2] - ViF[k-1][MX2])*dtdl + dt*ViS[k][MX2];
       du[MX3] += (ViF[k][MX3] - ViF[k-1][MX3])*dtdl + dt*ViS[k][MX3];
       #if HAVE_ENERGY
        du[ENG] += (ViF[k][ENG] - ViF[k-1][ENG])*dtdl;
       #endif
      #endif/* -- VISCOSITY -- */

    /* ----------------------------------------------
        RESISTIVITY: build rhs in the x3 direction
       ---------------------------------------------- */

      #if ADD_RESISTIVITY
       du[BX1] -= (res_flx[k][BX1] - res_flx[k-1][BX1])*dtdl;
       du[BX2] -= (res_flx[k][BX2] - res_flx[k-1][BX2])*dtdl;
       #if HAVE_ENERGY
        du[ENG] -= (res_flx[k][ENG] - res_flx[k-1][ENG])*dtdl; 
       #endif 
      #endif /* -- RESISTIVITY -- */
   
    /* ---------------------------------------------------
        THERMAL_CONDUCTION: build rhs in the x3 direction
       --------------------------------------------------- */

      #if ADD_TC
       du[ENG] += (tc_flx[k][ENG] - tc_flx[k-1][ENG])*dtdl;
      #endif
    }
  }}
#endif

/*  OLD VERSION 
D_EXPAND(th = max_inv_dtp[IDIR];            ,
         th = MAX(th, max_inv_dtp[JDIR]);  ,
         th = MAX(th, max_inv_dtp[KDIR]);)
return th; 
*/

/* --------------------------------------------------------------
     take the maximum of inverse dt over domain and zero 
     right hand side in the internal boundary zones.
   -------------------------------------------------------------- */

  th = 0.0;
  DOM_LOOP(k,j,i){
    #if ADD_VISCOSITY
     th = MAX(th, C_dtp[k][j][i][MX1]);   
    #endif
    #if ADD_RESISTIVITY
     EXPAND(th = MAX(th, C_dtp[k][j][i][BX1]);   ,  
            th = MAX(th, C_dtp[k][j][i][BX2]);   ,
            th = MAX(th, C_dtp[k][j][i][BX3]);)
    #endif
    #if ADD_TC
     th = MAX(th, C_dtp[k][j][i][ENG]);
    #endif
    #if INTERNAL_BOUNDARY == YES
     if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) {
       for (nv = NVAR; nv--;  ) dU[k][j][i][nv] = 0.0;
     }
    #endif
  }
  return th;
}
#undef ADD_VISCOSITY
#undef ADD_RESISTIVITY
#undef ADD_TC
