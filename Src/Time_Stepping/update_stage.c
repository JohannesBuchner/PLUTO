/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief   Single stage integration for RK time stepping.

  Advance the equations in conservative form by taking a single stage 
  in the form 
  \f[
      U \quad\Longrightarrow \quad  U + \Delta t R(V)
  \f]
  where \c U is a 3D array of conservative variables, \c V is a 3D array
  of primitive variables, \c R(V) is the right hand side containing 
  flux differences and source terms.
  Note that \c U and \c V may \e not necessarily be the map of 
  each other, i.e., \c U is \e not \c U(V).
  The right hand side can contain contributions from 
   
    - the direction set by the global variable ::g_dir, 
      when DIMENSIONAL_SPLITTING == YES;
    - all directions when DIMENSIONAL_SPLITTING == NO;
   
  When the integrator stage is the first one (predictor), this function 
  also computes the maximum of inverse time steps for hyperbolic and 
  parabolic terms (if the latters are included explicitly).
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
           T. Matsakos
  \date   March 02, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static void SaveAMRFluxes (const State_1D *, double **, int, int, Grid *);
static intList TimeStepIndexList();

/* ********************************************************************* */
void UpdateStage(const Data *d, Data_Arr UU, double **aflux,
                 Riemann_Solver *Riemann, double dt, Time_Step *Dts, 
                 Grid *grid)
/*!
 * 
 * \param [in,out]  d        pointer to PLUTO Data structure
 * \param [in,out]  UU       data array containing conservative variables
 *                           at the previous time step to be updated
 * \param [out]     aflux    interface fluxes needed for refluxing operations 
 *                           (only with AMR)
 * \param [in]      Riemann  pointer to a Riemann solver function
 * \param [in]      dt       the time step for the current update step
 * \param [in,out]  Dts      pointer to time step structure
 * \param [in]      grid     pointer to array of Grid structures
 *********************************************************************** */
{
  int  i, j, k;
  int  nv, dir, beg_dir, end_dir;
  int  *ip;
  double *inv_dl, dl2;
  static double ***T, ***C_dt[NVAR], **dcoeff;
  static State_1D state;
  Index indx;
  intList cdt_list;

  #if DIMENSIONAL_SPLITTING == YES
   beg_dir = end_dir = g_dir;
  #else
   beg_dir = 0;
   end_dir = DIMENSIONS-1;
  #endif

  cdt_list = TimeStepIndexList();

/* --------------------------------------------------------------
   1. Allocate memory
   -------------------------------------------------------------- */

  if (state.v == NULL){
    MakeState (&state);
    #if (PARABOLIC_FLUX & EXPLICIT)
     dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

/* --------------------------------------------------------------
   2. Reset arrays.
      C_dt is an array used to store the inverse time step for
      advection and diffusion.
      We use C_dt[RHO] for advection, 
             C_dt[MX1] for viscosity,
             C_dt[BX1...BX3] for resistivity and
             C_dt[ENG] for thermal conduction.
   --------------------------------------------------------------- */

  #if DIMENSIONAL_SPLITTING == NO
   if (C_dt[RHO] == NULL){
     FOR_EACH(nv, 0, (&cdt_list)) {
       C_dt[nv] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
     }
   }

   if (g_intStage == 1) KTOT_LOOP(k) JTOT_LOOP(j){
     FOR_EACH(nv, 0, (&cdt_list)) {
       memset ((void *)C_dt[nv][k][j],'\0', NX1_TOT*sizeof(double));
     }
   }
  #endif

/* ------------------------------------------------
   2a. Compute current arrays 
   ------------------------------------------------ */

  #if (RESISTIVITY == EXPLICIT) && (defined STAGGERED_MHD)
   GetCurrent (d, -1, grid);
  #endif

/* ------------------------------------------------
   2b. Compute Temperature array
   ------------------------------------------------ */

  #if THERMAL_CONDUCTION == EXPLICIT
   if (T == NULL) T = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
   TOT_LOOP(k,j,i) T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
  #endif

/* ----------------------------------------------------------------
   3. Main loop on directions
   ---------------------------------------------------------------- */

  for (dir = beg_dir; dir <= end_dir; dir++){

    g_dir = dir;  
    SetIndexes (&indx, grid);  /* -- set normal and transverse indices -- */
    ResetState (d, &state, grid);

    #if (RESISTIVITY == EXPLICIT) && !(defined STAGGERED_MHD)
     GetCurrent(d, dir, grid);
    #endif

    TRANSVERSE_LOOP(indx,ip,i,j,k){
      g_i = i;  g_j = j;  g_k = k;
      for ((*ip) = 0; (*ip) < indx.ntot; (*ip)++) {
        VAR_LOOP(nv) state.v[(*ip)][nv] = d->Vc[nv][k][j][i];
        state.flag[*ip] = d->flag[k][j][i];
        #ifdef STAGGERED_MHD
         state.bn[(*ip)] = d->Vs[g_dir][k][j][i];
        #endif
      }
      CheckNaN (state.v, 0, indx.ntot-1,0);
      States  (&state, indx.beg - 1, indx.end + 1, grid); 
      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
       CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
      #endif
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux(d->Vc, d->J, T, &state, dcoeff, indx.beg-1, indx.end, grid);
      #endif
      #if UPDATE_VECTOR_POTENTIAL == YES
       VectorPotentialUpdate (d, NULL, &state, grid);
      #endif
      #ifdef SHEARINGBOX
       SB_SaveFluxes (&state, grid);
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, dt, grid);

    /* -- update:  U = U + dt*R -- */

      #ifdef CHOMBO
       for ((*ip) = indx.beg; (*ip) <= indx.end; (*ip)++) { 
         VAR_LOOP(nv) UU[nv][k][j][i] += state.rhs[*ip][nv];
       }
       SaveAMRFluxes (&state, aflux, indx.beg-1, indx.end, grid);
      #else
       for ((*ip) = indx.beg; (*ip) <= indx.end; (*ip)++) { 
         VAR_LOOP(nv) UU[k][j][i][nv] += state.rhs[*ip][nv];
       }
      #endif

      if (g_intStage > 1) continue;

    /* -- compute inverse dt coefficients when g_intStage = 1 -- */

      inv_dl = GetInverse_dl(grid);
      for ((*ip) = indx.beg; (*ip) <= indx.end; (*ip)++) { 
        #if DIMENSIONAL_SPLITTING == NO

         #if !GET_MAX_DT
          C_dt[0][k][j][i] += 0.5*(  Dts->cmax[(*ip)-1] 
                                   + Dts->cmax[*ip])*inv_dl[*ip];
         #endif
         #if (PARABOLIC_FLUX & EXPLICIT)
          dl2 = 0.5*inv_dl[*ip]*inv_dl[*ip];
          FOR_EACH(nv, 1, (&cdt_list)) {  
            C_dt[nv][k][j][i] += (dcoeff[*ip][nv]+dcoeff[(*ip)-1][nv])*dl2;
          }
         #endif

        #elif DIMENSIONAL_SPLITTING == YES

         #if !GET_MAX_DT
          Dts->inv_dta = MAX(Dts->inv_dta, Dts->cmax[*ip]*inv_dl[*ip]);
         #endif
         #if (PARABOLIC_FLUX & EXPLICIT)
          dl2 = inv_dl[*ip]*inv_dl[*ip];
          FOR_EACH(nv, 1, (&cdt_list)) {
            Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[*ip][nv]*dl2);
          }
         #endif
        #endif 
      }
    }
  }

/* -------------------------------------------------------------------
   4. Additional terms here
   ------------------------------------------------------------------- */

  #if (ENTROPY_SWITCH)  && (RESISTIVITY == EXPLICIT)
   EntropyOhmicHeating(d, UU, dt, grid);
  #endif

  #ifdef SHEARINGBOX
   SB_CorrectFluxes (UU, 0.0, dt, grid);
  #endif

  #ifdef STAGGERED_MHD
   CT_Update(d, d->Vs, dt, grid);
  #endif

  #if DIMENSIONAL_SPLITTING == YES
   return;
  #endif

/* -------------------------------------------------------------------
   5. Reduce dt for dimensionally unsplit schemes.
   ------------------------------------------------------------------- */

  if (g_intStage > 1) return;

  for (k = KBEG; k <= KEND; k++){ g_k = k;
  for (j = JBEG; j <= JEND; j++){ g_j = j;
    for (i = IBEG; i <= IEND; i++){
      #if !GET_MAX_DT
       Dts->inv_dta = MAX(Dts->inv_dta, C_dt[0][k][j][i]);
      #endif
      #if (PARABOLIC_FLUX & EXPLICIT)
       FOR_EACH(nv, 1, (&cdt_list)) { 
         Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[nv][k][j][i]);
       }
      #endif
    }
  }}
  #if !GET_MAX_DT
   Dts->inv_dta /= (double)DIMENSIONS;
  #endif
  #if (PARABOLIC_FLUX & EXPLICIT)
   Dts->inv_dtp /= (double)DIMENSIONS;
  #endif

}

#ifdef CHOMBO
/* ********************************************************************* */
void SaveAMRFluxes (const State_1D *state, double **aflux, 
                    int beg, int end,  Grid *grid)
/*! 
 *  Rebuild fluxes in a way suitable for AMR operation
 *  by adding pressure and multiplying by area.
 *
 *********************************************************************** */
{
  int  i, j, k, nv, *in;
  int nxf, nyf, nzf;
  int nxb, nyb, nzb;
  long int indf;
  double wflux, r, area;

  for (i = beg; i <= end; i++) state->flux[i][MXn] += state->press[i];

#if (GEOMETRY == CARTESIAN) && (CH_SPACEDIM > 1)
  if ((g_dir == IDIR) && (g_stretch_fact != 1.)) {
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) state->flux[i][nv] *= g_stretch_fact;
    }
  }
  #if (CH_SPACEDIM == 3)
  if ((g_dir == JDIR) && (g_x3stretch != 1.)) {
    for (i = beg; i <= end; i++) {
     NVAR_LOOP(nv) state->flux[i][nv] *= g_x3stretch;
    }
  }
  if ((g_dir == KDIR) && (g_x2stretch != 1.)) {
    for (i = beg; i <= end; i++) {
     NVAR_LOOP(nv) state->flux[i][nv] *= g_x2stretch;
    }
  }  
  #endif
#endif

#if GEOMETRY == CYLINDRICAL
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= grid[IDIR].A[i];
        #if CH_SPACEDIM > 1
        state->flux[i][nv] *= g_x2stretch;
        #endif
      }
    }
  }else{
    area = fabs(grid[IDIR].x[g_i]);
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) state->flux[i][nv] *= area;
    }
  }
#endif

#if GEOMETRY == SPHERICAL
  if (g_dir == IDIR){
  #if CH_SPACEDIM > 1
    area = grid[JDIR].dV[g_j]/g_level_dx;
    #if CH_SPACEDIM == 3
    area *= g_x3stretch;
    #endif
  #endif 
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= grid[IDIR].A[i];
        #if CH_SPACEDIM > 1
        state->flux[i][nv] *= area;
        #endif
      }
    #if (COMPONENTS == 3) && CHOMBO_CONS_AM
      state->flux[i][iMPHI] *= grid[IDIR].xr[i]*sin(grid[JDIR].x[g_j]);
    #endif
    }
  }
  if (g_dir == JDIR){
    area = fabs(grid[IDIR].x[g_i]);
  #if CHOMBO_LOGR == YES
    area *= grid[IDIR].dx[g_i]/g_level_dx;
  #endif
  #if CH_SPACEDIM == 3
    area *= g_x3stretch;
  #endif
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= grid[JDIR].A[i]*area;
      }
  #if (COMPONENTS == 3) && CHOMBO_CONS_AM
      state->flux[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].xr[i]);
  #endif
    }
  }
  if (g_dir == KDIR){
    for (i = beg; i <= end; i++) {
      area = g_x2stretch* fabs(grid[IDIR].x[g_i]);
    #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/g_level_dx;
    #endif
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= area;
      }
    #if (COMPONENTS == 3) && CHOMBO_CONS_AM
      state->flux[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].x[g_j]);
     #endif
    }
  }
#endif

#if GEOMETRY == POLAR
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= grid[IDIR].A[i];
        #if CH_SPACEDIM > 1
        state->flux[i][nv] *= g_x2stretch;
        #endif
        #if CH_SPACEDIM == 3
        state->flux[i][nv] *= g_x3stretch; 
        #endif
      }
      #if (COMPONENTS > 1) && CHOMBO_CONS_AM
      state->flux[i][iMPHI] *= grid[IDIR].xr[i];
      #endif
    }
  }
  if (g_dir == JDIR){
    area = g_x3stretch;
    #if CHOMBO_LOGR == YES
    area *= grid[IDIR].dx[g_i]/g_level_dx;
    #endif
    #if CH_SPACEDIM == 3
    area *= g_x3stretch;
    #endif 
    if (area != 1.) {
      for (i = beg; i <= end; i++) {
        NVAR_LOOP(nv) {
          state->flux[i][nv] *= area;
        }
      }
    }
    #if (COMPONENTS > 1) && CHOMBO_CONS_AM
    for (i = beg; i <= end; i++) state->flux[i][iMPHI] *= grid[IDIR].x[g_i];
    #endif
  }
  if (g_dir == KDIR){
    for (i = beg; i <= end; i++) {
      area = g_x2stretch*fabs(grid[IDIR].x[g_i]);
      #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/g_level_dx;
      #endif
      NVAR_LOOP(nv) {
        state->flux[i][nv] *= area;
      }
      #if (COMPONENTS > 1) && CHOMBO_CONS_AM
      state->flux[i][iMPHI] *= grid[IDIR].x[g_i];
      #endif
    }
  }
#endif
  
/* -- store fluxes for re-fluxing operation -- */

  nxf = grid[IDIR].np_int + (g_dir == IDIR);
  nyf = grid[JDIR].np_int + (g_dir == JDIR);
  nzf = grid[KDIR].np_int + (g_dir == KDIR);

  nxb = grid[IDIR].lbeg - (g_dir == IDIR);
  nyb = grid[JDIR].lbeg - (g_dir == JDIR);
  nzb = grid[KDIR].lbeg - (g_dir == KDIR);

#if TIME_STEPPING == RK2 
  wflux = 0.5;
#else
  wflux = 1.0;
#endif
 
  i = g_i; j = g_j; k = g_k;
  if (g_dir == IDIR) in = &i;
  if (g_dir == JDIR) in = &j;
  if (g_dir == KDIR) in = &k;
  
  for ((*in) = beg; (*in) <= end; (*in)++) {
    #if HAVE_ENERGY && ENTROPY_SWITCH
     state->flux[*in][ENTR] = 0.0;
    #endif
    NVAR_LOOP(nv) {
      indf = nv*nzf*nyf*nxf + (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);
      aflux[g_dir][indf] = wflux*state->flux[(*in)][nv];
    }
  }
 
}
#endif

/* ********************************************************************* */
intList TimeStepIndexList()
/*!
 * Return the intList of inverse time step indices.
 *
 *********************************************************************** */
{
  int i = 0;
  intList cdt;
  
  cdt.indx[i++] = RHO;
  #if VISCOSITY == EXPLICIT
   cdt.indx[i++] = MX1;
  #endif
  #if RESISTIVITY == EXPLICIT
   EXPAND(cdt.indx[i++] = BX1;  ,
          cdt.indx[i++] = BX2;  ,
          cdt.indx[i++] = BX3;)
  #endif
  #if THERMAL_CONDUCTION == EXPLICIT
   cdt.indx[i++] = ENG;
  #endif
  cdt.nvar = i;
  
  return cdt;
}
