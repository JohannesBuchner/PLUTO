#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::advanceStep(FArrayBox&       a_U,
                             FArrayBox&       a_Utmp,
                             const FArrayBox& a_dV,
                             FArrayBox&  split_tags,
                             BaseFab<unsigned char>& a_Flags,
                             FluxBox&         a_F,
                             Time_Step        *Dts,
                             const Box&       UBox, 
                             Grid *grid)
/*
 *
 *
 *
 *
 *********************************************************************** */
{
  CH_assert(isDefined());
  CH_assert(UBox == m_currentBox);

  int nv, *in;
  int nxf, nyf, nzf, indf;
  int nxb, nyb, nzb;
  int i, j, k;

  int    errp, errm, errh;
  double ***UU[NVAR], *du;
 #ifdef SKIP_SPLIT_CELLS 
  double ***splitcells;
 #endif
  double inv_dtp, *inv_dl;
  static Data d;
  static Data_Arr UH, dU;
  static Data_Arr UP[DIMENSIONS], UM[DIMENSIONS];
  static double **u, ***T;
  #if (PARABOLIC_FLUX & EXPLICIT)
   static double **dcoeff;
  #endif   
  RBox *rbox = GetRBox(DOM, CENTER);
  Index indx;
  static State_1D state;
  static unsigned char *flagp, *flagm;  // these should go inside state !!

  Riemann_Solver *Riemann = rsolver;

/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL)
   print1 ("! CTU only works in cartesian or cylindrical coordinates\n");
   QUIT_PLUTO(1);
  #endif     

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("! advanceStep(): need to re-allocate matrix\n");
    QUIT_PLUTO(1);
  }

/* -----------------------------------------------------------------
                          Allocate memory
   ----------------------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   for (nv = 0; nv < NVAR; nv++) a_U.divide(a_dV,0,nv);
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     Box curBox = a_U.box();
     for(BoxIterator bit(curBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_U(iv,iMPHI) /= a_dV(iv,1);
       a_U(iv,iMPHI) -= a_U(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
     }
    #else
     a_U.divide(a_dV,1,iMPHI);
    #endif
   #endif
  #else
   if (g_stretch_fact != 1.) a_U /= g_stretch_fact;
  #endif

  for (nv = 0; nv < NVAR; nv++){
    UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, a_U.dataPtr(nv));
  }
  #ifdef SKIP_SPLIT_CELLS
   splitcells = ArrayBoxMap(KBEG, KEND, JBEG, JEND, IBEG, IEND, 
                            split_tags.dataPtr(0));
  #endif
  #if RESISTIVITY != NO
   if (d.J == NULL) d.J = ARRAY_4D(3,NX3_MAX, NX2_MAX, NX1_MAX, double);
  #endif

/* -----------------------------------------------------------
         Allocate static memory areas
   -----------------------------------------------------------  */

  if (state.flux == NULL){

    MakeState (&state);

    nxf = nyf = nzf = 1;
    D_EXPAND(nxf = NMAX_POINT;  ,
             nyf = NMAX_POINT;  ,
             nzf = NMAX_POINT;)

    d.Vc   = ARRAY_4D(NVAR, nzf, nyf, nxf, double);
    d.flag = ARRAY_3D(nzf, nyf, nxf, unsigned char);
 
    flagp = ARRAY_1D(NMAX_POINT, unsigned char);
    flagm = ARRAY_1D(NMAX_POINT, unsigned char);
    u     = ARRAY_2D(NMAX_POINT, NVAR, double);
    
    UH = ARRAY_4D(nzf, nyf, nxf, NVAR, double);
    dU = ARRAY_4D(nzf, nyf, nxf, NVAR, double);

    D_EXPAND(UM[IDIR]  = ARRAY_4D(nzf, nyf, nxf, NVAR, double);
             UP[IDIR]  = ARRAY_4D(nzf, nyf, nxf, NVAR, double);  ,
  
             UM[JDIR]  = ARRAY_4D(nzf, nxf, nyf, NVAR, double);
             UP[JDIR]  = ARRAY_4D(nzf, nxf, nyf, NVAR, double);  ,

             UM[KDIR]  = ARRAY_4D(nyf, nxf, nzf, NVAR, double);
             UP[KDIR]  = ARRAY_4D(nyf, nxf, nzf, NVAR, double);)

    #if (PARABOLIC_FLUX & EXPLICIT)
     dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
    #if THERMAL_CONDUCTION == EXPLICIT
     T = ARRAY_3D(nzf, nyf, nxf, double);
    #endif
  }

  g_intStage = 1;
  TOT_LOOP(k,j,i) d.flag[k][j][i] = 0;

  #ifdef SKIP_SPLIT_CELLS
   DOM_LOOP(k,j,i){
     if (splitcells[k][j][i] < 0.5){
       d.flag[k][j][i] |= FLAG_SPLIT_CELL;
     }
   }
  #endif
  getPrimitiveVars (UU, &d, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH)
   FlagShock (&d, grid);
  #endif
  #if THERMAL_CONDUCTION == EXPLICIT
   TOT_LOOP(k,j,i) T[k][j][i] = d.Vc[PRS][k][j][i]/d.Vc[RHO][k][j][i];
  #endif

/* ----------------------------------------------------
      1. Normal predictors 
   ---------------------------------------------------- */

  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);
    ResetState (&d, &state, grid);
    #if (RESISTIVITY == EXPLICIT)
     GetCurrent(&d, g_dir, grid);
    #endif
    TRANSVERSE_LOOP(indx,in,i,j,k){  
      g_i = i; g_j = j; g_k = k;

      state.up = UP[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uL = state.up;
      state.um = UM[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uR = state.um + 1;

      for ((*in) = 0; (*in) < indx.ntot; (*in)++) {
        NVAR_LOOP(nv) state.v[*in][nv] = d.Vc[nv][k][j][i];
        state.flag[*in] = d.flag[k][j][i];
      }

      CheckNaN (state.v, indx.beg-1, indx.end+1, 0);
      PrimToCons (state.v, u, 0, indx.ntot-1);

#if !(PARABOLIC_FLUX & EXPLICIT)
      States  (&state, indx.beg-1, indx.end+1, grid);
      Riemann (&state, indx.beg-1, indx.end, Dts->cmax, grid);
      RightHandSide (&state, Dts, indx.beg, indx.end, 0.5*g_dt, grid);

      if (g_dir == IDIR){  /* -- initialize UU, UH and dU -- */
        for (nv = NVAR; nv--; ){
          state.rhs[indx.beg-1][nv] = 0.0;
          state.rhs[indx.end+1][nv] = 0.0;
        }
        for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++){
        for (nv = NVAR; nv--; ){
          UU[nv][k][j][i] = u[*in][nv];
          UH[k][j][i][nv] = u[*in][nv] + state.rhs[*in][nv];
          dU[k][j][i][nv] = state.rhs[*in][nv];
          state.up[*in][nv]  -= state.rhs[*in][nv];
          state.um[*in][nv]  -= state.rhs[*in][nv];
        }}
      }else{
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++){
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] += state.rhs[*in][nv];
          UH[k][j][i][nv] += state.rhs[*in][nv];
          state.up[*in][nv] -= state.rhs[*in][nv];
          state.um[*in][nv] -= state.rhs[*in][nv];
        }}
      }
#else
      for ((*in) = 0; (*in) < indx.ntot; (*in)++) {
      for (nv = NVAR; nv--;  ) {
        state.vp[*in][nv] = state.vm[*in][nv] = state.vh[*in][nv] = state.v[*in][nv];
      }}
      PrimToCons(state.vm, state.um, 0, indx.ntot-1);
      PrimToCons(state.vp, state.up, 0, indx.ntot-1);
      
      Riemann (&state, indx.beg-1, indx.end, Dts->cmax, grid);

  /* -----------------------------------------------------------
          compute rhs using the hyperbolic fluxes only
     ----------------------------------------------------------- */

      #if (VISCOSITY == EXPLICIT)
       for ((*in) = 0; (*in) < indx.ntot; (*in)++) for (nv = NVAR; nv--;  )
         state.par_src[*in][nv] = 0.0;
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, 0.5*g_dt, grid);
      ParabolicFlux (d.Vc, d.J, T, &state, dcoeff, indx.beg-1, indx.end, grid);

  /* ----------------------------------------------------------
      compute LR states and subtract normal (hyperbolic) 
      rhs contribution.
      NOTE: states are computed from IBEG - 1 (= indx.beg)
            up to IEND + 1 (= indx.end) since EMF has already
            been evaluated and stored using 1st order states
            above.
     ---------------------------------------------------------- */

      States (&state, indx.beg, indx.end, grid);
      for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
      for (nv = NVAR; nv--; ){
        state.up[*in][nv] -= state.rhs[*in][nv];
        state.um[*in][nv] -= state.rhs[*in][nv];
      }}

  /* -----------------------------------------------------------
       re-compute the full rhs using the total (hyp+par) rhs
     ----------------------------------------------------------- */

      RightHandSide (&state, Dts, indx.beg, indx.end, 0.5*g_dt, grid);
      if (g_dir == IDIR){
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] = state.rhs[*in][nv];
          UH[k][j][i][nv] = u[*in][nv] + state.rhs[*in][nv];
        }}
      }else{
        for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
        for (nv = NVAR; nv--; ){
          dU[k][j][i][nv] += state.rhs[*in][nv];
          UH[k][j][i][nv] += state.rhs[*in][nv];
        }}
      }
#endif

    }
  }

/* ------------------------------------------------
    2. compute time and cell centered state. 
       Useful for source terms like gravity, 
       curvilinear terms and Powell's 8wave. 
   ------------------------------------------------ */
   
  g_dir = IDIR;
  SetIndexes (&indx, grid);
  TRANSVERSE_LOOP(indx,in,i,j,k){
    g_i = i; g_j = j; g_k = k;
    errp = ConsToPrim(UH[k][j], state.v, indx.beg, indx.end, d.flag[k][j]);
    WARNING(
      if (errp != 0)  print("! PatchUnsplit: error recovering U^{n+1/2}\n");
    )
    for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
      NVAR_LOOP(nv) d.Vc[nv][k][j][i] = state.v[*in][nv];
    }
    #if THERMAL_CONDUCTION == EXPLICIT
     for ((*in) = indx.beg; (*in) <= indx.end; (*in)++){
       T[k][j][i] = d.Vc[PRS][k][j][i]/d.Vc[RHO][k][j][i];
     }
    #endif
  }

/* ----------------------------------------------------
           2. Final Conservative Update
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

  g_intStage = 2;
  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);

    nxf = grid[IDIR].np_int + (g_dir == IDIR);
    nyf = grid[JDIR].np_int + (g_dir == JDIR);
    nzf = grid[KDIR].np_int + (g_dir == KDIR);

    nxb = grid[IDIR].lbeg - (g_dir == IDIR);
    nyb = grid[JDIR].lbeg - (g_dir == JDIR);
    nzb = grid[KDIR].lbeg - (g_dir == KDIR);

    ResetState (&d, &state, grid);
    #if (RESISTIVITY == EXPLICIT)
     GetCurrent(&d, g_dir, grid);
    #endif    
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      g_i = i;  g_j = j;  g_k = k;

      state.up = UP[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uL = state.up;
      state.um = UM[g_dir][*(indx.pt2)][*(indx.pt1)]; state.uR = state.um + 1;

    /* --------------------------------------------
        Correct normal predictors with transverse 
        fluxes to obtain corner coupled states. 
       -------------------------------------------- */

      for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++){
        for (nv = NVAR; nv--;  ) {
          state.up[*in][nv] += dU[k][j][i][nv];
          state.um[*in][nv] += dU[k][j][i][nv];
          state.vh[*in][nv]  = d.Vc[nv][k][j][i];
        }
        flagp[*in] = flagm[*in] = state.flag[*in] = d.flag[k][j][i];
      }

    /* ------------------------------------------------
        convert time and cell centered state to 
        conservative vars and corner-coupled states to 
        primitive.
      ------------------------------------------------ */
       
      errp = ConsToPrim (state.up, state.vp, indx.beg-1, indx.end+1, flagp);
      errm = ConsToPrim (state.um, state.vm, indx.beg-1, indx.end+1, flagm);

  /* ---- check admissibility of corner coupled states ---- */

      if (errm || errp){
        WARNING(
          print ("! Corner coupled states not physical: reverting to 1st order (level=%d)\n",
                m_level);  
          showPatch(grid);
        )
        for ((*in) = indx.beg-1; (*in) <= indx.end+1; (*in)++){
          if (    (flagp[*in] & FLAG_CONS2PRIM_FAIL)
               || (flagm[*in] & FLAG_CONS2PRIM_FAIL)){
            for (nv = 0; nv < NVAR; nv++) state.v[*in][nv] = d.Vc[nv][k][j][i];

            for (nv = 0; nv < NVAR; nv++) {
              state.vm[*in][nv] = state.vp[*in][nv] = state.vh[*in][nv] = state.v[*in][nv];
            }
          }
        }
      }

    /* -------------------------------------------------------
           compute hyperbolic and parabolic fluxes 
       ------------------------------------------------------- */

      Riemann (&state, indx.beg-1, indx.end, Dts->cmax, grid);
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux (d.Vc, d.J, T, &state, dcoeff, indx.beg - 1, indx.end, grid);
       inv_dl  = GetInverse_dl(grid);
       for ((*in) = indx.beg-1; (*in) <= indx.end; (*in)++) {
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
      #endif

      RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);
      saveFluxes (&state, indx.beg-1, indx.end, grid);

      for ((*in) = indx.beg; (*in) <= indx.end; (*in)++) {
      for (nv = 0; nv < NVAR; nv++) {
        UU[nv][k][j][i] += state.rhs[*in][nv];
      }}   
        
// Put fluxes in the FarrayBox a_F to be passed to Chombo

      for ((*in) = indx.beg-1; (*in) <= indx.end; (*in)++) {
        #if HAVE_ENERGY && ENTROPY_SWITCH
         state.flux[*in][ENG] = 0.0;
         state.flux[*in][ENTR] = 0.0;
        #endif
        for (nv = 0; nv < NVAR; nv++) {
          indf = nv*nzf*nyf*nxf + (k - nzb)*nyf*nxf 
                                + (j - nyb)*nxf 
                                + (i - nxb);
          a_F[g_dir].dataPtr(0)[indf] = state.flux[*in][nv];
        }
      }
    }
  }

  #ifdef GLM_MHD
   glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->inv_dta*m_dx);
//   glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->inv_dta); /* If subcycling is turned off */
   double dtdx = g_dt/g_coeff_dl_min/m_dx;
//    double dtdx = g_dt/g_coeff_dl_min; /* If subcycling is turned off */
   GLM_Source (UU, dtdx, grid);
  #endif

/* ----------------------------------------------
    Source terms included via operator splitting
   ---------------------------------------------- */

#if COOLING != NO
  ConsToPrim3D(UU, d.Vc, d.flag, rbox);
  SplitSource (&d, g_dt, Dts, grid);
  PrimToCons3D(d.Vc, UU, rbox);
#endif

#if ENTROPY_SWITCH
/* -------------------------------------------------------------------
    At this stage we have U^(n+1) that contains both total energy (E)
    and entropy (sigma_c) although they have evolved differently.
    To synchronize them we convert UU to primitive and then again
    to conservative. This will ensure that in every cell
    E and sigma_c can be mapped one into another consistently.
    This step is *essential* when, at then next step,
    primitive variables will be computed in every zone from entropy
    rather than selectively from energy and entropy.
   ------------------------------------------------------------------- */
   
  ConsToPrim3D(UU, d.Vc, d.flag, rbox);
  PrimToCons3D(d.Vc, UU, rbox);
#endif

/* ---------------------------------------------------------------
    We pass U*dV/m_dx^3 back to Chombo rather than U.
   --------------------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     for(BoxIterator bit(curBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_U(iv,iMPHI) += a_U(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
       a_U(iv,iMPHI) *= a_dV(iv,1);
     }
    #else
     a_U.mult(a_dV,1,iMPHI);
    #endif
   #endif
   for (nv = 0; nv < NVAR; nv++) a_U.mult(a_dV,0,nv);
  #else
   if (g_stretch_fact != 1.) a_U *= g_stretch_fact;
  #endif

/* -------------------------------------------------
               Free memory 
   ------------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);

  #ifdef SKIP_SPLIT_CELLS
   FreeArrayBoxMap (splitcells, KBEG, KEND, JBEG, JEND, IBEG, IEND);
  #endif
}
