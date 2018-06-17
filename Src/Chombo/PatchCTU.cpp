#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

static void StatesFlat (const Sweep *sweep, int beg, int end);

/* ********************************************************************* */
void PatchPluto::advanceStep(FArrayBox&       a_U,
                             FArrayBox&       a_Utmp,
                             const FArrayBox& a_dV,
                             FArrayBox&  split_tags,
                             BaseFab<unsigned char>& a_Flags,
                             FluxBox&         a_F,
                             timeStep        *Dts,
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
  int nbeg, nend, ntot; 
  int    errp, errm, errh;

  static unsigned char *flagp, *flagm;  // these should go inside sweep !!

  static Sweep sweep;

  State *stateC = &(sweep.stateC);
  State *stateL = &(sweep.stateL);
  State *stateR = &(sweep.stateR);

  double ***UU[NVAR], *du;
 #ifdef SKIP_SPLIT_CELLS 
  double ***splitcells;
 #endif
  double invDt_par, *inv_dl;
  double **up, **um, **vp, **vm;
  static Data d;
  static Data_Arr Uh, Vc0;
  static Data_Arr rhs, Up[DIMENSIONS], Um[DIMENSIONS];
  static double **u;

  Riemann_Solver *Riemann = rsolver;
  RBox   sweepBox;

/* ---------------------------------------------------------
   0. Check algorithm compatibilities
   --------------------------------------------------------- */

#if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL)
  print1 ("! CTU only works in cartesian or cylindrical coordinates\n");
  QUIT_PLUTO(1);
#endif     

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("! advanceStep(): need to re-allocate matrix\n");
    QUIT_PLUTO(1);
  }

/* ---------------------------------------------------------
   1. Map Chombo data structure.
   --------------------------------------------------------- */

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

/* ---------------------------------------------------------
   2. Allocate static memory areas
   ---------------------------------------------------------  */

  if (sweep.flux == NULL){

    MakeState (&sweep);

    d.Vc   = ARRAY_4D(NVAR, NX3_MAX, NX2_MAX, NX1_MAX, double);
    d.Uc   = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
    d.flag = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, unsigned char);
    #if RESISTIVITY
    d.J    = ARRAY_4D(3,NX3_MAX, NX2_MAX, NX1_MAX, double);
    #endif
    #if THERMAL_CONDUCTION
    d.Tc   = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #endif

    flagp = ARRAY_1D(NMAX_POINT, unsigned char);
    flagm = ARRAY_1D(NMAX_POINT, unsigned char);
    u     = ARRAY_2D(NMAX_POINT, NVAR, double);
    
    Uh  = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
    Vc0 = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);

    D_EXPAND(Um[IDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
             Up[IDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);  ,
  
             Um[JDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
             Up[JDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);  ,

             Um[KDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
             Up[KDIR] = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);)

    rhs = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
  }

  up = stateL->u;   
  um = stateR->u-1;
  vp = stateL->v;   
  vm = stateR->v-1;

  g_intStage = 1;
  TOT_LOOP(k,j,i) {
    d.flag[k][j][i] = 0;
    NVAR_LOOP(nv) rhs[k][j][i][nv] = 0.0;
  }

// Transpose array to have nv as fastest running index

  TOT_LOOP(k,j,i) NVAR_LOOP(nv) d.Uc[k][j][i][nv] = UU[nv][k][j][i];

#ifdef SKIP_SPLIT_CELLS
  DOM_LOOP(k,j,i){
    if (splitcells[k][j][i] < 0.5){
      d.flag[k][j][i] |= FLAG_SPLIT_CELL;
    }
  }
#endif
  getPrimitiveVars (d.Uc, &d, grid);
#if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH)
  FlagShock (&d, grid);
#endif

// Copy solution values at t^n as safe replacement values of
// corner-coupled-sweeps.

  TOT_LOOP(k,j,i) NVAR_LOOP(nv) {
    Vc0[k][j][i][nv] = d.Vc[nv][k][j][i];
  }

/* ----------------------------------------------------
   3. Normal predictors 
   ---------------------------------------------------- */

  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

  /* -- Set integration box for predictor step -- */

    RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
    RBoxSetDirections (&sweepBox, g_dir);
    SetVectorIndices (g_dir);

    D_EXPAND(                                      ;  ,
             (*sweepBox.tbeg)--; (*sweepBox.tend)++;  ,
             (*sweepBox.bbeg)--; (*sweepBox.bend)++;)
  #if (defined STAGGERED_MHD)
    D_EXPAND((*sweepBox.nbeg)--; (*sweepBox.nend)++;  ,
             (*sweepBox.tbeg)--; (*sweepBox.tend)++;  ,
             (*sweepBox.bbeg)--; (*sweepBox.bend)++;)
  #else
    #if (PARABOLIC_FLUX & EXPLICIT)
    (*sweepBox.nbeg)--; (*sweepBox.nend)++;
    #endif
  #endif
    ResetState(&d, &sweep, grid);

    nbeg = *sweepBox.nbeg;
    nend = *sweepBox.nend;
    ntot = grid->np_tot[g_dir];
    
    NVAR_LOOP(nv) sweep.rhs[nbeg-1][nv] = sweep.rhs[nend+1][nv] = 0.0;

    BOX_TRANSVERSE_LOOP(&sweepBox, k,j,i){  
      in  = sweepBox.n;
      g_i = i; g_j = j; g_k = k;

    /* ---- Get 1-D arrays of primitive quantities ---- */

      for ((*in) = 0; (*in) < ntot; (*in)++) {
        NVAR_LOOP(nv) stateC->v[*in][nv] = d.Vc[nv][k][j][i];
        #ifdef GLM_MHD
        NVAR_LOOP(nv) sweep.vn[*in][nv] = d.Vc[nv][k][j][i];
        #endif
        sweep.flag[*in] = d.flag[k][j][i];
      }

      CheckNaN (stateC->v, nbeg-1, nend+1, 0);
      PrimToCons (stateC->v, u, 0, ntot-1);

#if !(PARABOLIC_FLUX & EXPLICIT)
      States  (&sweep, nbeg-1, nend+1, grid);

 /* -- Copy sweeps before Riemann solver changes them (for GLM) -- */

      #ifdef GLM_MHD
      for (*in = nbeg-1; *in <= nend+1; (*in)++) {
        NVAR_LOOP(nv){
          Up[g_dir][k][j][i][nv] = up[*in][nv];
          Um[g_dir][k][j][i][nv] = um[*in][nv];
        }
      }
      #endif

      Riemann (&sweep, nbeg-1, nend, Dts->cmax, grid);
      RightHandSide (&sweep, Dts, nbeg, nend, 0.5*g_dt, grid);

      for ((*in) = nbeg-1; (*in) <= nend+1; (*in)++){
        NVAR_LOOP(nv){
          d.Uc[k][j][i][nv] = u[*in][nv];
          #ifdef GLM_MHD
          Up[g_dir][k][j][i][nv] -= sweep.rhs[*in][nv];
          Um[g_dir][k][j][i][nv] -= sweep.rhs[*in][nv];
          #else
          Up[g_dir][k][j][i][nv] = up[*in][nv] - sweep.rhs[*in][nv];
          Um[g_dir][k][j][i][nv] = um[*in][nv] - sweep.rhs[*in][nv];
          #endif
          rhs[k][j][i][nv]       += sweep.rhs[*in][nv];
      }}
#else
      StatesFlat (&sweep, 0, ntot-1);            
      Riemann (&sweep, nbeg-1, nend, Dts->cmax, grid);
      RightHandSide (&sweep, Dts, nbeg, nend, 0.5*g_dt, grid);
      for (*in = nbeg; *in <= nend; (*in)++) {
        NVAR_LOOP(nv) rhs[k][j][i][nv] += sweep.rhs[*in][nv];
      }

      States  (&sweep, nbeg, nend, grid);

      for (*in = nbeg; *in <= nend; (*in)++) {
        NVAR_LOOP(nv){
          Up[g_dir][k][j][i][nv] = up[*in][nv] - sweep.rhs[*in][nv];
          Um[g_dir][k][j][i][nv] = um[*in][nv] - sweep.rhs[*in][nv];
        }
      }
#endif

    /*  -- Compute inverse time step -- */

      inv_dl = GetInverse_dl(grid);
      for (*in = nbeg; *in <= nend; (*in)++) { 
        Dts->invDt_hyp = MAX(Dts->invDt_hyp, Dts->cmax[*in]*inv_dl[*in]);
      }

    }
  }

/* ----------------------------------------------------------
   4. Advance conservative variables by dt/2.
   ---------------------------------------------------------- */

  RBoxDefine (IBEG, IEND,JBEG, JEND,KBEG, KEND, CENTER, &sweepBox);  
#if (PARABOLIC_FLUX & EXPLICIT) || (defined STAGGERED_MHD)
  D_EXPAND(sweepBox.ibeg--; sweepBox.iend++;  ,
           sweepBox.jbeg--; sweepBox.jend++;  ,
           sweepBox.kbeg--; sweepBox.kend++;);
#endif

  BOX_LOOP(&sweepBox, k,j,i){
    NVAR_LOOP(nv) Uh[k][j][i][nv] = d.Uc[k][j][i][nv] + rhs[k][j][i][nv];
  }

#if (PARABOLIC_FLUX & EXPLICIT)
  ParabolicUpdate (&d, Uh, &sweepBox, NULL, 0.5*g_dt, Dts, grid);
  D_EXPAND(ParabolicUpdate (NULL, Up[IDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);
           ParabolicUpdate (NULL, Um[IDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);  ,
           
           ParabolicUpdate (NULL, Up[JDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);
           ParabolicUpdate (NULL, Um[JDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);  ,
           
           ParabolicUpdate (NULL, Up[KDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);
           ParabolicUpdate (NULL, Um[KDIR], &sweepBox, NULL, 0.5*g_dt, Dts, grid);)
#endif

/* ---------------------------------------------------------------
   5. Convert time-centered conserved array to primitive array 
   --------------------------------------------------------------- */

  g_dir = IDIR;
  for (k = sweepBox.kbeg; k <= sweepBox.kend; k++){
  for (j = sweepBox.jbeg; j <= sweepBox.jend; j++){
    g_j = j; g_k = k;
    errp = ConsToPrim(Uh[k][j], stateC->v, sweepBox.ibeg, sweepBox.iend, d.flag[k][j]);
    for (i = sweepBox.ibeg; i <= sweepBox.iend; i++) {
      NVAR_LOOP(nv) d.Vc[nv][k][j][i] = stateC->v[i][nv];
    }
  }}

/* ----------------------------------------------------
   6.  Final Conservative Update
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

  static double *aflux[3];
  for (i = 0; i < DIMENSIONS; i++) aflux[i] = a_F[i].dataPtr(0);

  g_intStage = 2;
  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

  /* -- Set integration box for corrector step -- */

    RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
    RBoxSetDirections (&sweepBox, g_dir);
    SetVectorIndices (g_dir);
    ResetState(&d, &sweep, grid);
     
    nbeg = *sweepBox.nbeg;
    nend = *sweepBox.nend;
     
    BOX_TRANSVERSE_LOOP(&sweepBox,k,j,i){  
      in  = sweepBox.n;
      g_i = i;  g_j = j;  g_k = k;

    /* --------------------------------------------
        Correct normal predictors with transverse 
        fluxes to obtain corner coupled sweeps. 
       -------------------------------------------- */

      for ((*in) = nbeg-1; (*in) <= nend+1; (*in)++){
        NVAR_LOOP(nv){
          up[*in][nv] = Up[g_dir][k][j][i][nv] + rhs[k][j][i][nv];
          um[*in][nv] = Um[g_dir][k][j][i][nv] + rhs[k][j][i][nv];
          stateC->v[*in][nv] = d.Vc[nv][k][j][i];
        }
        flagp[*in] = flagm[*in] = d.flag[k][j][i];
        sweep.flag[*in] = d.flag[k][j][i];
      }

    /* ------------------------------------------------
        convert time and cell centered sweep to 
        conservative vars and corner-coupled sweeps to 
        primitive.
      ------------------------------------------------ */

      errp = ConsToPrim (up, vp, nbeg-1, nend+1, flagp);
      errm = ConsToPrim (um, vm, nbeg-1, nend+1, flagm);

  /* ---- check admissibility of corner coupled sweeps ---- */

      if (errm || errp){
        WARNING( pout() << "! PatchPluto::advanceStep():"
                        << " corner coupled sweeps not physical: "
                        << "reverting to 1st order (level= " <<  m_level
                        << ")" << endl;)

        for ((*in) = nbeg-1; (*in) <= nend+1; (*in)++){
          if (    (flagp[*in] & FLAG_CONS2PRIM_FAIL)
               || (flagm[*in] & FLAG_CONS2PRIM_FAIL)){
            NVAR_LOOP(nv) stateC->v[*in][nv] = Vc0[k][j][i][nv];

            for (nv = 0; nv < NVAR; nv++) {
              stateL->v[*in][nv] = stateR->v[*in-1][nv] = stateC->v[*in][nv];
            }
          }
        }
      }

    /* -------------------------------------------------------
           compute hyperbolic and parabolic fluxes 
       ------------------------------------------------------- */

      Riemann (&sweep, nbeg-1, nend, Dts->cmax, grid);
      RightHandSide (&sweep, Dts, nbeg, nend, g_dt, grid);

    /* -- Store fluxes for regridding operations -- */

      for ((*in) = nbeg-1; (*in) <= nend; (*in)++){
        sweep.flux[*in][MXn] += sweep.press[*in];
        #if HAVE_ENERGY && ENTROPY_SWITCH
        sweep.flux[*in][ENG] = 0.0;
        sweep.flux[*in][ENTR] = 0.0;
        #endif
      }   
      StoreAMRFlux (sweep.flux, aflux, 0, 0, NVAR-1, nbeg-1, nend, grid);

      for ((*in) = nbeg; (*in) <= nend; (*in)++) {
      for (nv = 0; nv < NVAR; nv++) {
        d.Uc[k][j][i][nv] += sweep.rhs[*in][nv];
      }}

    /*  -- Compute inverse time step -- */

      inv_dl = GetInverse_dl(grid);
      for (*in = nbeg; *in <= nend; (*in)++) { 
        Dts->invDt_hyp = MAX(Dts->invDt_hyp, Dts->cmax[*in]*inv_dl[*in]);
      }
    }
  }

//pout() << "g_stepNumber = " << g_stepNumber << endl;

/* ---------------------------------------------------------------
   7. Update solution array with parabolic (diffusion) terms.
   --------------------------------------------------------------- */

  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);

#if (PARABOLIC_FLUX & EXPLICIT)
  ParabolicUpdate (&d, d.Uc, &sweepBox, aflux, g_dt, Dts, grid);
#endif

  saveFluxes (aflux, grid);

  #ifdef GLM_MHD
   glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->invDt_hyp*m_dx);
//   glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->inv_dta); /* If subcycling is turned off */
   double dtdx = g_dt/g_coeff_dl_min/m_dx;
//    double dtdx = g_dt/g_coeff_dl_min; /* If subcycling is turned off */
   GLM_Source (d.Uc, dtdx, grid);
  #endif

/* ------------------------------------------------
   8. Source terms included via operator splitting
   ------------------------------------------------ */

#if COOLING != NO
  ConsToPrim3D(d.Uc, d.Vc, d.flag, &sweepBox);
  SplitSource (&d, g_dt, Dts, grid);
  PrimToCons3D(d.Vc, d.Uc, &sweepBox);
#endif

#if ENTROPY_SWITCH
/* ---------------------------------------------------------------------
   9. At this stage we have U^(n+1) that contains both total energy (E)
      and entropy (sigma_c) although they have evolved differently.
      To synchronize them we convert UU to primitive and then again
      to conservative. This will ensure that in every cell
      E and sigma_c can be mapped one into another consistently.
      This step is *essential* when, at then next step,
      primitive variables will be computed in every zone from entropy
      rather than selectively from energy and entropy.
   --------------------------------------------------------------------- */
   
  ConsToPrim3D(d.Uc, d.Vc, d.flag, &sweepBox);
  PrimToCons3D(d.Vc, d.Uc, &sweepBox);
#endif

/* ---------------------------------------------------------------
    We pass U*dV/m_dx^3 back to Chombo rather than U.
   --------------------------------------------------------------- */

  TOT_LOOP(k,j,i) NVAR_LOOP(nv)  UU[nv][k][j][i] = d.Uc[k][j][i][nv];

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

/* ********************************************************************* */
void StatesFlat (const Sweep *sweep, int beg, int end)
/*
 *  Compute first order sweeps.
 *********************************************************************** */
{
  int i, nv;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);
  double **vp = stateL->v;
  double **vm = stateR->v - 1;
  double **up = stateL->u;
  double **um = stateR->u - 1;

  for (i = beg; i <= end; i++) {
  for (nv = NVAR; nv--;  ) {
    vp[i][nv] = vm[i][nv] = stateC->v[i][nv];
  }}
#ifdef STAGGERED_MHD
  for (i = beg; i <= end-1; i++) {
    stateL->v[i][BXn] = stateR->v[i][BXn] = sweep->bn[i];
  }
#endif
   PrimToCons(vm, um, beg, end);
   PrimToCons(vp, up, beg, end);
}


