#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::advanceStep(FArrayBox&       a_U,
                             FArrayBox&       a_Utmp,
                             const FArrayBox& a_dV,
                             FArrayBox&       split_tags,
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

  int nv, in;
  int nxf, nyf, nzf, indf;
  int nxb, nyb, nzb;
  int i, j, k;

  double ***UU[NVAR];
  double *inv_dl, dl2, cylr;
  static Data d;
#ifdef SKIP_SPLIT_CELLS
  double ***splitcells;
#endif
#if (PARABOLIC_FLUX & EXPLICIT)
  static double **dcoeff;
#endif   
  Index indx;
  static State_1D state;

  Riemann_Solver *Riemann;
  Riemann = rsolver;
#if TIME_STEPPING == RK2 
  double wflux = 0.5;
#else
  double wflux = 1.;
#endif
  RBox *rbox = GetRBox(DOM, CENTER);
  
/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

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
#if (TIME_STEPPING == RK2)
  d.flag = ArrayCharMap(NX3_TOT, NX2_TOT, NX1_TOT,a_Flags.dataPtr(0));
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

    d.Vc = ARRAY_4D(NVAR, nzf, nyf, nxf, double);
#if (TIME_STEPPING != RK2)
    d.flag = ARRAY_3D(nzf, nyf, nxf, unsigned char);
#endif 
#if (PARABOLIC_FLUX & EXPLICIT)
    dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
#endif
  }

  if (g_intStage == 1) TOT_LOOP(k,j,i) d.flag[k][j][i] = 0;      

  getPrimitiveVars (UU, &d, grid);
#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
  if (g_intStage == 1) FlagShock (&d, grid); /* Do it only at predictor */
#endif

#ifdef SKIP_SPLIT_CELLS
  if (g_intStage == 1) {
    DOM_LOOP(k,j,i){
      if (splitcells[k][j][i] < 0.5) d.flag[k][j][i] |= FLAG_SPLIT_CELL;
    }
  }
#endif

/* ----------------------------------------------------
    Reset arrays
   ---------------------------------------------------- */

#if (TIME_STEPPING == RK2)
  if (g_intStage == 1) a_Utmp.copy(a_U); //Temporary copy of old conserved variables
#endif

/* ----------------------------------------------------
    Loop on directions
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

  static double *aflux[3];
  for (in = 0; in < DIMENSIONS; in++) aflux[in] = a_F[in].dataPtr(0);

  UpdateStage(&d, UU, aflux, Riemann, g_dt, Dts, grid);

// Compute advective/diffusive timestep (predictor only)

#ifdef GLM_MHD
  if (g_intStage == 1) glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->inv_dta*m_dx);
//  if (g_intStage == 1) glm_ch_max_loc = MAX(glm_ch_max_loc, Dts->inv_dta); /* If subcycling is turned off */
#endif

// Final update: average old and new conservative variables

#if (TIME_STEPPING == RK2)
  if (g_intStage == 2) {
    a_U.plus(a_Utmp);
    a_U *= 0.5;
#endif

/* ----------------------------------------------
    Source terms included via operator splitting
   ---------------------------------------------- */

#ifdef GLM_MHD
   double dtdx = g_dt/g_coeff_dl_min/m_dx;
//    double dtdx = g_dt/g_coeff_dl_min; /* If subcycling is turned off */
   GLM_Source (UU, dtdx, grid);
#endif


#if COOLING != NO
    ConsToPrim3D(UU, d.Vc, d.flag, rbox);
    SplitSource (&d, g_dt, Dts, grid);
    PrimToCons3D(d.Vc, UU, rbox);
#endif

#if (TIME_STEPPING == RK2)
  }
#endif

#if ENTROPY_SWITCH

/* -------------------------------------------------------------------
    At this point we have U^(n+1) that contains both total energy (E)
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
 
  #if (TIME_STEPPING == RK2)
   FreeArrayCharMap(d.flag);
  #endif
}
