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

  \date   Sep 07, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void UpdateStage(Data *d, Data_Arr UU, double **aflux,
                 Riemann_Solver *Riemann, double dt, timeStep *Dts, 
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
 * \param [in]      grid     pointer to Grid structure
 *********************************************************************** */
{
  int  i, j, k;
  int  nv, dir, beg_dir, end_dir;
  int  *ip;

  static Sweep sweep;
  State *stateC = &(sweep.stateC);
  State *stateL = &(sweep.stateL);

  double *inv_dl, dl2;
  static double ***C_dt;
  RBox  sweepBox;

#if DIMENSIONAL_SPLITTING == YES
  beg_dir = end_dir = g_dir;
#else
  beg_dir = 0;
  end_dir = DIMENSIONS-1;
#endif

/* --------------------------------------------------------
   0. Allocate memory & reset arrays.
      C_dt is an array used to store the inverse time
      step for the hyperbolic solve.
   -------------------------------------------------------- */

  if (stateC->v == NULL){
    MakeState (&sweep);
    #if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
    C_dt = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #endif
  }

#if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
  if (g_intStage == 1){
    KTOT_LOOP(k) JTOT_LOOP(j){
      memset ((void *)C_dt[k][j],'\0', NX1_TOT*sizeof(double));
    }
  }
#endif

/* --------------------------------------------------------
   1. Compute Fcr force only at predictor step
      (g_intStage == 1).
   -------------------------------------------------------- */

#ifdef PARTICLES
  #if PARTICLES_TYPE == COSMIC_RAYS
  if (g_intStage == 1) Particles_CR_ComputeForce (d->Vc, d, grid);
  #elif PARTICLES_TYPE == DUST
  if (g_intStage == 1) Particles_Dust_ComputeForce (d->Vc, d, grid);
  #endif
#endif

#if FORCED_TURB == YES
  ForcedTurb *Ft;
  Ft = d->Ft;

/* Force only at every St_Decay Time interval */

  if(g_stepNumber%Ft->StirFreq == 0 ? 1:0){
    ForcedTurb_OUNoiseUpdate(Ft->OUPhases, 6*(Ft->NModes), Ft->OUVar,
                             Ft->StirFreq*dt, Ft->StirDecay);
    ForcedTurb_CalcPhases(Ft);
    ForcedTurb_ComputeAcceleration(Ft, grid);
  }
#endif

/* --------------------------------------------------------
   2. Update conservative solution array with hyperbolic 
      terms only.
   -------------------------------------------------------- */

  /* -- 2a. Compute current for Hall MHD -- */
  
  #if (HALL_MHD == EXPLICIT)
  #ifdef STAGGERED_MHD
    #error HALL_MHD not compatible with CT scheme
  #endif
  GetCurrent (d, grid);
  #endif

  for (dir = beg_dir; dir <= end_dir; dir++){

    g_dir = dir;  

  /* -- 2b. Set integration box for current update -- */

    RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
    RBoxSetDirections (&sweepBox, g_dir);
    SetVectorIndices (g_dir);

    #if (defined STAGGERED_MHD)
    D_EXPAND(                                      ;  ,
             (*sweepBox.tbeg)--; (*sweepBox.tend)++;  ,
             (*sweepBox.bbeg)--; (*sweepBox.bend)++;)
    #endif
    ResetState(d, &sweep, grid);

    int ntot = grid->np_tot[g_dir];
    int nbeg = *sweepBox.nbeg;
    int nend = *sweepBox.nend;
  
    BOX_TRANSVERSE_LOOP(&sweepBox, k,j,i){
      ip  = sweepBox.n;
      g_i = i;  g_j = j;  g_k = k;
      for ((*ip) = 0; (*ip) < ntot; (*ip)++) {
        NVAR_LOOP(nv) stateC->v[*ip][nv] = d->Vc[nv][k][j][i];
        sweep.flag[*ip] = d->flag[k][j][i];
        #ifdef STAGGERED_MHD
        sweep.bn[*ip] = d->Vs[g_dir][k][j][i];
        #endif
      }

      #if (HALL_MHD == EXPLICIT)
      double ***Jx = d->J[IDIR];
      double ***Jy = d->J[JDIR];
      double ***Jz = d->J[KDIR];
      for ((*ip) = 0; (*ip) < ntot-1; (*ip)++) {

        if (g_dir == IDIR){  
          stateL->J[*ip][IDIR] = AVERAGE_XYZ(Jx,k-1,j-1,i);
          stateL->J[*ip][JDIR] = AVERAGE_Z(Jy,k-1,j,i);
          stateL->J[*ip][KDIR] = AVERAGE_Y(Jz,k,j-1,i);
        }else if (g_dir == JDIR){  
          stateL->J[*ip][IDIR] = AVERAGE_Z(Jx,k-1,j,i);
          stateL->J[*ip][JDIR] = AVERAGE_XYZ(Jy,k-1,j,i-1);
          stateL->J[*ip][KDIR] = AVERAGE_X(Jz,k,j,i-1);
        }else if (g_dir == KDIR){  
          stateL->J[*ip][IDIR] = AVERAGE_Y(Jx,k,j-1,i);
          stateL->J[*ip][JDIR] = AVERAGE_X(Jy,k,j,i-1);
          stateL->J[*ip][KDIR] = AVERAGE_XYZ(Jz,k,j-1,i-1);
        }
      }
      #endif

      #if (defined PARTICLES) && (PARTICLES_TYPE == COSMIC_RAYS)
      Particles_CR_States1DCopy(d, &sweep, 1, ntot-2);
      #endif
      
      CheckNaN (stateC->v, 0, ntot-1,0);
      States  (&sweep, nbeg - 1, nend + 1, grid);
      Riemann (&sweep, nbeg - 1, nend, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
      CT_StoreUpwindEMF (&sweep, d->emf, nbeg - 1, nend, grid);
      #endif

      #if UPDATE_VECTOR_POTENTIAL == YES
      VectorPotentialUpdate (d, NULL, &sweep, grid);
      #endif
      #ifdef SHEARINGBOX
      SB_SaveFluxes (&sweep, grid);
      #endif
      RightHandSide (&sweep, Dts, nbeg, nend, dt, grid);
      #if FORCED_TURB == YES
      if (g_stepNumber%Ft->StirFreq == 0 ? 1:0){
        ForcedTurb_CorrectRHS(d, &sweep, nbeg, nend, dt,  grid);
      }  
      #endif

    /* -- Update:  U = U + dt*R -- */

      for ((*ip) = nbeg; (*ip) <= nend; (*ip)++) { 
        NVAR_LOOP(nv) UU[k][j][i][nv] += sweep.rhs[*ip][nv];
      }
      #ifdef CHOMBO
      for ((*ip) = nbeg-1; (*ip) <= nend; (*ip)++){
        sweep.flux[*ip][MXn] += sweep.press[*ip];
        #if HAVE_ENERGY && ENTROPY_SWITCH
        sweep.flux[*ip][ENTR] = 0.0;
        #endif
      }   
      StoreAMRFlux (sweep.flux, aflux, 0, 0, NVAR-1, nbeg-1, nend, grid);
      #endif 

    /* -- Compute inverse hyperbolic time step - */

      #if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
      if (g_intStage == 1){
        inv_dl = GetInverse_dl(grid);
        for ((*ip) = nbeg; (*ip) <= nend; (*ip)++) { 
          C_dt[k][j][i] += 0.5*(Dts->cmax[(*ip)-1] + Dts->cmax[*ip])*inv_dl[*ip];
        }
      }
      #else
      inv_dl = GetInverse_dl(grid);
      for ((*ip) = nbeg-1; (*ip) <= nend; (*ip)++) { 
        Dts->invDt_hyp = MAX(Dts->invDt_hyp, Dts->cmax[*ip]*inv_dl[*ip]);
      }
      #endif
    }
  }

/* --------------------------------------------------------
   3. Compute (hyperbolic) emf
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD
  CT_ComputeEMF(d,grid);
#endif

/* --------------------------------------------------------
   4. Correct fluxes for shearingbox 
   -------------------------------------------------------- */

#ifdef SHEARINGBOX
  SB_CorrectFluxes (UU, 0.0, dt, grid);
  #ifdef STAGGERED_MHD
  SB_CorrectEMF(d->emf, d->Vs, grid);
  #endif
#endif

/* ----------------------------------------------------------
   5. Update solution array with parabolic (diffusion) terms.
   ---------------------------------------------------------- */

#if (PARABOLIC_FLUX & EXPLICIT)
  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
  ParabolicUpdate (d, UU, &sweepBox, aflux, dt, Dts, grid);
  #if (defined STAGGERED_MHD) && (RESISTIVITY == EXPLICIT)
  CT_ResistiveEMF (d, 1, grid);
  #endif
#endif

/* --------------------------------------------------------
   6. Update staggered magnetic field with total (hyp+par)
      emf. Note that this call cannot be moved before to
      ParabolicUpdate() which may compute J based on the
      current value of d->Vs.
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD
  CT_Update(d, d->Vs, dt, grid);
#endif

/* ---------------------------------------------------------------
   7. Update solution array with particle feedback
      (note that, at corrector, d->Fcr is computed from
       total momentum variation in the particle pusher).
   --------------------------------------------------------------- */

#ifdef PARTICLES
  #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
  RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
  Particles_CR_ConservativeFeedback (UU, d->Fcr, dt, &sweepBox);
  #endif
  #if (PARTICLES_TYPE == DUST) && (PARTICLES_DUST_FEEDBACK == YES)
  RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
  Particles_Dust_ConservativeFeedback (UU, d->Fdust, dt, &sweepBox);
  #endif
#endif

/* -------------------------------------------------------------------
   8. Reduce dt for dimensionally unsplit schemes.
   ------------------------------------------------------------------- */

#if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
  if (g_intStage == 1){
    DOM_LOOP(k,j,i) Dts->invDt_hyp = MAX(Dts->invDt_hyp, C_dt[k][j][i]);
    Dts->invDt_hyp /= (double)DIMENSIONS;
  }
#endif
}
