/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Store or retrieve the Electromotive Force (EMF).              

  This file provides a database functionality for storing or 
  retrieving EMF components and related information at different 
  points and times in the code.

  The CT_StoreEMF() function is called immediately after a 1D Riemann 
  solver during the hydro sweeps in order to save Fluxes and 
  characteristic signal velocities into the emf structure for later 
  reuse. 
  The fluxes coming from different sweeps are the different components
  of the advective part (-v X B) part of the electric field.

  The function CT_GetEMF() is used to obtain the edge-centered electric 
  field by properly averaging the EMF components previously stored
  at the zone faces during the 1D sweeps.

  \author  A. Mignone (mignone@ph.unito.it)
  \date    Aug 27, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static EMF emf;

#define eps_UCT_CONTACT   1.e-6
#define EX(k,j,i)  (vz[k][j][i]*By[k][j][i] - vy[k][j][i]*Bz[k][j][i])
#define EY(k,j,i)  (vx[k][j][i]*Bz[k][j][i] - vz[k][j][i]*Bx[k][j][i])
#define EZ(k,j,i)  (vy[k][j][i]*Bx[k][j][i] - vx[k][j][i]*By[k][j][i])

/* ********************************************************************* */
void CT_StoreEMF (const State_1D *state, int beg, int end, Grid *grid)
/*!
 * Store EMF components and related information available 
 * during 1D sweeps.
 *
 * \param [in]      state pointer to State_1D structure
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 * \param [in]      grid  pointer to Grid structure;
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int i, j, k, s;

/* ----------------------------------------------------
     Allocate memory for EMF structure and 
     check for incompatible combinations of algorithms 
   ---------------------------------------------------- */

  if (emf.ez == NULL){

    emf.ibeg = emf.iend = 0;
    emf.jbeg = emf.jend = 0;
    emf.kbeg = emf.kend = 0;

  /* -- memory allocation -- */

    D_EXPAND(                                          ;  ,
      emf.ez = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      emf.ex = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
      emf.ey = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )

    #if CT_EMF_AVERAGE == UCT_CONTACT
     D_EXPAND(
       emf.svx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
       emf.svy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
       emf.svz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);
     )
    #endif

     D_EXPAND(                                           ;  ,
       emf.ezi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); 
       emf.ezj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.exj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.exk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

       emf.eyi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.eyk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     )

    #if CT_EMF_AVERAGE == UCT_HLL

     D_EXPAND(
       emf.SxL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SxR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.SyL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SyR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.SzL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SzR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     )
    #endif
  }

/* ------------------------------------------------------
     Store emf components or other necessary 1-D data
   ------------------------------------------------------ */

  if (g_dir == IDIR){

    emf.ibeg = beg; emf.iend = end;
    for (i = beg; i <= end; i++) { 

      D_EXPAND(emf.ezi[g_k][g_j][i] = -state->flux[i][BX2];  ,
                                                              ,
               emf.eyi[g_k][g_j][i] =  state->flux[i][BX3]; ) 

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[i][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[i][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;

       emf.svx[g_k][g_j][i] = s;
      #endif 

      #if CT_EMF_AVERAGE == UCT_HLL
       emf.SxL[g_k][g_j][i] = MAX(0.0, -state->SL[i]); 
       emf.SxR[g_k][g_j][i] = MAX(0.0,  state->SR[i]); 
      #endif

      #if CT_EMF_AVERAGE == RIEMANN_2D
       emf.ezi[g_k][g_j][i] = -2.0*(state->pnt_flx[i][BX2] + 
                                    state->dff_flx[i][BX2]); 
      #endif

    }

  }else if (g_dir == JDIR){

    emf.jbeg = beg; emf.jend = end;
    for (j = beg; j <= end; j++) {

      D_EXPAND(                                           ;   ,
               emf.ezj[g_k][j][g_i] =  state->flux[j][BX1];   ,
               emf.exj[g_k][j][g_i] = -state->flux[j][BX3]; )

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[j][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[j][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;
       emf.svy[g_k][j][g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
       emf.SyL[g_k][j][g_i] = MAX(0.0, -state->SL[j]); 
       emf.SyR[g_k][j][g_i] = MAX(0.0,  state->SR[j]); 
      #endif

      #if CT_EMF_AVERAGE == RIEMANN_2D 
//       emf.ezj[g_k][j][g_i] += state->dff_flx[j][BX1];  
       emf.ezj[g_k][j][g_i]   = 2.0*(state->pnt_flx[j][BX1] + 
                                     state->dff_flx[j][BX1]); 
      #endif
    }

  }else if (g_dir == KDIR){

    emf.kbeg = beg; emf.kend = end;
    for (k = beg; k <= end; k++) {

      emf.eyk[k][g_j][g_i] = -state->flux[k][BX1]; 
      emf.exk[k][g_j][g_i] =  state->flux[k][BX2]; 

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[k][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[k][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;
       emf.svz[k][g_j][g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
       emf.SzL[k][g_j][g_i] = MAX(0.0, -state->SL[k]); 
       emf.SzR[k][g_j][g_i] = MAX(0.0,  state->SR[k]); 
      #endif

    }
  }

/* ------------------------------------------------------
         Store velocity slopes if necessary 
   ------------------------------------------------------ */

  #if CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU 
    if (g_intStage == 2) return;    
   #endif

   /* -- "end+1" needed to save dvx_dx -- */

   CT_StoreVelSlopes (&emf, state, beg, end + 1); 

  #endif

}

/* ********************************************************************* */
EMF *CT_GetEMF (const Data *d, Grid *grid)
/*!
 * Retrieve EMF by suitable average of 1D face-centered fluxes.
 * 
 * \param [in] d
 * \param [in] grid
 *
 * \return a pointer to an edge-centered EMF.
 *********************************************************************** */
{
  int    i, j, k;
  double ***vx, ***vy, ***vz;
  double ***Bx, ***By, ***Bz;

/* -----------------------------------------------------
    Return only the resistive part of emf when using 
    super time stepping. 
   ----------------------------------------------------- */

  #if RESISTIVITY == SUPER_TIME_STEPPING
   if (g_operatorStep == PARABOLIC_STEP){
     TOT_LOOP(k,j,i) {
       #if DIMENSIONS == 3
        emf.ex[k][j][i] = emf.ey[k][j][i] = 0.0;
       #endif
       emf.ez[k][j][i] = 0.0;
     }
     CT_AddResistiveEMF(d, grid); 
     return (&emf);
   }
  #endif

/* -------------------------------------
       set boundary conditions on 
       face-centered electric fields
   ------------------------------------- */
/*
  #ifdef CTU
   if (g_intStage == 2)
  #endif
  EMF_BOUNDARY (&emf, grid);
*/

/* ------------------------------------------------------
       Compute slopes of staggered magnetic fields 
   ------------------------------------------------------ */

  #if CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU
    if (g_intStage == 1)
   #endif
   CT_GetStagSlopes(d->Vs, &emf, grid);
  #endif

/* -----------------------------------------------------
                 Select average 
   ----------------------------------------------------- */

  #if CT_EMF_AVERAGE == ARITHMETIC || CT_EMF_AVERAGE == RIEMANN_2D

   CT_EMF_ArithmeticAverage (&emf, 0.25);

  #elif CT_EMF_AVERAGE == UCT_CONTACT

   CT_EMF_ArithmeticAverage (&emf, 1.0);
   CT_EMF_IntegrateToCorner (d, &emf, grid);
   for (k = emf.kbeg; k <= emf.kend; k++){
   for (j = emf.jbeg; j <= emf.jend; j++){
   for (i = emf.ibeg; i <= emf.iend; i++){      
     #if DIMENSIONS == 3
      emf.ex[k][j][i] *= 0.25;
      emf.ey[k][j][i] *= 0.25;
     #endif
     emf.ez[k][j][i] *= 0.25;
   }}}

  #elif CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU
    if (g_intStage == 1) CT_EMF_CMUSCL_Average (d, &emf, grid);
    else   
   #endif
   CT_EMF_HLL_Solver (d, &emf, grid);

  #elif CT_EMF_AVERAGE == UCT0

/* -- Subtract cell-centered contribution -- */

   EXPAND(vx = d->Vc[VX1]; Bx = d->Vc[BX1];  ,
          vy = d->Vc[VX2]; By = d->Vc[BX2];  ,
          vz = d->Vc[VX3]; Bz = d->Vc[BX3];)

   for (k = emf.kbeg; k <= emf.kend + KOFFSET; k++){
   for (j = emf.jbeg; j <= emf.jend + JOFFSET; j++){
   for (i = emf.ibeg; i <= emf.iend + IOFFSET; i++){       
     #if DIMENSIONS == 3
      emf.exj[k][j][i] *= 2.0;
      emf.exk[k][j][i] *= 2.0;
      emf.eyi[k][j][i] *= 2.0;
      emf.eyk[k][j][i] *= 2.0;

      emf.exj[k][j][i] -= 0.5*(EX(k,j,i) + EX(k,j+1,i));
      emf.exk[k][j][i] -= 0.5*(EX(k,j,i) + EX(k+1,j,i));

      emf.eyi[k][j][i] -= 0.5*(EY(k,j,i) + EY(k,j,i+1));
      emf.eyk[k][j][i] -= 0.5*(EY(k,j,i) + EY(k+1,j,i));
     #endif
     emf.ezi[k][j][i] *= 2.0;
     emf.ezj[k][j][i] *= 2.0;
     emf.ezi[k][j][i] -= 0.5*(EZ(k,j,i) + EZ(k,j,i+1));
     emf.ezj[k][j][i] -= 0.5*(EZ(k,j,i) + EZ(k,j+1,i));
   }}}

   CT_EMF_ArithmeticAverage (&emf, 0.25);

  #else
    print1 ("! CT_GetEMF: unknown EMF average.\n");
    QUIT_PLUTO(1);
  #endif


/* ------------------------------------------------------
    Add contributions from resistive terms with explicit
    time stepping.
   ------------------------------------------------------ */

  #if RESISTIVITY == EXPLICIT
   CT_AddResistiveEMF(d, grid); 
  #endif
   
/* -------------------------------------------------------------
    Fine Tuning: EMF_USERDEF_BOUNDARY can be used to directly
    set the edge-centered electric field
   ------------------------------------------------------------- */
/*
  #ifdef CTU
   if (step == 2)
  #endif
  {
    int lside[3] = {X1_BEG, X2_BEG, X3_BEG};
    int rside[3] = {X1_END, X2_END, X3_END};
    int dir;

    for (dir = 0; dir < DIMENSIONS; dir++){
      if (grid[dir].lbound == USERDEF)
        EMF_USERDEF_BOUNDARY (&emf, lside[dir], EDGE_EMF, grid);  
      if (grid[dir].rbound == USERDEF)
        EMF_USERDEF_BOUNDARY (&emf, rside[dir], EDGE_EMF, grid);  
    }   
  }
*/

  return (&emf);
}

#if RESISTIVITY != NO
/* ********************************************************************* */
void CT_AddResistiveEMF (const Data *d, Grid *grid)
/*!
 * Add resistive terms to EMF.
 * 
 * \param [in] d
 * \param [in] grid
 *
 * \return a pointer to an edge-centered EMF.
 *********************************************************************** */
{
  int    i, j, k;
  double ***vx, ***vy, ***vz;
  double ***Bx, ***By, ***Bz;
  Data_Arr eta;
/*
  #ifdef CTU
    if (g_intStage == 2) return (&emf);  
  #endif
*/
  eta = GetStaggeredEta();

  for (k = emf.kbeg; k <= emf.kend; k++){
  for (j = emf.jbeg; j <= emf.jend; j++){
  for (i = emf.ibeg; i <= emf.iend; i++){

/*  -------------------------------------------------------------
     Simple implementation in Cartesian coordinates by computing 
     currents directly here (bypass GetCurrent). 
    -------------------------------------------------------------
    
double ***Bx = d->Vs[BX1s], Jx;
double ***By = d->Vs[BX2s], Jy;
double ***Bz = d->Vs[BX3s], Jz;
double dx = grid[IDIR].dx[i];
double dy = grid[JDIR].dx[j];
double dz = grid[KDIR].dx[k];

#if DIMENSIONS == 3
Jx = (Bz[k][j+1][i] - Bz[k][j][i])/dy - (By[k+1][j][i] - By[k][j][i])/dz;
Jy = (Bx[k+1][j][i] - Bx[k][j][i])/dz - (Bz[k][j][i+1] - Bz[k][j][i])/dx;
emf.ex[k][j][i] += Jx*g_inputParam[ETAX];
emf.ey[k][j][i] += Jy*g_inputParam[ETAY];
#endif
Jz = (By[k][j][i+1] - By[k][j][i])/dx - (Bx[k][j+1][i] - Bx[k][j][i])/dy;
emf.ez[k][j][i] += Jz*g_inputParam[ETAZ];
*/
     #if DIMENSIONS == 3
      emf.ex[k][j][i] += eta[IDIR][k][j][i]*d->J[IDIR][k][j][i];
      emf.ey[k][j][i] += eta[JDIR][k][j][i]*d->J[JDIR][k][j][i];
     #endif 
     emf.ez[k][j][i] += eta[KDIR][k][j][i]*d->J[KDIR][k][j][i];

  }}}
}
#endif /* RESISTIVITY != NO */

#undef EX
#undef EY
#undef EZ
#undef eps_UCT_CONTACT
