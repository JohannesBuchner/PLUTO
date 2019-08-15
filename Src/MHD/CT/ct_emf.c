/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Store or retrieve the Electromotive Force (EMF).              

  This file provides a database functionality for storing or 
  retrieving EMF components and related information at different 
  points and times in the code.

  The CT_Allocate() allocates memory for EMF structure members.
  The CT_StoreEMF() function is called immediately after a 1D Riemann 
  solver during the hydro sweeps in order to save Fluxes at cell
  interfaces and characteristic signal velocities into the emf
  structure for later reuse. 
  The fluxes coming from different sweeps are the different components
  of the advective part (-v X B) part of the electric field.

  The function CT_GetEMF() is used to obtain the edge-centered electric 
  field by properly averaging the EMF components previously stored
  at the zone faces during the 1D sweeps.

  \author  A. Mignone (mignone@ph.unito.it)
  \date    Sep 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_Allocate (EMF *emf0)
/*!
 *  Allocate memory for EMF structure and 
 *  check for incompatible combinations of algorithms 
 *
 *
 *********************************************************************** */
{
  emf0->ibeg = emf0->iend = 0;
  emf0->jbeg = emf0->jend = 0;
  emf0->kbeg = emf0->kend = 0;

  /* -- memory allocation -- */

  D_EXPAND(                                                      ;  ,
           emf0->ez = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
           emf0->ex = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->ey = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);)

#if CT_EMF_AVERAGE == UCT_CONTACT
  D_EXPAND(emf0->svx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
           emf0->svy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
           emf0->svz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);)
#endif
  D_EXPAND(                                                       ;  ,
           emf0->ezi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); 
           emf0->ezj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

           emf0->exj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->exk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

           emf0->eyi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->eyk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);)

#if CT_EMF_AVERAGE == UCT_HLL
  D_EXPAND(emf0->SxL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->SxR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

           emf0->SyL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->SyR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

           emf0->SzL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
           emf0->SzR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);)

             
  emf0->dvx_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dvx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  emf0->dvy_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dvy_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  #if DIMENSIONS == 3
  emf0->dvx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dvy_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  emf0->dvz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dvz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dvz_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif
    
  emf0->dbx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dby_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  #if DIMENSIONS == 3
  emf0->dbx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dby_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dbz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  emf0->dbz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

#endif

}

#define eps_UCT_CONTACT   1.e-6
/* ********************************************************************* */
void CT_StoreUpwindEMF (const Sweep *sweep, EMF *emf, int beg, int end,
                        Grid *grid)
/*!
 * Store EMF components and related information available 
 * during 1D sweeps.
 *
 * \param [in]  sweep  pointer to Sweep structure
 * \param [in]  beg    initial index of computation 
 * \param [in]  end    final   index of computation
 * \param [in]  grid   pointer to Grid structure;
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int i, j, k, s;

#if (HALL_MHD) && (CT_EMF_AVERAGE != ARITHMETIC)
  #error HALL_MHD should be used with ARITHMETIC average.
#endif

/* ------------------------------------------------------
     Store emf components or other necessary 1-D data
   ------------------------------------------------------ */

  if (g_dir == IDIR){

    emf->ibeg = beg; emf->iend = end;
    for (i = beg; i <= end; i++) { 
      D_EXPAND(emf->ezi[g_k][g_j][i] = -sweep->flux[i][BX2];  ,
                                                              ,
               emf->eyi[g_k][g_j][i] =  sweep->flux[i][BX3]; ) 

      #if CT_EMF_AVERAGE == UCT_CONTACT
      if      (sweep->flux[i][RHO] >  eps_UCT_CONTACT) s = 1;
      else if (sweep->flux[i][RHO] < -eps_UCT_CONTACT) s = -1;
      else s = 0;

      emf->svx[g_k][g_j][i] = s;
      #endif 

      #if CT_EMF_AVERAGE == UCT_HLL
      emf->SxL[g_k][g_j][i] = MAX(0.0, -sweep->SL[i]); 
      emf->SxR[g_k][g_j][i] = MAX(0.0,  sweep->SR[i]); 
      #endif
    }

  }else if (g_dir == JDIR){

    emf->jbeg = beg; emf->jend = end;
    for (j = beg; j <= end; j++) {
      D_EXPAND(                                           ;   ,
               emf->ezj[g_k][j][g_i] =  sweep->flux[j][BX1];   ,
               emf->exj[g_k][j][g_i] = -sweep->flux[j][BX3]; )

      #if CT_EMF_AVERAGE == UCT_CONTACT
      if      (sweep->flux[j][RHO] >  eps_UCT_CONTACT) s = 1;
      else if (sweep->flux[j][RHO] < -eps_UCT_CONTACT) s = -1;
      else s = 0;
      emf->svy[g_k][j][g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
      emf->SyL[g_k][j][g_i] = MAX(0.0, -sweep->SL[j]); 
      emf->SyR[g_k][j][g_i] = MAX(0.0,  sweep->SR[j]); 
      #endif
    }

  }else if (g_dir == KDIR){

    emf->kbeg = beg; emf->kend = end;
    for (k = beg; k <= end; k++) {
      emf->eyk[k][g_j][g_i] = -sweep->flux[k][BX1]; 
      emf->exk[k][g_j][g_i] =  sweep->flux[k][BX2]; 

      #if CT_EMF_AVERAGE == UCT_CONTACT
      if      (sweep->flux[k][RHO] >  eps_UCT_CONTACT) s = 1;
      else if (sweep->flux[k][RHO] < -eps_UCT_CONTACT) s = -1;
      else s = 0;
      emf->svz[k][g_j][g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
      emf->SzL[k][g_j][g_i] = MAX(0.0, -sweep->SL[k]); 
      emf->SzR[k][g_j][g_i] = MAX(0.0,  sweep->SR[k]); 
      #endif
    }
  }

/* ------------------------------------------------------
         Store velocity slopes if necessary 
   ------------------------------------------------------ */

#if CT_EMF_AVERAGE == UCT_HLL
  #ifdef CTU 
  print ("! UCT_HLL average not compatible with CTU schemes (stencil too small)\n");
  QUIT_PLUTO(1);
  #endif

  /* -- "end+1" needed to save dvx_dx -- */

  CT_StoreVelSlopes (emf, sweep, beg, end + 1); 
#endif

}

/* ********************************************************************* */
void CT_ComputeEMF (const Data *data, Grid *grid)
/*!
 *  Compute the total electromotive force (EMF) using the 
 *  hyperbolic fluxes stored during the most recent upwind step.
 * 
 * \param [in] d
 * \param [in] grid
 *
 *********************************************************************** */
{
  int  i, j, k;
  EMF  *emf = data->emf;
  
/* --------------------------------------------------------
   0. Check if cell-centered emf is needed.
      This may depend on the scheme.
   -------------------------------------------------------- */
  
#if    CT_EMF_AVERAGE == UCT_CONTACT \
    || CT_EMF_AVERAGE == UCT0 \
    || ((defined PARTICLES) && (PARTICLES_TYPE == COSMIC_RAYS))
  CT_ComputeCenterEMF (data);
#endif
  
/* --------------------------------------------------------
   1. Select averaging scheme
   -------------------------------------------------------- */

#if CT_EMF_AVERAGE == ARITHMETIC

  CT_EMF_ArithmeticAverage (emf, 0.25);

#elif CT_EMF_AVERAGE == UCT_CONTACT

  CT_EMF_ArithmeticAverage (emf, 1.0);
  CT_EMF_IntegrateToCorner (data, emf, grid);
  for (k = emf->kbeg; k <= emf->kend; k++){
  for (j = emf->jbeg; j <= emf->jend; j++){
  for (i = emf->ibeg; i <= emf->iend; i++){      
    #if DIMENSIONS == 3
    emf->ex[k][j][i] *= 0.25;
    emf->ey[k][j][i] *= 0.25;
    #endif
    emf->ez[k][j][i] *= 0.25;
  }}}

#elif CT_EMF_AVERAGE == UCT_HLL

  CT_GetStagSlopes(data->Vs, emf, grid);
  CT_EMF_HLL_Solver (data, emf, grid);

#elif CT_EMF_AVERAGE == UCT0

  for (k = emf->kbeg; k <= emf->kend + KOFFSET; k++){
  for (j = emf->jbeg; j <= emf->jend + JOFFSET; j++){
  for (i = emf->ibeg; i <= emf->iend + IOFFSET; i++){       
  #if DIMENSIONS == 3
    emf->exj[k][j][i] *= 2.0;
    emf->exk[k][j][i] *= 2.0;
    emf->eyi[k][j][i] *= 2.0;
    emf->eyk[k][j][i] *= 2.0;

    emf->exj[k][j][i] -= 0.5*(data->Ex1[k][j][i] + data->Ex1[k][j+1][i]);
    emf->exk[k][j][i] -= 0.5*(data->Ex1[k][j][i] + data->Ex1[k+1][j][i]);

    emf->eyi[k][j][i] -= 0.5*(data->Ex2[k][j][i] + data->Ex2[k][j][i+1]);
    emf->eyk[k][j][i] -= 0.5*(data->Ex2[k][j][i] + data->Ex2[k+1][j][i]);
   #endif
    emf->ezi[k][j][i] *= 2.0;
    emf->ezj[k][j][i] *= 2.0;
    emf->ezi[k][j][i] -= 0.5*(data->Ex3[k][j][i] + data->Ex3[k][j][i+1]);
    emf->ezj[k][j][i] -= 0.5*(data->Ex3[k][j][i] + data->Ex3[k][j+1][i]);
  }}}

  CT_EMF_ArithmeticAverage (emf, 0.25);

#else
  print ("! CT_ComputeEMF: unknown EMF average.\n");
  QUIT_PLUTO(1);
#endif
  
}

#if RESISTIVITY != NO
/* ********************************************************************* */
void CT_ResistiveEMF (const Data *data, int op, Grid *grid)
/*!
 * Compute/Add resistive terms to EMF (used during parabolic update
 * or STS operator splitting)
 * 
 * \param [in] d
 * \param [in] op     operation: op = 0 means initialize,
 *                               op = 1 means add.                                
 * \param [in] grid
 *
 * \return a pointer to an edge-centered EMF.
 *********************************************************************** */
{
  int    i, j, k;
  double ***vx, ***vy, ***vz;
  double ***Bx, ***By, ***Bz;
  Data_Arr eta;
  EMF  *emf = data->emf; 

  eta = GetStaggeredEta();

  if (op == 0){
    for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){
      #if DIMENSIONS == 3
      emf->ex[k][j][i] = eta[IDIR][k][j][i]*data->J[IDIR][k][j][i];
      emf->ey[k][j][i] = eta[JDIR][k][j][i]*data->J[JDIR][k][j][i];
      #endif 
      emf->ez[k][j][i] = eta[KDIR][k][j][i]*data->J[KDIR][k][j][i];

    }}}
  }else if (op == 1){
    for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){
      #if DIMENSIONS == 3
      emf->ex[k][j][i] += eta[IDIR][k][j][i]*data->J[IDIR][k][j][i];
      emf->ey[k][j][i] += eta[JDIR][k][j][i]*data->J[JDIR][k][j][i];
      #endif 
      emf->ez[k][j][i] += eta[KDIR][k][j][i]*data->J[KDIR][k][j][i];
    }}}
   
  }else{
    print ("! CT_ResistiveEFM(): invalid op\n");
    QUIT_PLUTO(1); 
  }
}  
#endif /* RESISTIVITY != NO */

#undef eps_UCT_CONTACT

/* ********************************************************************* */
void CT_ComputeCenterEMF(const Data *data)
/*!
 *  Compute cell-center inductive electric field.
 *********************************************************************** */
{
  int i,j,k;
  double vx1, vx2, vx3;
  double Bx1, Bx2, Bx3;
  double qg, scrh;

  vx1 = vx2 = vx3 = 0.0;
  Bx1 = Bx2 = Bx3 = 0.0;

  TOT_LOOP(k,j,i){
    EXPAND(vx1 = data->Vc[VX1][k][j][i];  ,
           vx2 = data->Vc[VX2][k][j][i];  ,
           vx3 = data->Vc[VX3][k][j][i];)

    EXPAND(Bx1 = data->Vc[BX1][k][j][i];  ,
           Bx2 = data->Vc[BX2][k][j][i];  ,
           Bx3 = data->Vc[BX3][k][j][i];)

  /* -- Compute inductive electric field -- */

    data->Ex1[k][j][i] = (vx3*Bx2 - vx2*Bx3);
    data->Ex2[k][j][i] = (vx1*Bx3 - vx3*Bx1);
    data->Ex3[k][j][i] = (vx2*Bx1 - vx1*Bx2);

  /* -- Add CR Hall term  -- */

    #ifdef PARTICLES
    #if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_CR_FEEDBACK == YES)
    qg   = data->Vc[RHO][k][j][i]*PARTICLES_CR_E_MC_GAS;
    scrh = 1.0/qg;
    data->Ex1[k][j][i] -= data->Fcr[IDIR][k][j][i]*scrh;
    data->Ex2[k][j][i] -= data->Fcr[JDIR][k][j][i]*scrh;
    data->Ex3[k][j][i] -= data->Fcr[KDIR][k][j][i]*scrh;
    #endif
    #endif
  }
}
