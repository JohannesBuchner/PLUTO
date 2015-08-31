#include"pluto.h"

static int nstep0;
static Data_Arr V0, Vs0;
static double t0, dt0;

/* ************************************************************ */
void SAVE_SOL (Data *d, Grid *grid)
/*
 *
 * Save all necessary information to re-start the 
 * computation at runtime level.
 *
 ************************************************************** */
{
  int i, j, k;

  if (V0 == NULL){
    V0 = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #ifdef STAGGERED_MHD
     Vs0 = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

  for (nv = NVAR; nv--;  ){
  DOM_LOOP(k,j,i){
    V0[nv][k][j][i] = d->Vc[nv][k][j][i];
  }}

  #ifdef STAGGERED_MHD
   for (k = KBEG-1; k <= KEND; k++){
   for (j = JBEG-1; j <= JEND; j++){
   for (i = IBEG-1; i <= IEND; i++){
     D_EXPAND(
       Vs0[BX1s][k][j][i] = d->Vs[BX1s][k][j][i];  ,
       Vs0[BX2s][k][j][i] = d->Vs[BX2s][k][j][i];  ,
       Vs0[BX3s][k][j][i] = d->Vs[BX3s][k][j][i];)
   }}}
  #endif

  nstep0 = g_stepNumber; 
  dt0    = g_dt;
  t0     = g_time;
}

/* ************************************************************ */
void GET_SOL (Data *d, Grid *grid)
/*
 *
 *  Restore previously saved solution array
 *
 ************************************************************** */

{
  for (nv = NVAR; nv--;  ){
  DOM_LOOP(k,j,i){
    d->Vc[nv][k][j][i] = V0[nv][k][j][i];
  }}

  #ifdef STAGGERED_MHD
   for (k = KBEG-1; k <= KEND; k++){
   for (j = JBEG-1; j <= JEND; j++){
   for (i = IBEG-1; i <= IEND; i++){
     D_EXPAND(
       d->Vs[BX1s][k][j][i] = Vs0[BX1s][k][j][i];  ,
       d->Vs[BX2s][k][j][i] = Vs0[BX2s][k][j][i];  ,
       d->Vs[BX3s][k][j][i] = Vs0[BX3s][k][j][i];)
   }}}
  #endif
  g_stepNumber   = nstep0;
  g_dt = dt0;
  t0      = g_time;

}

