#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv;  
  double inv_dl2;
  double *inv_dl;
  double ***Ch_dt, ***Cp_dt;
  static double *cmax, **dcoeff;
  static double ***T;
  static State_1D state;
  Index indx;
  
  if (state.rhs == NULL) {
    MakeState(&state);
    cmax   = ARRAY_1D(NMAX_POINT, double);
    dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
    T      = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); 
  }
  
  Ch_dt = GetUserVar("Ch_dt");
  Cp_dt = GetUserVar("Cp_dt");
  
/* ------------------------------------------------------
    compute temperature array
   ------------------------------------------------------ */
   
  TOT_LOOP(k,j,i){
    T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
  }

/* ------------------------------------------------------
        X1  sweep over computational zones
   ------------------------------------------------------ */
    
  g_dir = IDIR;
  SetIndexes (&indx, grid);
  KDOM_LOOP(k){
  JDOM_LOOP(j){
    GetInverse_dl(grid);
    
    for (i = 0; i < NX1_TOT; i++){
      for (nv = NVAR; nv--;  ) state.v[i][nv] = d->Vc[nv][k][j][i];
    }
    States (&state, IBEG-1, IEND+1, grid);
    HLL_Solver(&state, IBEG-1, IEND, cmax, grid);
    TC_Flux (T, &state, dcoeff, IBEG-1, IEND, grid);  
    
    Ch_dt[k][j][i] = 0.5*(cmax[i] + cmax[i-1])*inv_dl[i];

    inv_dl2        = 0.5*inv_dl[i]*inv_dl[i];
    Cp_dt[k][j][i] = (dcoeff[i-1][ENG] + dcoeff[i][ENG])*inv_dl2;
  }}

/* ------------------------------------------------------
        X2  sweep over computational zones
   ------------------------------------------------------ */
    
  g_dir = JDIR;
  SetIndexes (&indx, grid);
  KDOM_LOOP(k){
  IDOM_LOOP(i){
    GetInverse_dl(grid);
    
    for (j = 0; j < NX2_TOT; j++){
      for (nv = NVAR; nv--;  ) state.v[j][nv] = d->Vc[nv][k][j][i];
    }
    States (&state, JBEG-1, JEND+1, grid);
    HLL_Solver(&state, JBEG-1, JEND, cmax, grid);
    TC_Flux (T, &state, dcoeff, JBEG-1, JEND, grid);  
    
    Ch_dt[k][j][i] += 0.5*(cmax[j] + cmax[j-1])*inv_dl[j];

    inv_dl2         = 0.5*inv_dl[j]*inv_dl[j];
    Cp_dt[k][j][i] += (dcoeff[j-1][ENG] + dcoeff[j][ENG])*inv_dl2;
  }}
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





