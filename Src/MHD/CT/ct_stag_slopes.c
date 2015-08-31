#include "pluto.h"
static double MC_LIM2 (double dp, double dm);

/* ********************************************************************* */
void CT_StoreVelSlopes (EMF *emf, const State_1D *state, int beg, int end)
/*!
 *
 *
 *
 *********************************************************************** */
{
  int i=g_i, j=g_j, k=g_k;

/* ----------------------------------------------------
             Allocate static memory areas 
   ---------------------------------------------------- */

  if (emf->dvx_dx == NULL){

    emf->dvx_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dvx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    emf->dvy_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dvy_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    #if DIMENSIONS == 3
     emf->dvx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvy_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

     emf->dvz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvz_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) {
      emf->dvx_dx[k][j][i] = state->vp[i][VX1] - state->vm[i][VX1];
      emf->dvy_dx[k][j][i] = state->vp[i][VX2] - state->vm[i][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dx[k][j][i] = state->vp[i][VX3] - state->vm[i][VX3];
      #endif
    }

  }else if (g_dir == JDIR){

    for (j = beg; j <= end; j++) {
      emf->dvx_dy[k][j][i] = state->vp[j][VX1] - state->vm[j][VX1];
      emf->dvy_dy[k][j][i] = state->vp[j][VX2] - state->vm[j][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dy[k][j][i] = state->vp[j][VX3] - state->vm[j][VX3];
      #endif
    }

  }else if (g_dir == KDIR){

    for (k = beg; k <= end; k++) {
      emf->dvx_dz[k][j][i] = state->vp[k][VX1] - state->vm[k][VX1];
      emf->dvy_dz[k][j][i] = state->vp[k][VX2] - state->vm[k][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dz[k][j][i] = state->vp[k][VX3] - state->vm[k][VX3];
      #endif
    }
  }
}
/* ********************************************************************* */
void CT_GetStagSlopes (const Data_Arr b, EMF *emf, Grid *grid)
/*!
 * Compute slopes of staggered magnetic fields components:
 * dBx/dy, dBx/dz, 
 * dBy/dx, dBy/dz,
 * dBz/dx, dBz/dy,
 *
 * Exclude normal derivatives (i.e. dbx_dx).
 *
 *
 *********************************************************************** */
{
  int    i,j,k;
  double ***bx, ***by, ***bz;

  D_EXPAND(bx = b[BX1s];  ,   
           by = b[BX2s];  ,
           bz = b[BX3s];)

  if (emf->dbx_dy == NULL){
    emf->dbx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dby_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #if DIMENSIONS == 3
     emf->dbx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dby_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dbz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dbz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }
 
  for (k = KOFFSET; k < NX3_TOT - KOFFSET; k++){
  for (j = JOFFSET; j < NX2_TOT - JOFFSET; j++){
  for (i = IOFFSET; i < NX1_TOT - IOFFSET; i++){

    emf->dbx_dy[k][j][i] = MC_LIM2(bx[k][j+1][i] - bx[k][j][i], 
                                   bx[k][j][i]   - bx[k][j-1][i]);
    emf->dby_dx[k][j][i] = MC_LIM2(by[k][j][i+1] - by[k][j][i], 
                                   by[k][j][i]   - by[k][j][i-1]);
    #if DIMENSIONS == 3
     emf->dbx_dz[k][j][i] = MC_LIM2(bx[k+1][j][i] - bx[k][j][i], 
                                    bx[k][j][i]   - bx[k-1][j][i]);
     emf->dby_dz[k][j][i] = MC_LIM2(by[k+1][j][i] - by[k][j][i], 
                                    by[k][j][i]   - by[k-1][j][i]);
     emf->dbz_dx[k][j][i] = MC_LIM2(bz[k][j][i+1] - bz[k][j][i], 
                                    bz[k][j][i]   - bz[k][j][i-1]);
     emf->dbz_dy[k][j][i] = MC_LIM2(bz[k][j+1][i] - bz[k][j][i], 
                                    bz[k][j][i]   - bz[k][j-1][i]);
    #endif
  }}}
}

/* ********************************************************* */
double MC_LIM2 (double dp, double dm)
/*
 *
 *
 *
 *********************************************************** */
{
  double dc, scrh;

  if (dp*dm < 0.0) return(0.0);

  dc   = 0.5*(dp + dm);
  scrh = 2.0*(fabs(dp) < fabs(dm) ? dp:dm);
  return (fabs(dc) < fabs(scrh) ? dc:scrh);
}


