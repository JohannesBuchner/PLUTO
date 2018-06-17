/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute rhs for thermal conduction 

  Compute the one-dimensional right hand side for the
  thermal conduction operator in the direction given by ::g_dir.

  \authors A. Mignone (mignone@ph.unito.it)\n

 \b References

  \date   May 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void TC_RHS (const Data *d, Data_Arr dU, double *dcoeff,
             double **aflux, double dt, int beg, int end, Grid *grid)
/*!
 * \param [in]   d           pointer to PLUTO Data structure
 * \param [out]  dU          a 4D array containing conservative variables
 *                           increment
 * \param [out]  dcoeff      1D array of diffusion coefficients   
 * \param [out]  aflux       pointer to 2D array for AMR re-fluxing
 *                           operations
 * \param [in]   dt          the current time-step                            
 * \param [in]   beg,end     initial and final interface indices
 * \param [in]   grid        pointer to Grid structure.
 *
 *********************************************************************** */
{
  int i = g_i;
  int j = g_j;
  int k = g_k;
  int nv;
  double dtdV, dtdx;
  static  double *fA;
  static  Sweep sweep;

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */
  
  if (sweep.vn == NULL) {
    MakeState (&sweep);
    fA = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
   1. Compute TC flux
   -------------------------------------------------------- */
 
  if (g_dir == IDIR) {
    ITOT_LOOP (i) NVAR_LOOP(nv) sweep.vn[i][nv] = d->Vc[nv][k][j][i];
  } else if (g_dir == JDIR) {
    JTOT_LOOP (j) NVAR_LOOP(nv) sweep.vn[j][nv] = d->Vc[nv][k][j][i];
  } else if (g_dir == KDIR) {
    KTOT_LOOP (k) NVAR_LOOP(nv) sweep.vn[k][nv] = d->Vc[nv][k][j][i];
  }
  TC_Flux (d->Tc, &sweep, dcoeff, beg-1, end, grid);

/* --------------------------------------------------------
   2. Multiply flux X area & compute rhs
   -------------------------------------------------------- */

  if (g_dir == IDIR){
    #if GEOMETRY != CARTESIAN
    for (i = beg-1; i <= end; i++){
      fA[i] = sweep.tc_flux[i][ENG]*grid->A[IDIR][k][j][i];
    }  
    #endif
    for (i = beg; i <= end; i++){
      #if GEOMETRY == CARTESIAN
      dtdx = dt/grid->dx[IDIR][i];
      dU[k][j][i][ENG] += dtdx*(sweep.tc_flux[i][ENG] - sweep.tc_flux[i-1][ENG]);
      #else
      dtdV = dt/grid->dV[k][j][i];
      dU[k][j][i][ENG] += dtdV*(fA[i] - fA[i-1]);
      #endif
    }    
  } else if (g_dir == JDIR){
    #if GEOMETRY != CARTESIAN
    for (j = beg-1; j <= end; j++){
      fA[j] = sweep.tc_flux[j][ENG]*grid->A[JDIR][k][j][i];
    }  
    #endif
    for (j = beg; j <= end; j++){
      #if GEOMETRY == CARTESIAN
      dtdx = dt/grid->dx[JDIR][j];
      dU[k][j][i][ENG] += dtdx*(sweep.tc_flux[j][ENG] - sweep.tc_flux[j-1][ENG]);
      #else
      dtdV = dt/grid->dV[k][j][i];
      dU[k][j][i][ENG] += dtdV*(fA[j] - fA[j-1]);
      #endif
    }    
  } else if (g_dir == KDIR){
    #if GEOMETRY != CARTESIAN
    for (k = beg-1; k <= end; k++){
      fA[k] = sweep.tc_flux[k][ENG]*grid->A[KDIR][k][j][i];
    }  
    #endif
    for (k = beg; k <= end; k++){
      #if GEOMETRY == CARTESIAN
      dtdx = dt/grid->dx[KDIR][k];
      dU[k][j][i][ENG] += dtdx*(sweep.tc_flux[k][ENG] - sweep.tc_flux[k-1][ENG]);
      #else
      dtdV = dt/grid->dV[k][j][i];
      dU[k][j][i][ENG] += dtdV*(fA[k] - fA[k-1]);
      #endif
    }    
  }   
    
  #ifdef CHOMBO
  StoreAMRFlux (sweep.tc_flux, aflux, -1, ENG, ENG, beg-1, end, grid);
  #endif
}
