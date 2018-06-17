/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Build the right hand side for the resistivity operator

  Compute the one-dimensional right hand side for the viscous
  operator in the direction given by ::g_dir.

  \authors A. Mignone (mignone@ph.unito.it)\n

 \b References

  \date   May 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ResistiveRHS (const Data *d, Data_Arr dU, double **dcoeff,
                  double **aflux, double dt, int beg, int end, Grid *grid)
/*!
 * \param [in]   d           pointer to PLUTO Data structure
 * \param [out]  dU          a 4D array containing conservative variables
 *                           increment
 * \param [out]  dcoeff      1D array of diffusion coefficients   
 * \param [out]  aflux       pointer to 2D array for AMR re-fluxing
 *                           operations
 * \param [in]   dt          the current time-step                            
 * \param [in]  beg,end      initial and final interface indices
 * \param [in]  grid         pointer to Grid structure.
 *********************************************************************** */
{
  int i, j, k, nv;
  double *x1 = grid->x[IDIR], *x1p = grid->xr[IDIR], *dx1 = grid->dx[IDIR];
  double *x2 = grid->x[JDIR], *x2p = grid->xr[JDIR], *dx2 = grid->dx[JDIR];
  double *x3 = grid->x[KDIR], *x3p = grid->xr[KDIR], *dx3 = grid->dx[KDIR];
  double *sp = grid->sp;
  double A, dtdV, dtdl, q, rhs[NVAR];  
  static double **res_flx, **fxA;
  intList var_list;
  #if HAVE_ENERGY 
  var_list.nvar = COMPONENTS+1;
  EXPAND(var_list.indx[i=0] = BX1;  ,
         var_list.indx[++i] = BX2;  ,
         var_list.indx[++i] = BX3;)
  var_list.indx[++i] = ENG;
  #else
  var_list.nvar = COMPONENTS;
  EXPAND(var_list.indx[i=0] = BX1;  ,
         var_list.indx[++i] = BX2;  ,
         var_list.indx[++i] = BX3;)
  #endif

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */
  
  if (res_flx == NULL){
    res_flx = ARRAY_2D(NMAX_POINT, NVAR, double);
    fxA     = ARRAY_2D(NMAX_POINT, NVAR, double);
  }
  
/* --------------------------------------------------------
   1. Add resistivity flux and source terms to
      total flux (sweep->flux) and total source terms.
   -------------------------------------------------------- */

  i = g_i; j = g_j; k = g_k;
  ResistiveFlux (d->Vc, d->J, res_flx, dcoeff, beg-1, end, grid);
  
  if (g_dir == IDIR){

  /* --------------------------------------------------------
     1a. Compute fluxes & sources in the X1 direction 
     -------------------------------------------------------- */

    for (i = beg-1; i <= end; i++){
      A = grid->A[IDIR][k][j][i];
      #if GEOMETRY != SPHERICAL
      FOR_EACH(nv, &var_list) fxA[i][nv] = A*res_flx[i][nv];
      #endif
      #if GEOMETRY == SPHERICAL
      EXPAND(fxA[i][iBR]   = A*res_flx[i][iBR];         ,
             fxA[i][iBTH]  = x1p[i]*res_flx[i][iBTH];   ,
             fxA[i][iBPHI] = x1p[i]*res_flx[i][iBPHI];)
      #if HAVE_ENERGY
      fxA[i][ENG] = A*res_flx[i][ENG];
      #endif
      #endif
    }

  /* -- 1b. Build rhs in the X1-direction -- */

    for (i = beg; i <= end; i++){
      dtdV = dt/grid->dV[k][j][i];
      dtdl = dt/dx1[i];
      #if GEOMETRY == SPHERICAL
      q = dtdl/x1[i];
      EXPAND(rhs[iBR]   = 0.0;                                  ,
             rhs[iBTH]  = q*(fxA[i][iBTH]  - fxA[i-1][iBTH]);   ,
             rhs[iBPHI] = q*(fxA[i][iBPHI] - fxA[i-1][iBPHI]);)
      #if HAVE_ENERGY 
      rhs[ENG] = dtdV*(fxA[i][ENG] - fxA[i-1][ENG]);
      #endif
      #else
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[i][nv] - fxA[i-1][nv]);
      #if (GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL) && (defined iBPHI)
      rhs[iBPHI] = dtdl*(res_flx[i][iBPHI] - res_flx[i-1][iBPHI]);
      #endif
      #endif

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == JDIR){

  /* --------------------------------------------------------
     2a. Compute fluxes & sources in the X2 direction 
     -------------------------------------------------------- */

    for (j = beg-1; j <= end; j++){
      A = grid->A[JDIR][k][j][i];
      FOR_EACH(nv, &var_list) fxA[j][nv] = A*res_flx[j][nv];
    }

  /* -- 2b. Build rhs in the X2-direction -- */

    double **dx_dl = grid->dx_dl[JDIR];
    for (j = beg; j <= end; j++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[j][nv] - fxA[j-1][nv]);
      #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
      dtdl = dt/dx2[j]*dx_dl[j][i];
      rhs[iBPHI] = dtdl*(res_flx[j][iBPHI] - res_flx[j-1][iBPHI]);
      #endif

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == KDIR){

  /* --------------------------------------------------------
     3a. Compute fluxes & sources in the X3 direction 
     -------------------------------------------------------- */

    for (k = beg-1; k <= end; k++){
      A = grid->A[KDIR][k][j][i];
      FOR_EACH(nv, &var_list)  fxA[k][nv] = A*res_flx[k][nv];
    }

  /* -- 3b. Build rhs in the X3-direction -- */

    for (k = beg; k <= end; k++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[k][nv] - fxA[k-1][nv]);
      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }
  }

  
  #ifdef CHOMBO
  StoreAMRFlux (res_flx, aflux,-1, BX1, BX1+COMPONENTS-1, beg-1, end,grid);
  #if HAVE_ENERGY
  StoreAMRFlux (res_flx, aflux,-1, ENG, ENG, beg-1, end, grid);
  #endif
  #endif
}
