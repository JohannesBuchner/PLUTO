/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute linear reconstruction coefficients.

  Reconstruction coefficients needed for linear interpolation are
  computed in the ::PLM_CoefficientsSet() each time a new grid 
  is created.
  There's a different set of 1D coefficients in each direction.
  
  The function ::PLM_CoefficientsGet() can be used to obtain
  a set of coefficients along a desired direction.
  
  \b References
     - "High-order conservative reconstruction schemes for finite
        volume methods in cylindrical and spherical coordinates"
        A. Mignone, JCP (2014), 270, 784.
        
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Feb 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double **cp3D, **cm3D;
static double **wp3D, **wm3D;
static double **dp3D, **dm3D;

/* ********************************************************************* */
void PLM_CoefficientsSet(Grid *grid)
/*!
 *  Compute interpolation coefficients for linear reconstruction.
 *
 * \param[in] grid    pointer to Grid structure
 *
 * \return  This function has no return value
 *********************************************************************** */
{
  int    i, d, beg, end;
  double *dx, *xr, *xgc;

  if (cp3D == NULL) {
    cp3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    cm3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    dp3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    dm3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    wp3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    wm3D = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
  }

/* -----------------------------------------------------
    Compute interpolation coefficients.
    This must be the first call to this function in
    order to store coefficients into memory.
   ----------------------------------------------------- */

  for (d = 0; d < DIMENSIONS; d++){
    dx  = grid->dx[d];
    xgc = grid->xgc[d];
    xr  = grid->xr[d];

  /* -- first and last zone are excluded -- */
    
    beg = 1;
    end = grid->np_tot[d] - 2;
    for (i = beg; i <= end; i++){
      wp3D[d][i] = dx[i]/(xgc[i+1] - xgc[i]); /* coeff. for dQF in Eq. [29] */
      wm3D[d][i] = dx[i]/(xgc[i] - xgc[i-1]); /* coeff. for dQB in Eq. [29] */

      cp3D[d][i] = (xgc[i+1] - xgc[i])/(xr[i] - xgc[i]);   /* Eq. [33], cF */
      cm3D[d][i] = (xgc[i] - xgc[i-1])/(xgc[i] - xr[i-1]); /* Eq. [33], cB */

      dp3D[d][i] = (xr[i] - xgc[i])/dx[i];     /* Eq. [30], plus sign */
      dm3D[d][i] = (xgc[i] - xr[i-1])/dx[i];   /* Eq. [30], minus sign */
    }
  }
}
/* ********************************************************************* */
void PLM_CoefficientsGet(PLM_Coeffs *plm_coeffs, int dir)
/*!
 *  Retrieve reconstruction coefficients in the PLM_Coeffs structure.
 *  This function can be called only if the previous one has been
 *  completed already.
 *
 * \param[out] plm_coeffs  a pointer to a PLM_Coeffs structure containing 
 *                         the 1D coefficients needed for reconstruction
 * \param[in] dir          the desired direction
 *
 *********************************************************************** */
{
  if (cp3D == NULL) {
    print ("! PLM_Coefficients: coefficients not set.\n");
    QUIT_PLUTO(1);
  }

  plm_coeffs->cp = cp3D[dir];
  plm_coeffs->cm = cm3D[dir];

  plm_coeffs->wp = wp3D[dir];
  plm_coeffs->wm = wm3D[dir];

  plm_coeffs->dp = dp3D[dir];
  plm_coeffs->dm = dm3D[dir];

}
