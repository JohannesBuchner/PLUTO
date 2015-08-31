/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the outermost wave speeds for HLL-based solvers.

  HLL_Speed() computes an estimate to the leftmost and rightmost
  wave signal speeds bounding the Riemann fan based on the input states
  ::vR and ::vL.
  Depending on the estimate, several variants are possible.
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 6, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLL_Speed (double **vL, double **vR, double *a2L, double *a2R,
                double *hL, double *hR, double *SL, double *SR, int beg, int end)
/*!
 * Compute leftmost (SL) and rightmost (SR) speed for the Riemann fan.
 * 
 * \param [in]  vL   left  state for the Riemann solver
 * \param [in]  vR   right state for the Riemann solver
 * \param [in] a2L   1-D array containing the square of the sound speed
 *                   for the left state
 * \param [in] a2R   1-D array containing the square of the sound speed
 *                   for the right state
 * \param [out] SL   the (estimated) leftmost speed of the Riemann fan
 * \param [out] SR   the (estimated) rightmost speed of the Riemann fan
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int    i, err;
  static real *sl_min, *sl_max;
  static real *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------
              DAVIS Estimate  
   ---------------------------------------------- */

  MaxSignalSpeed (vL, a2L, hL, sl_min, sl_max, beg, end);
  MaxSignalSpeed (vR, a2R, hR, sr_min, sr_max, beg, end);
/*
  err = MAX_CH_SPEED (vL, sl_min, sl_max, beg, end);
  if (err != 0) return err;
  err = MAX_CH_SPEED (vR, sr_min, sr_max, beg, end);
  if (err != 0) return err;
*/
  for (i = beg; i <= end; i++) {
    SL[i] = MIN(sl_min[i], sr_min[i]);
    SR[i] = MAX(sl_max[i], sr_max[i]);
  }
/*  return 0; */
}

