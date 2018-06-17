/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the outermost wave speeds for HLL-based solvers.

  HLL_Speed() computes an estimate to the leftmost and rightmost
  wave signal speeds bounding the Riemann fan based on the input states
  ::stateL->v, ::stateR->v.
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Oct 12, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLL_Speed (const State *stateL, const State *stateR,
                double *SL, double *SR, int beg, int end)
/*!
 * Compute leftmost (SL) and rightmost (SR) speeds for the Riemann fan.
 * 
 * \param [in]  stateL   pointer to a state structure for the left state
 * \param [in]  stateR   pointer to a state structure for the right state
 * \param [out] SL       the (estimated) leftmost speed of the Riemann fan
 * \param [out] SR       the (estimated) rightmost speed of the Riemann fan
 * \param [in]  beg      starting index of computation
 * \param [in]  end      final index of computation
 *
 *********************************************************************** */
{
  int    i;
  static double *sl_min, *sl_max;
  static double *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------
    use Davis estimate for the signal velocities
   ---------------------------------------------- */

  MaxSignalSpeed (stateL, sl_min, sl_max, beg, end);
  MaxSignalSpeed (stateR, sr_min, sr_max, beg, end);
  for (i = beg; i <= end; i++) {
    SL[i] = MIN(sl_min[i], sr_min[i]);
    SR[i] = MAX(sl_max[i], sr_max[i]);
  }
}
