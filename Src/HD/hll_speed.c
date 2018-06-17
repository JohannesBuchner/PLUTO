/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the outermost wave speeds for HLL-based solvers.

  HLL_Speed() computes an estimate to the leftmost and rightmost
  wave signal speeds bounding the Riemann fan based on the input states
  ::stateL->v, ::stateR->v.
  Depending on the estimate, several variants are possible.
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Feb 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/*   ****  switches, see below for description  ****  */

#define DAVIS_ESTIMATE     YES
#define EINFELDT_ESTIMATE  NO
#define ROE_ESTIMATE       NO

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
 * Switches:
 *
 *    ROE_ESTIMATE (YES/NO), DAVIS_ESTIMATE (YES/NO). TVD_ESTIMATE (YES/NO)
 *    JAN_HLL (YES/NO) 
 *
 *    These switches set how the wave speed estimates are
 *    going to be computed. Only one can be set to 'YES', and
 *    the rest of them must be set to 'NO'  
 *
 *     ROE_ESTIMATE:    b_m = \min(0, \min(u_R - c_R, u_L - c_L, u_{roe} - c_{roe}))     
 *                      b_m = \min(0, \min(u_R + c_R, u_L + c_L, u_{roe} + c_{roe}))
 * 
 *                      where u_{roe} and c_{roe} are computed using Roe averages.
 *  
 *     DAVIS_ESTIMATE:  b_m = \min(0, \min(u_R - c_R, u_L - c_L))     
 *                      b_m = \min(0, \min(u_R + c_R, u_L + c_L))  
 *
 *********************************************************************** */
{
  int    i;
  double scrh, s, c;
  double aL, dL;
  double aR, dR;
  double dvx, dvy, dvz;
  double a_av, du, vx;
  static double *sl_min, *sl_max;
  static double *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
              DAVIS Estimate  
   -------------------------------------------------------- */

#if DAVIS_ESTIMATE == YES
   MaxSignalSpeed (stateL, sl_min, sl_max, beg, end);
   MaxSignalSpeed (stateR, sr_min, sr_max, beg, end);
   for (i = beg; i <= end; i++) {

     SL[i] = MIN(sl_min[i], sr_min[i]);
     SR[i] = MAX(sl_max[i], sr_max[i]);

     aL = sqrt(stateL->a2[i]);
     aR = sqrt(stateR->a2[i]);

     scrh  = fabs(stateL->v[i][VXn]) + fabs(stateR->v[i][VXn]);    
     scrh /= aL + aR;

     g_maxMach = MAX(scrh, g_maxMach); 
   }
#endif /* DAVIS_ESTIMATE == YES */

/* --------------------------------------------------------
        Einfeldt Estimate for wave speed 
   -------------------------------------------------------- */

#if EINFELDT_ESTIMATE == YES
  for (i = beg; i <= end; i++) {

    aL = sqrt(a2L[i]);
    aR = sqrt(a2R[i]);

    dL   = sqrt(vL[i][RHO]);
    dR   = sqrt(vR[i][RHO]);
    a_av = 0.5*dL*dR/( (dL + dR)*(dL + dR));
    du   = vR[i][VXn] - vL[i][VXn];
    scrh = (dR*aL*aL + dR*aR*aR)/(dL + dR);
    scrh += a_av*du*du;
      
    SL[i] = (dl*ql[VXn] + dr*qr[VXn])/(dl + dr);
  
    bmin = MIN(0.0, um[VXn] - sqrt(scrh));
    bmax = MAX(0.0, um[VXn] + sqrt(scrh));      
  }
#endif /* EINFELDT_ESTIMATE == YES */

/* --------------------------------------------------------
              Roe-like Estimate  
   -------------------------------------------------------- */

#if ROE_ESTIMATE  == YES   
  for (i = beg; i <= end; i++) {

    aL = sqrt(a2L[i]);
    aR = sqrt(a2R[i]);

    s = 1.0/(1.0 + sqrt(vR[i][RHO]/vL[i][RHO]));
    c = 1.0 - s;

    vx = s*vL[i][VXn] + c*vR[i][VXn];

  /* ------------------------------------------------------
     The following definition of the averaged sound speed
     is equivalent to original formulation given by Roe;
     here we just simplify the expression without 
     explicitly computing the enthalpy:                                 
     ------------------------------------------------------ */
   
    EXPAND(dvx = vL[i][VX1] - vR[i][VX1];  ,
           dvy = vL[i][VX2] - vR[i][VX2];  ,
           dvz = vL[i][VX3] - vR[i][VX3];)

    scrh = EXPAND(dvx*dvx, + dvy*dvy, + dvz*dvz);
  
    a_av = sqrt(s*aL2 + c*aR2 + 0.5*s*c*gmm1*scrh);

    SL[i] = MIN(vL[i][VXn] - aL, vx - a_av);
    SR[i] = MAX(vR[i][VXn] + aR, vx + a_av);

    scrh = (fabs(vL[i][VXn]) + fabs(vR[i][VXn]))/(aL + aR);
    g_maxMach = MAX(scrh, g_maxMach);
  }
#endif /* ROE_ESTIMATE == YES */   
}
#undef DAVIS_ESTIMATE     
#undef EINFELDT_ESTIMATE  
#undef ROE_ESTIMATE       
