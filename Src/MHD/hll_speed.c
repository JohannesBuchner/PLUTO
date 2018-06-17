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

#define ROE_ESTIMATE    NO
#define DAVIS_ESTIMATE  YES   

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
 *    ROE_ESTIMATE (YES/NO), DAVIS_ESTIMATE (YES/NO)
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
  int  i;
  double scrh, s, c;
  static double *sl_min, *sl_max;
  static double *sr_min, *sr_max;
  static double *sm_min, *sm_max;
  static double **vm;

  if (sl_min == NULL){
    vm = ARRAY_2D(NMAX_POINT, NVAR, double);
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);

    sm_min = ARRAY_1D(NMAX_POINT, double);
    sm_max = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------------------
   1. Compute speeds using DAVIS Estimate  
   ------------------------------------------------------- */

#if DAVIS_ESTIMATE == YES
  MaxSignalSpeed (stateL, sl_min, sl_max, beg, end);
  MaxSignalSpeed (stateR, sr_min, sr_max, beg, end);

  for (i = beg; i <= end; i++) {

    SL[i] = MIN(sl_min[i], sr_min[i]);
    SR[i] = MAX(sl_max[i], sr_max[i]);

  /* -- define g_maxMach -- */

     scrh  = fabs(stateL->v[i][VXn]) + fabs(stateR->v[i][VXn]);    
     scrh /= sqrt(stateL->a2[i])     + sqrt(stateR->a2[i]);
/*
    #if EOS == IDEAL 
     scrh  = fabs(vL[i][VXn]) + fabs(vR[i][VXn]);    
     scrh = scrh/(sqrt(g_gamma*vL[i][PRS]/vL[i][RHO]) + sqrt(g_gamma*vR[i][PRS]/vR[i][RHO]));
    #elif EOS == ISOTHERMAL
     scrh  = 0.5*( fabs(vL[i][VXn]) + fabs(vR[i][VXn]) )/g_isoSoundSpeed;    
    #elif EOS == BAROTROPIC
     scrh  = 0.5*(fabs(vL[i][VXn]) + fabs(vR[i][VXn]));    
     s     = 0.5*(vL[i][RHO] + vR[i][RHO]);    
     c     = BAROTROPIC_PR(s);
     scrh /= sqrt(g_gamma*c/s); 
    #else
     print ("! HLL_SPEED: not defined for this EoS\n");
     QUIT_PLUTO(1);
    #endif
*/
    g_maxMach = MAX(scrh, g_maxMach);  

  }
#endif  /* DAVIS_ESTIMATE == YES */

/* --------------------------------------------------------
   2. Compute speeds using Roe-like estimates  
   ------------------------------------------------------- */

#if ROE_ESTIMATE == YES
  MaxSignalSpeed (vL, a2L, sl_min, sl_max, bgf, beg, end);
  MaxSignalSpeed (vR, a2R, sr_min, sr_max, bgf, beg, end);
  for (i = beg; i <= end; i++) {

    scrh = sqrt(vR[i][RHO]/vL[i][RHO]);
    s    = 1.0/(1.0 + scrh);
    c    = 1.0 - s;

    vm[i][RHO] = vL[i][RHO]*scrh;
    vm[i][VXn] = s*vL[i][VXn] + c*vR[i][VXn];

  /*  ----------------------------------
       the next definition is not the 
       same as the one given by Roe; it 
       is used here to define the sound 
       speed.
      ----------------------------------  */

  #if EOS == IDEAL
  vm[i][PRS] = s*gpl + c*gpr;  
    #endif
 
  /* ---------------------------------------
      this is the definitions given by
      Cargo and Gallice, see the Roe 
      Solver
     --------------------------------------- */

    EXPAND(vm[i][BXn] = c*vL[i][BXn] + s*vR[i][BXn];   ,
           vm[i][BXt] = c*vL[i][BXt] + s*vR[i][BXt];   ,
           vm[i][BXb] = c*vL[i][BXb] + s*vR[i][BXb];)

  }

  MAX_CH_SPEED(vm, sm_min, sm_max, bgf, beg, end);

  for (i = beg; i <= end; i++) {
    SL[i] = MIN(sl_min[i], sm_min[i]);
    SR[i] = MAX(sr_max[i], sm_max[i]);

  /* -- define g_maxMach -- */

  #if EOS == IDEAL
    scrh  = fabs(vm[i][VXn])/sqrt(vm[i][PRS]/vm[i][RHO]);
  #elif EOS == ISOTHERMAL
    scrh  = fabs(vm[i][VXn])/g_isoSoundSpeed;
  #elif EOS == BAROTROPIC
    print ("! HLL_SPEED: stop\n");
    QUIT_PLUTO(1);
  #endif   
    g_maxMach = MAX(scrh, g_maxMach); 
  }
#endif /* ROE_ESTIMATE == YES */

}
#undef DAVIS_ESTIMATE
#undef ROE_ESTIMATE
