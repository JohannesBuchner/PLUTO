#include"pluto.h"

/* ********************************************************************** */
void AUSMp_Solver (const State_1D *state, int beg, int end, 
            real *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the Euler equations using the AUSM+ 
 *     scheme given in 
 *
 *  "A Sequel to AUSM: AUSM+"
 *  Liou, M.S., JCP (1996), 129, 364
 *
 * LAST_MODIFIED
 *
 *   July 17th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 ************************************************************************ */
{
#if EOS == IDEAL
  int  i, nv;
  real aL, ML, MpL, PpL, asL2, asL, atL;
  real aR, MR, MmR, PmR, asR2, asR, atR;
  real a, m, p, mp, mm;
  real *vL, *vR, *uL, *uR;
  real alpha = 3.0/16.0, beta = 0.125;
  static real  **fl, **fr, **ul, **ur;

  if (fl == NULL){
    fl = ARRAY_2D(NMAX_POINT, NFLX, double);
    fr = ARRAY_2D(NMAX_POINT, NFLX, double);
    ul = ARRAY_2D(NMAX_POINT, NVAR, double);
    ur = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  PrimToCons (state->vL, ul, beg, end);
  PrimToCons (state->vR, ur, beg, end);

  for (i = beg; i <= end; i++)  {

    vL = state->vL[i];
    vR = state->vR[i];

    uL = ul[i];
    uR = ur[i];

    aL = sqrt(g_gamma*vL[PRS]/vL[RHO]);
    aR = sqrt(g_gamma*vR[PRS]/vR[RHO]);

    asL2  = EXPAND(vL[VX1]*vL[VX1], + vL[VX2]*vL[VX2], + vL[VX3]*vL[VX3]);
    asL2  = aL*aL/(g_gamma - 1.0) + 0.5*asL2;
    asL2 *= 2.0*(g_gamma - 1.0)/(g_gamma + 1.0);

    asR2  = EXPAND(vR[VX1]*vR[VX1], + vR[VX2]*vR[VX2], + vR[VX3]*vR[VX3]);
    asR2  = aR*aR/(g_gamma - 1.0) + 0.5*asR2;
    asR2 *= 2.0*(g_gamma - 1.0)/(g_gamma + 1.0);

    asL = sqrt(asL2);
    asR = sqrt(asR2);

    atL = asL2/MAX(asL, fabs(vL[VXn]));
    atR = asR2/MAX(asR, fabs(vR[VXn]));

    a = MIN(atL, atR);
/*
    a = 0.5*(aL + aR);
*/
  /* --------------------------------------------
          define split Mach numbers  
          define pressure terms
     -------------------------------------------- */

    ML = vL[VXn]/a;
    if ( fabs(ML) >= 1.0){
      MpL = 0.5*(ML  + fabs(ML));
      PpL = ML > 0.0 ? 1.0:0.0;
    }else{
      MpL = 0.25*(ML + 1.0)*(ML + 1.0) + beta*(ML*ML - 1.0)*(ML*ML - 1.0);
      PpL = 0.25*(ML + 1.0)*(ML + 1.0)*(2.0 - ML) + alpha*ML*(ML*ML - 1.0)*(ML*ML - 1.0);
    }

    MR = vR[VXn]/a;
    if ( fabs(MR) >= 1.0){
      MmR = 0.5*(MR  - fabs(MR));
      PmR = MR > 0.0 ? 0.0:1.0;
    }else{
      MmR = - 0.25*(MR - 1.0)*(MR - 1.0) - beta*(MR*MR - 1.0)*(MR*MR - 1.0);
      PmR =   0.25*(MR - 1.0)*(MR - 1.0)*(2.0 + MR) - alpha*MR*(MR*MR - 1.0)*(MR*MR - 1.0);
    }

    m = MpL + MmR;

    mp = 0.5*(m + fabs(m));
    mm = 0.5*(m - fabs(m));

  /* -------------------------------------------------------------
                     Compute fluxes 
     ------------------------------------------------------------- */

    state->press[i] = PpL*vL[PRS] + PmR*vR[PRS];
    state->flux[i][RHO]        = a*(mp*uL[RHO] + mm*uR[RHO]);
    EXPAND(state->flux[i][MX1] = a*(mp*uL[MX1] + mm*uR[MX1]); ,
           state->flux[i][MX2] = a*(mp*uL[MX2] + mm*uR[MX2]); ,
           state->flux[i][MX3] = a*(mp*uL[MX3] + mm*uR[MX3]); )
    state->flux[i][ENG]        = a*(mp*(uL[ENG] + vL[PRS]) + mm*(uR[ENG] + vR[PRS]));
 
  /*  ----  get max eigenvalue  ----  */

    cmax[i] = MAX(fabs(vL[VXn]) + aL, fabs(vR[VXn]) + aR);

    g_maxMach = MAX(fabs(ML), g_maxMach);
    g_maxMach = MAX(fabs(MR), g_maxMach);

  }
#else
  print1 ("! AUSMp_Solver: not defined for this EOS\n");
  QUIT_PLUTO(1);
#endif /* EOS == IDEAL */
}
