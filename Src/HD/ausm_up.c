#include"pluto.h"

/* ********************************************************************** */
void AUSMp_Solver (const State_1D *state, real *cmax, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   - Solve riemann problem for the Euler equations using the AUSM+ 
 *     scheme given in 
 *
 *  "A Sequel to AUSM: AUSM+"
 *  Liou, M.S., JCP (1996), 129, 364
 *
 * LAST_MODIFIED
 *
 *   July 17th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 ************************************************************************ */
{
  int  i, nv, beg, end;
  real aL, asL2, asL, atL,  phiL[NFLX];
  real aR, asR2, asR, atR,  phiR[NFLX];
  real ML, Mp1L, Mp2L, Mm2L, Mp4L, Pp5L;
  real MR, Mm1R, Mp2R, Mm2R, Mm4R, Pm5R;
  real a, m;
  real Ku, Kp, sigma, fa, Mm2, Mo2, Mo, M, scrh;
  real *vL, *vR, *uL, *uR;
  real alpha = 3.0/16.0, beta = 0.125;
  static real  **fl, **fr, **ul, **ur;


  beg = grid[g_dir].lbeg - 1;
  end = grid[g_dir].lend;

  if (fl == NULL){
    fl = array_2D(NMAX_POINT, NFLX);
    fr = array_2D(NMAX_POINT, NFLX);
    ul = array_2D(NMAX_POINT, NVAR);
    ur = array_2D(NMAX_POINT, NVAR);
  }

  PrimToCons (state->vL, ul, beg, end);
  PrimToCons (state->vR, ur, beg, end);

  Ku = 0.75; Kp = 0.25; sigma = 1.0;

  for (i = beg; i <= end; i++)  {

    vL = state->vL[i];
    vR = state->vR[i];

    uL = ul[i];
    uR = ur[i];

    phiL[RHO] = 1.0;
    EXPAND(phiL[MX1] = vL[VX1];  ,
           phiL[MX2] = vL[VX2];  ,
           phiL[MX3] = vL[VX3];)
    phiL[ENG] = (uL[ENG] - vL[PRS])/vL[RHO];

    phiR[RHO] = 1.0;
    EXPAND(phiR[MX1] = vR[VX1];  ,
           phiR[MX2] = vR[VX2];  ,
           phiR[MX3] = vR[VX3];)
    phiR[ENG] = (uR[ENG] - vR[PRS])/vR[RHO];


    aL = sqrt(g_gamma*vL[PRS]/vL[RHO]);
    aR = sqrt(g_gamma*vR[PRS]/vR[RHO]);

    asL2  = EXPAND(vL[VX1]*vL[VX1], + vL[VX2]*vL[VX2], + vL[VX3]*vL[VX3]);
    asL2  = aL*aL/(g_gamma - 1.0) + 0.5*asL2;
    asL2 *= 2.0*(g_gamma - 1.0)/(g_gamma + 1.0);

    asR2  = EXPAND(vR[VX1]*vR[VX1], + vR[VX2]*vR[VX2], + vR[VX3]*vR[VZ);
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
          Compute alpha and beta coeff
     -------------------------------------------- */

    Mm2 = 0.5*(vL[VXn]*vL[VXn] + vR[VXn]*vR[VXn])/a/a;
    Mo2 = MAX(Mm2, 0.4);
    Mo2 = MIN(1.0, Mo2);
    Mo  = sqrt(Mo2);
    fa  = Mo*(2.0 - Mo);

    alpha = 3.0/16.0*(-4.0 + 5.0*fa*fa);
    beta  = 1.0/8.0;

  /* --------------------------------------------
          define split Mach numbers  
          define pressure terms
     -------------------------------------------- */

    ML = vL[VXn]/a;
    MR = vR[VXn]/a;

    Mp1L = 0.5*(ML + fabs(ML));
    Mm1R = 0.5*(MR - fabs(MR));

    Mp2L =   0.25*(ML + 1.0)*(ML + 1.0);
    Mm2L = - 0.25*(ML - 1.0)*(ML - 1.0);

    Mp2R =   0.25*(MR + 1.0)*(MR + 1.0);
    Mm2R = - 0.25*(MR - 1.0)*(MR - 1.0);

    if (fabs(ML) >= 1.0){
      Mp4L = Mp1L;
      Pp5L = Mp1L/ML;
    }else{
      Mp4L = Mp2L*(1.0 - 16.0*beta*Mm2L);
      Pp5L = Mp2L*((2.0 - ML) - 16.0*alpha*ML*Mm2L);
    }

    if (fabs(MR) >= 1.0){
      Mm4R = Mm1R;
      Pm5R = Mm1R/MR;
    }else{
      Mm4R = Mm2R*(1.0 + 16.0*beta*Mp2R);
      Pm5R = Mm2R*(( - 2.0 - MR) + 16.0*alpha*MR*Mp2R);
    }

    scrh  = MAX(1.0 - sigma*Mm2, 0.0);

    M  = Mp4L + Mm4R - Kp/fa*scrh*(vR[PRS] - vL[PRS])/(vL[RHO] + vR[RHO])/a/a*2.0;

    m  = a*M;
    m *= M > 0.0 ? vL[RHO]: vR[RHO];

  /* -------------------------------------------------------------
                     Compute fluxes 
     ------------------------------------------------------------- */

    state->press[i]  = Pp5L*vL[PRS] + Pm5R*vR[PRS];
    state->press[i] -= Ku*Pp5L*Pm5R*(vL[RHO] + vR[RHO])*fa*a*(vR[VXn] - vL[VXn]);

    if (m > 0.0){
      for (nv = NFLX; nv--;   ) state->flux[i][nv] = m*phiL[nv];
    }else{
      for (nv = NFLX; nv--;   ) state->flux[i][nv] = m*phiR[nv];
    }
 
  /*  ----  get max eigenvalue  ----  */

    *cmax = MAX(*cmax, fabs(vL[VXn]) + aL);
    *cmax = MAX(*cmax, fabs(vR[VXn]) + aR);

    g_maxMach = MAX(fabs(ML), g_maxMach);
    g_maxMach = MAX(fabs(MR), g_maxMach);

  }
}
