#include"pluto.h"

/* ********************************************************************** */
void AUSMp_Solver (const Sweep *sweep, int beg, int end, 
                   double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the Euler equations using the AUSM+ 
 *     scheme given in 
 *
 *  "A Sequel to AUSM: AUSM+"
 *  Liou, M.S., JCP (1996), 129, 364
 *
 * LAST_MODIFIED
 *
 *   Oct 11, 2016 by Andrea Mignone  (mignone@to.astro.it)
 *
 ************************************************************************ */
{
#if EOS == IDEAL
  int  i, nv;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  double aL, ML, MpL, PpL, asL2, asL, atL;
  double aR, MR, MmR, PmR, asR2, asR, atR;
  double a, m, p, mp, mm;
  double *vL, *vR, *uL, *uR;
  double alpha = 3.0/16.0, beta = 0.125;

  for (i = beg; i <= end; i++)  {

    vL = stateL->v[i];
    vR = stateR->v[i];

    uL = stateL->u[i];
    uR = stateR->u[i];

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

    sweep->press[i] = PpL*vL[PRS] + PmR*vR[PRS];
    sweep->flux[i][RHO]        = a*(mp*uL[RHO] + mm*uR[RHO]);
    EXPAND(sweep->flux[i][MX1] = a*(mp*uL[MX1] + mm*uR[MX1]); ,
           sweep->flux[i][MX2] = a*(mp*uL[MX2] + mm*uR[MX2]); ,
           sweep->flux[i][MX3] = a*(mp*uL[MX3] + mm*uR[MX3]); )
    sweep->flux[i][ENG]        = a*(mp*(uL[ENG] + vL[PRS]) + mm*(uR[ENG] + vR[PRS]));
 
  /*  ----  get max eigenvalue  ----  */

    cmax[i] = MAX(fabs(vL[VXn]) + aL, fabs(vR[VXn]) + aR);

    g_maxMach = MAX(fabs(ML), g_maxMach);
    g_maxMach = MAX(fabs(MR), g_maxMach);

  }
#else
  print ("! AUSMp_Solver: not defined for this EOS\n");
  QUIT_PLUTO(1);
#endif /* EOS == IDEAL */
}
