#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, 
              double *kpar, double *knor, double *phi)
{
  double mu, nH, sqT, T, B2_cgs;
  double ctts_mu = 1.265060;

  mu   = ctts_mu/2.0;
  T    = v[PRS]/v[RHO]*(mu*CONST_mp*UNIT_VELOCITY*UNIT_VELOCITY/CONST_kB);
  sqT  = sqrt(T);

  *kpar = 5.6e-7*T*T*sqT;

  #if PHYSICS == MHD
   nH      = v[RHO]*UNIT_DENSITY/(mu*CONST_mp);
   B2_cgs  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]) + 1.e-12;
   B2_cgs *= 4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
   *knor   = 3.3e-16*nH*nH/sqT/B2_cgs;
  #else
   *knor = 0.0;
  #endif

  *kpar *= CONST_mp*mu/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);
  *knor *= CONST_mp*mu/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);

  *phi = 0.3;
}

