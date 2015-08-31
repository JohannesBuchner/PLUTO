#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, double *kappar,
              double *kapnor, double *phi)
{
  double a0, T0, a, ap, scrh;

  ap = 4.412;
  T0 = 6.05736872274111e7;
  a0 = 1.e16/pow(T0,2.5);
  a = ap/a0;

  *phi = 1.0e14;

  scrh = v[RHO] * a * pow(v[PRS]/v[RHO], 2.5) / (g_gamma - 1.0);

  *kappar = scrh;
  *kapnor = 0.0;

}
