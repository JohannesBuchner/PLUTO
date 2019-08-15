#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, double *kappar,
              double *kapnor, double *phi)
{
  double n=2.5, a, at;
  double kpar = g_inputParam[KPAR];
  double nref = g_inputParam[NREF];
  double T0   = TREF;
  
  a       = (g_gamma - 1.0)*kpar/(2.0*nref*CONST_kB);  /* in c.g.s units */
  at      = a/(UNIT_LENGTH*UNIT_LENGTH)*pow(T0,n);     /* in code units  */
  *kappar = at*pow(v[PRS]/v[RHO],n)/(g_gamma-1.0);
  *kapnor = 0.0;

}
