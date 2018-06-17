#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, double *kappar,
              double *kapnor, double *phi)
{
  double a0, T0, a, ap, scrh;

  ap = 4.412;
T0 = 0.5*CONST_mp/CONST_kB*UNIT_LENGTH*UNIT_LENGTH;

  a0 = UNIT_LENGTH*UNIT_LENGTH/pow(T0,2.5);
  a = ap/a0;

  *phi = 1.0e14;

  scrh = v[RHO] * a * pow(v[PRS]/v[RHO], 2.5) / (g_gamma - 1.0);

//print ("k(titos) = %12.6e\n",2.0*1.e10*CONST_kB*ap/(g_gamma-1.0)); exit(1);

double n=2.5, at;
double kpar = g_inputParam[KPAR];
double nref = g_inputParam[NREF];
T0   = TREF;
a    = (g_gamma - 1.0)*kpar/(2.0*nref*CONST_kB);
at   = a/(UNIT_LENGTH*UNIT_LENGTH)*pow(T0,n);
scrh = at*pow(v[PRS]/v[RHO],n)/(g_gamma-1.0);

  *kappar = scrh;
  *kapnor = 0.0;

}
