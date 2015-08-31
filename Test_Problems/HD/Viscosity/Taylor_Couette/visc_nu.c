#include "pluto.h"

/*---------------------------------------------------------------------------*/
/*---- Specification of explicit first and second viscosity coefficients ----*/
/*---------------------------------------------------------------------------*/

void Visc_nu(double *v, double x1, double x2, double x3, 
             double *nu1, double *nu2)
{
  double rmin = g_domBeg[IDIR];
  double rmax = g_domEnd[IDIR];
  double Re   = g_inputParam[REYN]; 
	
 *nu1 = g_inputParam[OMEGA]*rmin*(rmax - rmin)*v[RHO]/Re;
 *nu2 = 0.0;
}
