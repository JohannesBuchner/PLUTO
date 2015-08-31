#include "pluto.h"

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  eta[IDIR] = g_inputParam[ETAX];
  eta[JDIR] = g_inputParam[ETAY];
  eta[KDIR] = g_inputParam[ETAZ];
}
