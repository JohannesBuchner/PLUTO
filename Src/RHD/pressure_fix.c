/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RHD using a pressure fix.

  Fix p to a small value, solve for the square of velocity.
  This step involved re-computing W at each step of the iteration.
  Once the root has been found, we recompute total energy E.
  Return 0 if succesful, 1 otherwise.

  \authors A. Mignone \n
           C. Zanni

  \date Oct 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define MAX_ITER 20
static double VelocitySquareFunc(double, Map_param *);
/* ********************************************************************* */
int PressureFix(Map_param *par)
/*!
 * Fix p to a small value, solve for the square of velocity by using
 * secant algorithm applied to Eq (9) of Mignone, Plewa \& Bodo (2005).
 * This step involved re-computing W at each step of the iteration.
 * Once the root has been found, we recompute total energy E.
 * Return 0 if succesful, 1 otherwise.
 *
 *********************************************************************** */
{
  int    k, done=0;
  double v2, v2c, fc, f, dW, S2_W2;
  double fmin, fmax, v2min, v2max;
  
  par->prs = g_smallPressure; 

  v2max = 1.0-1.e-8;
  v2c = 0.95;
  fc  = VelocitySquareFunc(v2c, par);
  v2  = 0.96;
  for (k = 1; k < MAX_ITER; k++){
    f   = VelocitySquareFunc(v2, par);
    if (done == 1) break;
    dW  = (v2 - v2c)/(f - fc)*f;
    v2c = v2; fc = f;
    v2 -= dW;
    v2 = MIN(v2max,v2);
    v2 = MAX(v2, 0.0);
    if (fabs(f) < 1.e-9) done = 1;
  }
  if (v2 >= 1.0 || k >= MAX_ITER) {
    print ("! PressureFix: too many iter while fixing p , v^2 = %f\n", v2);
    return (1);
  }

/* -- Redefine energy, density and entropy -- */
  
  par->E   = par->W - par->prs;
  par->rho = par->D/par->lor;

#if ENTROPY_SWITCH
{
  double rho = par->rho;
  double th  = par->prs/rho; 
  #if EOS == IDEAL
  par->sigma_c = par->prs*par->lor/pow(rho,g_gamma-1);
  #elif EOS == TAUB
  th = par->prs/rho;  
  par->sigma_c = par->prs*par->lor/pow(rho,2.0/3.0)*(1.5*th + sqrt(2.25*th*th + 1.0);
  #endif
}
#endif

  return(0);  /* -- success -- */
} 

/* ****************************************************************** */
double VelocitySquareFunc (double v2, Map_param *par)
/*!
 * 
 * Implement Eq (A3).
 * 
 * 
 ******************************************************************** */
{
  double lor2, pg;

  lor2     = 1.0/(1.0 - v2);
  par->lor = sqrt(lor2);
  pg       = par->prs*par->lor;
  #if EOS == IDEAL
   par->W = (par->D + pg*g_gamma/(g_gamma - 1.0))*par->lor;
  #elif EOS == TAUB
   par->W = (2.5*pg + sqrt(2.25*pg*pg + par->D*par->D))*par->lor;
  #endif

  return par->m2/(par->W*par->W) - v2;
}
#undef MAX_ITER
