/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RMHD using a pressure fix.

  Fix p to a small value, solve for the square of velocity by using
  secant or bisection algorithm applied to Eq (A3).
  This step involved re-computing W at each step of the iteration.
  Once the root has been found, we recompute total energy E.
  Return 0 if succesful, 1 otherwise.

  \authors A. Mignone \n
           C. Zanni

  \date  June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER 50
static double VelocitySquareFunc(double, void *);
/* ********************************************************************* */
int PressureFix(Map_param *par)
/*!
 *
 *********************************************************************** */
{
  int    k, done=0;
  double v2, f, dW, S2_W2;
  double fmin, fmax, v2min, v2max;
  
  par->prs = g_smallPressure; 

/* ------------------------------------------------
    Use secant algorithm
   ------------------------------------------------ */

/*
  v2max = 1.0-1.e-8;
  v2c = 0.95;
  fc  = VelocitySquareFunc(v2c, par);
  v2  = 0.96;
  for (k = 1; k < MAX_ITER; k++){
    f   = VelocitySquareFunc(v2, par);
print ("%d,  v2 = %12.6e  f = %12.6e\n",k,v2,f);
    if (done == 1) break;
    dW  = (v2 - v2c)/(f - fc)*f;
    v2c = v2; fc = f;
    v2 -= dW;
    v2 = MIN(v2max,v2);
    v2 = MAX(v2, 0.0);
    if (fabs(f) < 1.e-9) done = 1;
  }

  if (v2 >= 1.0 || k >= MAX_ITER) {
    FILE *fp;
    fp = fopen("pressure_fix.dat","w");
    for (v2 = 0.99; v2 < 1.0; v2 += 1.e-5){
      f = VelocitySquareFunc (v2, par);
      fprintf (fp, "%12.6e  %12.6e\n", v2, f);
    }
    fclose(fp);
    print ("! PressureFix: too many iter while fixing p , v^2 = %f\n", v2);
     
    return (1);
  }
*/

/* ------------------------------------------------
    Use bisection algorithm
   ------------------------------------------------ */

  v2min = 0.0;       fmin = VelocitySquareFunc(v2min, par);
  v2max = 1.0-1.e-9; fmax = VelocitySquareFunc(v2max, par);

  for (k = 1; k < MAX_ITER; k++){

    v2 = 0.5*(v2min + v2max);
    f  = VelocitySquareFunc(v2, par); 

    if (f*fmin > 0.0){
      v2min = v2; fmin = f;
    }else{
      v2max = v2; fmax = f;
    }
    if (fabs(f) < 1.e-9) break;

    if (v2 >= 1.0 || k >= MAX_ITER) {
/*
      FILE *fp;
      fp = fopen("pressure_fix.dat","w");
      for (v2 = 0.99; v2 < 1.0; v2 += 1.e-5){
        f = VelocitySquareFunc (v2, par);
        fprintf (fp, "%12.6e  %12.6e\n", v2, f);
      }
      fclose(fp);
*/
      print ("! PressureFix(): too many iterations while fixing p , v^2 = %f\n", v2);
      return (1);
    }
  }
  
/* -- Redefine energy, proper density and entropy -- */
    
  S2_W2    = par->S2/(par->W*par->W);
#if RMHD_REDUCED_ENERGY == YES
   par->E = par->W - par->D - par->prs + 0.5*(1.0 + v2)*par->B2 - 0.5*S2_W2;
#else
   par->E = par->W - par->prs + 0.5*(1.0 + v2)*par->B2 - 0.5*S2_W2;
#endif
  par->rho = par->D/par->lor;

/* -- Recompute entropy consistently -- */

#if ENTROPY_SWITCH
{
  double rho = par->rho;
  double th  = par->prs/rho;  
  #if EOS == IDEAL
  par->sigma_c = par->prs*par->lor/pow(rho,g_gamma-1);
  #elif EOS == TAUB
  par->sigma_c = par->prs*par->lor/pow(rho,2.0/3.0)*(1.5*th + sqrt(2.25*th*th + 1.0));
  #endif
}
#endif

  return(0);  /* -- success -- */
} 

/* ****************************************************************** */
double VelocitySquareFunc (double v2, void *par)
/*!
 * 
 * Implement Eq (A3) of Mignone \& McKinney (2007).
 * 
 ******************************************************************** */
{
  double lor2, W2, f, pg;
  Map_param *p = (Map_param *) par; 

  lor2   = 1.0/(1.0 - v2);
  p->lor = sqrt(lor2);

  pg     = p->prs*p->lor;
  #if EOS == IDEAL
   p->W = (p->D + pg*g_gamma/(g_gamma - 1.0))*p->lor;
  #elif EOS == TAUB
   p->W = (2.5*pg + sqrt(2.25*pg*pg + p->D*p->D))*p->lor;
  #endif
  
  W2  = p->W*p->W;

  f  =  p->S2*(2.0*p->W + p->B2) + p->m2*W2;
  f /= (p->W + p->B2)*(p->W + p->B2)*W2;
  f -= v2;
  return f;
}
#undef MAX_ITER
