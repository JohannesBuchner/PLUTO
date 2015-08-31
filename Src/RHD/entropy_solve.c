/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RHD using entropy.

  Convert the conservative variables u=[D, m, sigma_c]
  (where sigma_c = D*sigma is the conserved entropy) to
  primitive variable using a Newton-Raphson/Bisection scheme.

  \authors C. Zanni \n
           A. Mignone

  \date Oct 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

#define MAX_ITER 20
/* ********************************************************************* */
int EntropySolve (Map_param *par)
/*!
 * Convert the conservative variables {D, m, sigma_c}
 * (where sigma_c = D*sigma is the conserved entropy) to 
 * primitive variable using a Newton-Raphson/Bisection scheme.
 *
 * LAST MODIFIED
 *
 *  February 5th (2012) by C. Zanni & A. Mignone 
 *  (zanni@oato.inaf.it, mignone@oato.inaf.it)
 *
 *********************************************************************** */
{
  int    k, done;
  double v2, v2max, v2min, acc = 1.e-13;
  double h, W, W2, W3;
  double scrh, x, temp;
  double lor, lor2, rho, sigma, th;
  double dW_dlor;
  double fv2, dfv2, dv2, dv2old;
  double m2, D;

  D  = par->D;
  m2 = par->m2;

  sigma = par->sigma_c/D;
  v2max = 1.0 - 1.e-8;
  v2min = 0.0;

  done = 0;

  v2     = 0.5*(v2max+v2min);
  dv2old = v2max-v2min;
  dv2    = dv2old;

  lor2 = 1.0/(1.0 - v2);
  lor  = sqrt(lor2);
  rho  = D/lor;

  #if EOS == IDEAL
   x  = pow(rho,g_gamma - 1.0);
   h  = 1.0 + g_gamma/(g_gamma - 1.0)*sigma*x;
   W  = D*h*lor;
   dW_dlor = (2 - g_gamma)*W/lor + (g_gamma - 1.0)*D;
  #elif EOS == TAUB
   x  = pow(rho,2.0/3.0);
   th = sigma*x/sqrt(1.0 + 3.0*sigma*x);
   h  = 2.5*th + sqrt(2.25*th*th + 1.0);
   W  = D*h*lor;

   scrh  = -2.0*x/(3.0*lor)*sigma*(2.0*sigma*x - 3.0*th*th);
   scrh /=  2.0*th*(1.0 + 3.0*sigma*x);
   dW_dlor = W/lor + D*lor*(5.0*h - 8.0*th)/(2.0*h - 5.0*th)*scrh;
  #endif

  W2    = W*W;
  W3    = W2*W;

  fv2  = v2 - m2/W2;
  dfv2 = 1.0 + m2/W3*dW_dlor*lor2*lor;

  for (k = 1; k < MAX_ITER; k++) {
    if ((((v2-v2max)*dfv2-fv2)*((v2-v2min)*dfv2-fv2) > 0.0)
        || (fabs(2.*fv2) > fabs(dv2old*dfv2))) {
       dv2old = dv2;
       dv2 = 0.5*(v2max-v2min);
       v2 = v2min+dv2;
       if (v2min == v2) done = 1;
    } else { 
       dv2old = dv2;
       dv2 = fv2/dfv2;
       temp = v2;
       v2 -= dv2;
       if (temp == v2) done = 1;
    }

    if (fabs(dv2) < acc*v2 || fabs(fv2) < acc) done = 1;
    
    lor2 = 1.0/(1.0 - v2);
    lor  = sqrt(lor2);

    rho = D/lor;

    #if EOS == IDEAL
     x  = pow(rho,g_gamma - 1.0);
     h  = 1.0 + g_gamma/(g_gamma - 1.0)*sigma*x;
     W  = D*h*lor;
     dW_dlor = (2 - g_gamma)*W/lor + (g_gamma - 1.0)*D;
    #elif EOS == TAUB
     x  = pow(rho,2.0/3.0);
     th = sigma*x/sqrt(1.0 + 3.0*sigma*x);
     h  = 2.5*th + sqrt(2.25*th*th + 1.0);
     W  = D*h*lor;

     scrh  = -2.0*x/(3.0*lor)*sigma*(2.0*sigma*x - 3.0*th*th);
     scrh /=  2.0*th*(1.0 + 3.0*sigma*x);
     dW_dlor = W/lor + D*lor*(5.0*h-8.0*th)/(2.0*h-5.0*th)*scrh;
    #endif

    W2    = W*W;

    if (done) break;

    W3    = W2*W;
    fv2 = v2 - m2/W2;

    if (fv2 < 0.0) v2min = v2;
    else           v2max = v2;

    dfv2 = 1.0 + m2/W3*dW_dlor*lor2*lor;
  }

  #if EOS == IDEAL
   par->prs = sigma*x*rho;
  #elif EOS == TAUB
   par->prs = th*rho;
  #endif

  if (k == MAX_ITER) {
    WARNING(
      print ("! EntropySolve: too many iterations,%d, ",k);
    )
    par->prs = g_smallPressure;
    par->E   = W - par->prs;
    return(1);
  }

  if (par->prs < 0.0 || sigma < 0.0 || W < 0.0) {
    WARNING(
      print ("! EntropySolve: negative pressure, p = %12.6e\n",par->prs);
    )
    par->prs = g_smallPressure;
    par->E   = W - par->prs;
    return(1);
  }

  par->rho = par->D/lor;
  par->lor = lor;
  par->E   = W - par->prs;  /* redefine energy */
  return(0); /* -- success -- */
}
#undef MAX_ITER
