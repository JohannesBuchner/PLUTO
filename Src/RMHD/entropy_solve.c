/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RMHD using entropy.

  Convert the conservative variables u=[D, m, B, sigma_c] 
  (where sigma_c = D*sigma is the conserved entropy) to 
  primitive variable v=[rho,v,B,p] using a Newton-Raphson/Bisection scheme.
 
  \authors C. Zanni \n
           A. Mignone

  \date  June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER 50
/* ********************************************************************* */
int EntropySolve (Map_param *par)
/*!
 *
 *********************************************************************** */
{
  int    k, done;
  double v2, v2max, v2min, acc = 1.e-13;
  double h, W, W2, W3, S2_W2;
  double scrh, x, temp;
  double lor, lor2, rho, sigma, th;
  double dv2_dW, dW_dlor;
  double fv2, dfv2, dv2, dv2old;
  double m2, B2, S, S2, D;
  
  D  = par->D;
  m2 = par->m2;
  B2 = par->B2;
  S  = par->S;
  S2 = S*S;

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
  S2_W2 = S2/W2;
  scrh  = W + B2;

  fv2 = v2 - (m2 + S2_W2*(scrh + W))/(scrh*scrh);

  dv2_dW  = S2*(3.0*W*scrh + B2*B2) + m2*W3;
  dv2_dW *= -2.0/(W3*scrh*scrh*scrh);

  dfv2 = 1.0 - 0.5*dv2_dW*dW_dlor*lor2*lor;

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
    S2_W2 = S2/W2;

    if (done) break;

    W3    = W2*W;
    scrh  = W + B2;

    fv2 = v2 - (m2 + S2_W2*(scrh + W))/(scrh*scrh);

    if (fv2 < 0.0) v2min = v2;
              else v2max = v2;

    dv2_dW  = S2*(3.0*W*scrh + B2*B2) + m2*W3;
    dv2_dW *= -2.0/(W3*scrh*scrh*scrh);

    dfv2 = 1.-0.5*dv2_dW*dW_dlor*lor2*lor;
  }

/* -- Compure pressure */

  #if EOS == IDEAL
   par->prs = sigma*x*rho;
  #elif EOS == TAUB
   par->prs = th*rho;
  #endif

  if (k == MAX_ITER) {
    WARNING(
      print ("! EntropySolve: too many iterations,%d, ",k);
    )
    return(1);
  }

  if (par->prs < 0.0 || sigma < 0.0 || W < 0.0) {
    WARNING(
      print ("! EntropySolve: negative pressure, p = %12.6e\n",par->prs);
    )
    return(1);
  }

  par->W   = W;
  par->lor = lor;
  par->rho = rho;

/* -- Redefine energy -- */

  #if RMHD_REDUCED_ENERGY == YES
   #if EOS == IDEAL
    par->E = par->prs*(g_gamma/(g_gamma-1.0)*lor2 - 1.0) + D*lor2*v2/(1.0 + lor)
            + 0.5*(1.0 + v2)*B2 - 0.5*S2_W2;
   #elif EOS == TAUB
    scrh = 2.25*th*th + v2;
    par->E = par->prs*(2.5*lor2 - 1.0) 
            + D*lor2/(lor*sqrt(2.25*th*th + 1.0) + 1.0)*scrh
            + 0.5*(1.0 + v2)*B2 - 0.5*S2_W2;
   #endif
  #else
   par->E = W - par->prs + 0.5*(1.0 + v2)*B2 - 0.5*S2_W2;
  #endif
  return(0); /* -- success -- */
}
#undef MAX_ITER
