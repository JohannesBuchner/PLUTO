/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RHD using total energy density.

  Try to recover gas pressure from conserved variables {D, m, E}
  using the algorithm outlined in Section 2 of Mignone, Plewa Bodo (2005).
  Specifically, we solve Eq. (12) using a Newton-Raphson root finder.

  \b References
     - "The Piecewise Parabolic Method for Multidimensional Relativistic
        Fluid Dynamics"\n
        Mignone, Plewa \& Bodo, ApJS (2005) 160, 199.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER 20 
/* ********************************************************************* */
int EnergySolve (Map_param *par)
/*!
 *  Use Newton algorithm to recover pressure p0 from 
 *  conservative quantities  D, m and E 
 *
 *********************************************************************** */
{
  int iter;
  double p, D_1, alpha, alpha2, lor2, lor, m, tol=1.e-14;
  double tau, theta, h, dh_dp, dh_dtau, gmmr;
  double yp, dyp, dp, scrh;
  double D, E, m2;

  D    = par->D;
  E    = par->E;
  m2   = par->m2;

  #if EOS == IDEAL
   gmmr = g_gamma/(g_gamma - 1.0);
  #endif
  m    = sqrt(m2);
  p    = m - E;
  p    = MAX(p, 1.e-18);
  D_1  = 1.0/D;
  for (iter = 0; iter < MAX_ITER; iter++) {

    alpha  = E + p;
    alpha2 = alpha*alpha;
    lor2   = 1.0 - m2/alpha2;

    lor2 = MAX(lor2,1.e-9);
    lor2 = 1.0/lor2;
    lor  = sqrt(lor2);

    tau   = lor*D_1;
    theta = p*tau;

    #if EOS == IDEAL
     h       = 1.0 + gmmr*theta;
     dh_dp   = gmmr*tau;
     dh_dtau = gmmr*p;
    #elif EOS == TAUB
     h       = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
     scrh    = (5.0*h - 8.0*theta)/(2.0*h - 5.0*theta);
     dh_dp   = tau*scrh;
     dh_dtau = p*scrh;
    #endif

    yp  = D*h*lor - E - p;
    dyp = D*lor*dh_dp - m2*lor2*lor/(alpha2*alpha)*(lor*dh_dtau + D*h) - 1.0;
    dp  = yp/dyp;
    p  -= dp;
    if (fabs (dp) < tol*p) break;
  }

/* ----------------------------------------------------------
            check if solution is consistent
   ---------------------------------------------------------- */

  if (p != p) {
    WARNING(
      print("! EnergySolve: NaN found while recovering pressure, ");
    )
    return 1;
  } 

  if (p < 0.0){
    WARNING(  
      print("! EnergySolve: negative pressure (%8.2e), ", p);
    )
    return 1;
  }

  par->rho = par->D/lor;
  par->W   = D*h*lor;
  par->prs = p;
  par->lor = lor;

/* -- Recompute entropy consistently -- */

#if ENTROPY_SWITCH
  #if EOS == IDEAL
  par->sigma_c = p*lor/pow(par->rho,g_gamma-1);
  #elif EOS == TAUB
  theta = p/rho;  
  par->sigma_c = p*lor/pow(par->rho,2.0/3.0)*(1.5*theta + sqrt(2.25*theta*theta + 1.0);
  #endif
#endif

  return 0; /* -- success -- */
}
#undef MAX_ITER
