/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RMHD using total energy density.

  Try to recover gas pressure from conserved variables {D, m, E, B}
  using the algorithm outlined in Section A3 of Mignone \& McKinney (2007). 
  Specifically, we solve Eq. (A4) or (A6) (depending on the value of
  \c ::RMHD_REDUCED_ENERGY) using a Newton-Raphson scheme.
  Here W = rho*h*lorentz^2, E, D, etc... have the same meaning
  as in the mentioned paper.

  \author A. Mignone (mignone@ph.unito.it)

  \date   June 25, 2015

  \b References
     - "Equation of state in relativistic magnetohydrodynamics: variable versus
        constant adiabatic index"\n
        Mignone \& Mc Kinney, MNRAS (2007) 378, 1118. 

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER 50
/* ********************************************************************* */
int EnergySolve (Map_param *par)
/*!
 *
 * Solve f(W) = 0, where f(W) is Equation (A4) or (A6).
 * \return Return (0) is successful, (1) otherwise.
 *
 *********************************************************************** */
{
  int  k;
  int  done;
  double Y1, Y2, scrh, th;
  double W, W2, S2_W2, fW, dfW, dW;
  double dv2_dW, chi, chi_p_rho;
  double rho, p, lor, lor2;
  double dp, dp_drho, dp_dchi;
  double dchi_dW, drho_dW;
  double acc = 1.e-11;
  double one_m_v2, vel2, Wp;
  double m2, B2, S, S2;
  double D, E;

/* -- copy structure members to local variables -- */

  D  = par->D;
  E  = par->E;
  S  = par->S;
  m2 = par->m2;
  B2 = par->B2;
  S2 = par->S2;

/* -------------------------------------------
    Provide initial guess by taking the
    positive root of Eq. (A27).
   ------------------------------------------- */

  #if RMHD_REDUCED_ENERGY == YES
   Y1 = -4.0*(E + D - B2);
   Y2 = m2 - 2.0*(E + D)*B2 + B2*B2;
  #else
   Y1 = -4.0*(E - B2);
   Y2 = m2 - 2.0*E*B2 + B2*B2;
  #endif
  chi = Y1*Y1 - 12.0*Y2;
  chi = MAX(0.0,chi);
  W = ( - Y1 + sqrt(chi))/6.0;
  W = MAX(D, W);

  #if RMHD_REDUCED_ENERGY == YES
   Wp = W - D;
  #endif
 
  done = 0; p = -1.0;
  for (k = 1; k < MAX_ITER; k++) {

    #if RMHD_REDUCED_ENERGY == YES
     W = Wp + D;
    #endif
  
    W2    = W*W;
    S2_W2 = S2/W2;
    Y1    = 1.0/(W + B2);
    Y2    = Y1*Y1;

    vel2 = S2_W2*Y1*(Y1*W + 1.0) + m2*Y2;       /* Eq (A3) */
    one_m_v2 = 1.0 - vel2;
    lor2 = 1.0/one_m_v2;
    lor  = sqrt(lor2);
    if (vel2 > 1.0){
      WARNING(
        print ("! EnergySolve: |v| = %f > 1 during iter # %d, ",
                 vel2,k);
      )
      return (1);
    }
    #if RMHD_REDUCED_ENERGY == YES
     chi = Wp/lor2 - D*vel2/(lor + 1.0);
    #else
     chi = (W - D*lor)*one_m_v2;
    #endif

    dv2_dW  = -2.0*Y2*(3.0*S2_W2 + Y1*(S2_W2*B2*B2/W + m2)); /* Eq (A16) */

   /* -- if chi < 0 we let it converge anyway -- */

    rho = D/lor;

   /* -- kinematical terms -- */

    dchi_dW =  one_m_v2 - 0.5*lor*(D + 2.0*chi*lor)*dv2_dW;
    drho_dW = -0.5*D*lor*dv2_dW;

    #if EOS == IDEAL

     dp_dchi = (g_gamma - 1.0)/g_gamma;    /* Eq. (A 18) */
     dp_drho = 0.0;
     p       = chi*dp_dchi;
 
    #elif EOS == TAUB

     chi_p_rho = chi + rho;
     scrh = sqrt(9.0*chi*chi + 18.0*rho*chi + 25.0*rho*rho);
     p    = 2.0*chi*(chi_p_rho + rho)/(5.0*chi_p_rho + scrh);  /* Eq (A22) */

     scrh    = 1.0/(5.0*chi_p_rho - 8.0*p);
     dp_dchi = (2.0*chi_p_rho - 5.0*p)*scrh;    /* Eq (A20) */
     dp_drho = (2.0*chi - 5.0*p)*scrh;          /* Eq (A21) */

    #endif
    
    if (done) break;

    dp  = dp_dchi*dchi_dW + dp_drho*drho_dW;
    #if RMHD_REDUCED_ENERGY == YES
     fW  = Wp + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p);   /* Eq. (A25) */
     dfW = 1.0 - dp - (B2*m2 - S2)*Y2*Y1;               /* Eq. (A8)  */
     dW  = fW/dfW;
     Wp -= dW;
     if (fabs(dW) < acc*Wp || fabs(fW) < acc) done = 1;
    #else
     fW  = W + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p);  /* Eq. (A25) */
     dfW = 1.0 - dp - (B2*m2 - S2)*Y2*Y1;             /* Eq. (A8) */
     dW  = fW/dfW;
     W -= dW;
     if (fabs(dW) < acc*W || fabs(fW) < acc) done = 1;
    #endif

  }

  if (k == MAX_ITER) {
    WARNING(
      print ("! EnergySolve(): too many iterations, %12.6e %12.6e %12.6e, ", 
              W, dW, fW); 
    )
    return(1);
  }
  if (p < 0.0) {
    WARNING(
      print ("! EnergySolve(): negative pressure, p = %12.6e, ",p);
    )
    return(1);
  }

/* -- set output parameters -- */

  #if RMHD_REDUCED_ENERGY == YES
   par->W = Wp + D;
  #else
   par->W = W;
  #endif

  par->rho = rho;
  par->lor = lor;
  par->prs = p;

/* -- Recompute entropy consistently -- */

#if ENTROPY_SWITCH
  #if EOS == IDEAL
  par->sigma_c = par->prs*lor/pow(rho,g_gamma-1);
  #elif EOS == TAUB
  th = p/rho;  
  par->sigma_c = par->prs*lor/pow(rho,2.0/3.0)*(1.5*th + sqrt(2.25*th*th + 1.0));
  #endif
#endif

  return(0);  /* -- normal exit -- */
}
#undef MAX_ITER
