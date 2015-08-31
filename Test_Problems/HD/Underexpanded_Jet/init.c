/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Underexpanded jet.

  The domain is initialized with a static medium, whereas in the bottom
  boundary the following conditions are applied (see Section 5.2 of [Mig07])
  \f[
    P_{\rm jet}   = \frac{P_{\rm ratio}}{\Gamma}\left(\frac{2}{\Gamma + 1}\right)^{\Gamma/(\Gamma - 1)}
    \,,\quad
    \rho_{\rm jet} = \rho_{\rm ratio}\left(\frac{2}{\Gamma + 1}\right)^{1/(\Gamma - 1)}
    \,,\quad
    v_{\rm jet}   = \sqrt{\frac{\Gamma P_{\rm jet}}{\rho_{\rm jet}}}\,.
  \f]

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[DN_RATIO]</tt>: sets the value of \f$\rho_{\rm ratio}\f$;
  - <tt>g_inputParam[PR_RATIO]</tt>: sets the value of \f$P_{\rm ratio}\f$;

  Configurations:

  - #01: Second order accuracy;
  - #02: Third order accuracy;

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 09, 2014 

  \b References: 
     - [Mig07]: "PLUTO: a numerical code for computational astrophysics",
       Mignone et al., ApJS (2007), 170, 228

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  g_gamma = 5./3.;

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = 1.0/g_gamma;
  us[TRC] = 0.0;

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  *R;
  real    pjet, dnjet, vjet;
  real    scrh;

  scrh = 1.0/(g_gamma - 1.0);

  R     = grid[IDIR].xgc;
  pjet  = g_inputParam[PR_RATIO]*pow(2.0/(g_gamma + 1.0),g_gamma*scrh)/g_gamma;
  dnjet = g_inputParam[DN_RATIO]*pow(2.0/(g_gamma + 1.0),scrh);
  vjet  = sqrt(g_gamma*pjet/dnjet);

  if (side == X2_BEG){

    X2_BEG_LOOP(k,j,i){

      if (R[i] <= 1.) {
        d->Vc[RHO][k][j][i] = dnjet;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = vjet;
        d->Vc[PRS][k][j][i] = pjet;
      } else {
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][2*JBEG - j - 1][i];
      }
    }
  } 
}

