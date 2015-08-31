/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Thermal conduction (TC) module header file.

  Contains prototypes for the thermal conduction module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos  
  \date   Sep 13, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

void   TC_Flux (double ***, const State_1D *, double **, int, int, Grid *);
void   TC_kappa (double *, double, double, double, double *, double *, double *);
void   GetGradient (double ***, double **, int, int, Grid *);
