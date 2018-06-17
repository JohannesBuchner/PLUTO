/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Thermal conduction (TC) module header file.

  Contains prototypes for the thermal conduction module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos  
  \date    May 11, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */

void   TC_Flux (double ***, const Sweep *, double *, int, int, Grid *);
void   TC_kappa (double *, double, double, double, double *, double *, double *);
void   TC_RHS (const Data *, Data_Arr, double *,
               double **, double, int, int, Grid *);

void   GetGradient (double ***, double **, int, int, Grid *);
