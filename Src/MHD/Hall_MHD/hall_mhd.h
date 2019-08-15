/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Header file for Hall MHD.

  Provides macros, function prototypes and structure definitions for 
  the Hall-MHD module.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Feb 23, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
  

void   HallMHD_WhistlerSpeed (const State *, int, int, Grid *);
double HallMHD_ne(double *);
void   HallMHD_Flux (const Sweep*, Data_Arr, double **, double **, int, int, Grid *);
