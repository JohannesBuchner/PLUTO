#include "pluto.h"
#include "cooling_defs.h"

/* ********************************************************************* */
double GetMaxRate (double *v0, double *k1, double T0)
/*!
 * Return an estimate of the maximum rate (dimension 1/time) in the  
 * chemical network. 
 * This will serve as a "stiffness" detector in the main ode integrator.
 *   
 * For integration to be carried explicitly all the time,
 * return a small value (1.e-12).
 *
 *********************************************************************** */
{
  int  nv;
  double   ev1[NIONS+1];
  double   scrh, maxrate;
  static double **J1;

  if (J1 == NULL)  {
    J1 = ARRAY_2D(NIONS+1, NIONS+1, double);
  }

  maxrate = 0.0;
  for (nv = 0; nv < NIONS - Fe_IONS; nv++) {  /* -- exclude Iron -- */

  /* ---------------------------------------------- 
       if the initial concentration is close to 0,
       do not take Crate as a stiffness indicator
     ---------------------------------------------- */

    if (v0[NFLX + nv] < 1.e-6) continue;

    maxrate = MAX(maxrate, CoolCoeffs.Crate[nv]);
  } 
  maxrate *= UNIT_LENGTH/UNIT_VELOCITY;

/*
  Numerical_Jacobian    (v0, J1);
  lmax = Decompose(J1[0], NIONS + 1, ev1);
  Radiat(v0, k1);
  scrh = maxrate/fabs(lmax) - 1.0;
  if (fabs(scrh) > 5.e-2){
     printf ("Crate are not the max eig, %12.6e  %12.6e\n", 
              maxrate, fabs(lmax));
     QUIT_PLUTO(1);
  }
  return(fabs(lmax));
*/

  return (maxrate);
}
