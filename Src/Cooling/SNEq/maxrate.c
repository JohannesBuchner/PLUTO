#include "pluto.h"

/* ********************************************************** */
double GetMaxRate (double *v0, double *k1, double T0)
/*
 *
 *  PURPOSE:
 *
 *    return an estimate of the maximum rate (dimension 1/time) 
 *    in the chemical network. This will serve as a
 *    "stiffness" detector in the main ode integrator.
 *   
 *    For integration to be carried explicitly all the time,
 *    return a small value (1.e-12).
 *
 ************************************************************ */
{
  return (1.e-12);
}

