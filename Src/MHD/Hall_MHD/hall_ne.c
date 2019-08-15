/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the 'ne' term in Hall-MHD.

  
  \authors A. Mignone \n
           E. Striani \n
           B. Vaidya

  \date    April 04, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
double HallMHD_ne(double *v)
/*!
 * Compute the ne term for Hall MHD
 * variables, coordinates and currents.
 *
 * \param [in]  v    array of primitive variables
 *
 *********************************************************************** */
{
  double e_ne = 1.0;
  return (e_ne*v[RHO]);
}
