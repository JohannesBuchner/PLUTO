/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the components of the diagonal resistive tensor. 

  Use this function to supply the resistivity in the three directions
  \f$ \eta_{x1}\f$, \f$ \eta_{x2}\f$ and \f$ \eta_{x3}\f$.
  
  \authors T. Matsakos \n
           A. Mignone (mignone@ph.unito.it)\n
  \date    March 22, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Resistive_eta(double *v, double x1, double x2, double x3,
                   double *J, double *eta)
/*!
 * Compute the resistive tensor components as function of the primitive
 * variables, coordinates and currents.
 *
 * \param [in]  v    array of primitive variables
 * \param [in]  x1   coordinate in the X1 direction
 * \param [in]  x2   coordinate in the X2 direction
 * \param [in]  x3   coordinate in the X3 direction
 * \param [in]  J    current components, J[IDIR], J[JDIR], J[KDIR]
 * \param [out] eta  an array containing the three components of
 *                   \f$ \tens{\eta}\f$.
 *
 *********************************************************************** */
{
  eta[IDIR] = 1.0;
  eta[JDIR] = 1.0;
  eta[KDIR] = 1.0;
}
