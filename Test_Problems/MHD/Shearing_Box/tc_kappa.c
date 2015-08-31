/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the thermal conduction coefficients.

  Use this function to supply the thermal conduction coefficients \f$ 
  \kappa_\| \f$ and \f$ \kappa_\bot \f$ along and across magnetic 
  field lines and the \f$ \phi \f$ parameter used to control the magnitude 
  of the saturated flux \f$ F_{\rm sat} = 5\phi\rho c_{\rm iso}^3 \f$.
  To exclude saturation, simply set \f$ \phi \f$ to a very large number.
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos  
  \date    Oct 31, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void TC_kappa(double *v, double x1, double x2, double x3,
              double *kpar, double *knor, double *phi)
/*!
 * Compute thermal conduction coefficients.
 *
 * \param [in] v   array of primitive variables
 * \param [in] x1  coordinate in the X1 direction
 * \param [in] x2  coordinate in the X2 direction
 * \param [in] x3  coordinate in the X3 direction
 * \param [out] kpar pointer to the conduction coefficient 
 *                   \f$ \kappa_\parallel \f$ in the direction of magnetic
 *                     field
 * \param [out] knor pointer to the conduction coefficient 
 *                   \f$ \kappa_\perp \f$ perpendicular to magnetic field 
 * \param [out] phi    pointer to the parameter \f$ \phi \f$ controlling the
 *                     magnitude of the saturated flux.
 *    
 *********************************************************************** */
{
  *kpar = g_inputParam[KAPPA]*v[RHO];
  *knor = *kpar;
  *phi  = 1.e40;
}
