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
  double mu, T, B2_cgs, nH;

  mu   = 0.5;
  T    = v[PRS]/v[RHO]*mu*KELVIN;

  *kpar = 5.6e-7*T*T*sqrt(T);
  #if PHYSICS == MHD
   B2_cgs  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]) + 1.e-12;
   B2_cgs *= 4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
   nH    = v[RHO]*UNIT_DENSITY/CONST_mp;
   *knor = 3.3e-16*nH*nH/(sqrt(T)*B2_cgs);
  #endif
  
  *phi = 0.3;
}
