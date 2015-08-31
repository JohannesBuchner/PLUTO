/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the hydro (HD) flux.                                             

  Compute the flux of the conservative HD equations in the direction 
  given by ::g_dir.
  This function defines the component of the hyperbolic flux tensor 
  of the standard HD equations.\n
  In what follows:
  - \c VXn, \c MXn are the velocity, momentum components in the direction 
    given by ::g_dir (normal, \c "n")
  - \c VXt, \c MXt and \c VXb, \c MXb are the transverse components
    (tangent \c "t" and bi-tangent \c "b").

 \author A. Mignone (mignone@ph.unito.it)
 \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (double **u, double **w, double *a2, double **fx, double *p, 
           int beg, int end)
/*!
 * \param [in]    u    1D array of conserved quantities
 * \param [in]    w    1D array of primitive quantities
 * \param [in]   a2    1D array of sound speeds
 * \param [out]  fx    1D array of fluxes (total pressure excluded)
 * \param [out]   p    1D array of pressure values
 * \param [in]   beg   initial index of computation 
 * \param [in]   end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   nv, ii;

  for (ii = beg; ii <= end; ii++) {
    fx[ii][RHO] = u[ii][MXn];
    EXPAND(fx[ii][MX1] = u[ii][MX1]*w[ii][VXn]; ,
           fx[ii][MX2] = u[ii][MX2]*w[ii][VXn]; ,
           fx[ii][MX3] = u[ii][MX3]*w[ii][VXn];)
#if HAVE_ENERGY
    p[ii] = w[ii][PRS];
    fx[ii][ENG] = (u[ii][ENG] + w[ii][PRS])*w[ii][VXn];
#elif EOS == ISOTHERMAL
    p[ii] = a2[ii]*w[ii][RHO];
#endif
/*
#if DUST == YES
    fx[ii][RHO_D] = u[ii][MXn_D];
    EXPAND(fx[ii][MX1_D] = u[ii][MX1_D]*w[ii][VXn_D]; ,
           fx[ii][MX2_D] = u[ii][MX2_D]*w[ii][VXn_D]; ,
           fx[ii][MX3_D] = u[ii][MX3_D]*w[ii][VXn_D];)
#endif
*/
  }
}
