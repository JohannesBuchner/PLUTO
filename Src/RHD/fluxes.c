/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the flux for the relativistic hydro (RHD) equations.

  Compute the flux of the conservative RHD equations in the direction 
  given by ::g_dir.
  This function defines the component of the hyperbolic flux tensor 
  of the standard RHD equations.\n
  In what follows:
  - \c VXn, \c MXn are the velocity, momentum components in the direction 
    given by ::g_dir (normal, \c "n")
  - \c VXt, \c MXt and \c VXb, \c MXb are the transverse components
    (tangent \c "t" and bi-tangent \c "b").

 \author A. Mignone (mignone@ph.unito.it)
 \date   Oct 12, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (const State *state, int beg, int end)
/*
 * \param [in,out]  state   Pointer to a state structure
 * \param [in]      beg     initial index of computation 
 * \param [in]      end     final   index of computation
 *
 *********************************************************************** */
{
  int    nv, i;
  double *u, *v, *fx, vn;

  for (i = beg ; i <= end; i++) {
    v  = state->v[i];
    u  = state->u[i];
    fx = state->flux[i];

    vn = v[VXn];

    fx[RHO]  = u[RHO]*vn;
    EXPAND(fx[MX1] = u[MX1]*vn;  ,
           fx[MX2] = u[MX2]*vn;  ,
           fx[MX3] = u[MX3]*vn;)
    fx[ENG] = u[MXn];
    state->prs[i] = v[PRS];
  }
}
