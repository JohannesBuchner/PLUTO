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
 \date   Oct 12, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (const State *state, int beg, int end)
/*!
 * \param [in,out]  state   Pointer to a state structure
 * \param [in]      beg     initial index of computation 
 * \param [in]      end     final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   nv, i;
  double *u, *v, *flux;

  for (i = beg; i <= end; i++) {
    u    = state->u[i];
    v    = state->v[i];
    flux = state->flux[i];

    flux[RHO] = u[MXn];
    EXPAND(flux[MX1] = u[MX1]*v[VXn]; ,
           flux[MX2] = u[MX2]*v[VXn]; ,
           flux[MX3] = u[MX3]*v[VXn];)
#if HAVE_ENERGY
    state->prs[i] =  v[PRS];
    flux[ENG]     = (u[ENG] + v[PRS])*v[VXn];
#elif EOS == ISOTHERMAL
    state->prs[i] = state->a2[i]*v[RHO];
#endif
/*
#if DUST_FLUID == YES
    fx[ii][RHO_D] = u[ii][MXn_D];
    EXPAND(fx[ii][MX1_D] = u[ii][MX1_D]*w[ii][VXn_D]; ,
           fx[ii][MX2_D] = u[ii][MX2_D]*w[ii][VXn_D]; ,
           fx[ii][MX3_D] = u[ii][MX3_D]*w[ii][VXn_D];)
#endif
*/
  }
}
