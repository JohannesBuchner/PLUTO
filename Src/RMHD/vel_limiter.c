/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Velocity limiter for relativistic hydro or MHD.
  
  Flatten reconstruction when either the left or right 
  reconstructed velocity values exceeds one.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 11, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void VelocityLimiter (double *v, double *vp, double *vm)
/*!
 * Check whether the total reconstructed velocity is > 1
 * If a superluminal value occurs, flatten distribution.
 *
 *********************************************************************** */
{
#if RECONSTRUCT_4VEL == NO && COMPONENTS > 1
   int    nv;
   double v2m, v2p;

   v2m = EXPAND(vm[VX1]*vm[VX1],  + vm[VX2]*vm[VX2], + vm[VX3]*vm[VX3]);
   v2p = EXPAND(vp[VX1]*vp[VX1],  + vp[VX2]*vp[VX2], + vp[VX3]*vp[VX3]);
   if (v2m >= 1.0 || v2p >= 1.0){
     for (nv = NVAR; nv--;  ) vm[nv] = vp[nv] = v[nv];  
   }
#endif
}
