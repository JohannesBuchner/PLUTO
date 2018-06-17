/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Check if primitive states are physically admissible.

  This function is called at the end of the Hancock or Characteristic
  Tracing functions to verify that the L/R states at the half-step
  are physically admissible.
  
  \authors A. Mignone (mignone@ph.unito.it)

  \date   Sept 7, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* *********************************************************************** */
void CheckPrimStates(double **vM, double **vP, double **v0,  int beg, int end)
/*!
 *
 *
 ************************************************************************** */
{
#if PHYSICS != ADVECTION
  int    i, nv, switch_to_1st;
  double scrhm, scrhp;
  double *ac, *ap, *am;

  for (i = beg; i <= end; i++){
  
    switch_to_1st = 0;

    ac = v0[i];
    am = vM[i];
    ap = vP[i];

  /*  ----  Prevent unphysical states by reverting to first
            order in time and space,  i.e.  set dw = 0      ----  */

#if HAVE_ENERGY
    switch_to_1st = (ap[PRS] < 0.0) || (am[PRS] < 0.0) ;
#endif
    switch_to_1st = switch_to_1st || 
                    (ap[RHO] < 0.0) || (am[RHO] < 0.0) ;

    /*  ----  Check for superluminal velocities  ---- */

#if (PHYSICS == RHD) || (PHYSICS == RMHD)
    #if RECONSTRUCT_4VEL == NO
    scrhm = EXPAND(am[VX1]*am[VX1], + am[VX2]*am[VX2], + am[VX3]*am[VX3]);
    scrhp = EXPAND(ap[VX1]*ap[VX1], + ap[VX2]*ap[VX2], + ap[VX3]*ap[VX3]);
    switch_to_1st = switch_to_1st || (scrhm >= 1.0);
    switch_to_1st = switch_to_1st || (scrhp >= 1.0);
    #endif
#endif

    if (switch_to_1st){
#ifdef STAGGERED_MHD 
       scrhp = ap[BXn];
       scrhm = am[BXn];
#endif

      NVAR_LOOP(nv) am[nv] = ap[nv] = ac[nv];
 
#ifdef STAGGERED_MHD
       ap[BXn] = scrhp;
       am[BXn] = scrhm;
#endif
      
    }
  }
#endif
}


/*
#if GEOMETRY == CYLINDRICAL
if (NSWEEP == 1) {
 for (nv = 0; nv < NVAR; nv++){
   if (nv == VX) scrh = fabs(v1[3][nv] + v1[4][nv]);
   else          scrh = fabs(v1[3][nv] - v1[4][nv]);
   
   #if PHYSICS == RMHD 
    if (nv == BZ) scrh = fabs(v1[3][nv] + v1[4][nv]);
   #endif
   
   if (scrh > 1.e-8){
     printf ("symmetry violated, z : %d,  var: %d\n", *nyp, nv);
     SHOW(rhs,4);
     SHOW(rhs,3);
     printf (" --- centered:\n");
     SHOW(vv,0); SHOW(vv,1); SHOW(vv,2); SHOW(vv,3); 
     SHOW(vv,4); SHOW(vv,5); SHOW(vv,6); SHOW(vv,7);
     printf (" --- edge\n");
     SHOW(vl,3); SHOW(vr,3);
     SHOW(vr,2); SHOW(vl,4);

     printf ("Source: \n");
     SHOW(src,3);SHOW(src,4);
     scrhp =  fl[4][MXn]*GG->A[4]/GG->dV[4] - src[4][MXn];
     scrhm = -fr[2][MXn]*GG->A[2]/GG->dV[3] - src[3][MXn];
     printf ("%12.6e  %12.6e\n",scrhp, scrhm);
     exit(1);
   }
 }
}
#endif
*/


