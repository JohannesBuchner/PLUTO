/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief 2D Blast wave problem with thermal conduction.

  Sets the initial conditions for a 2D Blast wave problem with thermal conduction.
  For the HD case:
  \f[
    \rho = \rho_{\rm out} + \frac{\rho_{\rm in} - \rho_{\rm out}}{\cosh\left[10(r/r_0)^{10}\right]}\,,
    \quad\quad\quad
    T = \left\{\begin{array}{ll}
          T_{\rm in}  & \quad{\rm for} \quad r \le r_0 \\ \noalign{\medskip}
          T_{\rm out} & \quad{\rm for} \quad r > r_0
     \end{array}\right.\,,
  \f]
  where \f$r_0\f$ is the cloud radius.

  Input Parameters are that are read from \c pluto.ini are 
  - <tt>g_inputParam[T_IN], g_inputParam[T_OUT]</tt>:     Temperature inside and outside the circle (in K);
  - <tt>g_inputParam[RHO_IN], g_inputParam[RHO_OUT]</tt>: Density inside and outside (in dimensionless units);
  - <tt>g_inputParam[BMAG]</tt>:                          Magnetic field strength (in Gauss);
  - <tt>g_inputParam[THETA]</tt>:                         Orientation of the field (in degrees);

  Configurations:
  - #01-04 have the same initial condition and are done either with
    an explicit time stepping or STS, HD and MHD and do not show
    evidence for any numerical artifact.
  - #05-06, on the other hand, show that STS suffers from some kind of
    unstable behavior due to the flux limiter switching from classical
    to saturated regimes. Only small CFL (0.1 or less) or larger values
    of STS_NU (e.g 0.05) mitigate the problem.
    Future improvement (RKC/RKL ?) should address this issue.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 27, 2015 

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  static int first_call=1;
  double r, r0, mu, T, prs_ref, prof;

  mu  = 1.26;
  g_gamma = 5.0/3.0;        

  prs_ref = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
  
/* ----------------------------------------------
             Use c.g.s units 
   ---------------------------------------------- */

  r    = sqrt(EXPAND(x1*x1, + x2*x2, + x3*x3));
  r0   = 1.0;   /* -- cloud radius -- */
  prof = 1.0/cosh(10.0*pow(r/r0,10));

  us[VX1] = us[VX2] = 0.0;
  T       = g_inputParam[T_OUT]   + (g_inputParam[T_IN]   - g_inputParam[T_OUT])*(r <= r0);
  us[RHO] = g_inputParam[RHO_OUT] + (g_inputParam[RHO_IN] - g_inputParam[RHO_OUT])*prof;

  us[PRS] = T*us[RHO]/KELVIN;

  #if PHYSICS == MHD
   us[BX1] = g_inputParam[BMAG]*cos(g_inputParam[THETA]*CONST_PI/180.0);
   us[BX2] = g_inputParam[BMAG]*sin(g_inputParam[THETA]*CONST_PI/180.0);
   us[BX3] = 0.0;

   us[BX1] /= sqrt(prs_ref*4.0*CONST_PI);
   us[BX2] /= sqrt(prs_ref*4.0*CONST_PI);
   us[BX3] /= sqrt(prs_ref*4.0*CONST_PI);

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = x2*us[BX1] - x1*us[BX2];
  #endif

  #ifdef GLM_MHD 
   us[PSI_GLM] = 0.0;
  #endif

}   
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{ }

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  static int  first_call = 1;
  int   i, j, k, nv;
  static double vin[256];

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (first_call){
      Init(vin, 0.0, -10.0, 0.0);
      first_call = 0;
    }
    if (box->vpos == CENTER){
      for (nv = 0; nv < NVAR; nv++ ) BOX_LOOP(box,k,j,i){
        d->Vc[nv][k][j][i] = vin[nv];
      }
    }
  }
}

