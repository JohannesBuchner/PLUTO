/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic magnetized Kelvin-Helmholtz problem.

  Implements the initial condition for a relativistic KH flow in
  2D Cartesian coordinates as in section 6.6 of Mignone et al. (ApJS, 2012).
  Density and pressure are initially constant and equal to 
  \f$\rho = 1\f$ and \f$p = p_0\f$ where the latter value is retrieved 
  from the definition of the Mach number \f$M = v/c_s\f$.
  The sound speed is \f$c_s=\sqrt{\Gamma p/w}\f$ for the \c IDEAL equation 
  of state and \f$c_s = \sqrt{\frac{p}{3w}\frac{5w - 8p}{w - p}}\f$ for the 
  \c TAUB EoS (see Eq. [17], [31] and [32] of Mignone \& McKinney 
  [MNRAS, 2007]), where \f$ w = \rho h\f$ is the gas enthalpy.
  
  The velocity is parallel to the \c x axis with a hyperbolic tangent 
  shear profile.
  Magnetic field has the form 
  \f[
     \vec{B} = \left(\sqrt{2\sigma_{\rm pol}p_0},\, 0,\, 
                     \sqrt{2\sigma_{\rm tor}p_0}\right)
  \f]
  where \f$\sigma_{\rm pol}\f$ and \f$\sigma_{\rm tor}\f$ control the strength
  of the poloidal and toroidal component, respectively.
  
  Periodic boundary conditions are assumed in the horizontal direction while
  zero-gradient (outflow) conditions hold at the lower and upper 
  boundaries in the \c y direction.
  
  The presence of a magnetic field gives rise to a filamentary structure
  and to the formation of a number of vortices arranged symmetrically 
  with respect to the central point. 
  The figure below shows density for configuration #03 using AMR
  with 4 levels of refinement at \c t = 5 

  The control parameters for this problem are:
  
  -# <tt>g_inputParam[SIGMA_TOR]</tt>:  control the strength of the toroidal 
     (\c z) component of magnetic field;
  -# <tt>g_inputParam[SIGMA_POL]</tt>:  control the strength of the poloidal
     (\c x) component of magnetic field;
  -# <tt>g_inputParam[VEL0]</tt>: the flow velocity (in units of the speed of 
     light);
  -# <tt>g_inputParam[MACH]</tt>: the flow Mach number used to recover the
     pressure value;

  
  \image html rmhd_kh.03.jpg "...".

  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics", Mignone et al., ApJS (2012) 198, 7.
     - "Equation of state in relativistic magnetohydrodynamics: 
        variable versus constant adiabatic index",
        Mignone \& McKinney, MNRAS (2007) 378, 1118.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 22, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double x,y, arg, eps, kx;
  double alpha, beta;
  static double pr0;

  if (first_call == 1){
    arg = g_inputParam[MACH]/g_inputParam[VEL0];
    arg = arg*arg;
    #if EOS == IDEAL
     g_gamma = 4./3.;
     pr0 = (g_gamma - 1.0)/((g_gamma - 1.0)*arg - 1.0)/g_gamma;
    #else
    { 
      double a,b,c;  /* See Eq. 32 of Mignone \& McKinney (MNRAS 2007) */
      a = 15.0 - 6.0*arg; b = 24.0 - 10.0*arg; c = 9.0;
      arg = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;
      pr0 = 2.0/3.0*sqrt(arg*arg/(1.0 - arg*arg));
    }
    #endif
    first_call = 0;
  }

  x = x1;
  y = x2;

  kx    = 2.0*CONST_PI*x1;
  alpha = 1.0/100.0;
  beta  = 1.0/10.0;

  arg = y*y/(beta*beta);
  eps = 0.01*g_inputParam[VEL0];

  v[RHO] = 1.0;
  v[VX1] = -g_inputParam[VEL0]*tanh(y/alpha);
  v[VX2] =  eps*0.5*(sin(kx) - sin(-kx))*exp(-arg);
  v[VX3] = 0.0;

/* --------------------------------------------
    find pressure values given the Mach number
   -------------------------------------------- */

  v[PRS] = pr0;

/* v[PRS] = 20.0;  */

  v[TRC] = (y < 0.0 ? 1.0:-1.0);

  #if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = sqrt(2.0*v[PRS]*g_inputParam[SIGMA_POL]);
   v[BX2] = 0.0;
   v[BX3] = sqrt(2.0*v[PRS]*g_inputParam[SIGMA_TOR]);

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = y*v[BX1];

  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *
 *********************************************************************** */
{ 
}
