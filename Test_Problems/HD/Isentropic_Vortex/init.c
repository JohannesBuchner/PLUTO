/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Isentropic vortex problem.

  The isentropic vortex is a smooth exact solution of the 2D Euler 
  equations.
  It consists of a single vortex centered at \c (5,5) in pressure 
  equilibrium:
  \f[
    \Big(\delta v_x\,, \delta v_y\Big) 
    = - (y_c,\,x_c)\frac{\epsilon}{2\pi}
                   \exp\left(\frac{1-r_c^2}{2}\right) \,,\quad
    T = \frac{p}{\rho} 
      = 1 - \frac{(\Gamma-1)\epsilon^2}{8\Gamma\pi^2}\exp(1-r_c^2)\,,\quad
    s = \frac{p}{\rho^\Gamma}    
  \f] 
  where \f$ (x_c,\,y_c) = (x-5,\,y-5),\, r_c=\sqrt{x_c^2+y_c^2},\,\epsilon=5\f$
  while \f$\Gamma = 1.4\f$ is the adiabatic index.
  A constant entropy is used with value \c s=1.
  The vortex shifts along the main diagonal of the computational domain 
  with uniform velocity \c (1,1) and returns in its original position 
  after \c t=10.
  The domain is assumed to be periodic and its size must be taken large 
  enough to ensure there is no interaction between off-domain vortices.

  The vortex problem is often used as computational benchmark to test 
  the accuracy and dissipation of a numerical scheme in reproducing 
  the vortex structure after several revolutions.
  
  There're no input parameters for this problem.

  \image html hd_isentropic_vortex.13.jpg "Density cut at y=5 after 10 revolutions using PPM and Finite difference WENOZ and PPM (conf #01 and #03)."

  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 3, 2014

  \b References
     - C.W. Shu, "High-order Finite Difference and Finite Volume WENO 
       Schemes and Discontinuous Galerkin Methods for CFD",
       ICASE Report No.2001-11,  NASA/CR-2001-210865
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double s, T, r2, xt, yt, eps;     
  double k0;

  g_gamma = 1.4;
  
  xt = x1 - 5.0;
  yt = x2 - 5.0;
   
  r2  = xt*xt + yt*yt;
  
  eps = 5.0;
  
  s   = 1.0;
  T   = 1.0 - (g_gamma - 1.0)*eps*eps/(8.0*g_gamma*CONST_PI*CONST_PI)*exp(1.0 - r2);
  
  k0     = eps/(2.0*CONST_PI)*exp(0.5*(1.0 - r2));
  us[RHO] = pow(T/s, 1.0/(g_gamma - 1.0));
  us[VX1] = 1.0 - k0*yt;
  us[VX2] = 1.0 + k0*xt;
  us[PRS] = T*us[RHO];
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *********************************************************************** */
{}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{ }
