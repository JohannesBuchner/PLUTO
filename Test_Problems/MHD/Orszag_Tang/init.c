/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Orszag-Tang MHD vortex.

  The Orszag Tang vortex system describes a doubly periodic fluid 
  configuration leading to two-dimensional supersonic MHD turbulence.
  Although an analytical solution is not known, its simple and reproducible 
  set of initial conditions has made it a widespread benchmark for 
   inter-scheme comparison.

  The computational domain is the periodic box \f$[0,2\pi]^D\f$ where 
  \c D is the number of spatial dimensions.
  In 2D, the initial condition is given by
  \f[
     \vec{v} = \left(-\sin y,\, \sin x, 0\right) \,,\qquad
     \vec{B} = \left(-\sin y,\, \sin 2x, 0\right) \,,\qquad
     \rho = 25/9\,,\qquad
     p    = 5/3
  \f]
  
  This test problem does not have any input parameter.

  A snapshot of the solution on a \c 512x512 grid is shown below.

  \image html mhd_ot.02.jpg "Density at t=3.1 (configuration #02)."

  \author A. Mignone (mignone@ph.unito.it)
  \date   April 13, 2014

  \b References
     - "Comparison of some Flux Corrected Transport and TVD numerical 
        schemes for hydrodynamic and magnetohydrodynamic problems",
        Toth & Odstrcil, JCP (1996) 128, 82
     - "High-order conservative finite difference GLM-MHD schemes for 
        cell-centered MHD", Mignone, Tzeferacos & Bodo, JCP (2010) 229, 5896.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double x=x1, y=x2 ,z=x3;

  us[VX1] = - sin(y);
  us[VX2] =   sin(x);
  us[VX3] = 0.0;
  us[BX1] = - sin(y);
  us[BX2] =   sin(2.0*x);
  us[BX3] = 0.0;
  us[RHO] = 25./9.;
  #if EOS != ISOTHERMAL && EOS != BAROTROPIC
   us[PRS] = 5.0/3.0;
  #endif
  us[TRC] = (y>CONST_PI)*1.0;

  us[AX1] = 0.0;
  us[AX2] = 0.0;
  us[AX3] = cos(y) + 0.5*cos(2.0*x);

 #if DIMENSIONS == 3
 {
   double c0 = 0.8;
   us[VX2] = - sin(z);
   us[VX3] =   sin(y);
   us[VX1] =   0.0;
  
   us[BX2] = c0*( - 2.0*sin(2.0*z) + sin(x));
   us[BX3] = c0*(       sin(x) + sin(y));
   us[BX1] = c0*(       sin(y) + sin(z));
   us[RHO] = 25./9.;
   us[PRS] = 5.0/3.0;

   us[AX2] = c0*( cos(z) - cos(x));
   us[AX3] = c0*(-cos(y) + cos(x));
   us[AX1] = c0*( cos(y) + cos(2.0*z));

  }
  #endif

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
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
