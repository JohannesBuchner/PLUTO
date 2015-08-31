/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Sedov-Taylor blast wave.

  Set the initial condition for a Sedov-Taylor blast wave problem.
  Ambient density and pressure are constant everywhere and equal to 
  \c rho0 and \c p0. 
  While \c rho0 is an input parameter, \c p0=1.e-5.
  A large amount of energy is deposited in a small volume (a
  line, circle or sphere depending on geometry and dimensions) 
  consisting of 2 or 3 computational zones.

  The input parameters read from pluto.ini are labeled as:
  - <tt>g_inputParam[ENRG0]</tt>: initial energy of the blast wave
    (\e not energy density);
  - <tt>g_inputParam[DNST0]</tt>: initial constant density;
  - <tt>g_inputParam[GAMMA]</tt>: specific heat ratio.

  The available configurations refer to:
  -# Cartesian (1D)
  -# Cylindrical (1D)
  -# Spherical (1D)
  -# Cartesian (2D, equivalent to cylindrical)
  -# Cartesian (3D, equivalent to spherical)

  \image html hd_sedov.02.jpg "Density plot at t=0.5 for configuration #02 (cylindrical blast)"

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 08, 2014
  
  \b References:
     -  L.I. Sedov, "Similarity and Dimensional Methods in Mechanics", 
        Academic Press, New York, 1959.
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
  double dr, vol, r;

  g_gamma = g_inputParam[GAMMA];

/* --------------------------------------------------
    dr is the size of the initial energy deposition 
    region: 2 ghost zones.
   -------------------------------------------------- */

  #if DIMENSIONS == 1
   dr = 2.0*(g_domEnd[IDIR]-g_domBeg[IDIR])/(double)NX1;
  #else
   dr = 3.5*(g_domEnd[IDIR]-g_domBeg[IDIR])/(double)NX1;
  #endif

/* ---------------------------------------
     compute region volume 
   --------------------------------------- */

  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 1)
   vol = 2.0*dr;
  #elif (GEOMETRY == CYLINDRICAL && DIMENSIONS == 1)|| \
        (GEOMETRY == CARTESIAN   && DIMENSIONS == 2)
   vol = CONST_PI*dr*dr;
  #elif (GEOMETRY == SPHERICAL   && DIMENSIONS == 1)|| \
        (GEOMETRY == CYLINDRICAL && DIMENSIONS == 2)|| \
        (GEOMETRY == CARTESIAN   && DIMENSIONS == 3)
   vol = 4.0/3.0*CONST_PI*dr*dr*dr;
  #else
   print1 ("! Init: geometrical configuration not allowed\n");
   QUIT_PLUTO(1);
  #endif

  r = EXPAND(x1*x1, + x2*x2, +x3*x3);
  r = sqrt(r);

  us[RHO] = g_inputParam[DNST0];
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;

  if (r <= dr)  us[PRS] =  (g_gamma - 1.0)*g_inputParam[ENRG0]/vol;
  else          us[PRS] = 1.e-5;

  us[TRC]  = 0.0;
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
{ }

