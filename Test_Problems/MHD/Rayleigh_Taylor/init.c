/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Rayleigh-Taylor instability setup for hydro or MHD.

  Sets the initial condition for a Rayleigh-Taylor instability 
  problem in 2 or 3 dimensions.
  The initial condition consists of an interface separating two 
  fluids with different densities in hydrostatic balance:
  \f[
     \rho(y) = \left\{\begin{array}{ll}
          1    & \quad{\rm for} \quad y \le 0 \\ \noalign{\medskip}
          \eta & \quad{\rm for} \quad y > 0
     \end{array}\right.
     \,,\qquad 
     p = p_0 + \rho y g
  \f]
  where \f$ \eta \f$ is the density of the fluid on top, \c g is the 
  (constant) gravity pointing in the negative \c y direction.
  The value of the pressure at the interface (\c p0) is chosen in such 
  a way that the sound speed in the light fluid is 1.
  Likewise, density is normalized to the value in the light fluid.
  The horizontal extent of the computational domain defines the unit length:
  \c Lx=1.

  For magnetized setups, the magnetic field is purely horizontal:
  \f[
     \vec{B} =  \chi B_c\hvec{i} \,,\qquad
     B_c \equiv \sqrt{(\rho_{\rm hi} - \rho_{\rm lo})L|g|}\quad \to\quad
              \sqrt{(\eta - 1.0)|g|} 
  \f]
  where \c Bc is the critical magnetic field above which perturbations
  parallel to the magnetic field are suppressed (see Boyd 
  & Sanderson, page 99) while the value of \f$\chi\f$ is user-supplied
  (note that a factor \f$ 1/\sqrt{4\pi} \f$ needs to be incorporated 
   when initializing \c B in code units).

  The system is destabilized by perturbing the vertical velocity in
  proximity of the interface using a single mode (in 2D) or a Gaussian
  perturbation (in 3D).
  Aleternatively, a random perturbation can be used by setting
  \c USE_RANDOM_PERTURBATION to \c YES in your \c definitions.h.
  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[ETA]</tt>:  sets the density of the upper fluid;
  - <tt>g_inputParam[GRAV]</tt>: sets gravity (must be < 0);
  - <tt>g_inputParam[CHI]</tt>:  sets the magnetic field strength in unit 
                                 of the critical magnetic field \c Bc;

  \note Increasing the value of \c p0 by a factor \c q^2 is equivalent
        to increase gravity of the same factor and decrease velocity  
        and time by a factor \c q.

  The Rayleigh-Taylor setup has been tested with the following 
  configurations:

  <CENTER>
  Conf.|PHYSICS| GEOMETRY |DIM| T. STEP.|INTERP.  |divB|  Notes
  -----|-------|----------|---|---------| --------|----|----------
   #01 |  HD   |CARTESIAN | 2 | HANCOCK |LINEAR   | -  |  -
   #02 |  HD   |CARTESIAN | 2 |  RK3    |MP5_FD   | -  |  -
   #03 |  HD   |CARTESIAN | 2 |  RK3    |PARABOLIC| -  |  -
   #04 |  HD   |CARTESIAN | 2 |  RK2    |LINEAR   | -  |  (*)   
   #05 |  MHD  |CARTESIAN | 2 | HANCOCK |LINEAR   | CT |  -
   #06 |  MHD  |CARTESIAN | 2 |  RK3    |MP5_FD   |GLM |  - 
   #07 |  MHD  |CARTESIAN | 2 |  RK3    |PARABOLIC| CT |  - 
   #08 |  MHD  |CARTESIAN | 2 | ChTr    |PARABOLIC|GLM |  (*)     
   #09 |  MHD  |CARTESIAN | 3 |  RK2    |LINEAR   |GLM |  - 
   #10 |  MHD  |CARTESIAN | 3 | HANCOCK |LINEAR   |GLM | [Mig12](*)   
  </CENTER>
 
  (*) Setups for AMR-Chombo

  \image html mhd_rt.01.jpg "Density plot at the end of configuration #01"

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 06, 2014 

  \b References: \n
     - [Mig12]: Section 5.5 of Mignone et al., ApJS (2012)  198:7 \n
     - [SG07]: Stone & Gardiner, Physics of Fluids (2007), 19
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double x=x1, y=x2, z=x3;
  double rnd, Bc, g;
  
  g   = g_inputParam[GRAV];

  if (y < 0.0) v[RHO] = 1.0;
  else         v[RHO] = g_inputParam[ETA];

  v[PRS] = 1.0/g_gamma + v[RHO]*g*y;  /* Hydrostatic balance */
  v[VX1] = v[VX2] = v[VX3] = 0.0;

/* -----------------------------------
    Add perturbation
   ----------------------------------- */

  #if USE_RANDOM_PERTURBATION == YES
   rnd    = RandomNumber();
   v[VX2] = 1.e-2*rnd*exp(-y*y*200.0);
  #else
   #if DIMENSIONS == 2
    v[VX2] = -1.e-2*(1.0 + cos(2.0*CONST_PI*x))*exp(-y*y*200.0);
   #elif DIMENSIONS == 3
    rnd = sqrt(x*x + z*z);
    v[VX2] = -1.e-2*exp(-rnd*rnd/(0.2*0.2))/cosh(y*y/(0.1*0.1)); 
   #endif
  #endif

/* -----------------------------------
    Set magnetic field in units of Bc
   ----------------------------------- */

  #if PHYSICS == MHD
   Bc = sqrt((g_inputParam[ETA] - 1.0)*fabs(g)); /* Critical field */
   v[BX1] = g_inputParam[CHI]*Bc/sqrt(4.0*CONST_PI);
   v[BX2] = 0.0;
   v[BX3] = 0.0;
   v[AX1] = v[AX2] = 0.0;
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

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*  
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = g_inputParam[GRAV];
  g[KDIR] = 0.0;
  #ifdef CTU 
   if      (x2 > g_domEnd[JDIR]) g[JDIR] *= -1.0;
   else if (x2 < g_domBeg[JDIR]) g[JDIR] *= -1.0;
  #endif
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *************************************************************************** */
{
  #ifdef CTU
   if      (x2 < g_domBeg[JDIR]) x2 = 2.0*g_domBeg[JDIR] - x2;
   else if (x2 > g_domEnd[JDIR]) x2 = 2.0*g_domEnd[JDIR] - x2;
  #endif

  return -g_inputParam[GRAV]*x2;
}
#endif
