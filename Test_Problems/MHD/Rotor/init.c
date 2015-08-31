/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief MHD Rotor test problem.

  The rotor problem consists of a rapidly spinning cylinder embedded in 
  a static background medium with uniform density and pressure and
  a constant magnetic field along the x direction 
  \f$B_x = 5/\sqrt{4\pi}\f$.
  The cylinder rotates uniformly with constant angular velocity 
  \f$ \omega \f$ and has larger density:
  \f[
     (\rho, v_\phi) = \left\{\begin{array}{ll}
     (10, \omega r) & \qquad {\rm for}\quad r < r_0 
     \\ \noalign{\medskip}
     (1+9f,f\omega r_0) & \qquad {\rm for}\quad r_0 \le r \le r_1
     \\ \noalign{\medskip}
     (1,0) & \qquad {\rm otherwise}
   \end{array}\right.
  \f]
  Here \f$ r_0=0.1\f$ and \f$ = (r_1-r)/(r_1-r_0)\f$ is a taper function.
  The ideal equation of state with \f$\Gamma = 1.4\f$ is used.
  As the disk rotates, strong torsional Alfven waves form and propagate 
  outward carrying angular momentum from the disk to the ambient.

  A list of tested configurations is given in the following table:

  <CENTER>
  Conf.| GEOMETRY| divB  |BCK_FIELD| AMR
  -----|---------| ------|---------|-----
  #01  |CARTESIAN| CT    | NO      | NO
  #02  |POLAR    | CT    | NO      | NO
  #03  |POLAR    | 8W    | NO      | NO
  #04  |POLAR    | CT    | YES     | NO
  #05  |POLAR    | GLM   | NO      | NO
  #06  |CARTESIAN| GLM   | NO      | NO
  #07  |CARTESIAN| GLM   | NO      | YES
  #08  |CARTESIAN| 8W    | NO      | YES
  #09  |POLAR    | CT    | YES     | NO
  #10  |POLAR    | CT    | YES     | NO
  #11  |POLAR    | GLM   | YES     | YES
  </CENTER>
  
  A snapshot of the solution using static and AMR grid is given below.

  \image html mhd_rotor.01.jpg "Density map (in log scale) at t = 0.15 in Cartesian coordinates using 400 x 400 grid zones (configuration #01). Field lines are super-imposed."
  \image html mhd_rotor.08.jpg "Density map (in log scale) at t = 0.15 in Cartesian coordinates using a base grid of 64x64 and 4 levels of refinement (conf. #08)."

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 08, 2014

  \b Reference: 
     - "On the divergence-free condition in Godunov-type schemes 
        for ideal MHD: the upwind constrained transport method"
        P. Londrillo, L. Del Zanna,  JCP (2004), 195, 17
     - "PLUTO: A numerical code for computational astrophysics"
        Mignone et al., ApJS (2007) 170, 228 (see Sect. 5.4)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  double r, r0, r1, Bx, f, omega;

  g_gamma   = 1.4;
  omega = 20.0;
  Bx    = 5.0/sqrt(4.0*CONST_PI);
  r0    = 0.1;
  r1    = 0.115;

  #if GEOMETRY == CARTESIAN
   r = sqrt(x1*x1 + x2*x2);

   v[PRS] = 1.0;
   v[BX1] = Bx;
   v[BX2] = 0.0;

   f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX1] = -omega*x2;
     v[VX2] =  omega*x1;
   }else if (r < r1) {
     v[RHO] = 1.0 + 9.0*f;
     v[VX1] = -f*omega*x2*r0/r;
     v[VX2] =  f*omega*x1*r0/r;
   } else {
     v[RHO] = 1.0;
     v[VX1] = 0.0;
     v[VX2] = 0.0;
   }

   v[AX1] = v[AX2] = 0.0;
   v[AX3] = Bx*x2;
  #elif GEOMETRY == POLAR
   r = x1;

   v[BX1] =  Bx*cos(x2);
   v[BX2] = -Bx*sin(x2);
   v[VX1] = 0.0;
   v[PRS] = 1.0;

  f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX2] = omega*r;
   }else if (r < r1) {
     v[RHO] = 1.0 + 9.0*f;
     v[VX2] = f*omega*r0;
   } else {
     v[RHO] = 1.0;
     v[VX2] = 0.0;
   }
  
   v[AX1] = v[AX2] = 0.0;
   v[AX3] = - r*sin(x2)*Bx;
  #endif

  #if BACKGROUND_FIELD == YES  /* With background splitting, the initial 
                                  condition should contain deviation only */
   v[BX1] = v[BX2] = v[BX3] =
   v[AX1] = v[AX2] = v[AX3] = 0.0;
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

#if PHYSICS == MHD
/* ************************************************************** */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* 
 * Sets the initial curl-free magnetic field.
 **************************************************************** */
{
  double Bx;
  Bx = 5.0/sqrt(4.0*CONST_PI);

  #if GEOMETRY == CARTESIAN
   B0[0] = Bx;
   B0[1] = 0.0;
   B0[2] = 0.0;
  #elif GEOMETRY == POLAR
   B0[0] =   Bx*cos(x2);
   B0[1] = - Bx*sin(x2);
   B0[2] = 0.0;
  #endif
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 * Provide inner radial boundary condition in polar geometry.
 * Zero gradient is prescribed on density, pressure and magnetic field.
 * For the velocity, zero gradient is imposed on v/r (v = vr, vphi).
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *r, slp;

  if (side == X1_BEG){
    r = grid[IDIR].x;
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        slp = r[i]/r[IBEG];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[VX1][k][j][i] = slp*d->Vc[VX1][k][j][IBEG];
        d->Vc[VX2][k][j][i] = slp*d->Vc[VX2][k][j][IBEG]; 
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
      #endif
    }
  }
}
