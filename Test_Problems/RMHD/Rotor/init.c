/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief RMHD rotor test problem in 2 or 3 dimensions.

  The 2D rotor problem [dZBL03][Mig12] consists of a rapidly spinning disk
  embedded in a uniform background medium treaded by a costant magnetic field.
  The initial conditin reads as
  \f[
     \left(\rho, p, v_x, v_y, B_x\right) = \left\{\begin{array}{ll}
      \left(10,  1,-\omega y,\,\omega x,\, 1\right) & \quad\mathrm{for}\quad r < 0.1
      \\ \noalign{\medskip}
      \left(1,\, 1,\, 0,\, 0,\, 1\right) & \quad\mathrm{otherwise}
    \end{array}\right.
  \f]
  where \f$\omega = 9.95\f$ is the angular frequency of rotation.
  The computational domain is the unit square and outflow boundary conditions
  are imposed everywhere.

  As the disk rotates, strong torsional Alfven waves form and propagate 
  outward carrying angular momentum from the disk to the ambient.
  The emerging flow structure is enclosed by a circular fast forward shock
  traveling into the surrounding medium.
  An inward fast shock bounds the innermost oval region where density has
  been depleted to lower values.
  The presence of the magnetic field slows down the rotor, and the maximum
  Lorentz factor decreases from the nominal value of 10 to 2.2 (approx).

  A list of tested configurations is given in the following table:

  <CENTER>
  Conf.| GEOMETRY| DIM |T. STEPPING|RECONSTRUCTION| divB | AMR
  -----|---------|-----|-----------|--------------|------|-----
  #01  |CARTESIAN|  2  |  RK3      |  LINEAR      | CT   |  NO
  #02  |CARTESIAN|  2  |  RK2      |  LINEAR      | CT   |  NO
  #03  |CARTESIAN|  2  |  HANCOCK  |  LINEAR      | 8W   |  NO
  #04  |CARTESIAN|  2  |  HANCOCK  |  LINEAR      | CT   |  NO
  #05  |CARTESIAN|  3  |  RK2      |  LINEAR      | GLM  |  NO
  #06  |CARTESIAN|  3  |  RK2      |  LINEAR      | GLM  |  NO
  #07  |CARTESIAN|  3  |  HANCOCK  |  LINEAR      | CT   |  NO
  #08  |CARTESIAN|  3  |  RK3      |  LimO3       | GLM  |  NO
  </CENTER>

  The final solutin at t = 0.4 on a grid of 400x400 grid points is shown below.

  \image html rmhd_rotor.03.jpg "Density map (in log scale) at t = 0.4 in Cartesian coordinates using 400 x 400 grid zones (configuration #03). Field lines are super-imposed."

  The three dimensional version extends the current configuration to a
  spinning sphere and it is described in [MUB09].

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 6, 2015

  \b Reference: 
     - [dZBL03] "An efficient shock-capturing central-type scheme for multidimensional
        relativistic flows"
        Del Zanna, Bucciantini, Londrillo, A&A (2003) 400, 397
     - [MUB09] "A five-wave Harten-Lax-van Leer Riemann solver for relativistic
        magnetohydrodynamics"
        Mignone, Ugliano \& Bodo, MNRAS (2009) 393, 1141.
     - [Mig12] "The PLUTO Code for Adaptive Mesh Computations in Astrophysical Fluid Dynamics"
        Mignone et al., ApJS (2012) 198, 7 (see Sect. 6.3)
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
  double r, r0, Bx, omega;

  g_gamma   = 5./3.;
  omega = 0.995;
  Bx    = 1.0;
  r0    = 0.1;

  #if GEOMETRY == CARTESIAN

   r = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
   r = sqrt(r);

   v[RHO] = 1.0;
   v[VX1] = 0.0;
   v[VX2] = 0.0;
   v[PRS] = 1.0;
   v[BX1] = Bx;
   v[BX2] = 0.0;

   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX1] = -omega*x2/r0;
     v[VX2] =  omega*x1/r0;
   }

   v[AX1] = v[AX2] = 0.0;
   v[AX3] = v[BX1]*x2;

  #elif GEOMETRY == POLAR

   r = x1;

   v[RHO] = 1.0;
   v[VX2] = 0.0;
   v[VX1] = 0.0;
   v[BX1] =  Bx*cos(x2);
   v[BX2] = -Bx*sin(x2);
   v[PRS] = 1.0;
  
   if (r <= r0) {
     v[RHO] = 10.0;
     v[VX2] = omega*r/r0;
   }

   v[AX1] = v[AX2] = 0.0;
   v[AX3] = - r*sin(x2)*Bx;

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
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *r, slp;

  if (side == X1_BEG){
    r = grid[IDIR].x;
    if (box->vpos == CENTER){
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

