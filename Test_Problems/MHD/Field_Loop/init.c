/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advection of a magnetic field loop.

  This problem consists of a weak magnetic field loop being advected in a
  uniform velocity field. Since the total pressure is dominated by the
  thermal contribution, the magnetic field is essentially transported as
  a passive scalar.
  The preservation of the initial circular shape tests the scheme
  dissipative properties and the correct discretization balance of
  multidimensional terms. 

  Following [GS05][MT10][MTB10] (see also references therein),
  the computational box is defined by \f$ x\in[-1,1],\, y\in[-0.5,0.5] \f$
  discretized on \f$ 2N_y\times N_y\f$ grid cells (Ny=64).  
  Density and pressure are initially constant and equal to 1. 
  The velocity of the flow is given by
  \f[  
    \vec{v} = V_0(\cos\alpha, \sin\alpha)
  \f]
  with \f$V_0 = \sqrt{5},\,\sin \alpha = 1/\sqrt{5},\, \cos \alpha = 2/\sqrt{5}\f$.
  The magnetic field is defined through its  magnetic vector potential as 
  \f[   
    A_z = \left\{ \begin{array}{ll}
      A_0(R-r) & \textrm{if} \quad R_1 < r \leq R \,, \\ \noalign{\medskip}
      0   & \textrm{if} \quad r > R \,,
  \end{array} \right.
  \f]
  with \f$ A_0 = 10^{-3},\, R = 0.3,\, r = \sqrt{x^2+y^2}\f$.
  A slightly different variant is used for the finite difference schemes
  as explained in [MTB10]:
  \f[   
    A_z = \left\{ \begin{array}{ll}
    a_0 + a_2r^2 & \textrm{if} \quad 0 \leq r \leq R_1 \,, \\ \noalign{\medskip}  
    A_0(R-r) & \textrm{if} \quad R_1 < r \leq R \,, \\ \noalign{\medskip}
    0   & \textrm{if} \quad r > R \,,
  \end{array} \right.
  \f]     
  where \f$R_1=0.2 R,\, a_2 = -0.5A_0/R_1,\, a_0 = A_0(R-R_1) - a_2R_1^2\f$.

  Double periodic boundary conditions are imposed.

  A snapshot of the solution on a \c 128x64 grid at t=0.2 is shown below.

  \image html mhd_fl.03.jpg "Magetic pressure at t=0.2 (configuration #03)."

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 6, 2015

  \b References
     - [GS05] "An unsplit Godunov method for ideal MHD via constrained transport",
        Gardiner \& Stone JCP (2005) 205, 509
     - [MT10] "A second-order unsplit Godunov scheme for cell-centered MHD:
               The CTU-GLM scheme", Mignone \& Tzeferacos, JCP (2010) 229, 2117.

     - [MTB10] "High-order conservative finite difference GLM-MHD schemes for 
        cell-centered MHD", Mignone, Tzeferacos & Bodo, JCP (2010) 229, 5896.
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
  double c, s, v0=sqrt(5.0);
  double A=1.e-3, R=0.3, R1, r, Bphi;

  r = sqrt(x1*x1 + x2*x2);
  c = 2.0/sqrt(5.0);
  s = 1.0/sqrt(5.0);
/*
  c = 1.0/sqrt(2.0);
  s = 1.0/sqrt(2.0);
*/
  v[RHO] = 1.0;
  v[VX1] = v0*c;
  v[VX2] = v0*s;
  v[VX3] = 0.0;
  v[PRS] = 1.0;
  v[TRC] = 0.0;

#if PHYSICS == MHD || PHYSICS == RMHD

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = A*(R - r)*(r <= R);
  
  #ifdef FINITE_DIFFERENCE
  R1 = 0.2*R;  
  if (r < R1) {
    double a0, a2; 
    a2 = -A/(2.0*R1);
    a0 = A*(R-R1) - a2*R1*R1;
    v[AX3] = a0 + a2*r*r;
    A      = -2.0*a2*r;
  } else if (r <= R){
    v[AX3] = A*(R - r);
  } else {
    v[AX3] = 0.0;
    A      = 0.0;
  }
  #endif

  v[BX1] = -A*x2/r*(r <= R);  /*  =   dAz/dy  */
  v[BX2] =  A*x1/r*(r <= R);  /*  = - dAz/dx  */
  v[BX3] = 0.0;
  
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
 *********************************************************************** */
{ }

