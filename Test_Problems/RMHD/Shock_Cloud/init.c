/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic shock-cloud interaction.

  Sets initial condition for the relativistic shock-cloud interaction
  problem in 2D or 3D. 
  A high-density spherical cloud with \f$ \rho = 10 \f$ is initially
  located at \f$ (0.8,0) \f$ adjacent to a shock wave located at
  \f$ x=0 \f$.
  The post-shock region (\f$ x < 0 \f$) is at rest with
  \f[
     \begin{array}{lcll}
      (\rho,\, v_z,\, B_z,\, p) &=& (42.5942,\,0,\, -2.12971,\, 127.9483)
      & \qquad\mathrm{Ideal}\quad\mathrm{EOS, [MB06]}
     \\ \noalign{\medskip}
      (\rho,\, v_z,\, B_z,\, p) &=& (39.5052,\,0,\, 1.9753,\, 129.72386)
      & \qquad\mathrm{TM}\quad\mathrm{EOS, [Mig12]}
    \end{array}
  \f]
  while the upstream region is characterized by
  \f[
      (\rho,\, v_z,\, B_z,\, p) = (1, \, -\sqrt{0.99}, \,0.5,\, 10^{-3})
  \f]
  This test has no input parameters

  The first two configurations employ adaptive mesh refinement.

  \image html rmhd_shock_cloud.01.jpg "Density map (in log scale) for configuration #01"

  \b References
     - [MB06]: "An HLLC Riemann solver for relativistic flows -- 
       II. Magnetohydrodynamics", Mignone & Bodo MNRAS(2006) 368,1040
     - [Mig12]: "The PLUTO Code for Adaptive Mesh Computations in Astrophysical
       Fluid dynamics", Mignone et al., ApJS (2012), 198:7
       
  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Aug 31, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ************************************************************** */
void Init (double *us, double x1, double x2, double x3)
/* 
 *
 * 
 *
 *
 **************************************************************** */
{
  double x,y,z,r;
  double c,s,alpha;

  #if EOS == IDEAL
   g_gamma = 4./3.;
  #endif

  alpha = 0.0*CONST_PI/4.0;

  x = x1; y = x2; z = x3;

  if (x < 0.6){

    #if EOS == IDEAL
     us[RHO] = 42.5942;
     us[PRS] = 127.9483;
     us[BX3] = -2.12971;
    #elif EOS == TAUB
     us[RHO] = 39.5052;
     us[PRS] = 129.72386;
     us[BX3] = 1.9753;
    #endif
    us[VX1] = 0.0;
    us[VX2] = 0.0;
    us[VX3] = 0.0;
    us[BX1] = 0.0;
    us[BX2] = 0.0;

  }else{

    us[RHO] = 1.0;
    us[PRS] = 1.e-3;
    us[VX1] = -sqrt(1.0 - 1.0/100.);
    us[VX2] = 0.0;
    us[VX3] = 0.0;
    us[BX1] = 0.0;
    us[BX2] = 0.0;
    us[BX3] = 0.5;

  }
/*
  #if DIMENSIONS == 3
   us[BX3] /= sqrt(2.0);
   us[BX2]  = us[BX3];
  #endif
*/

  r = D_EXPAND( (x-0.8)*(x-0.8), + y*y, + z*z);
  r = sqrt(r);
  if (r < 0.15) us[RHO] = 10.0;

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

/* ************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *
 **************************************************************** */
{

}

/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 * 
 *
 **************************************************************** */
{
  int   i, j, k, nv;
  real  x1, x2, x3;

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    X2_BEG_LOOP(k,j,i){}
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    X2_END_LOOP(k,j,i){}
  }
}


