/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Reconnection test (Harris sheet) in 2D.
 
  In this setup we reproduce a 2D Harris current sheet in the Hall MHD framework.
  The setup, where we make use of an Isothermal EoS, is the same as the 2D Harris current sheet provided in 
  /Test_Problems/MHD/Resistive_MHD/Current_Sheet/ (see also Section 5.3 of Mignone et al., ApJS (2012) 198:7).
   
  Configuration #01 includes Hall MHD while configuration #02 use classical MHD without resistivity (ideal case).
  Fig. 1 shows current density \f$ J_z\f$ and magnetic field lines for configuration #01 at \f$ t = 57 \f$. 
  A comparison between Hall and classical MHD is shown in fig. 2, where we plot the reconnected flux 
  \f$ \int_{0}^{L_x/2} B_y dx \f$ as a function of time.



  \image html hall_cs_fig1.png "Fig. 1 - The reconnected magnetic flux as a function of time for Hall and classical MHD for a Harris current sheet"
  
  \image html hall_cs_fig2.png "Fig. 2 - Pseudocolor rendering of the current density (color) with magnetic field lines at  t = 57 for the Hall MHD simulation of the Harris current sheet."

  \authors E. Striani (edoardo.striani@iaps.inaf.it)\n
           A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya (bhargav.vaidya@unito.it) 

 \date    May 09, 2017
  
  \b Reference
     - "The PLUTO Code for Adaptive Mesh Computations in
        Astrophysical FLuid Dynamics" Mignone et al., ApJS (2012) 198,7
     - "Geospace Environmental Modeling (GEM) magnetic
        reconnection challenge" Birn et al., JGR (2001) 106, 3715
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double x3)
/*
 *********************************************************************** */
{
  double cs2 = 0.5, b0 = 1.0, l, Psi0;
  double Lx, Ly, kx, ky;

  l = g_inputParam[WIDTH];
  v[RHO] = 0.2 + 1.0/(cosh(y/l)*(cosh(y/l)));

#if HAVE_ENERGY
  v[PRS] = cs2*v[RHO];  /* v[PRS] = 0.5 in the original PLUTO (2012) paper. */
#else 
  g_isoSoundSpeed = sqrt(cs2);
#endif

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

#if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = b0*tanh(y/l);
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR]; kx = CONST_PI/Lx;
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR]; ky = CONST_PI/Ly;

  Psi0    = g_inputParam[PSI0];
  v[BX1] += -Psi0*ky*sin(ky*y)*cos(2.0*kx*x);
  v[BX2] +=  Psi0*2.0*kx*sin(2.0*kx*x)*cos(ky*y);
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
/* *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* *********************************************************************** */
{
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/* *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* *********************************************************************** */
{
  return 0.0;
}
#endif


